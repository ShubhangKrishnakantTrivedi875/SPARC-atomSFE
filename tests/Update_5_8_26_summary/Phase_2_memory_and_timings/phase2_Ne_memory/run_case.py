"""
Run ONE phase-1 case and write its log, result and memory trace.

    python run_case.py --family sparc   --index 0
    python run_case.py --family oneatom --index 3
    python run_case.py --cid new__Z07__ae__RPA__loose

Output, all under results/<code>/Z<zz>/<ae|psp>/<xc>/<setting>/ :

    params.json     exactly what was asked of the code
    run.log         the code's own terminal output, verbatim
    result.json     energies, eigenvalues, timings, memory summary
    rss_trace.csv   RSS of the whole process tree vs wall time

Memory is sampled from the parent, so it covers the child (and, for
OneAtomFEM, every MPI rank) without instrumenting either code.  The trace is
what localises the peak: for XC = RPA the maximum sits in the RPA phase, and
rss_trace.csv shows when.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import resource
import subprocess
import sys
import threading
import time

import psutil

import cases as C

HERE       = C.HERE

# Interpreter used for the child runs.  Pinned explicitly because OneAtomFEM's
# Pulay step takes an inverse of a 0x0 matrix on the first iteration: scipy
# >= 1.11 returns an empty array, scipy 1.10 raises through MKL and kills the
# run.  Set PHASE1_PYTHON to override; otherwise the harness's own interpreter
# is used, which under the sbatch is the activated conda env.
PYTHON = os.environ.get("PHASE1_PYTHON", sys.executable)
SAMPLE_DT  = 0.005         # s between RSS samples
DISCOVER_DT = 0.5          # s between full rescans for new/reparented processes
KB, MB     = 1024.0, 1024.0 ** 2


# ---------------------------------------------------------------------------
#  RSS sampling over the whole child process tree
# ---------------------------------------------------------------------------
class RssSampler(threading.Thread):
    """
    Peak RSS of everything belonging to one case.

    The process tree alone is not enough: MVAPICH2 reparents the MPI ranks, so
    they are not descendants of the mpirun we spawned and a tree-only sampler
    reports just the launcher (~16 MB).  Any process of this user whose command
    line contains `match` is therefore counted too.
    """

    def __init__(self, pid, match=None):
        super().__init__(daemon=True)
        self.pid    = pid
        self.match  = match
        self.uid    = os.getuid()
        self.stop   = threading.Event()
        self.trace  = []                     # (t_since_start, total_mb, n_procs)
        self.peak_total_mb  = 0.0
        self.peak_single_mb = 0.0

    def _discover(self, root):
        """Tree walk plus a command-line scan, for ranks the tree does not reach."""
        try:
            procs = {p.pid: p for p in [root] + root.children(recursive=True)}
        except psutil.Error:
            procs = {}
        if self.match:
            for p in psutil.process_iter(["pid", "uids", "cmdline"]):
                try:
                    if p.info["uids"].real != self.uid or p.pid in procs:
                        continue
                    if any(self.match in (a or "") for a in (p.info["cmdline"] or [])):
                        procs[p.pid] = p
                except (psutil.Error, TypeError, AttributeError):
                    continue
        return procs

    def run(self):
        t0 = time.time()
        try:
            root = psutil.Process(self.pid)
        except psutil.Error:
            return
        # Discovery is the expensive part -- a full process_iter scan costs far more
        # than the sample interval.  So it runs every DISCOVER_DT and the pid list is
        # cached in between; sampling then only reads /proc/<pid>/statm per process.
        procs, t_discover = self._discover(root), time.time()
        while not self.stop.is_set():
            if time.time() - t_discover >= DISCOVER_DT:
                procs, t_discover = self._discover(root), time.time()
            total, single, n, dead = 0.0, 0.0, 0, None
            for pid, p in procs.items():
                try:
                    rss = p.memory_info().rss / MB
                except psutil.Error:
                    (dead := dead if dead is not None else []).append(pid)
                    continue
                total += rss
                single = max(single, rss)
                n += 1
            if dead:
                for pid in dead:
                    procs.pop(pid, None)
            if n:
                self.trace.append((time.time() - t0, total, n))
                self.peak_total_mb  = max(self.peak_total_mb, total)
                self.peak_single_mb = max(self.peak_single_mb, single)
            self.stop.wait(SAMPLE_DT)

    def write_trace(self, path):
        with open(path, "w") as fh:
            fh.write("t_s,rss_total_mb,n_procs\n")
            for t, m, n in self.trace:
                fh.write("%.2f,%.1f,%d\n" % (t, m, n))


# ---------------------------------------------------------------------------
#  Command construction
# ---------------------------------------------------------------------------
def sparc_cmd(case, params, outdir):
    params_path = os.path.join(outdir, "params.json")
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2, sort_keys=True)
    return [PYTHON, "-u", os.path.join(HERE, "_child_sparc.py"),
            C.CODE_ROOT[case["code"]], params_path], dict(os.environ)


def oneatom_cmd(case, params, outdir):
    """OneAtomFEM is driven through its NU_* environment overrides."""
    env = dict(os.environ)
    env.update(
        NU_Z      = str(case["Z"]),
        NU_XC     = params["xc_functional"],
        NU_AE     = "1" if params["all_electron_flag"] else "0",
        NU_NE     = str(params["finite_element_number"]),
        NU_P      = str(params["polynomial_order"]),
        NU_Q      = str(params["quadrature_point_number"]),
        NU_DOMAIN = str(params["domain_size"]),
        NU_MESH   = str(params["mesh_type"]),
        NU_CONC   = str(params["mesh_concentration"]),
        NU_SCFTOL = str(params["scf_tolerance"]),
    )
    if "frequency_quadrature_point_number" in params:
        env["NU_QOMEGA"] = str(params["frequency_quadrature_point_number"])
        env["NU_LMAX"]   = str(params["angular_momentum_cutoff"])
    with open(os.path.join(outdir, "params.json"), "w") as fh:
        json.dump({"params": params,
                   "env": {k: v for k, v in env.items() if k.startswith("NU_")}},
                  fh, indent=2, sort_keys=True)
    # ranks are only useful for the OEP functionals; LDA/HF collapse to one
    n_ranks = int(os.environ.get("SLURM_NTASKS", "1")) if C.USES_OEP[case["xc"]] else 1
    mpirun  = os.environ.get("MPIRUN", "mpirun")
    # The NU_* vars go through an explicit `env` wrapper rather than the parent
    # environment: mpirun launchers do not reliably forward arbitrary variables
    # to the ranks, and a silently unforwarded NU_Z means the run quietly uses
    # the defaults baked into the source file instead of this case.
    nu = ["%s=%s" % (k, v) for k, v in sorted(env.items()) if k.startswith("NU_")]
    return ([mpirun, "-n", str(n_ranks), "env"] + nu +
            [PYTHON, "-u", "Onescrpa_KS_solver_with_inversion.py"], env)


# ---------------------------------------------------------------------------
#  Output parsing
# ---------------------------------------------------------------------------
def parse_sparc(text):
    m = re.search(r"RESULT_JSON_BEGIN\s*\n(.*?)\nRESULT_JSON_END", text, re.S)
    if not m:
        return {"ok": False, "error": "no RESULT_JSON block in child output"}
    return json.loads(m.group(1))


_ONEATOM_EIG = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s*$")


def parse_oneatom(text):
    out = {"ok": False, "error": None}

    def grab(label):
        m = re.search(r"^\s*" + re.escape(label) + r"\s*:\s*(-?\d+\.\d+)\s*$",
                      text, re.M)
        return float(m.group(1)) if m else None

    E = grab("Total energy")
    if E is None:
        out["error"] = "no 'Total energy' line in output"
        return out
    out["ok"] = True
    out["total_energy"] = E
    for key, label in (("kinetic", "Total kinetic energy"),
                       ("external", "  External potential"),
                       ("hartree", "  Hartree"),
                       ("exact_exchange", "Exact exchange (HF)"),
                       ("rpa_correlation", "RPA correlation")):
        v = grab(label)
        if v is not None:
            out.setdefault("energy_components", {})[key] = v

    block = text.split("OCCUPIED EIGENVALUES (Hartree)")
    if len(block) > 1:
        eig, nv, lv = [], [], []
        for line in block[-1].splitlines():
            if line.startswith("====") and eig:
                break
            m = _ONEATOM_EIG.match(line)
            if m:
                nv.append(int(m.group(2)))
                lv.append(int(m.group(3)))
                eig.append(float(m.group(4)))
        out.update(eigenvalues=eig, n_values=nv, l_values=lv, n_states=len(eig))
    return out


# ---------------------------------------------------------------------------
#  Driver
# ---------------------------------------------------------------------------
def run(case, timeout_s=None):
    outdir = case["outdir"]
    os.makedirs(outdir, exist_ok=True)
    params = C.params_for(case)

    if case["code"] in C.SPARC_CODES:
        cmd, env = sparc_cmd(case, params, outdir)
        cwd, parse = C.CODE_ROOT[case["code"]], parse_sparc
        match = "_child_sparc.py"
    else:
        cmd, env = oneatom_cmd(case, params, outdir)
        cwd, parse = C.CODE_ROOT[case["code"]], parse_oneatom
        match = "Onescrpa_KS_solver_with_inversion.py"

    log_path = os.path.join(outdir, "run.log")
    header = (
        "=" * 78 + "\n"
        "  cid        : %s\n"
        "  code root  : %s\n"
        "  command    : %s\n"
        "  host       : %s\n"
        "  slurm job  : %s  task %s\n"
        "  started    : %s\n"
        "  params     : %s\n" % (
            case["cid"], cwd, " ".join(cmd), os.uname()[1],
            os.environ.get("SLURM_JOB_ID", "-"),
            os.environ.get("SLURM_ARRAY_TASK_ID", "-"),
            time.strftime("%Y-%m-%d %H:%M:%S"),
            json.dumps(params, sort_keys=True))
        + "=" * 78 + "\n\n")

    rusage_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    chunks, t0 = [], time.time()
    with open(log_path, "w") as log:
        log.write(header)
        log.flush()
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        sampler = RssSampler(proc.pid, match=match)
        sampler.start()
        try:
            for line in proc.stdout:
                chunks.append(line)
                log.write(line)
                log.flush()
            proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            proc.kill()
            chunks.append("\n*** killed: exceeded timeout of %s s ***\n" % timeout_s)
        finally:
            sampler.stop.set()
            sampler.join(timeout=2.0)
        wall = time.time() - t0
        log.write("\n" + "=" * 78 + "\n")
        log.write("  exit code  : %s\n  wall       : %.2f s\n" % (proc.returncode, wall))
        log.write("  peak RSS   : %.1f MB (tree)   %.1f MB (single proc)\n"
                  % (sampler.peak_total_mb, sampler.peak_single_mb))
        log.write("=" * 78 + "\n")

    text = "".join(chunks)
    res  = parse(text)

    rusage_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    res.update(
        cid          = case["cid"],
        code         = case["code"],
        Z            = case["Z"],
        electrons    = case["electrons"],
        xc           = case["xc"],
        setting      = case["setting"],
        params       = params,
        exit_code    = proc.returncode,
        wall_s       = wall,
        memory       = dict(
            peak_rss_tree_mb    = round(sampler.peak_total_mb, 1),
            peak_rss_single_mb  = round(sampler.peak_single_mb, 1),
            ru_maxrss_child_mb  = round(max(rusage_after - rusage_before, rusage_after) / KB, 1),
            n_samples           = len(sampler.trace),
            sample_dt_s         = SAMPLE_DT,
            discover_dt_s       = DISCOVER_DT,
        ),
        slurm        = dict(
            job_id     = os.environ.get("SLURM_JOB_ID"),
            array_task = os.environ.get("SLURM_ARRAY_TASK_ID"),
            ntasks     = os.environ.get("SLURM_NTASKS"),
            cpus       = os.environ.get("SLURM_CPUS_PER_TASK"),
        ),
        python       = PYTHON,
        threads      = {k: os.environ.get(k) for k in
                        ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                         "RPA_N_WORKERS", "RPA_BLAS_THREADS")},
        finished_at  = time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    with open(os.path.join(outdir, "result.json"), "w") as fh:
        json.dump(res, fh, indent=2, sort_keys=True)
    sampler.write_trace(os.path.join(outdir, "rss_trace.csv"))

    print("[%s] exit=%s wall=%.1fs peakRSS=%.0fMB E=%s"
          % (case["cid"], proc.returncode, wall, sampler.peak_total_mb,
             res.get("total_energy", "-")))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=("sparc", "oneatom"))
    ap.add_argument("--index", type=int)
    ap.add_argument("--cid")
    ap.add_argument("--timeout", type=float, default=None,
                    help="seconds; kill the case and record the partial log")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.cid:
        case = C.by_cid(a.cid)
    else:
        if a.family is None or a.index is None:
            ap.error("give --cid, or both --family and --index")
        pool = C.sparc_cases() if a.family == "sparc" else C.oneatom_cases()
        if a.list:
            for i, c in enumerate(pool):
                print(i, c["cid"])
            return 0
        if a.index >= len(pool):
            print("index %d beyond %d cases -- nothing to do" % (a.index, len(pool)))
            return 0
        case = pool[a.index]

    res = run(case, timeout_s=a.timeout)
    return 0 if res.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
