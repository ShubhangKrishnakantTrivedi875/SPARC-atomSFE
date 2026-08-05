"""
Runs ONE SPARC-atomSFE case and emits a JSON block on stdout.

Invoked as a child process by run_case.py so that a crash, an OOM or a
pinned-thread setting in one case cannot affect the harness or the next case.
The parent samples this process's RSS, so nothing is measured here except the
wall time of solve() itself.

    python -u _child_sparc.py <code_root> <params.json>

Everything before RESULT_JSON_BEGIN is the code's own terminal output and is
kept verbatim in the log.
"""

import json
import os
import sys
import time

code_root   = sys.argv[1]
params_path = sys.argv[2]

with open(params_path) as fh:
    params = json.load(fh)

# psps/ is resolved relative to the code root, so run from there
os.chdir(code_root)
sys.path.insert(0, code_root)

from src.solver import AtomicDFTSolver          # noqa: E402  (after chdir/path)

kwargs = dict(params)
kwargs.setdefault("verbose", True)
kwargs.setdefault("use_preconditioner", False)
kwargs["enable_parallelization"] = True
if not kwargs.get("all_electron_flag", True):
    kwargs.setdefault("psp_dir_path",  os.path.join(code_root, "psps"))
    kwargs.setdefault("psp_file_name", "%02d.psp8" % kwargs["atomic_number"])

result = {"ok": False, "error": None}
t_construct = t_solve = float("nan")
try:
    t0 = time.time()
    solver = AtomicDFTSolver(**kwargs)
    t_construct = time.time() - t0

    t0 = time.time()
    out = solver.solve(save_full_spectrum=False)
    t_solve = time.time() - t0

    occ = out["occupation_info"]
    eps = out["eigen_energies"]
    result.update(
        ok               = True,
        total_energy     = float(out["energy"]),
        n_states         = int(occ.n_states),
        eigenvalues      = [float(v) for v in eps],
        occupations      = [float(v) for v in occ.occupations],
        n_values         = [int(v)   for v in occ.n_values],
        l_values         = [int(v)   for v in occ.l_values],
        converged        = bool(out.get("converged", False)),
        scf_iterations   = int(out.get("iterations", -1)),
        outer_iterations = int(out.get("outer_iterations", -1) or -1),
    )
    comps = out.get("energy_components", None)
    if comps is not None:
        result["energy_components"] = {
            k: float(v) for k, v in vars(comps).items()
            if isinstance(v, (int, float))
        }
except BaseException as exc:                      # noqa: BLE001 - report, never mask
    import traceback
    result["error"] = "%s: %s" % (type(exc).__name__, exc)
    result["traceback"] = traceback.format_exc()

result["t_construct_s"] = t_construct
result["t_solve_s"]     = t_solve

print("\nRESULT_JSON_BEGIN")
print(json.dumps(result))
print("RESULT_JSON_END")
sys.stdout.flush()
sys.exit(0 if result["ok"] else 1)
