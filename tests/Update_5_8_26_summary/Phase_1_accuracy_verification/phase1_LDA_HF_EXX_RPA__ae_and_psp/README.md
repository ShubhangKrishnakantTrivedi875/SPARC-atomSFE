# Phase-1 validation harness

Compares three codes on the same physics, to check that the open-shell and
pseudopotential corrections in `SPARC-atomSFE` are right, and to quantify the
memory reduction in the RPA path.

| code key | root |
|---|---|
| `new`     | `SPARC-atomSFE` (corrected) |
| `orig`    | `SPARC-atomSFE_original/SPARC-atomSFE` |
| `oneatom` | `OneAtomFEM` (reference, MPI) |

**240 cases** = 3 codes x 5 atoms (Z = 2, 7, 25, 42, 79) x 2 electron treatments
(all-electron, pseudopotential) x 4 functionals (LDA, HF, EXX, RPA) x 2 settings.

| setting | N_fe | p | q | domain | omega | Lmax |
|---|---|---|---|---|---|---|
| `loose` | 15 | 20 | 45 | 13 Bohr | 8 | 4 |
| `tight` | 30 | 20 | 45 | 13 Bohr | 8 | 25 |

Polynomial mesh, concentration 2, SCF tolerance 1e-8 throughout.

## Layout

```
results/<code>/Z<zz>/<ae|psp>/<xc>/<setting>/
    params.json     exactly what was asked of the code
    run.log         that code's own terminal output, verbatim, plus a header
    result.json     energy, eigenvalues, timings, memory, exit status
    rss_trace.csv   RSS of the whole run vs wall time  (t_s, rss_total_mb, n_procs)
summary/            comparison tables written by collect.py
slurm_out/          SLURM stdout + one sacct dump per array task
```

## Running

```bash
sbatch submit_sparc.sbatch                  # 160 cases, ntasks=1  cpus-per-task=8
sbatch submit_oneatom.sbatch                #  80 cases, ntasks=8  cpus-per-task=1
sbatch --array=0-39 submit_sparc.sbatch     # cheap end only
```

Cases are ordered **cheapest first** by `cases.py`, so a partial array still
gives a usable table. `python cases.py sparc` lists index -> cid.

One case at a time, for debugging:

```bash
python run_case.py --cid new__Z07__ae__RPA__loose
python run_case.py --family sparc --index 12 --timeout 3600
```

Then:

```bash
python collect.py                       # all tables + terminal digest
python collect.py --xc RPA --setting loose
```

## Why the two sbatch files differ

The parallel axis is the same in both codes -- the imaginary-frequency loop --
but the mechanism is not. `SPARC-atomSFE` uses a thread pool, so its cores go
into `--cpus-per-task` and `RPA_N_WORKERS`. OneAtomFEM uses MPI ranks, so its
cores go into `--ntasks`. Both are set to 8, matching `omega = 8`, so one worker
or rank handles one frequency.

BLAS is pinned to one thread in both. That is not a limitation being worked
around: at these matrix sizes (n_quad = 675 loose, 1350 tight) multi-threaded
MKL measured **4-6x slower** for a full SCF, because `solve` and small `matmul`
lose more to thread barriers than they gain.

## Memory

`peak_rss_tree_mb` in `result.json` is the maximum over the run of the summed
RSS of every process belonging to the case, sampled at 10 Hz from the parent.
The MPI ranks are found by command-line match, not by process tree -- MVAPICH2
reparents them, and a tree-only sampler reports just the launcher (~16 MB).

`rss_trace.csv` is what localises the peak. For `xc = RPA` the maximum sits in
the RPA phase; the trace shows when, so the RPA contribution can be separated
from SCF setup rather than assumed.

`collect.py` also pulls `MaxRSS` from the per-task `sacct` dump as an
independent check on the sampled figure.

## Two gotchas worth knowing

**Interpreter.** OneAtomFEM's Pulay step inverts a 0x0 matrix on the first
iteration. scipy >= 1.11 returns an empty array; scipy 1.10 raises through MKL
and kills the run. The harness uses its own interpreter for children (override
with `PHASE1_PYTHON`), so submit from an environment where `python` is `krr_env`.

**Environment forwarding.** The `NU_*` overrides reach OneAtomFEM's ranks through
an explicit `env` wrapper on the `mpirun` line. Launchers do not reliably forward
arbitrary variables, and an unforwarded `NU_Z` means the run silently uses the
defaults compiled into the source file — it looks like a successful run with the
wrong system. If OneAtomFEM energies come back identical across different Z,
that is the cause.

## Not comparable without care

OneAtomFEM sets `double_hybrid_flag = 1` for `XC = 'RPA'`, so its RPA is run as a
double hybrid with a GGA_PBE sub-functional, whereas `SPARC-atomSFE` hardcodes
`double_hybrid_flag = False`. `alpha_x = alpha_c = 1` makes them agree in
principle; confirm against the energy components in `run.log` before treating an
RPA discrepancy as a bug.
