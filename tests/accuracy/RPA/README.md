# RPA accuracy / timing / memory tests

Test suite for RPA: self-consistent RPA (`RPA-OEP`, both all-electron and
pseudopotential) and non-self-consistent RPA (`RPA@GGA_PBE`, all-electron only),
across 5 elements (He, Be, Al, Ar, Au) -- 15 cases total.

## Running

```bash
python -m pytest tests/accuracy/RPA -m rpa_accuracy -v
```

Run from the repo root (SPARC-atomSFE folder, outside the tests directory),
in an environment with this package installed (`pip install -e .`) 
and pytest available.

**Strongly recommended: run with a 24-core allocation.** The stored reference
data (`reference.json` and each case's `result.json`) was generated at
`omega=24` (24 RPA frequency workers, one per core). Running with fewer cores
does not affect the hard accuracy checks -- total energy, RPA correlation
energy, exact exchange energy, and ionization potential (RPA-OEP only) 
are core-count-independent -- but it will make the timing/memory
checks to generate warnings, since wall time and peak memory both
scale with the number of RPA workers actually available. `conftest.py` warns
automatically if fewer than 24 cores are detected.
If the memory and timings are off by 20% of the reference results memory and timings,
the test will generate a warning. 
Only accuracy checks are used to determine if the tests passes or fails.

To run just a subset, select by element or by full case ID:

```bash
python -m pytest tests/accuracy/RPA -m rpa_accuracy -k He -v
python -m pytest "tests/accuracy/RPA/test_rpa_accuracy.py::test_rpa_case[RPA-OEP/AE/He]" -v
```

List every valid case ID:

```bash
python -m pytest tests/accuracy/RPA --collect-only -q
```

## What's checked

Each case is re-run fresh, in its own subprocess (`_run_one.py`, invoked via
`cases.run_subprocess_case`), and compared two ways:

- **Hard** (fails the test, 1e-5 absolute tolerance): total energy, RPA
  correlation energy, exact exchange energy, and -- self-consistent RPA
  only -- ionization potential (`-HOMO` eigenvalue). These four are the only
  quantities stored in `reference.json`.
- **Soft** (warns, never fails, 20% relative tolerance): wall time and peak
  RSS, read from each case's own `result.json` (not from `reference.json`).

A case with no `reference.json` entry, or no baseline `result.json`, is
skipped, not failed -- it means the data hasn't been produced yet.

## Files in this folder

- `cases.py` -- the 15 case definitions, plus the shared subprocess runner
  used by both the tests and the reference-data generator
- `_run_one.py` -- runs one case in isolation; reports total energy, RPA
  correlation energy, exact exchange energy, HOMO, wall time, peak RSS
- `conftest.py` -- warns once per session if fewer than 24 cores are available
- `test_rpa_accuracy.py` -- the actual pytest tests
- `reference.json` -- the four physics quantities per case; the only thing the
  hard checks compare against
- `RPA-OEP/<AE|psp>/<element>/result.json` and
  `RPA@GGA_PBE/AE/<element>/result.json` -- per-case baseline (also the source
  of the timing/memory soft-check baseline)

`cases.py` and `_run_one.py` are required to run the tests -- they are not
generation-only tools. The reference-data generator itself
(`generate_reference_data.py`), its SLURM launcher
(`generate_reference.sbatch`), a pytest SLURM launcher
(`run_rpa_accuracy_tests.sbatch`), and the raw `run.log` files for each case
are kept outside this folder (not needed to run the tests, only to regenerate
the reference data) -- ask whoever maintains this suite for their current
location if you need to add elements or regenerate after a code change.
