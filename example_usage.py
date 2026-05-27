from pathlib import Path

from src.solver import AtomicDFTSolver
import numpy as np
import pickle
from atomic_number_to_and_fro_symbol import element_names

import pandas as pd

SUBFOLDERS    = ['PBE', 'rSCAN', 'PBE0_0.25', 'PBE0_0.5', 'PBE0_0.75', 'PBE0_1', 'RPA']
#SUBSUBFOLDERS = ['All_electron', 'PseudoDojo_standard', 'PseudoDojo_stringent', 'SG15', 'SPMS']

SUBFOLDERS     = ['PBE0_1']
SUBSUBFOLDERS  = ['SPMS']

Z_all = np.r_[np.arange(21, 22)]
#Z_all = np.array([2])
#Z_all = np.array([22,23,24,25,26,27,28])

meta_folder = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/PSEUDOPOTENTIAL_ACCURACY'

# ─── once, near the top of your script ─────────────────────────────────────
eigenvalue_csv_dir = meta_folder      # ← edit
#os.makedirs(eigenvalue_csv_dir, exist_ok=True)

# ─── before the atomic_number loop (inside the psp_type loop) ──────────────
comparison_rows = []



run_SCF = 1
write_SCF_results_to_file = 1
reference_data_from_psps = 0
write_eigs_and_reference_eigs_to_csv = 0


"""Parse <INPUT> section of ONCVPSP .psp8 pseudopotential files."""
import pickle
from pathlib import Path

import numpy as np


_float = lambda s: float(s.upper().replace('D', 'E'))           # handles Fortran 'D' exponent


def parse_psp8(filepath):
    """Parse the <INPUT> section of an ONCVPSP .psp8 file.

    Returns dict with keys:
        'iexc'                  : int          (XC functional code)
        'lmax'                  : int          (max angular momentum)
        'ref_config'            : list[(n, l)] (reference electron configuration)
        'l', 'rc', 'ep'         : np.ndarray   (one per channel, length lmax+1)
        'n'                     : np.ndarray   (principal qn paired with each ep,
                                                = max n in ref_config sharing that l)
        'epsh1','epsh2','depsh' : float        (log-derivative analysis params)
    """
    text = Path(filepath).read_text()
    if '<INPUT>' not in text:
        raise ValueError(f'{filepath}: no <INPUT> section')
    lines = text.split('<INPUT>', 1)[1].split('</INPUT>', 1)[0].splitlines()
    out   = {}

    def next_data(start):
        """First non-blank, non-comment line .strip()ed, or None."""
        for j in range(start, len(lines)):
            sj = lines[j].strip()
            if sj and not sj.startswith('#'):
                return sj
        return None

    for i, ln in enumerate(lines):
        s = ln.strip()
        if not s.startswith('#'):
            continue
        toks = s.lstrip('#').replace(',', ' ').split()           # header tokens, comma-tolerant

        # ── iexc — fifth column after "# atsym z nc nv iexc psfile" ─────────────
        # in the iexc block, also store nc and nv:
        if 'atsym' in s and 'iexc' in s:
            d = next_data(i + 1)
            if d:
                f = d.split()
                out['nc']   = int(f[2])
                out['nv']   = int(f[3])
                out['iexc'] = int(f[4])

        # ── reference configuration — FIRST "# n l f …" block only ──────────────
        elif toks[:3] == ['n', 'l', 'f'] and 'ref_config' not in out:
            ref_cfg = []
            for j in range(i + 1, len(lines)):
                sj = lines[j].strip()
                if not sj or sj.startswith('#'):
                    break
                fj = sj.split()
                ref_cfg.append((int(fj[0]), int(fj[1])))
            out['ref_config'] = ref_cfg

        # ── lmax — first int after "# lmax" ─────────────────────────────────────
        elif s == '# lmax':
            d = next_data(i + 1)
            if d:
                out['lmax'] = int(d.split()[0])

        # ── l, rc, ep — read exactly (lmax+1) data lines after the header ───────
        elif toks[:3] == ['l', 'rc', 'ep']:
            n_chan = out['lmax'] + 1 if 'lmax' in out else None
            l_, rc_, ep_ = [], [], []
            for j in range(i + 1, len(lines)):
                sj = lines[j].strip()
                if not sj:
                    continue
                if sj.startswith('#'):
                    break
                f = sj.split()
                l_.append(int(f[0])); rc_.append(_float(f[1])); ep_.append(_float(f[2]))
                if n_chan is not None and len(l_) == n_chan:
                    break
            out['l']  = np.array(l_)
            out['rc'] = np.array(rc_)
            out['ep'] = np.array(ep_)

        # ── epsh1, epsh2, depsh — three floats after "# epsh1 epsh2 depsh" ──────
        elif 'epsh1' in s:
            d = next_data(i + 1)
            if d:
                f = d.split()
                out['epsh1'], out['epsh2'], out['depsh'] = _float(f[0]), _float(f[1]), _float(f[2])

    # ── derive n per ep channel: max n in ref_config sharing that l ─────────────
    # replace the bottom "derive n" block with this:
    if 'ref_config' in out and 'l' in out and 'nv' in out:
        valence = out['ref_config'][-out['nv']:]                  # last nv entries
        def _n_for_l(li):
            v_ns = [n for n, lc in valence if lc == li]
            if v_ns:
                return min(v_ns)                                  # deepest pseudized state at this l
            all_ns = [n for n, lc in out['ref_config'] if lc == li]
            return min(all_ns) if all_ns else -1                  # fallback for unbound channels
        out['n'] = np.array([_n_for_l(li) for li in out['l']])
    # ── sanity check ────────────────────────────────────────────────────────────
    if 'lmax' in out and 'l' in out and len(out['l']) != out['lmax'] + 1:
        print(f'WARNING  {filepath}: parsed {len(out["l"])} channels but lmax+1 = {out["lmax"]+1}')

    return out


# def dump_folder(psp_file, psps_type, out_dir):
#     """Parse every '<Z>.psp8' in `folder` and pickle each dict as
#     `<out_dir>/<Z>_<psps_type>.pkl`."""
#     folder, out_dir = Path(folder), Path(out_dir)
#     #out_dir.mkdir(parents=True, exist_ok=True)
#     for psp_file in sorted(folder.glob('*.psp8')):
#         atomic_number = psp_file.stem                              # '06', '79', …
        
#         with (out_dir / f'{atomic_number}_{psps_type}.pkl').open('wb') as f:
#             pickle.dump(data, f)

def _decimal_tol(x, cap=6):
    """tol = 10^(-n) where n = # decimals in str(x), capped at `cap`."""
    s = f'{abs(x):.{cap}f}'.rstrip('0').rstrip('.')
    return 10.0 ** -(len(s.split('.')[1]) if '.' in s else 0)


for i in range(0,len(Z_all)):
    atomic_number = int(Z_all[i])
    element = element_names(atomic_number)
    
    for j in range(len(SUBFOLDERS)):
        XC = SUBFOLDERS[j]
        
        hybrid_mixing_parameter = 0
        if XC == 'PBE':
            xc_functional = 'GGA_PBE'
        
        elif XC == 'PBE0_0.25':
            xc_functional = 'PBE0'
            hybrid_mixing_parameter = 0.25
            
        elif XC == 'rSCAN':
            xc_functional = 'RSCAN'
            
        elif XC == 'PBE0_0.5':
            xc_functional = 'PBE0'
            hybrid_mixing_parameter = 0.5
            
        elif XC == 'PBE0_0.75':
            xc_functional = 'PBE0'
            hybrid_mixing_parameter = 0.75
            
        elif XC == 'PBE0_1':
            xc_functional = 'PBE0'
            hybrid_mixing_parameter = 1
            
        ## NEED TO LOOK UP HOW TO RUN PBE0 with variable alpha parameter
        
        for k in range(len(SUBSUBFOLDERS)):
            psp_type = SUBSUBFOLDERS[k]
            
            if psp_type=='All_electron':
                all_electron_flag, finite_element_number, mesh_concentration = (True, 12, 101)
                psp_dir_path = None
                psp_file_name = None
                
            else:
                all_electron_flag, finite_element_number, mesh_concentration  = (False, 10, 101)
                if psp_type == 'SPMS':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/SPMS_PBE'
                elif psp_type == 'SG15':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/SG15'
                elif psp_type == 'PseudoDojo_standard':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/PseudoDojo/standard/nc-sr-05_pbe_standard_psp8'
    
                elif psp_type == 'PseudoDojo_stringent':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/PseudoDojo/stringent/nc-sr-05_pbe_stringent_psp8'
                
                
                if atomic_number < 10:
                    psp_file_name = "0"+str(atomic_number)+".psp8"
                else:
                    psp_file_name = str(atomic_number)+".psp8"
            
            
            
            if run_SCF == 1:
                solver = AtomicDFTSolver(
                    atomic_number             = atomic_number,
                    xc_functional             = xc_functional,
                    domain_size               = 40.0,
                    finite_element_number     = finite_element_number,
                    polynomial_order          = 20,
                    quadrature_point_number   = 60,
                    mesh_type                 = "exponential",
                    mesh_concentration        = mesh_concentration,
                    scf_tolerance             = 1e-11,
                    mesh_spacing              = 0.001,
                    hybrid_mixing_parameter   = float(hybrid_mixing_parameter),
                    verbose                   = True, 
                    all_electron_flag         = all_electron_flag,
                    use_oep                   = False,
                    use_preconditioner        = True,
                    pulay_mixing_frequency    = 2,
                    pulay_mixing_history      = 11,
                    psp_dir_path              = psp_dir_path,
                    psp_file_name             = psp_file_name,
                )

                results = solver.solve()
            #print("eigen_energies:", results['eigen_energies'])
            write_to_path = meta_folder +"/"+ f'{atomic_number:02d}_{element}_data'+"/"+XC+"/"+psp_type+"/"
                
            if write_SCF_results_to_file == 1:
                filename = f'{atomic_number:02d}_{XC}_{psp_type}_{finite_element_number}_{mesh_concentration}.pkl'
                with open(write_to_path + filename,"wb") as file:
                    pickle.dump(results, file)

            print(results['energy'])
            if reference_data_from_psps == 1:
                if not all_electron_flag and XC == 'PBE':
                    reference_psps_data = parse_psp8(psp_dir_path + "/" + psp_file_name)
                    print('ref_eigen:', reference_psps_data['ep'])
                    print('run_eigen:', results['eigen_energies'])

                    # ── match run eigenvalues to reference by (n, l); NaN where missing ────
                    occ_n_psps = results['occupation_info'].occ_n
                    occ_l_psps = results['occupation_info'].occ_l
                    run_nl_map = {(int(n), int(l)): i
                                for i, (n, l) in enumerate(zip(occ_n_psps, occ_l_psps))}

                    ref_n  = reference_psps_data['n'].astype(int)
                    ref_l  = reference_psps_data['l'].astype(int)
                    ref_nl = list(zip(ref_n, ref_l))
                    ref_ep = reference_psps_data['ep']

                    run_ep = np.array([results['eigen_energies'][run_nl_map[nl]] if nl in run_nl_map else np.nan
                                    for nl in ref_nl])
                    eigen_diff = run_ep - ref_ep

                    missing = [nl for nl in ref_nl if nl not in run_nl_map]
                    if missing:
                        print(f"WARNING  {atomic_number}: ref (n,l) pairs {missing} not in run "
                            f"(run has {sorted(run_nl_map.keys())}) — recorded as NaN")

                    if run_SCF == 1:
                        valid = ~np.isnan(eigen_diff)
                        tols  = np.array([_decimal_tol(v) for v in ref_ep])
                        matched_nl = [nl for nl, v in zip(ref_nl, valid) if v]
                        max_err    = np.nanmax(np.abs(eigen_diff)) if valid.any() else np.nan
                        print(f"{atomic_number}, matched nl = {matched_nl}")
                        print(f"{atomic_number}, eigenvalues_error = {max_err:.2e}, "
                            f"per-channel tol = {tols}")
                        assert np.all(np.abs(eigen_diff[valid]) < tols[valid]), \
                            f"|diff| = {np.abs(eigen_diff)} exceeds per-channel tol = {tols}"

                    save_reference_to_filename = f'{atomic_number:02d}_{psp_type}_reference.pkl'
                    with open(write_to_path + "/" + save_reference_to_filename, "wb") as f:
                        pickle.dump(reference_psps_data, f)

                    # ── accumulate one row for the CSV using (n,l)-matched arrays ──────────
                    N_L = 4
                    pad = lambda a: list(a) + [np.nan] * (N_L - len(a))
                    comparison_rows.append({
                        'Z':            atomic_number,
                        **{f'n_l{i}':    v for i, v in enumerate(pad(ref_n))},
                        **{f'ref_l{i}':  v for i, v in enumerate(pad(ref_ep))},
                        **{f'run_l{i}':  v for i, v in enumerate(pad(run_ep))},
                        **{f'diff_l{i}': v for i, v in enumerate(pad(eigen_diff))},
                        'max_abs_diff': float(np.nanmax(np.abs(eigen_diff))) if np.any(~np.isnan(eigen_diff))
                                        else np.nan,
                    })

# ─── after the atomic_number loop completes (still inside psp_type / XC) ───
if write_eigs_and_reference_eigs_to_csv == 1:
    csv_out = f'{eigenvalue_csv_dir}/{XC}_{psp_type}_eigenvalues_comparison.csv'
    pd.DataFrame(comparison_rows).sort_values('Z').to_csv(csv_out, index=False)
    print(f'Wrote {csv_out}')
