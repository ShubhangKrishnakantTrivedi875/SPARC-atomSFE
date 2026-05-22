import pickle
import sys
import numpy as np


atomic_number = 6
xc_functional = 'GGA_PBE'

psp_type1 = 'SPMS'
psp_type2 = 'all_electron'

load_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/first_check/'

psps_load = f'{atomic_number}_{xc_functional}_{psp_type1}.pkl'
ae_load = f'{atomic_number}_{xc_functional}_{psp_type2}.pkl'

with open(load_path + psps_load, "rb") as file_psps:
    psps_results = pickle.load(file_psps)
    
with open(load_path + ae_load, "rb") as file_ae:
    ae_results = pickle.load(file_ae)
    

## Get all necessary variables


occ_n_psps = psps_results['occupation_info'].occ_n
occ_l_psps = psps_results['occupation_info'].occ_l

occ_n_ae   = ae_results['occupation_info'].occ_n
occ_l_ae   = ae_results['occupation_info'].occ_l

occ_e_psps = psps_results['occupation_info'].occ_spin_up_plus_spin_down
occ_e_ae   = ae_results['occupation_info'].occ_spin_up_plus_spin_down

psps_pairs = list(zip(occ_n_psps, occ_l_psps))
ae_pairs   = list(zip(occ_n_ae  , occ_l_ae))

ae_valence_idx = [ae_pairs.index(p) for p in psps_pairs]

rc_Vloc = psps_results['pseudopotential_info'].r_grid_local[-1]
rc_projectors_list = psps_results['pseudopotential_info'].r_cutoff_max_per_l

uniform_grid = ae_results['uniform_grid'][:,None]
uniform_grid[0] = 1e-11

## Check 1 (compare valence eigenvalues)

ae_valence_eigenvalues = ae_results['eigen_energies'][ae_valence_idx]
psps_eigenvalues       = psps_results['eigen_energies']

diff_eigenvalues = np.zeros(len(ae_valence_idx))
for i in range(len(ae_valence_idx)):
    diff_eigenvalues[i] = psps_eigenvalues[i] - ae_valence_eigenvalues[i]
    print(f'(n,l) = ({occ_n_psps[i], occ_l_psps[i]}), eigenvalue_ae = {ae_valence_eigenvalues[i]}, eigenvalue_psps = {psps_eigenvalues[i]}, diff_eigenvalue = {diff_eigenvalues[i]}')
    

## Check 2 compare valence orbitals for each index l beyond rc_l: (psp_orbital_{nl}(r))|r>r_c_l = (ae_orbital_{nl}(r))|r>r_c_l

ae_valence_orbitals = ae_results['orbitals_on_uniform_grid'][:,ae_valence_idx]
psps_orbitals       = psps_results['orbitals_on_uniform_grid']


diff_orbitals = np.zeros(len(ae_valence_idx))
for i in range(len(ae_valence_idx)):
    l_index = occ_l_psps[i]
    rc_l = rc_projectors_list[l_index]
    print(f"rc_l for l = {l_index}, is {rc_l}")
    psps_orbital_current = psps_orbitals[:,i]
    ae_orbital_current = ae_valence_orbitals[:,i]
    diff_orbitals[i] = np.max(np.abs(psps_orbital_current[uniform_grid[:,0]>rc_l] + ae_orbital_current[uniform_grid[:,0] > rc_l]))
    
    w_trapz = np.full(len(uniform_grid[uniform_grid[:,0]<rc_l]), uniform_grid[2,0] - uniform_grid[1,0])
    w_trapz[0] = w_trapz[-1] = 0.5*w_trapz[1]
    psps_integral_norm = np.sum(w_trapz * psps_orbital_current[uniform_grid[:,0]<rc_l] ** 2)
    ae_integral_norm   = np.sum(w_trapz * ae_orbital_current[uniform_grid[:,0]<rc_l] ** 2)
    
    print(f'(n,l) = ({occ_n_psps[i], occ_l_psps[i]}), abs_diff_orbitals = {diff_orbitals[i]}')
    print(f'(n,l) = ({occ_n_psps[i], occ_l_psps[i]}), psps_integral_norm = {psps_integral_norm}, ae_integral_norm = {ae_integral_norm}, norm_diff = {ae_integral_norm - psps_integral_norm}')

## Check 3 compare valence electron density beyond rc

ae_valence_density = (1 / (4 * np.pi)) * np.sum(occ_e_ae[ae_valence_idx] * (ae_valence_orbitals / uniform_grid) ** 2, axis = 1)
psps_density       = (1 / (4 * np.pi)) * np.sum(occ_e_psps * (psps_orbitals / uniform_grid) ** 2, axis = 1)


## Check 4 Check Orbital norm conservation for each index l within rc_l: int_{r=0}^{r = r_c_l} (psp_orbital_{nl}(r) * psp_orbital_{n'l}(r)) = integral (ae_orbital_{nl}(r) * ae_orbital_{n'l}(r))

unique_l = sorted(set(occ_l_psps))

for l in unique_l:
    
    # indices (in the psps arrays) of all valence shells with this l
    shells_in_l = [i for i, l_i in enumerate(occ_l_psps) if l_i == l]
    
    rc_l = rc_projectors_list[l]
    mask = uniform_grid[:, 0] <= rc_l
    r_inside = uniform_grid[mask, 0]
    
    shell_labels = [(occ_n_psps[i], occ_l_psps[i]) for i in shells_in_l]
    print(f"\n--- l = {l}, rc_l = {rc_l:.4f}, shells in this channel: {shell_labels} ---")
    
    for i in shells_in_l:
        for j in shells_in_l:
            if j < i:
                continue   # symmetric, only need upper triangle
            
            psps_int = np.trapezoid(psps_orbitals[mask, i] * psps_orbitals[mask, j],     r_inside)
            ae_int   = np.trapezoid(ae_valence_orbitals[mask, i] * ae_valence_orbitals[mask, j], r_inside)
            
            ni, li = occ_n_psps[i], occ_l_psps[i]
            nj, lj = occ_n_psps[j], occ_l_psps[j]
            print(f"  <u_{ni}{li} | u_{nj}{lj}>  psps = {psps_int:+.6f},  "
                  f"ae = {ae_int:+.6f},  diff = {ae_int - psps_int:+.2e}")


import matplotlib.pyplot as plt
# %%
for i in range(len(ae_valence_idx)):
    n, l = occ_n_psps[i], occ_l_psps[i]
    rc_l = rc_projectors_list[l]
    
    plt.figure()
    plt.plot(uniform_grid[:, 0], ae_valence_orbitals[:, i], label='all-electron')
    plt.plot(uniform_grid[:, 0], psps_orbitals[:, i], '--', label='pseudo')
    plt.axvline(rc_l, color='k', linestyle=':', label=f'$r_c$ = {rc_l:.3f}')
    plt.xlabel('r (Bohr)')
    plt.ylabel(f'$u_{{{n}{l}}}(r)$')
    plt.title(f'(n, l) = ({n}, {l})')
    plt.legend()
    plt.tight_layout()
    plt.show()
# %%
