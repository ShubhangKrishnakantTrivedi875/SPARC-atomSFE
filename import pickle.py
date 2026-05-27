

import pickle
import sys
import json
import os
from types import SimpleNamespace
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from atomic_number_to_and_fro_symbol import element_names
from scipy.integrate import simpson

common_path = "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/PSEUDOPOTENTIAL_ACCURACY"
suffix1     = "psp_accuracy"
suffix2     = "SPMS_SG15"
json_path   = common_path + "/" + suffix1 + "/" + suffix2 + "/" + suffix1 +"_" + suffix2 +".json"   
savefig     = 1

Z_all = np.r_[np.arange(1,58), np.arange(72, 84)]
Z_all = np.array([6])

xc_functional_psp_type = [('PBE', 'SPMS'), ('rSCAN', 'SG15'),  ('PBE0_0.25', 'SPMS'), 
                          ('PBE0_0.5', 'SPMS'), ('PBE0_0.75', 'SPMS'), ('PBE0_1', 'SPMS')]

#xc_functional_psp_type = [('PBE', 'SPMS')]
REF = 'All_electron'
        
S = SimpleNamespace(atm = {})

# error messages for class LagrangeShapeFunctions
X_NODE_NOT_NUMPY_ARRAY_ERROR = \
    "x_node must be a numpy array, get {} instead."
X_NODE_NOT_AT_LEAST_2_ERROR = \
    "x_node must have length >= 2, get {} instead."
X_NODE_NOT_1D_OR_2D_ARRAY_ERROR = \
    "x_node must be a 1D or 2D numpy array, get dimension {} instead."
X_EVAL_NOT_NUMPY_ARRAY_ERROR = \
    "x_eval must be a numpy array, get {} instead."
X_EVAL_NOT_1D_OR_2D_ARRAY_ERROR = \
    "x_eval must be a 1D or 2D numpy array, get dimension {} instead."
X_EVAL_NOT_AT_LEAST_1_ERROR = \
    "x_eval must have length >= 1, get {} instead."
X_NODE_AND_X_EVAL_NOT_THE_SAME_SHAPE_ERROR = \
    "x_node and x_eval must have the same shape, get dimension {} and {} instead."


'''Function for evaluating the derivative of the orbitals (consistent with FEM approximation used within any element, not on boundaries)'''
def lagrange_basis_and_derivatives(x_node: np.ndarray, x_eval: np.ndarray, atol: float = 1e-20):
    ''' Adopted from SPARC-atomSFE code''' #https://github.com/SPARC-X/SPARC-atomSFE
    # basic checks
    assert isinstance(x_node, np.ndarray), \
        X_NODE_NOT_NUMPY_ARRAY_ERROR.format(type(x_node))
    assert isinstance(x_eval, np.ndarray), \
        X_EVAL_NOT_NUMPY_ARRAY_ERROR.format(type(x_eval))
    if x_node.ndim == 1:
        x_node = x_node[None, :]
    if x_eval.ndim == 1:
        x_eval = x_eval[None, :]
    assert x_node.ndim == 2, \
        X_NODE_NOT_1D_OR_2D_ARRAY_ERROR.format(x_node.ndim)
    assert x_eval.ndim == 2, \
        X_EVAL_NOT_1D_OR_2D_ARRAY_ERROR.format(x_eval.ndim)
    assert x_node.shape[0] == x_eval.shape[0], \
        X_NODE_AND_X_EVAL_NOT_THE_SAME_SHAPE_ERROR.format(x_node.shape[0], x_eval.shape[0])
    assert x_node.shape[1] >= 2, \
        X_NODE_NOT_AT_LEAST_2_ERROR.format(x_node.shape[1])
    assert x_eval.shape[1] >= 1, \
        X_EVAL_NOT_AT_LEAST_1_ERROR.format(x_eval.shape[1])

    # shapes
    n_elem, n_node = x_node.shape
    n_elem, n_eval = x_eval.shape

    # 1) Precompute barycentric weights: omega = 1 / c, where
    #    c_k = prod_{t≠k} (x_k - x_t), done per element.

    # diffs_nodes[i, k, t] = x_k - x_t (on element i)
    diffs_nodes = x_node[:, :, None] - x_node[:, None, :]                 # (n_elem, n_node, n_node)
    mask_offdiag = ~np.eye(n_node, dtype=bool)                            # (n_node, n_node) True for a≠b

    # c[i, k] = prod_{t≠k} (x_k - x_t): set diagonal to 1 so the product ignores it
    c = np.where(mask_offdiag[None, :, :], diffs_nodes, 1.0).prod(axis=2) # (n_elem, n_node)
    omega = 1.0 / c                                                       # (n_elem, n_node)

    # 2) Distances from evaluation points to nodes: dx = x - x_k
    dx = x_eval[:, :, None] - x_node[:, None, :]                          # (n_elem, n_eval, n_node)

    # 3) Identify nodal evaluations (x == x_k) using absolute tolerance
    is_nodal = np.isclose(x_eval[:, :, None], x_node[:, None, :], atol=atol, rtol=0.0)  # (n_elem, n_eval, n_node)
    has_nodal_row = is_nodal.any(axis=2)                                      # (n_elem, n_eval)

    # 4) Basis values at non-nodal points via barycentric formula
    #    L_k(x) = (omega_k /(x - x_k)) / sum_t (omega_t /(x - x_t))
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_dx = 1.0 / dx                                               # (n_elem, n_point, n_node)
        numer = omega[:, None, :] * inv_dx                              # (n_elem, n_point, n_node)
        denom = np.sum(numer, axis=2, keepdims=True)                    # (n_elem, n_point, 1)
        L = numer / denom                                               # (n_elem, n_point, n_node)

    # Enforce nodal rows to be exact one-hot (interpolation property)
    # First zero out rows that contain nodal matches; then set matching entries to 1.
    L = np.where(np.isfinite(L), L, 0.0)
    L = np.where(has_nodal_row[:, :, None], 0.0, L)
    L = np.where(is_nodal, 1.0, L)

    # 5) Derivatives at non-nodal points:
    #    L'_k(x) = L_k(x) * ( sum_t 1/(x - x_t) - 1/(x - x_k) )
    with np.errstate(divide='ignore', invalid='ignore'):
        harmonic_sum = np.sum(inv_dx, axis=2, keepdims=True)            # (n_elem, n_point, 1)
        dLdx = L * (harmonic_sum - inv_dx)                              # (n_elem, n_point, n_node)

    # 6) Nodal differentiation matrix D on the nodes {x_k} (per element)
    #    Off-diagonal: D[a,b] = (c_a / c_b) / (x_a - x_b)
    #    Diagonal:     D[a,a] = -sum_{b≠a} D[a,b]
    with np.errstate(divide='ignore', invalid='ignore'):
        D = (c[:, :, None] / c[:, None, :]) / diffs_nodes               # (n_elem, n_node, n_node)
    D = np.where(mask_offdiag[None, :, :], D, 0.0)
    D[..., np.arange(n_node), np.arange(n_node)] = -np.sum(D, axis=2)

    # Overwrite derivative rows where x coincides with a node:
    # For each (i, j) with a nodal match, find which node index k* and set dLdx[i, j, :] = D[i, k*, :]
    if has_nodal_row.any():
        nodal_index = np.argmax(is_nodal, axis=2)                       # (n_elem, n_point), argmax is safe since rows with no nodal are masked by has_nodal_row
        i_idx, j_idx = np.nonzero(has_nodal_row)                        # indices of rows to overwrite
        k_idx = nodal_index[i_idx, j_idx]
        dLdx[i_idx, j_idx, :] = D[i_idx, k_idx, :]

    # 7) Cleanup any residual non-finite entries (should be rare)
    L = np.where(np.isfinite(L), L, 0.0)
    dLdx = np.where(np.isfinite(dLdx), dLdx, 0.0)

    return L, dLdx


# Check in which finite-element does each rc_l exist
def fetch_ae_and_psps_orbital_coeffs_crucial(element_boundaries, FEM_nodes, interpolation_nodes, orbital_coefficients, ae_flag):
    
    Flag = 2
    
    if element_boundaries[0] == 0: 
        Flag = 0
       
    if ae_flag == 1:
        if element_boundaries[-1] == element_boundaries_ae[-1]:
            Flag = 1
    
    elif ae_flag == 0:
        if element_boundaries[-1] == element_boundaries_psps[-1]:
            Flag = 1
    
    else:
        raise ValueError("ae flag must be either 0 indicating pseudopotential variables or 1 indicating all electron variables")
    
    FEM_nodes_current_element_idx = np.argwhere((FEM_nodes >= element_boundaries[0]) & (FEM_nodes <= element_boundaries[-1]))[:,0]
    interpolation_nodes_current_element_idx = np.argwhere((interpolation_nodes >= element_boundaries[0]) & (interpolation_nodes <= element_boundaries[-1]))[:,0]
    
    FEM_nodes_current_element = FEM_nodes[FEM_nodes_current_element_idx]
    interpolation_nodes_current_element = interpolation_nodes[interpolation_nodes_current_element_idx]
    
    orbital_coefficients_current_element = orbital_coefficients[FEM_nodes_current_element_idx]
    
    basis_func, grad_basis_func = lagrange_basis_and_derivatives(FEM_nodes_current_element, interpolation_nodes_current_element)
    
    orbital_values = (basis_func@orbital_coefficients_current_element)[0,:]
    grad_orbital_values = (grad_basis_func@orbital_coefficients_current_element)[0,:]
    
    if Flag == 0:  
        log_derivative = grad_orbital_values[1:] / orbital_values[1:]
    elif Flag == 1:
        log_derivative = grad_orbital_values[:-1] / orbital_values[:-1]
    else:
        log_derivative = grad_orbital_values / orbital_values

    return interpolation_nodes_current_element_idx, orbital_values, grad_orbital_values, log_derivative
    

for atom in range(len(Z_all)):
    atomic_number = int(Z_all[atom])
    element = element_names(atomic_number)
    S.atm[atomic_number] = SimpleNamespace(element = element, xc_psp={})
    
    for xc_psp_idx in range(len(xc_functional_psp_type)):
        xc_functional, psp_type = xc_functional_psp_type[xc_psp_idx]
        S.atm[atomic_number].xc_psp[(xc_functional, psp_type)] = SimpleNamespace(nl={})
        
        D = S.atm[atomic_number].xc_psp[(xc_functional, psp_type)]
        
        ## READ/LOAD DATA
        load_path_psps = f'{common_path}/{atomic_number:02d}_{element}_data/{xc_functional}/{psp_type}/'
        load_path_ae = f'{common_path}/{atomic_number:02d}_{element}_data/{xc_functional}/{REF}/'
        
        psps_load = f'{atomic_number:02d}_{xc_functional}_{psp_type}.pkl'
        ae_load = f'{atomic_number:02d}_{xc_functional}_{REF}.pkl'
        
        with open(load_path_psps + psps_load, "rb") as file_psps:
            psps_results = pickle.load(file_psps)
            
        with open(load_path_ae + ae_load, "rb") as file_ae:
            ae_results = pickle.load(file_ae)


        ## Get all necessary variables
        '''Occupational information'''
        occ_n_psps = psps_results['occupation_info'].occ_n
        occ_l_psps = psps_results['occupation_info'].occ_l
        
        occ_n_ae   = ae_results['occupation_info'].occ_n
        occ_l_ae   = ae_results['occupation_info'].occ_l
        
        occ_e_psps = psps_results['occupation_info'].occ_spin_up_plus_spin_down
        occ_e_ae   = ae_results['occupation_info'].occ_spin_up_plus_spin_down
        
        psps_pairs = list(zip(occ_n_psps, occ_l_psps))
        ae_pairs   = list(zip(occ_n_ae  , occ_l_ae))
        
        ae_valence_idx = [ae_pairs.index(p) for p in psps_pairs]
        for i in range(len(psps_pairs)):
            D.nl[(psps_pairs[i][0], psps_pairs[i][1])] = SimpleNamespace(n = psps_pairs[i][0], l = psps_pairs[i][1])
        
        
        '''Pseudopotential related information'''
        rc_Vloc            = psps_results['pseudopotential_info'].r_grid_local[-1]
        rc_projectors_list = psps_results['pseudopotential_info'].r_cutoff_max_per_l
        
        unique_l1 = np.unique(occ_l_psps)
        for i in range(len(unique_l1)):
            n1 = occ_n_psps[np.argwhere(occ_l_psps == unique_l1[i])][:,0]
            for j in range(len(n1)):
                D.nl[n1[j], unique_l1[i]].rc_l = rc_projectors_list[unique_l1[i]]
        
       
        '''Grid related information'''
        uniform_grid = ae_results['uniform_grid'][:,None]
        r = uniform_grid
        FEM_nodes_ae   = ae_results['grid_data'].physical_nodes
        FEM_nodes_psps = psps_results['grid_data'].physical_nodes
        
        N_fe_ae   = ae_results['grid_data'].finite_element_number
        N_fe_psps = psps_results['grid_data'].finite_element_number
        
        Polynomial_order_ae   = (len(FEM_nodes_ae) - 1) / N_fe_ae
        Polynomial_order_psps = (len(FEM_nodes_psps) - 1) / N_fe_psps
        
        assert ( Polynomial_order_ae - int(Polynomial_order_ae) == 0)
        assert ( Polynomial_order_psps - int(Polynomial_order_psps) == 0)
        
        Polynomial_order_ae   = int(Polynomial_order_ae)
        Polynomial_order_psps = int(Polynomial_order_psps)
        
        element_boundaries_ae   = FEM_nodes_ae[::int(Polynomial_order_ae)]
        element_boundaries_psps = FEM_nodes_psps[::int(Polynomial_order_psps)]
        
        
        ''' ae valence and psps orbitals and eigenvalues'''
        ae_valence_eigenvalues = ae_results['eigen_energies'][ae_valence_idx]
        psps_eigenvalues       = psps_results['eigen_energies']
        
        ae_valence_orbitals_uniform_grid = ae_results['orbitals_on_uniform_grid'][:, ae_valence_idx]
        psps_orbitals_uniform_grid       = psps_results['orbitals_on_uniform_grid']
        
        #ae_valence_orbitals   = ae_results['orbitals'][:, ae_valence_idx] 
        #psps_valence_orbitals = psps_results['orbitals']
        
        ae_valence_orbital_coefficients = np.pad(ae_results['orbital_coefficients'][:, ae_valence_idx], ((1,1), (0,0)))
        psps_orbital_coefficients       = np.pad(psps_results['orbital_coefficients'], ((1,1), (0,0)))


        '''Obtain FEM nodes for the element in which rc_l resides'''    
        ae_element_idx_crucial    = np.zeros(len(rc_projectors_list), dtype = np.int64)
        psps_element_idx_crucial  = np.zeros(len(rc_projectors_list), dtype = np.int64)
        ae_FEM_node_idx_crucial   = np.zeros((len(rc_projectors_list), Polynomial_order_ae + 1), dtype = np.int64)
        psps_FEM_node_idx_crucial = np.zeros((len(rc_projectors_list), Polynomial_order_psps + 1), dtype = np.int64)
        ae_FEM_node_crucial       = np.zeros((len(rc_projectors_list), Polynomial_order_ae + 1))
        psps_FEM_node_crucial     = np.zeros((len(rc_projectors_list), Polynomial_order_psps + 1))
        for i in range(len(rc_projectors_list)):
            rc_l = rc_projectors_list[i]
            ae_element_idx_crucial[i]       = np.argwhere(element_boundaries_ae < rc_l)[-1,0]
            psps_element_idx_crucial[i]     = np.argwhere(element_boundaries_psps < rc_l)[-1,0]
            ae_FEM_node_idx_crucial[i, :]   = np.arange(ae_element_idx_crucial[i] * Polynomial_order_ae, (ae_element_idx_crucial[i] + 1) * Polynomial_order_ae + 1, dtype = np.int64)
            psps_FEM_node_idx_crucial[i, :] = np.arange(psps_element_idx_crucial[i] * Polynomial_order_psps, (psps_element_idx_crucial[i] + 1) * Polynomial_order_psps + 1, dtype = np.int64)
            
            ae_FEM_node_crucial[i, :]       = FEM_nodes_ae[ae_FEM_node_idx_crucial[i, :]]
            psps_FEM_node_crucial[i, :]     = FEM_nodes_psps[psps_FEM_node_idx_crucial[i, :]] 
    


        ## Check 1 (compare valence eigenvalues)
        diff_eigenvalues = np.zeros(len(ae_valence_idx))
        # Header
        print(f"{'n':>3} {'l':>3} {'eig_AE':>14} {'eig_PSP':>14} {'|diff|':>12}")
        print("-" * 51)
        # Rows
        for i in range(len(ae_valence_idx)):
            n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
            D.nl[n, l].psps_eigenvalues = psps_eigenvalues[i]
            D.nl[n, l].ae_eigenvalues = ae_valence_eigenvalues[i]
            
            diff_eigenvalues[i] = np.abs(psps_eigenvalues[i] - ae_valence_eigenvalues[i])
            D.nl[n, l].diff_eigenvalues = diff_eigenvalues[i]
            
            print(
                f"{occ_n_psps[i]:>3d} "
                f"{occ_l_psps[i]:>3d} "
                f"{ae_valence_eigenvalues[i]:>14.6f} "
                f"{psps_eigenvalues[i]:>14.6f} "
                f"{diff_eigenvalues[i]:>12.2e}"
            )
        print(f'Check 1: eigenvalues: max abs error = {np.max(np.abs(diff_eigenvalues)):.3e}')
        print(" ")
        print(" ")
        D.Check1_max_diff_eigenvalues = np.max(diff_eigenvalues)
        

        ## Check 2 compare valence orbitals for each index l beyond rc_l: (psp_orbital_{nl}(r))|r>r_c_l = (ae_orbital_{nl}(r))|r>r_c_l
        max_diff_orbitals  = np.zeros(len(ae_valence_idx))
        norm_diff_orbitals = np.zeros(len(ae_valence_idx))
        print("Check 2: valence orbitals beyond r_c")
        print(f"{'n':>3} {'l':>3} {'rc_l':>10} {'max|diff|':>14} {'rel_L2_err':>14}")
        print("-" * 48)
        
        for i in range(len(ae_valence_idx)):
            n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
            psps_orbital_current_uniform_grid = psps_orbitals_uniform_grid[uniform_grid[:,0]>rc_l,i]
            ae_orbital_current_uniform_grid   = ae_valence_orbitals_uniform_grid[uniform_grid[:,0]>rc_l,i]
            diff_minus = psps_orbital_current_uniform_grid - ae_orbital_current_uniform_grid
            diff_add   = psps_orbital_current_uniform_grid + ae_orbital_current_uniform_grid
            max_minus  = np.max(np.abs(diff_minus))
            max_add    = np.max(np.abs(diff_add))
            if max_minus <= max_add:
                max_diff_orbitals[i]  = max_minus
                norm_diff_orbitals[i] = np.linalg.norm(diff_minus, 2) / np.linalg.norm(ae_orbital_current_uniform_grid, 2)    
            else:
                max_diff_orbitals[i]  = max_add
                norm_diff_orbitals[i] = np.linalg.norm(diff_add, 2) / np.linalg.norm(ae_orbital_current_uniform_grid, 2)
            D.nl[n, l].max_diff_orbitals  = max_diff_orbitals[i]
            D.nl[n, l].norm_diff_orbitals = norm_diff_orbitals[i]
        
            print(f"{occ_n_psps[i]:>3d} {l:>3d} {rc_l:>10.5f} {max_diff_orbitals[i]:>14.2e} {norm_diff_orbitals[i]:>14.2e}")
        
        print("-" * 48)
        print(f"Check 2: orbitals : max (max (abs diff))   = {max(max_diff_orbitals):.2e}")
        print(f"Check 2: orbitals : max rel_L2_err         = {max(norm_diff_orbitals):.2e}")
        print(" ")
        print(" ")
        D.Check2_max_diff_orbitals  = np.max(max_diff_orbitals)
        D.Check2_norm_diff_orbitals = np.max(norm_diff_orbitals)


        ## Check 3 compare valence electron density beyond rc
        rc_dens = float(rc_projectors_list.max()) 
        mask_d = r[:,0] > rc_dens
        ae_valence_density_uniform_grid = (1 / (4 * np.pi)) * np.sum(occ_e_ae[ae_valence_idx] * (ae_valence_orbitals_uniform_grid[mask_d,:] / uniform_grid[mask_d,:]) ** 2, axis = 1)
        psps_density_uniform_grid       = (1 / (4 * np.pi)) * np.sum(occ_e_psps * (psps_orbitals_uniform_grid[mask_d,:] / uniform_grid[mask_d,:]) ** 2, axis = 1)
              # rc for density is set as max of rc_l, i.e. max out of all rc for pseudopotential orbitals
        max_diff_density  = float(np.max(np.abs(psps_density_uniform_grid - ae_valence_density_uniform_grid))) if mask_d.any() else np.nan
        norm_diff_density = np.linalg.norm(psps_density_uniform_grid - ae_valence_density_uniform_grid) / np.linalg.norm(ae_valence_density_uniform_grid)
        print(f"Check 3: density  : max abs error          = {max_diff_density:.3e}")
        print(f"Check 3: density  : relative L2-norm error = {norm_diff_density:.3e}")
        print(" ")
        print(" ")
        D.Check3_max_diff_density  = np.max(max_diff_density)
        D.Check3_norm_diff_density = np.max(norm_diff_density)
        

        ## Check 4 Check Orbital norm conservation for each index l within rc_l: int_{r=0}^{r = r_c_l} (psp_orbital_{nl}(r) * psp_orbital_{n'l}(r)) = integral (ae_orbital_{nl}(r) * ae_orbital_{n'l}(r))
        norm_check_rows = []
        occ_l_arr = np.asarray(occ_l_psps)
        unique_l = np.unique(occ_l_arr)
        max_relative_error_orbital_norm_conservation = 0.0
        for il in range(len(unique_l)):
                idx = np.argwhere(occ_l_arr == unique_l[il])[:, 0]
                n_shells = np.asarray(occ_n_psps)[idx]
                l = unique_l[il]
                psp_orbs_current_l = psps_orbitals_uniform_grid[:, idx]
                ae_valence_orbs_current_l = ae_valence_orbitals_uniform_grid[:, idx]
                
                rc_l = rc_projectors_list[l]
                mask = r[:,0] <= rc_l
                r_in = r[mask,0]
                for i in range(len(n_shells)):
                    D.nl[n_shells[i], l].psp_int = []
                    D.nl[n_shells[i], l].ae_int  = []
                    D.nl[n_shells[i], l].norm_conservation = []
                    for j in range(len(n_shells)):
                        if j < i: continue
                        psp_int = np.trapezoid(psp_orbs_current_l[mask, i] * psp_orbs_current_l[mask, j], r_in)
                        ae_int  = np.trapezoid(ae_valence_orbs_current_l[mask, i] * ae_valence_orbs_current_l[mask, j], r_in)
                        
                        max_relative_error_orbital_norm_conservation = max(max_relative_error_orbital_norm_conservation, np.abs(float(ae_int - psp_int)) / np.abs(float(ae_int)))
                        
                        D.nl[n_shells[i], l].psp_int.append((n_shells[j], psp_int))
                        D.nl[n_shells[i], l].ae_int.append((n_shells[j], ae_int))
                        D.nl[n_shells[i], l].norm_conservation.append((n_shells[j], np.abs(float(ae_int - psp_int)) / np.abs(float(ae_int))))
                        
                        norm_check_rows.append({
                            "l": l,
                            "ni": n_shells[i], "nj": n_shells[j],
                            "psp_int": float(psp_int),
                            "ae_int":  float(ae_int),
                            "diff":    np.abs(float(ae_int - psp_int)) / np.abs(float(ae_int)),
                        })
        
        print("Check 4: norm conservation inside r_c")
        print(f"{'l':>3} {'ni':>4} {'nj':>4} {'psp_int':>14} {'ae_int':>14} {'rel_diff':>12}")
        print("-" * 54)
        for row in norm_check_rows:
            print(
                f"{row['l']:>3d} "
                f"{row['ni']:>4d} "
                f"{row['nj']:>4d} "
                f"{row['psp_int']:>14.8f} "
                f"{row['ae_int']:>14.8f} "
                f"{row['diff']:>12.2e}"
            )
        print(f"Check 4: Orbital norm conservation within rc: max relative integral norm error: {max_relative_error_orbital_norm_conservation}")
        print(" ")
        print(" ")
        D.Check4_orbital_norm_conservation_error = max_relative_error_orbital_norm_conservation
        

        # ── Check 5: log-derivative at r = rc_l: using cubic spline fitting ──
        log_der_rows_spline = []   # list of dicts {n, l, logder_psp, logder_ae, diff}
        max_log_derivative_error_spline = 0.0
        ae_all_log_derivatives_spline = np.zeros((len(uniform_grid)-2, len(ae_valence_idx)))
        psps_all_log_derivatives_spline = np.zeros((len(uniform_grid)-2, len(ae_valence_idx)))
        
        ae_spline_fit_error = np.zeros(len(ae_valence_idx))
        psps_spline_fit_error = np.zeros(len(ae_valence_idx))
        
        spline_order = 3;   D.spline_order = spline_order
        for i in range(len(ae_valence_idx)):
                   n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
                   spl_psps = UnivariateSpline(uniform_grid[:,0], psps_orbitals_uniform_grid[:, i], s=0, k = spline_order)
                   psps_continuity_values = np.asarray(spl_psps.derivatives(rc_l))
                   psps_log_derivative_spline = psps_continuity_values[1] / psps_continuity_values[0]
                   D.nl[n, l].psps_log_derivative_spline = psps_log_derivative_spline
                   
                   spl_ae = UnivariateSpline(uniform_grid[:,0], ae_valence_orbitals_uniform_grid[:, i], s = 0, k = spline_order)
                   ae_continuity_values = np.asarray(spl_ae.derivatives(rc_l))
                   ae_log_derivative_spline = ae_continuity_values[1] / ae_continuity_values[0]
                   D.nl[n, l].ae_log_derivative_spline = ae_log_derivative_spline
                   
                   psps_all_ders_spline = np.asarray(spl_psps([uniform_grid[:,0]], nu = 1))[0,:]
                   psps_all_vals_spline = np.asarray(spl_psps([uniform_grid[:,0]], nu = 0))[0,:]
                   psps_all_log_derivatives_spline[:, i] = psps_all_ders_spline[1:-1] / psps_all_vals_spline[1:-1]
                   psps_spline_fit_error[i] = np.max(np.abs(psps_orbitals_uniform_grid[:, i] - psps_all_vals_spline))
                   assert (psps_spline_fit_error[i] < 1e-15)
             
                   ae_all_ders_spline = np.asarray(spl_ae([uniform_grid[:,0]], nu = 1))[0,:]
                   ae_all_vals_spline = np.asarray(spl_ae([uniform_grid[:,0]], nu = 0))[0,:]
                   ae_all_log_derivatives_spline[:, i] = ae_all_ders_spline[1:-1] / ae_all_vals_spline[1:-1]
                   ae_spline_fit_error[i] = np.max(np.abs(ae_valence_orbitals_uniform_grid[:, i] - ae_all_vals_spline))
                   assert (ae_spline_fit_error[i] < 1e-15)
                   
                   max_log_derivative_error_spline = max(max_log_derivative_error_spline, np.abs(psps_log_derivative_spline - ae_log_derivative_spline))
                    
                   log_der_rows_spline.append({
                       "n": n,
                       "l": l,
                       "rc_l": rc_l,
                       "psps_log_der_spline": psps_log_derivative_spline,
                       "ae_log_der_spline": ae_log_derivative_spline,
                       "log_der_error_spline":np.abs(psps_log_derivative_spline - ae_log_derivative_spline)
                       })
                   
                   # plt.figure()
                   # plt.plot(uniform_grid[1:-1,0], np.abs(psps_all_log_derivatives_spline[:, i] - ae_all_log_derivatives_spline[:, i]), label = 'diff of log derivatives')
                   # plt.yscale('log')
                   # plt.title(f'SPLINE: log-derivative abs-error (n={n}, l={l})')
                   # plt.show()
                   
                   
        print("\nCheck 5: log-derivative match at r_c using spline fitting")
        print(f"{'n':>3} {'l':>3} {'rc_l':>8} {'psps_log_derivative':>20} {'ae_log_derivative':>20} {'|diff|':>12}")
        print("-" * 71)
        for row in log_der_rows_spline:
            print(
                f"{row['n']:>3d} "
                f"{row['l']:>3d} "
                f"{row['rc_l']:>8.4f} "
                f"{row['psps_log_der_spline']:>20.6f} "
                f"{row['ae_log_der_spline']:>20.6f} "
                f"{row['log_der_error_spline']:>12.2e}"
            )
        print(f"Check 5: Max_abs_error_log_derivative_error_SPLINE: {max_log_derivative_error_spline}")
        print(" ")
        print(" ")
        D.Check5_max_abs_error_log_derivative_spline = max_log_derivative_error_spline


        # ── Check 6: log-derivative at r = rc_l using finite element derivate matrix multiplying orbital coefficients──
        log_der_rows_FEM = []   # list of dicts {n, l, logder_psp, logder_ae, diff}
        max_log_derivative_error_FEM = 0.0
        for i in range(len(ae_valence_idx)):
                   n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
                   
                   ae_FEM_nodes_current_element = ae_FEM_node_crucial[l, :]
                   assert ((rc_l > ae_FEM_nodes_current_element[0]) & (rc_l < ae_FEM_nodes_current_element[-1]))
             
                   ae_orbital_coefficients_current_element = ae_valence_orbital_coefficients[ae_FEM_node_idx_crucial[l, :], i]
                   basis_func_ae, derivative_basis_func_ae = lagrange_basis_and_derivatives(ae_FEM_nodes_current_element, np.array([rc_l]))
                   
                   ae_orbital_value_FEM      = basis_func_ae@ae_orbital_coefficients_current_element
                   ae_grad_orbital_value_FEM = derivative_basis_func_ae@ae_orbital_coefficients_current_element
                   
                   ae_log_derivative_FEM = ae_grad_orbital_value_FEM / ae_orbital_value_FEM
                   D.nl[n, l].ae_log_derivative_FEM = ae_log_derivative_FEM[0][0]
                   
                   psps_FEM_nodes_current_element = psps_FEM_node_crucial[l, :]
                   assert ((rc_l > psps_FEM_nodes_current_element[0]) & (rc_l < psps_FEM_nodes_current_element[-1]))
                   psps_orbital_coefficients_current_element   = psps_orbital_coefficients[psps_FEM_node_idx_crucial[l, :], i]
                   basis_func_psps, derivative_basis_func_psps = lagrange_basis_and_derivatives(psps_FEM_nodes_current_element, np.array([rc_l]))
                   
                   psps_orbital_value_FEM      = basis_func_psps@psps_orbital_coefficients_current_element
                   psps_grad_orbital_value_FEM = derivative_basis_func_psps@psps_orbital_coefficients_current_element
                   
                   psps_log_derivative_FEM = psps_grad_orbital_value_FEM / psps_orbital_value_FEM
                   D.nl[n, l].psps_log_derivative_FEM = psps_log_derivative_FEM[0][0]
                   
                   max_log_derivative_error_FEM = max(max_log_derivative_error_FEM, np.abs(psps_log_derivative_FEM - ae_log_derivative_FEM)[0][0])
                   log_der_rows_FEM.append({
                       "n": n,
                       "l": l,
                       "rc_l": rc_l,
                       "psps_log_der_FEM": psps_log_derivative_FEM[0][0],
                       "ae_log_der_FEM": ae_log_derivative_FEM[0][0],
                       "log_der_error_FEM":np.abs(psps_log_derivative_FEM - ae_log_derivative_FEM)[0][0]
                       })
                    
        print("\nCheck 6: log-derivative match at r_c using FEM basis functions")
        print(f"{'n':>3} {'l':>3} {'rc_l':>8} {'psps_log_derivative':>20} {'ae_log_derivative':>20} {'|diff|':>12}")
        print("-" * 71)
        for row in log_der_rows_FEM:
            print(
                f"{row['n']:>3d} "
                f"{row['l']:>3d} "
                f"{row['rc_l']:>8.4f} "
                f"{(row['psps_log_der_FEM']):>20.6f} "
                f"{(row['ae_log_der_FEM']):>20.6f} "
                f"{(row['log_der_error_FEM']):>12.2e}"
            )
        print(f"Check 6: Max_abs_error_log_derivative_error_FEM: {max_log_derivative_error_FEM}")
        print(" ")
        print(" ")
        D.Check6_max_abs_error_log_derivative_FEM = max_log_derivative_error_FEM



        # ── Check 7: log-derivative at all points using finite element derivate matrix multiplying orbital coefficients──
        log_der_rows = []   # list of dicts {n, l, logder_psp, logder_ae, diff}
        ae_orbital_values = np.zeros((len(uniform_grid[:,0]), len(ae_valence_idx)))
        ae_grad_orbital_values = np.zeros((len(uniform_grid[:,0]), len(ae_valence_idx)))
        ae_all_log_derivatives_FEM = np.zeros((len(uniform_grid[:,0]), len(ae_valence_idx)))
        psps_orbital_values = np.zeros((len(uniform_grid[:,0]), len(ae_valence_idx)))
        psps_grad_orbital_values = np.zeros((len(uniform_grid[:,0]), len(ae_valence_idx)))
        psps_all_log_derivatives_FEM = np.zeros((len(uniform_grid[:,0]), len(ae_valence_idx)))
        
        ae_basis_func_fitting_error = np.zeros(len(ae_valence_idx))
        psps_basis_func_fitting_error = np.zeros(len(ae_valence_idx))
        
        max_error_all_ae_log_derivatives = 0.0
        max_error_all_psps_log_derivatives = 0.0
        for i in range(len(ae_valence_idx)):
                   n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
                   
                   for j in range(len(element_boundaries_ae)-1):
        
                      interpolation_node_idx_current_element, orbital_values_curret_element, grad_orbital_values_current_element, log_derivative_current_element = fetch_ae_and_psps_orbital_coeffs_crucial(np.array([element_boundaries_ae[j], element_boundaries_ae[j+1]]), FEM_nodes_ae , uniform_grid[:, 0], ae_valence_orbital_coefficients[:, i], ae_flag=1)
                      ae_orbital_values[interpolation_node_idx_current_element, i] = orbital_values_curret_element
                      ae_grad_orbital_values[interpolation_node_idx_current_element, i] = grad_orbital_values_current_element
                      if j == 0:
                          ae_all_log_derivatives_FEM[interpolation_node_idx_current_element[1:], i] = log_derivative_current_element
                      elif j == len(element_boundaries_ae) - 2:
                          ae_all_log_derivatives_FEM[interpolation_node_idx_current_element[:-1], i] = log_derivative_current_element
                      else:
                          ae_all_log_derivatives_FEM[interpolation_node_idx_current_element, i] = log_derivative_current_element
        
                   ae_basis_func_fitting_error[i] = np.max(np.abs(ae_orbital_values[:, i] - ae_valence_orbitals_uniform_grid[:, i]))
                   assert (ae_basis_func_fitting_error[i] < 1e-15) 
                   
                   for j in range(len(element_boundaries_psps) - 1):
            
                      interpolation_node_idx_current_element, orbital_values_curret_element, grad_orbital_values_current_element, log_derivative_current_element = fetch_ae_and_psps_orbital_coeffs_crucial(np.array([element_boundaries_psps[j], element_boundaries_psps[j+1]]), FEM_nodes_psps , uniform_grid[:, 0], psps_orbital_coefficients[:, i], ae_flag=0)
                      psps_orbital_values[interpolation_node_idx_current_element, i] = orbital_values_curret_element
                      psps_grad_orbital_values[interpolation_node_idx_current_element, i] = grad_orbital_values_current_element
                      if j == 0:
                          psps_all_log_derivatives_FEM[interpolation_node_idx_current_element[1:], i] = log_derivative_current_element
                      elif j == len(element_boundaries_psps) - 2:
                          psps_all_log_derivatives_FEM[interpolation_node_idx_current_element[:-1], i] = log_derivative_current_element
                      else:
                          psps_all_log_derivatives_FEM[interpolation_node_idx_current_element, i] = log_derivative_current_element
                        
                     
                   psps_basis_func_fitting_error[i] = np.max(np.abs(psps_orbital_values[:, i] - psps_orbitals_uniform_grid[:, i]))
                   assert (psps_basis_func_fitting_error[i] < 1e-15)  
                   
                   max_error_all_ae_log_derivatives = max(max_error_all_ae_log_derivatives, np.max(np.abs(ae_all_log_derivatives_FEM[1:-1, i] - ae_all_log_derivatives_spline[:, i])))
                   max_error_all_psps_log_derivatives = max(max_error_all_psps_log_derivatives, np.max(np.abs(psps_all_log_derivatives_FEM[1:-1, i] - psps_all_log_derivatives_spline[:, i])))
                   
                   
                   log_der_rows.append({
                       "n": n,
                       "l": l,
                       "rc_l": rc_l,
                       "psps_log_der": psps_all_log_derivatives_FEM[1:-1, i],
                       "ae_log_der": ae_all_log_derivatives_FEM[1:-1, i],
                       "log_der_error":np.abs(psps_all_log_derivatives_FEM[1:-1, i] - ae_all_log_derivatives_FEM[1:-1, i])
                       })
                   
                   
                   # plt.figure()
                   # plt.plot(uniform_grid[1:-1,0], np.abs(psps_all_log_derivatives_FEM-ae_all_log_derivatives_FEM)[1:-1, i], label = 'diff of log derivatives')
                   # plt.yscale('log')
                   # #plt.xscale('log')
                   # plt.title(f'FEM: log-derivative abs-error (n={n}, l={l})')
                   # plt.show()
        
         
        print(f"Check 7 (verification against spline fitting way) Max abs error: psps = {max_error_all_psps_log_derivatives}, ae = {max_error_all_ae_log_derivatives}")           
        print(" ")
        print(" ")
        D.Check7_FEM_spline_consistency_all_points_ae   = (max_error_all_ae_log_derivatives)
        D.Check7_FEM_spline_consistency_all_points_psps = (max_error_all_psps_log_derivatives)



        # ── Check 8: log-derivative at r = rc_l: using cubic spline fitting of log ──
        log_der_rows_splinec = []   # list of dicts {n, l, logder_psp, logder_ae, diff}
        max_log_derivative_error_spline = 0.0
        
        for i in range(len(ae_valence_idx)):
                   n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
                   spl_psps = UnivariateSpline((uniform_grid[1:-1,0]), np.log(np.abs(psps_orbitals_uniform_grid[1:-1, i])), s=0, k = spline_order)
                   psps_continuity_values = np.asarray(spl_psps.derivatives((rc_l)))
                   psps_log_derivative_spline = psps_continuity_values[1]
                   
                   spl_ae = UnivariateSpline((uniform_grid[1:-1,0]), np.log(np.abs(ae_valence_orbitals_uniform_grid[1:-1, i])), s = 0, k = spline_order)
                   ae_continuity_values = np.asarray(spl_ae.derivatives((rc_l)))
                   ae_log_derivative_spline = ae_continuity_values[1]
                   
                   max_log_derivative_error_spline = max(max_log_derivative_error_spline, np.abs(psps_log_derivative_spline - ae_log_derivative_spline))
                    
                   log_der_rows_splinec.append({
                       "n": n,
                       "l": l,
                       "rc_l": rc_l,
                       "psps_log_der_spline": psps_log_derivative_spline,
                       "ae_log_der_spline": ae_log_derivative_spline,
                       "log_der_error_spline":np.abs(psps_log_derivative_spline - ae_log_derivative_spline)
                       })
                   
                  
        print("\nCheck 8: log-derivative match at r_c using spline fitting")
        print(f"{'n':>3} {'l':>3} {'rc_l':>8} {'psps_log_derivative':>20} {'ae_log_derivative':>20} {'|diff|':>12}")
        print("-" * 71)
        for row in log_der_rows_splinec:
            print(
                f"{row['n']:>3d} "
                f"{row['l']:>3d} "
                f"{row['rc_l']:>8.4f} "
                f"{row['psps_log_der_spline']:>20.6f} "
                f"{row['ae_log_der_spline']:>20.6f} "
                f"{row['log_der_error_spline']:>12.2e}"
            )
        print(f"Check 8: Max_abs_error_log_derivative_error_SPLINE: {max_log_derivative_error_spline}")
        print(" ")
        print(" ")



        print(f"Check 9: {spline_order}^th-order Spline fitting and basis Func interpolation error: ")
        for i in range(len(ae_valence_idx)):
            n = occ_n_psps[i]; l = occ_l_psps[i];
            
            print(f"  spline fitting error          : (n = {n}, l = {l})  psps = {psps_spline_fit_error[i]:.2e},  ae = {ae_spline_fit_error[i]:.2e}")
            print(f"  basis functions fitting error : (n = {n}, l = {l})  psps = {psps_basis_func_fitting_error[i]:.2e}  ae = {ae_basis_func_fitting_error[i]:.2e}")
 
        print(" ")
        print(" ")
        
        print(f"Check 10: basis Func interpolation at all points vs at rc_l consistency: ")
        for i in range(len(ae_valence_idx)):
            n = occ_n_psps[i]; l = occ_l_psps[i]; rc_l = rc_projectors_list[l]
            
            idx_on_uniform_grid = int(np.ceil(rc_l / [uniform_grid[1, 0] - uniform_grid[0, 0]])[0])
            
            ae_basis_func_consistency = ae_all_log_derivatives_FEM[idx_on_uniform_grid, i] - log_der_rows_FEM[i]['ae_log_der_FEM']
            psps_basis_func_consistency = psps_all_log_derivatives_FEM[idx_on_uniform_grid, i] - log_der_rows_FEM[i]['psps_log_der_FEM']
            
            print(f"  basis_func_consistency : (n = {n}, l = {l})  psps = {psps_basis_func_consistency:.2e},  ae = {ae_basis_func_consistency:.2e}")
       
        print(" ")
        print(" ")



'''MAKING THE PLOTS AND SAVING IN THE FOLDER OF EACH ATOMIC NUMBER'''
for atomic_number, elem_data in S.atm.items():
    element = elem_data.element
    
    checks = [
              ("Check1_max_diff_eigenvalues",                "Check 1: Eigenvalues"),
              ("Check2_max_diff_orbitals",                   "Check 2: Orbitals agreement beyond rc"),
              ("Check3_max_diff_density",                    "Check 3: Density agreemenet beyond rc"),
              ("Check4_orbital_norm_conservation_error",     "Check 4: Norm conservation within rc"),
              ("Check5_max_abs_error_log_derivative_spline", "Check 5: Log-derivative at rc using Spline"),
             ]
    

    for check, Title_check_name in (checks):
        
        values = []
        categories = []
        for (xc, psp), D_info in elem_data.xc_psp.items():
        
            categories.append(f"{xc}\n{psp}")
            values.append(getattr(D_info, check, None))
         
        figure = plt.figure()
        bars = plt.bar(categories, values)
        plt.bar_label(bars, fmt='%.2e', padding=3)
        plt.title("Z = "+ str(atomic_number) + " " + element + ", " +Title_check_name)
        #plt.xlabel("XC and psp_types")
        plt.margins(y=0.1)
        if check == "Check4_orbital_norm_conservation_error":
            plt.ylabel("Relative integral norm error")
        else:
            plt.ylabel("Max absolute error")

        if savefig == 1:
            plt.savefig(f'{common_path}/{atomic_number:02d}_{element}_data/{suffix1}/{suffix2}/{atomic_number:02d}_{Title_check_name}.pdf', bbox_inches = 'tight', pad_inches=1, dpi = 600)
         
        plt.show()
        plt.close(figure)



'''PRINTING TABLE'''
# ── Per-element summary table (read straight from S.atm[atomic_number].xc_psp) ──
ALL_XC_PSP_FOR_PRINT = [('PBE', 'SPMS'), ('rSCAN', 'SG15'),
                        ('PBE0_0.25', 'SPMS'), ('PBE0_0.5', 'SPMS'),
                        ('PBE0_0.75', 'SPMS'), ('PBE0_1', 'SPMS')]

CHECK_FIELDS = [
    ('Check 1', 'Check1_max_diff_eigenvalues'),
    ('Check 2', 'Check2_max_diff_orbitals'),
    ('Check 3', 'Check3_max_diff_density'),
    ('Check 4', 'Check4_orbital_norm_conservation_error'),
    ('Check 5', 'Check5_max_abs_error_log_derivative_spline'),
]

for atomic_number in S.atm:
    elem_data = S.atm[atomic_number]
    print()
    print("=" * 92)
    print(f"  Summary table  :                  atomic_number = {atomic_number}  ({elem_data.element})")
    print("=" * 92)

    header = f"  {'XC (PSP):':<18}"
    for label, _ in CHECK_FIELDS:
        header += f" {label:>12}"
    print(header)
    print("  " + "-" * 90)

    for xc, psp in ALL_XC_PSP_FOR_PRINT:
        if (xc, psp) not in elem_data.xc_psp:
            continue
        D_combo = elem_data.xc_psp[(xc, psp)]
        row = f"  {xc + ' (' + psp + '):':<18}"
        for _, attr in CHECK_FIELDS:
            val = getattr(D_combo, attr, None)
            row += f" {val:>12.3e}" if val is not None else f" {'—':>12}"
        print(row)
    print()
    
    
    
'''SAVING ALL ATOMS DATA TO .json FILE'''
# ── Persist S to a single JSON file (human-readable, key-addressable, mergeable) ──
# Schema:
#   { "atoms": { "<Z>": { "element": str,
#                          "xc_psp": { "<XC>__<PSP>": { "checks": {...},
#                                                       "nl": { "<n>_<l>": {...} } } } } } }
SAVE_PATH = json_path
def _to_json_safe(x):
    """Convert numpy scalars/arrays and tuple keys to JSON-friendly Python types."""
    if isinstance(x, dict):
        return {(f"{k[0]}_{k[1]}" if isinstance(k, tuple) else str(k)): _to_json_safe(v)
                for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_to_json_safe(v) for v in x]
    if hasattr(x, "tolist"):           # numpy arrays
        return _to_json_safe(x.tolist())
    if hasattr(x, "item") and not isinstance(x, (str, bytes)):
        try: return x.item()           # numpy scalars
        except Exception: pass
    return x

# Load existing data if file already exists, otherwise start fresh
if os.path.exists(SAVE_PATH):
    with open(SAVE_PATH, "r") as f:
        out = json.load(f)
else:
    out = {"atoms": {}}

# Merge in the elements processed by this run (overwrites stale entries for the same Z)
for atomic_number, elem_data in S.atm.items():
    z_key = str(int(atomic_number))
    atom_entry = {"element": elem_data.element, "xc_psp": {}}

    for (xc, psp), D_combo in elem_data.xc_psp.items():
        combo_key = f"{xc}__{psp}"

        # Top-level Check values
        checks = {
            "Check1_max_diff_eigenvalues":              getattr(D_combo, "Check1_max_diff_eigenvalues", None),
            "Check2_max_diff_orbitals":                 getattr(D_combo, "Check2_max_diff_orbitals", None),
            "Check2_norm_diff_orbitals":                getattr(D_combo, "Check2_norm_diff_orbitals", None),
            "Check3_max_diff_density":                  getattr(D_combo, "Check3_max_diff_density", None),
            "Check3_norm_diff_density":                 getattr(D_combo, "Check3_norm_diff_density", None),
            "Check4_orbital_norm_conservation_error":   getattr(D_combo, "Check4_orbital_norm_conservation_error", None),
            "Check5_max_abs_error_log_derivative_spline": getattr(D_combo, "Check5_max_abs_error_log_derivative_spline", None),
            "Check6_max_abs_error_log_derivative_FEM":    getattr(D_combo, "Check6_max_abs_error_log_derivative_FEM",    None),
            "Check7_FEM_spline_consistency_all_points_ae":   getattr(D_combo, "Check7_FEM_spline_consistency_all_points_ae",   None),
            "Check7_FEM_spline_consistency_all_points_psps": getattr(D_combo, "Check7_FEM_spline_consistency_all_points_psps", None),
            "spline_order":                              getattr(D_combo, "spline_order", None),
        }

        # Per-(n, l) details — dumps every attribute attached to D.nl[(n,l)]
        nl_dict = {}
        for (n, l), nl_data in D_combo.nl.items():
            nl_dict[f"{int(n)}_{int(l)}"] = {
                k: _to_json_safe(v) for k, v in vars(nl_data).items()
            }

        atom_entry["xc_psp"][combo_key] = {
            "checks": _to_json_safe(checks),
            "nl":     nl_dict,
        }

    out["atoms"][z_key] = atom_entry

with open(SAVE_PATH, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nWrote summary to {SAVE_PATH}")