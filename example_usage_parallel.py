from pathlib import Path

from src.solver import AtomicDFTSolver
import numpy as np
import pickle
from atomic_number_to_and_fro_symbol import element_names
from multiprocessing import Pool

import pandas as pd

SUBFOLDERS    = ['PBE', 'rSCAN', 'PBE0_0.25', 'PBE0_0.5', 'PBE0_0.75', 'PBE0_1', 'RPA']
#SUBSUBFOLDERS = ['All_electron', 'PseudoDojo_standard', 'PseudoDojo_stringent', 'SG15', 'SPMS']

SUBFOLDERS     = ['PBE0_0.5']
SUBSUBFOLDERS  = ['SPMS']

Z_all = np.r_[np.arange(1,58), np.arange(72, 84)]
Z_all = np.array([40, 41, 42, 72, 73, 74])
#Z_all = np.array([22,23,24,25,26,27,28])

N_WORKERS = 8   # set to number of CPUs you want to use

meta_folder = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/PSEUDOPOTENTIAL_ACCURACY'

run_SCF = 1
write_SCF_results_to_file = 1


def run_one(atomic_number):
    atomic_number = int(atomic_number)
    element = element_names(atomic_number)

    for XC in SUBFOLDERS:
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

        for psp_type in SUBSUBFOLDERS:
            if psp_type == 'All_electron':
                all_electron_flag, finite_element_number, mesh_concentration = (True, 12, 101)
                psp_dir_path = None
                psp_file_name = None
            else:
                all_electron_flag, finite_element_number, mesh_concentration = (False, 10, 20)
                if psp_type == 'SPMS':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/SPMS_PBE'
                elif psp_type == 'SG15':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/SG15'
                elif psp_type == 'PseudoDojo_standard':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/PseudoDojo/standard/nc-sr-05_pbe_standard_psp8'
                elif psp_type == 'PseudoDojo_stringent':
                    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/PseudoDojo/stringent/nc-sr-05_pbe_stringent_psp8'

                psp_file_name = f"{atomic_number:02d}.psp8"

            if run_SCF == 1:
                solver = AtomicDFTSolver(
                    atomic_number           = atomic_number,
                    xc_functional           = xc_functional,
                    domain_size             = 40.0,
                    finite_element_number   = finite_element_number,
                    polynomial_order        = 20,
                    quadrature_point_number = 60,
                    mesh_type               = "exponential",
                    mesh_concentration      = mesh_concentration,
                    scf_tolerance           = 1e-11,
                    mesh_spacing            = 0.001,
                    hybrid_mixing_parameter = float(hybrid_mixing_parameter),
                    verbose                 = True,
                    all_electron_flag       = all_electron_flag,
                    use_oep                 = False,
                    use_preconditioner      = False,
                    psp_dir_path            = psp_dir_path,
                    psp_file_name           = psp_file_name,
                )
                results = solver.solve()

            write_to_path = meta_folder + "/" + f'{atomic_number:02d}_{element}_data' + "/" + XC + "/" + psp_type + "/"

            if write_SCF_results_to_file == 1:
                filename = f'{atomic_number:02d}_{XC}_{psp_type}_{finite_element_number}_{mesh_concentration}.pkl'
                with open(write_to_path + filename, "wb") as file:
                    pickle.dump(results, file)

            print(f"Z={atomic_number:3d}  energy={results['energy']:.10f}", flush=True)


if __name__ == '__main__':
    with Pool(N_WORKERS) as pool:
        pool.map(run_one, Z_all)
