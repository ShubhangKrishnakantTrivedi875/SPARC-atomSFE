
from src.solver import AtomicDFTSolver
import numpy as np
import pickle


atomic_number = 6
xc_functional = 'GGA_PBE'
all_electron_flag = True


finite_element_number = 12
mesh_concentration = 100.0
psp_type = "all_electron"
psp_dir_path = None
psp_file_name = None
if all_electron_flag == False:
    psp_type = 'SPMS'
    psp_dir_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/'
    if atomic_number < 10:
        psp_file_name = "0"+str(atomic_number)+".psp8"
    else:
        psp_file_name = str(atomic_number)+".psp8"
        
    finite_element_number = 10
    mesh_concentration = 20.0
        
save_path = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/first_check/'

solver = AtomicDFTSolver(
    atomic_number             = atomic_number,
    xc_functional             = xc_functional,
    domain_size               = 40.0,
    finite_element_number     = finite_element_number,
    polynomial_order          = 20,
    quadrature_point_number   = 60,
    mesh_type                 = "exponential",
    mesh_concentration        = mesh_concentration,
    scf_tolerance             = 1e-10,
    mesh_spacing              = 0.001,
    
    verbose                   = True, 
    all_electron_flag         = all_electron_flag,
    use_oep                   = False,
    use_preconditioner        = False,
    psp_dir_path              = psp_dir_path,
    psp_file_name             = psp_file_name,
    
   # psp_dir_path = "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE",
   # psp_file_name= "02_He_2_1.1_1.2_pbe_v1.0.psp8"
)

results = solver.solve()

filename = f'{atomic_number}_{xc_functional}_{psp_type}.pkl'
print("eigenvalues", (results['eigen_energies']))
with open(save_path + filename,"wb") as file:
    pickle.dump(results, file)
    

# print("\nSolver results:")
# print(f"\t total energy             = {results['energy']:.6f} (Ha)")
# print(f"\t density.shape            = {results['rho'].shape}")
# print(f"\t quadrature_nodes.shape   = {results['quadrature_nodes'].shape}")
# print(f"\t quadrature_weights.shape = {results['quadrature_weights'].shape}")