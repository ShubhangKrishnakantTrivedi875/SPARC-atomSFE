
from src.solver import AtomicDFTSolver
import numpy as np
import pickle

atomic_number = 2
xc_functional = 'GGA_PBE'
all_electron_flag = True


solver = AtomicDFTSolver(
    atomic_number             = 2,
    xc_functional             = 'GGA_PBE',
    domain_size               = 40.0,
    finite_element_number     = 4,
    polynomial_order          = 20,
    quadrature_point_number   = 60,
    mesh_type                 = "uniform",
    mesh_concentration        = 101.0,
    scf_tolerance             = 1e-9,
    mesh_spacing              = 0.01,
    
    verbose                   = True, 
    all_electron_flag         = True,
    use_oep                   = False,
    use_preconditioner        = False,
    
   # psp_dir_path = "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE",
   # psp_file_name= "02_He_2_1.1_1.2_pbe_v1.0.psp8"
)

results = solver.solve()


with open("struct1.pkl","wb") as file:
    pickle.dump(results, file)
    
    
# print("\nSolver results:")
# print(f"\t total energy             = {results['energy']:.6f} (Ha)")
# print(f"\t density.shape            = {results['rho'].shape}")
# print(f"\t quadrature_nodes.shape   = {results['quadrature_nodes'].shape}")
# print(f"\t quadrature_weights.shape = {results['quadrature_weights'].shape}")