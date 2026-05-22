import os
from atomic_number_to_and_fro_symbol import atomic_converter

folder = '/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE/psps/PseudoDojo/standard/nc-sr-05_pbe_standard_psp8'   # <- change this

for filename in os.listdir(folder):
    
    if not filename.endswith('.py'):
        continue
    
    # try 2-letter symbol first (He, Li, Fe, ...), then 1-letter (H, C, N, ...)
    
    symbol = filename[-3:]
    try:
        z = atomic_converter(symbol)
    except Exception:
        continue
        
    new_name = f"{z}{filename[-3:]}"
    os.rename(os.path.join(folder, filename),
                os.path.join(folder, new_name))
    print(f"{filename}  ->  {new_name}")
    break