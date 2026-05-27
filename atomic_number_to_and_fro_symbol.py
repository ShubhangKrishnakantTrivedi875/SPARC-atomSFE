def atomic_converter(x):
    
    elements = {
        1:  'H',   2:  'He',  3:  'Li',  4:  'Be',  5:  'B',
        6:  'C',   7:  'N',   8:  'O',   9:  'F',   10: 'Ne',
        11: 'Na',  12: 'Mg',  13: 'Al',  14: 'Si',  15: 'P',
        16: 'S',   17: 'Cl',  18: 'Ar',  19: 'K',   20: 'Ca',
        21: 'Sc',  22: 'Ti',  23: 'V',   24: 'Cr',  25: 'Mn',
        26: 'Fe',  27: 'Co',  28: 'Ni',  29: 'Cu',  30: 'Zn',
        31: 'Ga',  32: 'Ge',  33: 'As',  34: 'Se',  35: 'Br',
        36: 'Kr',  37: 'Rb',  38: 'Sr',  39: 'Y',   40: 'Zr',
        41: 'Nb',  42: 'Mo',  43: 'Tc',  44: 'Ru',  45: 'Rh',
        46: 'Pd',  47: 'Ag',  48: 'Cd',  49: 'In',  50: 'Sn',
        51: 'Sb',  52: 'Te',  53: 'I',   54: 'Xe',  55: 'Cs',
        56: 'Ba',  57: 'La',  58: 'Ce',  59: 'Pr',  60: 'Nd',
        61: 'Pm',  62: 'Sm',  63: 'Eu',  64: 'Gd',  65: 'Tb',
        66: 'Dy',  67: 'Ho',  68: 'Er',  69: 'Tm',  70: 'Yb',
        71: 'Lu',  72: 'Hf',  73: 'Ta',  74: 'W',   75: 'Re',
        76: 'Os',  77: 'Ir',  78: 'Pt',  79: 'Au',  80: 'Hg',
        81: 'Tl',  82: 'Pb',  83: 'Bi',  84: 'Po',  85: 'At',
        86: 'Rn',  87: 'Fr',  88: 'Ra',  89: 'Ac',  90: 'Th',
        91: 'Pa',  92: 'U',
    }
    
    numbers = {sym: z for z, sym in elements.items()}
    
    if isinstance(x, int):
        if x not in elements:
            raise Exception("Atomic number out of range: enter Z between 1 (H) and 92 (U).")
        return elements[x]
    
    elif isinstance(x, str):
        if x not in numbers:
            raise Exception("Symbol not recognized: must be an element between 1 (H) and 92 (U).")
        return numbers[x]
    
    else:
        raise Exception("Input must be an integer atomic number or a string atomic symbol.")
    
def element_names(x):
    elements = {
        1:  'Hydrogen',   2:  'Helium',  3:  'Lithium',  4:  'Beryllium',  5:  'Boron',
        6:  'Carbon',   7:  'Nitrogen',   8:  'Oxygen',   9:  'Fluorine',   10: 'Neon',
        11: 'Sodium',  12: 'Magnesium',  13: 'Aluminum',  14: 'Silicon',  15: 'Phosphorus',
        16: 'Sulfur',   17: 'Chlorine',  18: 'Argon',  19: 'Potassium',   20: 'Calcium',
        21: 'Scandium',  22: 'Titanium',  23: 'Vanadium',   24: 'Chromium',  25: 'Manganese',
        26: 'Iron',  27: 'Cobalt',  28: 'Nickel',  29: 'Copper',  30: 'Zinc',
        31: 'Gallium',  32: 'Germanium',  33: 'Arsenic',  34: 'Selenium',  35: 'Bromine',
        36: 'Krypton',  37: 'Rubidium',  38: 'Strontium',  39: 'Yttrium',   40: 'Zirconium',
        41: 'Niobium',  42: 'Molybdenum',  43: 'Technetium',  44: 'Ruthenium',  45: 'Rhodium',
        46: 'Palladium',  47: 'Silver',  48: 'Cadmium',  49: 'Indium',  50: 'Tin',
        51: 'Antimony',  52: 'Tellurium',  53: 'Iodine',   54: 'Xenon',  55: 'Caesium',
        56: 'Barium',  57: 'Lanthanum',  58: 'Cerium',  59: 'Praseodymium',  60: 'Neodymium',
        61: 'Promethium',  62: 'Samarium',  63: 'Europium',  64: 'Gadolinium',  65: 'Terbium',
        66: 'Dysprosium',  67: 'Holmium',  68: 'Erbium',  69: 'Thulium',  70: 'Ytterbium',
        71: 'Lutetium',  72: 'Hafnium',  73: 'Tantalum',  74: 'Tungsten',   75: 'Rhenium',
        76: 'Osmium',  77: 'Iridium',  78: 'Platinum',  79: 'Gold',  80: 'Mercury',
        81: 'Thallium',  82: 'Lead',  83: 'Bismuth',  84: 'Polonium',  85: 'Astatine',
        86: 'Radon',  87: 'Francium',  88: 'Radium',  89: 'Actinium',  90: 'Thorium',
        91: 'Protactinium',  92: 'Uranium',
    }
    
    numbers = {sym: z for z, sym in elements.items()}
    
    if isinstance(x, int):
        if x not in elements:
            raise Exception("Atomic number out of range: enter Z between 1 (Hydrogen) and 92 (Uranium).")
        return elements[x]
    
    elif isinstance(x, str):
        if x not in numbers:
            raise Exception("Element name not recognized: must be an element between 1 (Hydrogen) and 92 (Uranium).")
        return numbers[x]
    
    else:
        raise Exception("Input must be an integer atomic number or a string element.")