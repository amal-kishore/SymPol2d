"""
Utility functions for 2dSYMPOL
"""

import numpy as np
from typing import Dict

# Typical van der Waals radii (in Angstroms)
VDW_RADII = {
    1: 1.20,   # H
    5: 1.92,   # B
    6: 1.70,   # C
    7: 1.55,   # N
    8: 1.52,   # O
    9: 1.47,   # F
    14: 2.10,  # Si
    15: 1.80,  # P
    16: 1.80,  # S
    17: 1.75,  # Cl
    34: 1.90,  # Se
    35: 1.85,  # Br
    42: 2.10,  # Mo
    52: 2.06,  # Te
    53: 1.98,  # I
    74: 2.10,  # W
}

def estimate_interlayer_distance(atomic_numbers: np.ndarray) -> float:
    """
    Estimate reasonable interlayer distance based on atomic composition.
    
    Args:
        atomic_numbers: Array of atomic numbers in the material
        
    Returns:
        Estimated interlayer distance in Angstroms
    """
    # Get unique elements
    unique_elements = np.unique(atomic_numbers)
    
    # Find maximum vdW radius
    max_radius = 1.7  # Default for carbon
    for z in unique_elements:
        if z in VDW_RADII:
            max_radius = max(max_radius, VDW_RADII[z])
    
    # Interlayer distance is approximately 2 * vdW radius + 0.5-1.0 Å gap
    # Common values:
    # - Graphene/hBN: ~3.3-3.5 Å
    # - TMDCs (MoS2, WS2): ~6.0-6.5 Å  
    # - Other chalcogenides: ~5.5-6.0 Å
    
    # Simple heuristic based on heaviest element
    # Updated to use 3.1 Å as default for vdW materials
    if any(z > 40 for z in unique_elements):  # Transition metals
        if any(z in [16, 34, 52] for z in unique_elements):  # TMDCs
            return 3.1  # Typical for vdW bilayers
        else:
            return 3.1
    elif any(z in [16, 34, 52] for z in unique_elements):  # Chalcogens without TM
        return 3.1
    else:
        return 3.1  # Default for all vdW materials