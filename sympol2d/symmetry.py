"""
Core symmetry analysis module for SYMPOL2D
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


# 2x2 linear parts of in-plane operations (reduced lattice basis)
OPER = {
    'E':  np.array([[ 1,  0],
                    [ 0,  1]], float),

    'C2': np.array([[-1,  0],
                    [ 0, -1]], float),  # 180° about z

    'Mx': np.array([[ 1,  0],
                    [ 0, -1]], float),  # mirror about y-axis (reflects y, x fixed)

    'My': np.array([[-1,  0],
                    [ 0,  1]], float),  # mirror about x-axis (reflects x, y fixed)

    'C6': np.array([[ 0.5, -np.sqrt(3)/2],
                    [ np.sqrt(3)/2, 0.5]], float),

    'C3': np.array([[-0.5, -np.sqrt(3)/2],
                    [ np.sqrt(3)/2, -0.5]], float),

    'C4': np.array([[ 0, -1],
                    [ 1,  0]], float),
}

# Minimal layer-group → linear-ops dictionary for our test.
# Glides/centerings: same linear part as mirrors/rotations; translations are absorbed on RHS.
LAYER_GROUP_OPERATIONS = {
    # Rectangular family
    'p2mm': ['E', 'C2', 'Mx', 'My'],
    'pman': ['E', 'C2', 'Mx', 'My'],
    'pmmm': ['E', 'C2', 'Mx', 'My'],
    'cmmm': ['E', 'C2', 'Mx', 'My'],
    'pmm2': ['E', 'C2', 'Mx', 'My'],
    'cmm2': ['E', 'C2', 'Mx', 'My'],

    # Single-axis symmetry
    'p1': ['E'],
    'p-1': ['E', 'C2'],  # Triclinic with inversion (C2 = inversion in 2D)
    'p2': ['E', 'C2'],
    'pm': ['E', 'My'],
    'pg': ['E', 'My'],
    'cm': ['E', 'My'],
    'p2mg': ['E', 'C2', 'Mx', 'My'],
    'p2gg': ['E', 'C2'],
    'c2mm': ['E', 'C2', 'Mx', 'My'],

    # Hexagonal / triangular
    'p-6m2': ['E', 'C6', 'C3', 'C2', 'Mx', 'My'],
    'p6mm':  ['E', 'C6', 'C3', 'C2', 'Mx', 'My'],
    'p3m1':  ['E', 'C3', 'C3', 'Mx', 'My'],
    'p31m':  ['E', 'C3', 'C3', 'Mx', 'My'],
    'p3':    ['E', 'C3', 'C3'],
    'p-3':   ['E', 'C3', 'C3', 'C2'],
    'p-3m1': ['E', 'C3', 'C3', 'C2', 'Mx', 'My'],
    'p-31m': ['E', 'C3', 'C3', 'C2', 'Mx', 'My'],

    # Square
    'p4':     ['E', 'C4', 'C2'],
    'p4mm':   ['E', 'C4', 'C2', 'Mx', 'My'],
    'p-4m2':  ['E', 'C4', 'C2', 'Mx', 'My'],
    'p-4':    ['E', 'C4', 'C2'],

    # Hexagonal rotations only
    'p6':     ['E', 'C6', 'C3', 'C2'],
    'p-6':    ['E', 'C6', 'C3', 'C2'],
}

# Does AB↔BA (tau, 1-tau) predict opposite Pz under *sliding* for this group?
Z_SIGN_FLIP_EXPECTED = {
    # Hexagonal: AB↔BA related by inversion → Pz flips
    'p-6m2': True, 'p6mm': True, 'p3m1': True, 'p31m': True,
    'p-3m1': True, 'p-31m': True, 'p-6': True, 'p6': True,
    'p3': True, 'p-3': True,

    # Triclinic with inversion: AB↔BA related by inversion → Pz flips
    'p-1': True,

    # Rectangular: AB↔BA related by C2z (z-even) → Pz does NOT flip
    'p2mm': False, 'pman': False, 'pmmm': False, 'cmmm': False,
    'pmm2': False, 'cmm2': False, 'c2mm': False, 'p2mg': False,
    'p2': False, 'p2gg': False,

    # Square: AB↔BA related by C2z (z-even) → Pz does NOT flip
    'p4mm': False, 'p-4m2': False, 'p4': False, 'p-4': False,

    # Single mirror: depends on details, default False
    'pm': False, 'pg': False, 'cm': False,
    'p1': False,
}

def mod1(v):
    """Wrap to [0,1)"""
    return np.mod(v, 1.0)

def is_integer_vec(x, atol=1e-8):
    """Check if vector is integer (within tolerance)"""
    return np.all(np.abs(x - np.rint(x)) <= atol)

def survives(R, tau, atol=1e-8):
    """
    Test if operation R survives at stacking vector tau.
    Condition: (E + R)τ ∈ Z²
    """
    M = OPER[R]
    lhs = (np.eye(2) + M) @ tau
    return is_integer_vec(lhs, atol=atol)

def classify_z_allowed(survivors):
    """
    Classify if z-polarization is allowed based on surviving symmetries.

    Returns:
        (bool, str): (is_allowed, classification)
        - "z-only": C2 survives, both Mx and My broken (strong z-polar class)
        - "z-allowed": Pz not forbidden by surviving in-plane ops
    """
    ops = set(survivors)
    mx, my, c2 = ('Mx' in ops), ('My' in ops), ('C2' in ops)

    if c2 and (not mx) and (not my):
        return True, "z-only"
    return True, "z-allowed"


# Legacy compatibility classes (deprecated but kept for existing code)
@dataclass
class SymmetryOperation:
    """Represents a 2D symmetry operation (legacy)"""
    name: str
    matrix: np.ndarray
    type: str  # 'rotation', 'mirror', 'identity'

    def __repr__(self):
        return f"SymOp({self.name})"


class LayerGroupSymmetry:
    """Handles 2D layer group symmetry operations (legacy interface)"""

    # Keep old operation names for backwards compatibility
    OPERATIONS = OPER

    def __init__(self, layer_group: str):
        """Initialize with a layer group symbol"""
        self.layer_group = layer_group.lower()
        self.operations = self._generate_operations()

    def _generate_operations(self) -> List[SymmetryOperation]:
        """Generate symmetry operations for the layer group"""
        operations = []

        if self.layer_group not in LAYER_GROUP_OPERATIONS:
            # Default to p1 if unknown
            print(f"Warning: Unknown layer group '{self.layer_group}', defaulting to p1")
            op_names = ['E']
        else:
            op_names = LAYER_GROUP_OPERATIONS[self.layer_group]

        for op_name in op_names:
            if op_name in OPER:
                matrix = OPER[op_name]
            elif '^' in op_name:
                # Handle powers like C3^2
                base, power = op_name.split('^')
                base_matrix = OPER[base]
                matrix = np.linalg.matrix_power(base_matrix, int(power))
            else:
                print(f"Warning: Unknown operation '{op_name}'")
                continue

            # Determine type
            if 'M' in op_name:
                op_type = 'mirror'
            elif op_name == 'E':
                op_type = 'identity'
            else:
                op_type = 'rotation'

            operations.append(SymmetryOperation(op_name, matrix, op_type))

        return operations

    def test_symmetry_preservation(self, tau: np.ndarray, tolerance: float = 1e-6) -> Dict[str, bool]:
        """
        Test if a stacking vector tau preserves each symmetry operation.

        The condition is: (E + R)τ = n (integer vector)

        Args:
            tau: 2D stacking vector in fractional coordinates [0,1)^2
            tolerance: Numerical tolerance for integer check

        Returns:
            Dictionary mapping operation names to preservation status
        """
        results = {}

        for op in self.operations:
            results[op.name] = survives(op.name, tau, atol=tolerance)

        return results

    def classify_stacking(self, tau: np.ndarray) -> str:
        """
        Classify a stacking as AA (non-polar) or AB/BA (polar).

        Args:
            tau: 2D stacking vector in fractional coordinates

        Returns:
            'AA' if all symmetries preserved, 'polar' otherwise
        """
        preserved = self.test_symmetry_preservation(tau)

        # AA stacking preserves all symmetries
        if all(preserved.values()):
            return 'AA'

        # Check if any mirror or inversion symmetry is broken
        for op in self.operations:
            if op.type in ['mirror', 'inversion'] and not preserved[op.name]:
                return 'polar'

        # If only rotations are broken, might still be non-polar
        return 'AA'

    def get_broken_symmetries(self, tau: np.ndarray) -> List[str]:
        """Get list of symmetries broken by a stacking"""
        preserved = self.test_symmetry_preservation(tau)
        return [op for op, is_preserved in preserved.items() if not is_preserved]

    def get_preserved_symmetries(self, tau: np.ndarray) -> List[str]:
        """Get list of symmetries preserved by a stacking"""
        preserved = self.test_symmetry_preservation(tau)
        return [op for op, is_preserved in preserved.items() if is_preserved]
