"""
Core symmetry analysis module for SYMPOL2D
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


@dataclass
class SymmetryOperation:
    """Represents a 2D symmetry operation"""
    name: str
    matrix: np.ndarray
    type: str  # 'rotation', 'mirror', 'identity'
    
    def __repr__(self):
        return f"SymOp({self.name})"


class LayerGroupSymmetry:
    """Handles 2D layer group symmetry operations"""
    
    # Common 2D point group operations
    OPERATIONS = {
        'E': np.array([[1, 0], [0, 1]]),  # Identity
        'C2': np.array([[-1, 0], [0, -1]]),  # 180° rotation
        'C3': np.array([[-0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, -0.5]]),  # 120° rotation
        'C4': np.array([[0, -1], [1, 0]]),  # 90° rotation
        'C6': np.array([[0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, 0.5]]),  # 60° rotation
        'Mx': np.array([[-1, 0], [0, 1]]),  # Mirror x
        'My': np.array([[1, 0], [0, -1]]),  # Mirror y
        'Mxy': np.array([[0, 1], [1, 0]]),  # Mirror along x=y
        'Mxy-': np.array([[0, -1], [-1, 0]]),  # Mirror along x=-y
    }
    
    # Layer group to point group mapping (simplified, extend as needed)
    LAYER_GROUP_OPERATIONS = {
        'p1': ['E'],
        'p2': ['E', 'C2'],
        'pm': ['E', 'My'],
        'pg': ['E', 'My'],  # glide plane
        'cm': ['E', 'My'],
        'p2mm': ['E', 'C2', 'Mx', 'My'],
        'p2mg': ['E', 'C2', 'Mx', 'My'],
        'p2gg': ['E', 'C2'],
        'c2mm': ['E', 'C2', 'Mx', 'My'],
        'p3': ['E', 'C3', 'C3^2'],
        'p3m1': ['E', 'C3', 'C3^2', 'Mx', 'Mxy', 'Mxy-'],
        'p31m': ['E', 'C3', 'C3^2', 'My', 'Mxy', 'Mxy-'],
        'p4': ['E', 'C4', 'C2', 'C4^3'],
        'p4mm': ['E', 'C4', 'C2', 'C4^3', 'Mx', 'My', 'Mxy', 'Mxy-'],
        'p-4m2': ['E', 'C4', 'C2', 'C4^3', 'Mx', 'My', 'Mxy', 'Mxy-'],  # Example from MoS2
        'p6': ['E', 'C6', 'C3', 'C2', 'C3^2', 'C6^5'],
        'p6mm': ['E', 'C6', 'C3', 'C2', 'C3^2', 'C6^5', 'Mx', 'My', 'Mxy', 'Mxy-'],
        'p-6m2': ['E', 'C6', 'C3', 'C2', 'C3^2', 'C6^5', 'Mx', 'My', 'Mxy'],  # Common for TMDCs
    }
    
    def __init__(self, layer_group: str):
        """Initialize with a layer group symbol"""
        self.layer_group = layer_group.lower()
        self.operations = self._generate_operations()
    
    def _generate_operations(self) -> List[SymmetryOperation]:
        """Generate symmetry operations for the layer group"""
        operations = []
        
        if self.layer_group not in self.LAYER_GROUP_OPERATIONS:
            # Default to p1 if unknown
            print(f"Warning: Unknown layer group '{self.layer_group}', defaulting to p1")
            op_names = ['E']
        else:
            op_names = self.LAYER_GROUP_OPERATIONS[self.layer_group]
        
        for op_name in op_names:
            if op_name in self.OPERATIONS:
                matrix = self.OPERATIONS[op_name]
            elif '^' in op_name:
                # Handle powers like C3^2
                base, power = op_name.split('^')
                base_matrix = self.OPERATIONS[base]
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
        E = np.eye(2)
        
        for op in self.operations:
            # Calculate (E + R)τ
            result = (E + op.matrix) @ tau
            
            # Check if result is an integer vector
            is_integer = np.allclose(result, np.round(result), atol=tolerance)
            results[op.name] = is_integer
        
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
        # This is a simplified criterion - adjust based on specific needs
        return 'AA'
    
    def get_broken_symmetries(self, tau: np.ndarray) -> List[str]:
        """Get list of symmetries broken by a stacking"""
        preserved = self.test_symmetry_preservation(tau)
        return [op for op, is_preserved in preserved.items() if not is_preserved]