"""
Grid scanner module for finding polar stackings
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from .symmetry import LayerGroupSymmetry


@dataclass
class StackingConfiguration:
    """Represents a stacking configuration"""
    tau: np.ndarray
    type: str  # 'AA', 'AB', 'BA'
    preserved_symmetries: List[str]
    broken_symmetries: List[str]
    polar_direction: Optional[str] = None  # 'x', 'y', 'z', 'xy', etc.
    interlayer_distance: Optional[float] = None  # in Angstroms
    
    def __repr__(self):
        pol_str = f", {self.polar_direction}-polar" if self.polar_direction else ""
        d_str = f", d={self.interlayer_distance:.2f}Å" if self.interlayer_distance else ""
        return f"{self.type}(τ=[{self.tau[0]:.3f}, {self.tau[1]:.3f}]{pol_str}{d_str})"


class StackingScanner:
    """Scans grid of stacking vectors to find polar configurations"""
    
    def __init__(self, layer_group: str, grid_size: int = 50, interlayer_distance: float = 3.1):
        """
        Initialize scanner with layer group and grid size.
        
        Args:
            layer_group: 2D layer group symbol
            grid_size: Number of grid points along each direction
            interlayer_distance: Default interlayer distance in Angstroms (default: 3.1)
        """
        self.symmetry = LayerGroupSymmetry(layer_group)
        self.grid_size = grid_size
        self.interlayer_distance = interlayer_distance
        self.stackings = []
    
    def _determine_polar_direction(self, broken_symmetries: List[str], preserved_symmetries: List[str]) -> str:
        """
        Determine polarization direction based on broken and preserved symmetries.
        
        Args:
            broken_symmetries: List of broken symmetry operations
            preserved_symmetries: List of preserved symmetry operations
            
        Returns:
            Polarization direction: 'x', 'y', 'z', 'xy', or 'general'
        """
        # Parse symmetries
        broken_set = set(broken_symmetries)
        preserved_set = set(preserved_symmetries)
        
        # Check mirror planes
        mx_broken = 'Mx' in broken_set
        my_broken = 'My' in broken_set
        mxy_broken = 'Mxy' in broken_set
        mxy_minus_broken = 'Mxy-' in broken_set
        
        # Check if C2 (180° rotation) is preserved
        c2_preserved = 'C2' in preserved_set
        
        # Polarization logic:
        # - x-polar: My preserved, Mx broken
        # - y-polar: Mx preserved, My broken  
        # - z-polar: Both Mx and My broken, but C2 preserved (out-of-plane)
        # - xy-polar: Diagonal mirrors broken
        
        if mx_broken and not my_broken and 'My' in preserved_set:
            return 'x'
        elif my_broken and not mx_broken and 'Mx' in preserved_set:
            return 'y'
        elif mx_broken and my_broken and c2_preserved:
            # Both in-plane mirrors broken but 180° rotation preserved -> z-polar
            return 'z'
        elif (mxy_broken or mxy_minus_broken) and not (mx_broken and my_broken):
            return 'xy'
        else:
            return 'general'
    
    def scan_grid(self) -> Dict[str, List[StackingConfiguration]]:
        """
        Scan the full grid of stacking vectors.
        
        Returns:
            Dictionary with 'AA' and 'polar' stackings
        """
        results = {'AA': [], 'polar': []}
        
        # Generate grid points
        grid_1d = np.linspace(0, 1, self.grid_size, endpoint=False)
        
        for i, tau_x in enumerate(grid_1d):
            for j, tau_y in enumerate(grid_1d):
                tau = np.array([tau_x, tau_y])
                
                # Test symmetry preservation
                preserved = self.symmetry.test_symmetry_preservation(tau)
                preserved_ops = [op for op, is_preserved in preserved.items() if is_preserved]
                broken_ops = [op for op, is_preserved in preserved.items() if not is_preserved]
                
                # Classify stacking
                stacking_type = self.symmetry.classify_stacking(tau)
                
                # Determine polarization direction for polar stackings
                polar_dir = None
                if stacking_type == 'polar':
                    polar_dir = self._determine_polar_direction(broken_ops, preserved_ops)
                
                config = StackingConfiguration(
                    tau=tau,
                    type=stacking_type,
                    preserved_symmetries=preserved_ops,
                    broken_symmetries=broken_ops,
                    polar_direction=polar_dir,
                    interlayer_distance=self.interlayer_distance
                )
                
                if stacking_type == 'AA':
                    results['AA'].append(config)
                else:
                    results['polar'].append(config)
        
        return results
    
    def find_polar_pairs(self) -> List[Tuple[StackingConfiguration, StackingConfiguration]]:
        """
        Find AB-BA pairs among polar stackings, prioritizing high-symmetry configurations.
        
        AB and BA should be related by inversion: BA = -AB (mod 1) = 1 - AB
        Prioritizes standard TMDC stackings like (1/3, 1/3) and (2/3, 2/3).
        
        Returns:
            List of (AB, BA) configuration pairs, sorted by symmetry importance
        """
        results = self.scan_grid()
        polar_stackings = results['polar']
        
        if not polar_stackings:
            return []
        
        pairs = []
        used = set()
        
        # First, look for high-symmetry pairs (standard TMDC stackings)
        high_symmetry_fractions = [1/6, 1/4, 1/3, 1/2, 2/3, 3/4, 5/6]
        
        def is_high_symmetry(tau, tolerance=0.02):
            """Check if tau is close to high-symmetry positions"""
            for frac in high_symmetry_fractions:
                if (abs(tau[0] - frac) < tolerance and abs(tau[1] - frac) < tolerance):
                    return True, frac
            return False, None
        
        # Sort polar stackings by symmetry priority
        def symmetry_priority(config):
            is_hs, frac = is_high_symmetry(config.tau)
            if is_hs:
                # Prioritize 1/3 and 2/3 (standard AB/BA)
                if abs(frac - 1/3) < 0.01 or abs(frac - 2/3) < 0.01:
                    return 0  # Highest priority
                else:
                    return 1  # High symmetry but not 1/3, 2/3
            return 2  # Lower priority
        
        polar_stackings.sort(key=symmetry_priority)
        
        for i, config1 in enumerate(polar_stackings):
            if i in used:
                continue
                
            tau1 = config1.tau
            
            # Look for the inversion partner: tau2 = 1 - tau1
            target_tau = (1.0 - tau1) % 1.0
            
            best_match = None
            best_distance = float('inf')
            best_j = -1
            
            for j, config2 in enumerate(polar_stackings):
                if j <= i or j in used:
                    continue
                    
                tau2 = config2.tau
                
                # Check distance to inversion partner
                distance = np.linalg.norm(tau2 - target_tau)
                if distance < self.grid_size**(-1) * 2:  # Within 2 grid spacings
                    if distance < best_distance:
                        best_distance = distance
                        best_match = config2
                        best_j = j
            
            if best_match is not None:
                # Verify they break the same symmetries
                if set(config1.broken_symmetries) == set(best_match.broken_symmetries):
                    # Assign AB/BA based on which is smaller
                    if np.sum(config1.tau) < np.sum(best_match.tau):
                        config1.type = 'AB'
                        best_match.type = 'BA'
                        pairs.append((config1, best_match))
                    else:
                        config1.type = 'BA'
                        best_match.type = 'AB'
                        pairs.append((best_match, config1))
                    
                    used.add(i)
                    used.add(best_j)
        
        # Sort pairs by priority (high-symmetry first)
        def pair_priority(pair):
            ab, ba = pair
            ab_hs, ab_frac = is_high_symmetry(ab.tau)
            ba_hs, ba_frac = is_high_symmetry(ba.tau)
            
            if ab_hs and ba_hs:
                # Both high symmetry - prioritize 1/3, 2/3 pairs
                if (abs(ab_frac - 1/3) < 0.01 and abs(ba_frac - 2/3) < 0.01):
                    return 0  # Highest priority - standard TMDC
                else:
                    return 1  # High symmetry but not standard
            return 2  # Lower priority
        
        pairs.sort(key=pair_priority)
        return pairs
    
    def find_polar_pairs_by_direction(self, direction: str) -> List[Tuple[StackingConfiguration, StackingConfiguration]]:
        """
        Find AB-BA pairs with specific polarization direction.
        
        Args:
            direction: Polarization direction ('x', 'y', 'z', 'xy', or 'general')
            
        Returns:
            List of (AB, BA) pairs with the specified polarization
        """
        all_pairs = self.find_polar_pairs()
        filtered_pairs = []
        
        for ab, ba in all_pairs:
            if ab.polar_direction == direction and ba.polar_direction == direction:
                filtered_pairs.append((ab, ba))
        
        return filtered_pairs
    
    def get_representative_stackings(self) -> Dict[str, StackingConfiguration]:
        """
        Get representative AA, AB, and BA stackings, prioritizing high-symmetry configurations.
        
        Returns:
            Dictionary with representative stackings, with AB/BA being the most symmetric pairs
        """
        results = self.scan_grid()
        representatives = {}
        
        # AA: Find the one at origin [0, 0]
        if results['AA']:
            aa_stackings = results['AA']
            origin_dist = [np.linalg.norm(s.tau) for s in aa_stackings]
            representatives['AA'] = aa_stackings[np.argmin(origin_dist)]
        
        # Find AB-BA pairs (now correctly prioritized)
        pairs = self.find_polar_pairs()
        
        if pairs:
            # First pair is now the highest priority (e.g., 1/3, 2/3 for TMDCs)
            representatives['AB'] = pairs[0][0]
            representatives['BA'] = pairs[0][1]
            
            # Store additional pairs with descriptive names
            if len(pairs) > 1:
                for i, (ab, ba) in enumerate(pairs[1:], 2):
                    # Check if they're high-symmetry
                    ab_frac = self._get_symmetry_fraction(ab.tau)
                    ba_frac = self._get_symmetry_fraction(ba.tau)
                    
                    if ab_frac and ba_frac:
                        # Use fraction-based names for high-symmetry pairs
                        ab_name = f'AB_{ab_frac[0]}_{ab_frac[1]}'
                        ba_name = f'BA_{ba_frac[0]}_{ba_frac[1]}'
                        representatives[ab_name] = ab
                        representatives[ba_name] = ba
                    else:
                        # Use numeric names for lower-symmetry pairs
                        representatives[f'AB{i}'] = ab
                        representatives[f'BA{i}'] = ba
        
        return representatives
    
    def _get_symmetry_fraction(self, tau):
        """Get the closest simple fraction representation of tau"""
        common_fractions = [(1, 6), (1, 4), (1, 3), (1, 2), (2, 3), (3, 4), (5, 6)]
        
        for num, den in common_fractions:
            frac_val = num / den
            if abs(tau[0] - frac_val) < 0.02 and abs(tau[1] - frac_val) < 0.02:
                return (num, den)
        return None