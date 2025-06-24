#!/usr/bin/env python3
"""
Create CORRECT WS2 bilayer CIF files using proper AB/BA stackings
AB: τ = [1/3, 1/3] - Standard B-site over A-site  
BA: τ = [2/3, 2/3] - Standard A-site over B-site (opposite polarity)
"""

import numpy as np

def read_monolayer_cif(filename):
    """Read monolayer CIF and extract structure parameters"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    a = b = c = None
    alpha = beta = gamma = None
    atoms = []
    
    for line in lines:
        if '_cell_length_a' in line:
            a = float(line.split()[1])
        elif '_cell_length_b' in line:
            b = float(line.split()[1])
        elif '_cell_length_c' in line:
            c = float(line.split()[1])
        elif '_cell_angle_alpha' in line:
            alpha = float(line.split()[1])
        elif '_cell_angle_beta' in line:
            beta = float(line.split()[1])
        elif '_cell_angle_gamma' in line:
            gamma = float(line.split()[1])
        elif len(line.strip()) > 0 and not line.startswith('_') and not line.startswith('loop_') and not line.startswith('data_') and not line.startswith("'"):
            parts = line.strip().split()
            if len(parts) >= 6:
                element = parts[0]
                label = parts[1]
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])
                atoms.append((element, label, x, y, z))
    
    return {
        'a': a, 'b': b, 'c': c,
        'alpha': alpha, 'beta': beta, 'gamma': gamma,
        'atoms': atoms
    }

def create_bilayer_cif(monolayer_data, tau, interlayer_distance_angstrom, output_file, stacking_type, description):
    """Create bilayer CIF with exact interlayer distance"""
    
    # Extract actual layer thickness from atom positions
    z_coords = [atom[4] for atom in monolayer_data['atoms']]
    z_min, z_max = min(z_coords), max(z_coords)
    layer_thickness_frac = z_max - z_min
    layer_thickness_ang = layer_thickness_frac * monolayer_data['c']
    
    # Design bilayer structure
    vacuum = 15.0  # Vacuum on top and bottom
    new_c_angstrom = vacuum + layer_thickness_ang + interlayer_distance_angstrom + layer_thickness_ang + vacuum
    
    with open(output_file, 'w') as f:
        f.write(f"data_WS2_bilayer_{stacking_type}\n")
        f.write("_chemical_formula_structural       W2S4\n")
        f.write('_chemical_formula_sum              "W2 S4"\n')
        f.write(f"_cell_length_a       {monolayer_data['a']:.12f}\n")
        f.write(f"_cell_length_b       {monolayer_data['b']:.12f}\n")
        f.write(f"_cell_length_c       {new_c_angstrom:.12f}\n")
        f.write(f"_cell_angle_alpha    {monolayer_data['alpha']}\n")
        f.write(f"_cell_angle_beta     {monolayer_data['beta']}\n")
        f.write(f"_cell_angle_gamma    {monolayer_data['gamma']}\n")
        f.write("\n")
        f.write('_space_group_name_H-M_alt    "P 1"\n')
        f.write("_space_group_IT_number       1\n")
        f.write("\n")
        f.write("loop_\n")
        f.write("  _space_group_symop_operation_xyz\n")
        f.write("  'x, y, z'\n")
        f.write("\n")
        f.write("loop_\n")
        f.write("  _atom_site_type_symbol\n")
        f.write("  _atom_site_label\n")
        f.write("  _atom_site_symmetry_multiplicity\n")
        f.write("  _atom_site_fract_x\n")
        f.write("  _atom_site_fract_y\n")
        f.write("  _atom_site_fract_z\n")
        f.write("  _atom_site_occupancy\n")
        
        # Layer 1 (bottom) - normalize and center
        layer1_center_ang = vacuum + layer_thickness_ang/2
        for i, (element, label, x, y, z) in enumerate(monolayer_data['atoms']):
            z_normalized = (z - z_min) / layer_thickness_frac
            z_layer = (z_normalized - 0.5) * layer_thickness_ang
            new_z = (layer1_center_ang + z_layer) / new_c_angstrom
            f.write(f"  {element}   {label}_L1        1.0  {x:.12f}  {y:.12f}  {new_z:.12f}  1.0000\n")
        
        # Layer 2 (top) - apply tau shift and place above gap
        layer2_center_ang = vacuum + layer_thickness_ang + interlayer_distance_angstrom + layer_thickness_ang/2
        for i, (element, label, x, y, z) in enumerate(monolayer_data['atoms']):
            # Apply tau shift
            new_x = (x + tau[0]) % 1.0
            new_y = (y + tau[1]) % 1.0
            z_normalized = (z - z_min) / layer_thickness_frac
            z_layer = (z_normalized - 0.5) * layer_thickness_ang
            new_z = (layer2_center_ang + z_layer) / new_c_angstrom
            f.write(f"  {element}   {label}_L2        1.0  {new_x:.12f}  {new_y:.12f}  {new_z:.12f}  1.0000\n")
        
        f.write(f"\n# {description}\n")
        f.write(f"# tau = [{tau[0]:.6f}, {tau[1]:.6f}] (high-symmetry stacking)\n")
        f.write(f"# Interlayer distance = {interlayer_distance_angstrom:.1f} Å\n")
        f.write(f"# Layer thickness = {layer_thickness_ang:.2f} Å each\n")

def analyze_stacking_geometry(tau_AB, tau_BA):
    """Analyze the geometric relationship between AB and BA stackings"""
    
    print("Geometric analysis of AB vs BA:")
    print(f"AB: τ = [{tau_AB[0]:.6f}, {tau_AB[1]:.6f}]")
    print(f"BA: τ = [{tau_BA[0]:.6f}, {tau_BA[1]:.6f}]")
    print()
    
    # Check various relationships
    diff = tau_BA - tau_AB
    print(f"BA - AB = [{diff[0]:.6f}, {diff[1]:.6f}]")
    
    # Check if BA = -AB (mod 1)
    neg_AB = (-tau_AB) % 1.0
    is_inversion = np.allclose(tau_BA, neg_AB, atol=1e-6)
    print(f"-AB (mod 1) = [{neg_AB[0]:.6f}, {neg_AB[1]:.6f}]")
    print(f"Is BA = -AB? {is_inversion}")
    
    # Check if BA = 1 - AB
    one_minus_AB = 1.0 - tau_AB
    is_complement = np.allclose(tau_BA, one_minus_AB, atol=1e-6)
    print(f"1 - AB = [{one_minus_AB[0]:.6f}, {one_minus_AB[1]:.6f}]")
    print(f"Is BA = 1 - AB? {is_complement}")
    
    print()
    print("Physical meaning:")
    if tau_AB[0] == 1/3:
        print("AB (1/3, 1/3): Top layer B-sites over bottom layer A-sites")
        print("BA (2/3, 2/3): Top layer A-sites over bottom layer B-sites")
        print("These are the standard TMDC AB and BA stackings.")
        print("They represent opposite polarizations in 2H structure.")

def main():
    """Generate CORRECT AB/BA bilayer structures"""
    
    monolayer_data = read_monolayer_cif('1WS2-1.cif')
    print(f"WS2 monolayer structure:")
    print(f"  a = {monolayer_data['a']:.3f} Å")
    print(f"  atoms: {len(monolayer_data['atoms'])}")
    
    interlayer_distance = 3.1
    
    # CORRECT high-symmetry stackings
    stackings = {
        'AA': {
            'tau': np.array([0.0, 0.0]), 
            'description': 'Non-polar AA stacking (metal over metal)'
        },
        'AB': {
            'tau': np.array([1/3, 1/3]),  # Standard AB stacking
            'description': 'Polar AB stacking (B-sites over A-sites)'
        },
        'BA': {
            'tau': np.array([2/3, 2/3]),  # Standard BA stacking  
            'description': 'Polar BA stacking (A-sites over B-sites, opposite to AB)'
        }
    }
    
    print(f"\nGenerating CORRECT AB/BA structures:")
    
    # Analyze relationship between AB and BA
    analyze_stacking_geometry(stackings['AB']['tau'], stackings['BA']['tau'])
    print()
    
    for stacking_type, config in stackings.items():
        tau = config['tau']
        description = config['description']
        output_file = f"WS2_{stacking_type}_correct.cif"
        
        print(f"Creating {output_file}:")
        print(f"  τ = [{tau[0]:.6f}, {tau[1]:.6f}] - {description}")
        create_bilayer_cif(monolayer_data, tau, interlayer_distance, output_file, stacking_type, description)
    
    print(f"\n✅ CORRECT AB/BA bilayer structures created!")
    print(f"✅ Uses standard TMDC stackings (1/3, 1/3) and (2/3, 2/3)")
    print(f"✅ These represent true opposite polarities")
    print(f"✅ 3.1 Å interlayer distance maintained")

if __name__ == '__main__':
    main()