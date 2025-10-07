#!/usr/bin/env python3
"""
Create DFT-ready h-BN bilayer CIF files directly
"""

import numpy as np
from datetime import datetime

def write_cif_file(filename, stacking_type, tau, description):
    """
    Write a DFT-ready CIF file for h-BN bilayer.
    """

    # h-BN parameters
    a = 2.504  # lattice parameter in Angstroms
    c = 25.0   # c-axis with ~20 Å vacuum
    interlayer = 3.33  # interlayer distance in Angstroms

    # Calculate z positions for centering the bilayer in the cell
    z_center = 0.5
    z_offset = interlayer / (2 * c)
    z1 = z_center - z_offset  # First layer
    z2 = z_center + z_offset  # Second layer

    cif_content = f"""# DFT-ready h-BN bilayer structure
# Generated for DFT calculations with proper vacuum spacing
# Stacking: {stacking_type} - tau = [{tau[0]:.4f}, {tau[1]:.4f}]
# Interlayer distance: {interlayer:.3f} Angstroms
# Vacuum spacing: ~20 Angstroms
# Description: {description}
# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

data_{stacking_type}_hBN_DFT

_chemical_name_common     'h-BN {stacking_type} bilayer for DFT'

_cell_length_a            {a:.6f}
_cell_length_b            {a:.6f}
_cell_length_c            {c:.6f}
_cell_angle_alpha         90.0000
_cell_angle_beta          90.0000
_cell_angle_gamma         120.0000

_symmetry_space_group_name_H-M    'P 1'
_symmetry_Int_Tables_number       1

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
B1   B    0.000000   0.000000   {z1:.6f}   1.0
N1   N    0.333333   0.666667   {z1:.6f}   1.0
B2   B    {tau[0]:.6f}   {tau[1]:.6f}   {z2:.6f}   1.0
N2   N    {(0.333333 + tau[0]) % 1:.6f}   {(0.666667 + tau[1]) % 1:.6f}   {z2:.6f}   1.0
"""

    with open(filename, 'w') as f:
        f.write(cif_content)

    return a, c, interlayer


def create_vasp_poscar(filename, stacking_type, tau, description):
    """
    Create VASP POSCAR file for h-BN bilayer.
    """

    # h-BN parameters
    a = 2.504
    c = 25.0
    interlayer = 3.33

    # Calculate z positions
    z_center = 0.5
    z_offset = interlayer / (2 * c)
    z1 = z_center - z_offset
    z2 = z_center + z_offset

    poscar_content = f"""h-BN {stacking_type} bilayer - {description}
1.0
    {a:.6f}    0.000000    0.000000
   {-a/2:.6f}   {a*np.sqrt(3)/2:.6f}    0.000000
    0.000000    0.000000   {c:.6f}
B    N
2    2
Direct
    0.000000    0.000000   {z1:.6f}
    {tau[0]:.6f}    {tau[1]:.6f}   {z2:.6f}
    0.333333    0.666667   {z1:.6f}
    {(0.333333 + tau[0]) % 1:.6f}    {(0.666667 + tau[1]) % 1:.6f}   {z2:.6f}
"""

    with open(filename, 'w') as f:
        f.write(poscar_content)


def create_quantum_espresso_input(filename, stacking_type, tau, description):
    """
    Create Quantum ESPRESSO input file template.
    """

    a = 2.504
    c = 25.0
    interlayer = 3.33

    # Calculate z positions
    z_center = 0.5
    z_offset = interlayer / (2 * c)
    z1 = z_center - z_offset
    z2 = z_center + z_offset

    # Convert to Cartesian coordinates for QE
    positions = [
        (0.0, 0.0, z1 * c),
        (0.0, a/np.sqrt(3), z1 * c),
        (tau[0] * a - tau[1] * a/2, tau[1] * a * np.sqrt(3)/2, z2 * c),
        ((0.333333 + tau[0]) * a - (0.666667 + tau[1]) * a/2,
         (0.666667 + tau[1]) * a * np.sqrt(3)/2, z2 * c)
    ]

    qe_content = f"""&CONTROL
    calculation = 'scf'
    prefix = 'hBN_{stacking_type}'
    pseudo_dir = './pseudo/'
    outdir = './tmp/'
/

&SYSTEM
    ibrav = 0
    nat = 4
    ntyp = 2
    ecutwfc = 80.0
    ecutrho = 320.0
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.01
    vdw_corr = 'DFT-D3'
/

&ELECTRONS
    conv_thr = 1.0d-8
/

ATOMIC_SPECIES
    B  10.81  B.pbe-n-kjpaw_psl.1.0.0.UPF
    N  14.01  N.pbe-n-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS (angstrom)
    {a:.6f}    0.000000    0.000000
   {-a/2:.6f}   {a*np.sqrt(3)/2:.6f}    0.000000
    0.000000    0.000000   {c:.6f}

ATOMIC_POSITIONS (crystal)
    B    0.000000    0.000000   {z1:.6f}
    N    0.333333    0.666667   {z1:.6f}
    B    {tau[0]:.6f}    {tau[1]:.6f}   {z2:.6f}
    N    {(0.333333 + tau[0]) % 1:.6f}    {(0.666667 + tau[1]) % 1:.6f}   {z2:.6f}

K_POINTS (automatic)
    12 12 1 0 0 0
"""

    with open(filename, 'w') as f:
        f.write(qe_content)


def main():
    """Create all DFT-ready h-BN structures"""

    print("Creating DFT-ready h-BN bilayer structures...")
    print("=" * 60)

    # Define configurations
    configs = {
        'AA': {
            'tau': np.array([0.0, 0.0]),
            'description': 'AA\' (bulk) stacking - Non-polar'
        },
        'AB': {
            'tau': np.array([1/3, 2/3]),  # Standard AB stacking
            'description': 'AB stacking - Polar (downward, -P)'
        },
        'BA': {
            'tau': np.array([2/3, 1/3]),  # BA stacking (opposite of AB)
            'description': 'BA stacking - Polar (upward, +P)'
        }
    }

    for stacking, config in configs.items():
        print(f"\nCreating {stacking} configuration:")
        print("-" * 40)
        print(f"Description: {config['description']}")
        print(f"Stacking vector τ: [{config['tau'][0]:.4f}, {config['tau'][1]:.4f}]")

        # Create CIF file
        cif_file = f'hBN_{stacking}_DFT.cif'
        a, c, d = write_cif_file(cif_file, stacking, config['tau'], config['description'])
        print(f"  ✓ Created: {cif_file}")

        # Create VASP POSCAR
        poscar_file = f'POSCAR_{stacking}'
        create_vasp_poscar(poscar_file, stacking, config['tau'], config['description'])
        print(f"  ✓ Created: {poscar_file}")

        # Create QE input template
        qe_file = f'hBN_{stacking}.in'
        create_quantum_espresso_input(qe_file, stacking, config['tau'], config['description'])
        print(f"  ✓ Created: {qe_file}")

        print(f"\n  Structure parameters:")
        print(f"    Lattice a, b: {a:.3f} Å")
        print(f"    Lattice c: {c:.3f} Å (includes ~20 Å vacuum)")
        print(f"    Interlayer distance: {d:.3f} Å")

    print("\n" + "=" * 60)
    print("SUCCESS: DFT-ready files created!")
    print("\nRecommended DFT settings:")
    print("  • Functional: PBE-D3, optB88-vdW, or vdW-DF2")
    print("  • K-points: 12×12×1 (Gamma-centered)")
    print("  • Energy cutoff: 500 eV (VASP) or 80 Ry (QE)")
    print("  • Force threshold: 0.01 eV/Å")
    print("  • Electronic convergence: 10⁻⁸ eV")
    print("\nPolarization calculation:")
    print("  • Use Berry phase method")
    print("  • For VASP: LCALCPOL = .TRUE., DIPOL = 0.5 0.5 0.5")
    print("  • For QE: Use 'berry_phase' calculation")
    print("=" * 60)


if __name__ == "__main__":
    main()