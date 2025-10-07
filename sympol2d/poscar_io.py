"""
POSCAR file I/O for SYMPOL2D
"""

import re
import numpy as np


def load_poscar(path):
    """
    Read VASP POSCAR/CONTCAR file (VASP5 format).

    Args:
        path: Path to POSCAR file

    Returns:
        Tuple of (comment, lattice, symbols, counts, coords, formula)
        - comment: str, first line
        - lattice: 3×3 ndarray, lattice vectors in Angstroms
        - symbols: list of element symbols
        - counts: list of atom counts per element
        - coords: N×3 ndarray, fractional coordinates
        - formula: str, chemical formula (e.g., "MoS2")

    Note:
        - Assumes 'Direct' coordinates
        - Handles 'Selective dynamics' if present
    """
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    comment = lines[0]
    scale = float(lines[1])
    lattice = np.array([[float(x) for x in lines[2+i].split()] for i in range(3)], float) * scale

    # VASP5: symbols line then counts
    sym_line = lines[5].split()
    cnt_line = [int(x) for x in lines[6].split()]
    n = sum(cnt_line)

    # Check for Direct/Cartesian and selective dynamics
    mode_idx = 7
    mode = lines[mode_idx].lower()
    start = mode_idx + 1

    if mode.startswith("selective") or mode.startswith("s"):
        mode = lines[mode_idx + 1].lower()
        start = mode_idx + 2

    if not (mode.startswith("direct") or mode.startswith("d")):
        raise ValueError(f"Only 'Direct' coordinates supported, found: {mode}")

    # Parse coordinates (first 3 columns)
    coords = np.array([[float(x) for x in lines[start+i].split()[:3]] for i in range(n)], float)

    # Build formula string
    formula = "".join(s + (str(c) if c > 1 else "") for s, c in zip(sym_line, cnt_line))

    return comment, lattice, sym_line, cnt_line, coords, formula
