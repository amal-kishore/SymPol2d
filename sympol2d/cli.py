#!/usr/bin/env python3
"""
Command-line interface for SYMPOL2D
"""

import argparse
import json
import sys
import numpy as np
from pathlib import Path

from .scanner import scan_for_z, pick_best_pair, z_sign_flip_expected
from .symmetry import mod1
from .builder import make_bilayer_poscar
from .poscar_io import load_poscar
from . import c2db_interface as c2db


def main():
    """Main CLI entry point"""
    p = argparse.ArgumentParser(
        prog="SymPol2D",
        description="Sliding Ferroelectricity Prescreen (Out-of-Plane)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using monolayer POSCAR (recommended)
  python -m sympol2d.cli search --layer-group p-6m2 --poscar POSCAR --gap 3.1 --export

  # Using C2DB UID (requires c2db.db)
  python -m sympol2d.cli search --uid 1MoS2-1 --grid 60 --export

  # Rectangular system (no Pz flip expected)
  python -m sympol2d.cli search --layer-group pman --poscar BP_mono.vasp --export --allow-nonflipping
"""
    )

    s = p.add_subparsers(dest="cmd")

    sp = s.add_parser("search", help="Find ONE AB/BA pair for out-of-plane sliding ferroelectricity")
    sp.add_argument("--uid", type=str, help="C2DB UID (needs local c2db.db)")
    sp.add_argument("--layer-group", type=str, help="Layer group (e.g., p-6m2, p6mm, p2mm, pman)")
    sp.add_argument("--grid", type=int, default=60, help="tau-grid resolution (default: 60)")
    sp.add_argument("--poscar", type=str, help="Monolayer POSCAR (for structure export)")
    sp.add_argument("--gap", type=float, default=3.1, help="Target interlayer gap in Angstroms (default: 3.1)")
    sp.add_argument("--export", action="store_true", help="Write exactly one AB/BA pair as POSCARs")
    sp.add_argument("--out-prefix", type=str, default="sympol2d", help="Output file prefix (default: sympol2d)")
    sp.add_argument("--allow-nonflipping", action="store_true",
                    help="Export even if group is z-even (no opposite Pz under sliding)")
    sp.add_argument("--database", "-d", type=str, default="raw/c2db.db",
                    help="Path to c2db database (default: raw/c2db.db)")

    args = p.parse_args()

    if args.cmd != "search":
        p.print_help()
        return 0

    # Resolve layer group (UID → layer_group if DB present)
    lg = args.layer_group
    formula = "X2D"

    if args.uid and not lg:
        rec = c2db.fetch_by_uid(args.uid, args.database)
        if rec:
            lg, formula = rec["layer_group"], rec["formula"]
        else:
            print("Warning: c2db.db not available or UID not found; provide --layer-group.")

    if not lg:
        print("Error: provide --layer-group or a valid --uid with c2db.db.")
        return 2

    # 1) Scan τ-space
    print(f"Scanning {args.grid}×{args.grid} grid for layer group '{lg}'...")
    cand = scan_for_z(lg, ngrid=args.grid)

    if not cand:
        payload = {
            "ok": True,
            "layer_group": lg,
            "message": "No z-allowed stackings found."
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(f"Found {len(cand)} z-allowed candidates")

    # 2) Pick ONE best pair
    tau, tau_ba = pick_best_pair(lg, cand)

    if tau is None:
        payload = {
            "ok": True,
            "layer_group": lg,
            "message": "No candidate pair chosen."
        }
        print(json.dumps(payload, indent=2))
        return 0

    flip = z_sign_flip_expected(lg)

    payload = {
        "ok": True,
        "layer_group": lg,
        "z_sign_flip_expected": flip,
        "tau_AB": [float(tau[0]), float(tau[1])],
        "tau_BA": [float(tau_ba[0]), float(tau_ba[1])],
        "note": ("Opposite Pz under sliding is symmetry-expected; exporting exactly one AB/BA pair."
                 if flip else
                 "Opposite Pz under sliding is NOT symmetry-available for this group; see --allow-nonflipping."),
    }
    print(json.dumps(payload, indent=2))

    # 3) Optional export (POSCARs) — only when we have monolayer geometry
    if args.export:
        if not args.poscar:
            print("\nExport skipped: please supply --poscar MONOLAYER_POSCAR to build real bilayers.")
            return 0

        poscar_path = Path(args.poscar)
        if not poscar_path.exists():
            print(f"\nError: POSCAR file '{args.poscar}' not found.")
            return 1

        print(f"\nLoading monolayer structure from: {args.poscar}")
        comment, lattice, symbols, counts, mono_frac, formula_guess = load_poscar(args.poscar)
        formula = formula_guess or formula

        # Compute dz_frac = monolayer_thickness_frac + gap/c
        zmin, zmax = mono_frac[:, 2].min(), mono_frac[:, 2].max()
        mono_thick_frac = (zmax - zmin)  # fractional thickness of monolayer
        c = lattice[2, 2] if abs(lattice[2, 0]) < 1e-9 and abs(lattice[2, 1]) < 1e-9 else np.linalg.norm(lattice[2])
        dz_frac = mono_thick_frac + args.gap / float(c)

        print(f"Monolayer thickness: {mono_thick_frac*c:.3f} Å ({mono_thick_frac:.4f} fractional)")
        print(f"Interlayer gap: {args.gap:.3f} Å ({args.gap/c:.4f} fractional)")
        print(f"Total vertical offset: {dz_frac*c:.3f} Å ({dz_frac:.4f} fractional)")

        if flip or args.allow_nonflipping:
            poscar_ab = make_bilayer_poscar(
                formula=formula, lattice=lattice, frac_coords_mono=mono_frac,
                tau=np.array(tau), dz_frac=dz_frac,
                comment=f"{formula} bilayer AB (tau={tau[0]:.3f},{tau[1]:.3f})"
            )
            poscar_ba = make_bilayer_poscar(
                formula=formula, lattice=lattice, frac_coords_mono=mono_frac,
                tau=np.array(tau_ba), dz_frac=dz_frac,
                comment=f"{formula} bilayer BA (tau={tau_ba[0]:.3f},{tau_ba[1]:.3f})"
            )

            ab_file = f"{args.out_prefix}_AB.vasp"
            ba_file = f"{args.out_prefix}_BA.vasp"

            with open(ab_file, "w") as f:
                f.write(poscar_ab)
            with open(ba_file, "w") as f:
                f.write(poscar_ba)

            print(f"\nExported bilayer structures:")
            print(f"  AB: {ab_file}")
            print(f"  BA: {ba_file}")
        else:
            print("\nNo export: group is z-even; AB↔BA do not give opposite Pz.")
            print("Use --allow-nonflipping to export anyway.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
