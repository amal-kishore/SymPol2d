#!/usr/bin/env python3
"""
Command-line interface for SYMPOL2D
"""

import argparse
import sys
import json
from pathlib import Path
import numpy as np

from .c2db_interface import C2DBInterface
from .scanner import StackingScanner
from .symmetry import LayerGroupSymmetry
from .utils import estimate_interlayer_distance
from .cif_writer import save_all_stackings_cif


def search_material(args):
    """Search for polar stackings of a specific material"""
    
    # Check if database exists
    db_path = Path(args.database)
    if not db_path.exists():
        print(f"❌ Error: Database file '{args.database}' not found!")
        print("\nSuggestions:")
        print("  1. Check if the database path is correct")
        print("  2. Download the C2DB database from: https://cmr.fysik.dtu.dk/c2db/c2db.html")
        print("  3. Use --database option to specify the correct path")
        print(f"  4. For the example, the database is in: example/c2db.db")
        return 1
    
    print(f"✅ Database found: {args.database}")
    
    # Connect to c2db
    try:
        db = C2DBInterface(args.database)
    except Exception as e:
        print(f"❌ Error connecting to database: {e}")
        print("\nThe database file may be corrupted or in an unsupported format.")
        return 1
    
    # Get material
    if args.uid:
        material = db.get_material_by_uid(args.uid)
        if not material:
            print(f"Error: Material '{args.uid}' not found in database")
            return 1
    elif args.formula:
        # Search by formula
        results = db.search_materials(formula=args.formula, limit=10)
        if not results:
            print(f"Error: No materials found with formula '{args.formula}'")
            return 1
        
        if len(results) > 1:
            print(f"Found {len(results)} materials with formula '{args.formula}':")
            for i, mat in enumerate(results):
                print(f"  {i+1}. {mat['uid']} ({mat['layer_group']})")
            
            # Use first one or ask user to specify
            if not args.auto_select:
                print("\nPlease specify --uid with the exact material ID")
                return 1
            else:
                print(f"\nAuto-selecting first match: {results[0]['uid']}")
                material = db.get_material_by_uid(results[0]['uid'])
        else:
            material = db.get_material_by_uid(results[0]['uid'])
    else:
        print("Error: Must specify either --uid or --formula")
        return 1
    
    print(f"\nAnalyzing: {material.formula} ({material.uid})")
    print(f"Layer group: {material.layer_group}")
    print(f"Number of atoms: {material.natoms}")
    
    # Determine interlayer distance
    if args.interlayer_distance is None:
        # Auto-estimate based on material composition
        interlayer_d = estimate_interlayer_distance(material.numbers)
        print(f"Auto-estimated interlayer distance: {interlayer_d:.2f} Å")
    else:
        interlayer_d = args.interlayer_distance
        print(f"User-specified interlayer distance: {interlayer_d:.2f} Å")
    
    # Create scanner
    scanner = StackingScanner(material.layer_group, grid_size=args.grid, 
                              interlayer_distance=interlayer_d)
    
    print(f"\nScanning {args.grid}x{args.grid} grid for stacking configurations...")
    
    # Find representative stackings
    representatives = scanner.get_representative_stackings()
    
    if not representatives:
        print("Warning: No stackings found")
        return 1
    
    # Display results
    print("\n" + "="*60)
    print("STACKING CONFIGURATIONS FOUND:")
    print("="*60)
    
    # AA stacking
    if 'AA' in representatives:
        aa = representatives['AA']
        print(f"\nAA stacking (non-polar):")
        print(f"  τ = [{aa.tau[0]:.4f}, {aa.tau[1]:.4f}], d = {aa.interlayer_distance:.2f} Å")
        print(f"  Preserved symmetries: {', '.join(aa.preserved_symmetries)}")
    
    # Get all polar pairs
    if args.polar_direction == 'all':
        # Show all pairs grouped by direction
        all_pairs = scanner.find_polar_pairs()
        
        # Group by polarization direction
        pairs_by_direction = {}
        for ab, ba in all_pairs:
            direction = ab.polar_direction
            if direction not in pairs_by_direction:
                pairs_by_direction[direction] = []
            pairs_by_direction[direction].append((ab, ba))
        
        # Display by direction
        for direction in sorted(pairs_by_direction.keys()):
            pairs = pairs_by_direction[direction]
            print(f"\n{direction.upper()}-POLARIZED PAIRS ({len(pairs)} found):")
            
            # Show first few representative pairs
            for i, (ab, ba) in enumerate(pairs[:3]):  # Show first 3 of each type
                print(f"\n  Pair {i+1}:")
                print(f"    AB: τ = [{ab.tau[0]:.4f}, {ab.tau[1]:.4f}]")
                print(f"    BA: τ = [{ba.tau[0]:.4f}, {ba.tau[1]:.4f}]")
                print(f"    Broken symmetries: {', '.join(ab.broken_symmetries[:3])}...")
            
            if len(pairs) > 3:
                print(f"  ... and {len(pairs) - 3} more {direction}-polarized pairs")
    else:
        # Show only pairs with specific polarization
        pairs = scanner.find_polar_pairs_by_direction(args.polar_direction)
        
        if pairs:
            print(f"\n{args.polar_direction.upper()}-POLARIZED PAIRS ({len(pairs)} found):")
            
            for i, (ab, ba) in enumerate(pairs):
                print(f"\n  Pair {i+1}:")
                print(f"    AB: τ = [{ab.tau[0]:.4f}, {ab.tau[1]:.4f}]")
                print(f"    BA: τ = [{ba.tau[0]:.4f}, {ba.tau[1]:.4f}]")
                print(f"    Broken symmetries: {', '.join(ab.broken_symmetries)}")
        else:
            print(f"\nNo {args.polar_direction}-polarized pairs found.")
    
    # Save results if requested
    if args.output:
        output_path = Path(args.output)
        
        # Determine output format
        if args.output_format == 'json' or output_path.suffix == '.json':
            # Save as JSON
            output_data = {
                'material': {
                    'uid': material.uid,
                    'formula': material.formula,
                    'layer_group': material.layer_group
                },
                'grid_size': args.grid,
                'stackings': {}
            }
            
            for name, config in representatives.items():
                output_data['stackings'][name] = {
                    'tau': config.tau.tolist(),
                    'preserved_symmetries': config.preserved_symmetries,
                    'broken_symmetries': config.broken_symmetries,
                    'polar_direction': config.polar_direction
                }
            
            if output_path.suffix != '.json':
                output_path = output_path.with_suffix('.json')
            
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=2)
            
            print(f"\nResults saved to: {output_path}")
            
        else:  # Default to CIF format
            # Create output directory for CIF files
            if output_path.suffix in ['.cif', '']:
                # If user specified a .cif file or no extension, use parent as directory
                output_dir = output_path.parent / output_path.stem
            else:
                output_dir = output_path
            
            # Save all stackings as CIF files
            cif_files = save_all_stackings_cif(material, representatives, output_dir)
            
            print(f"\nCIF files saved to: {output_dir}/")
            for cif_file in cif_files:
                print(f"  - {cif_file.name}")
    
    # Also save CIF files by default if no output specified
    elif representatives:
        default_dir = Path(f"{material.formula}_stackings")
        cif_files = save_all_stackings_cif(material, representatives, default_dir)
        
        print(f"\nCIF files saved to: {default_dir}/")
        for cif_file in cif_files:
            print(f"  - {cif_file.name}")
    
    return 0


def list_materials(args):
    """List materials in the database"""
    
    # Check if database exists
    db_path = Path(args.database)
    if not db_path.exists():
        print(f"❌ Error: Database file '{args.database}' not found!")
        print("\nSuggestions:")
        print("  1. Check if the database path is correct")
        print("  2. Download the C2DB database from: https://cmr.fysik.dtu.dk/c2db/c2db.html")
        print("  3. Use --database option to specify the correct path")
        print(f"  4. For the example, the database is in: example/c2db.db")
        return 1
    
    print(f"✅ Database found: {args.database}")
    
    try:
        db = C2DBInterface(args.database)
    except Exception as e:
        print(f"❌ Error connecting to database: {e}")
        print("\nThe database file may be corrupted or in an unsupported format.")
        return 1
    
    if args.layer_groups:
        # List all layer groups
        print("\nLayer groups in database:")
        print("-" * 40)
        groups = db.get_all_layer_groups()
        for group, count in groups[:20]:  # Show top 20
            print(f"  {group:12s} : {count:5d} materials")
        if len(groups) > 20:
            print(f"  ... and {len(groups)-20} more groups")
    else:
        # List sample materials
        results = db.search_materials(
            formula=args.formula,
            layer_group=args.layer_group,
            limit=args.limit
        )
        
        print(f"\nFound {len(results)} materials:")
        print("-" * 60)
        print(f"{'UID':30s} {'Formula':15s} {'Layer Group':15s}")
        print("-" * 60)
        
        for mat in results:
            print(f"{mat['uid']:30s} {mat['formula']:15s} {mat['layer_group']:15s}")
    
    return 0


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='SYMPOL2D: SYMmetry-based prediction of POLarity in 2D bilayers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search by material UID
  sympol2d search --uid MoS2-165798ab5e18 --grid 50
  
  # Search by formula
  sympol2d search --formula MoS2 --grid 50 --auto-select
  
  # List materials
  sympol2d list --formula WS2
  
  # List layer groups
  sympol2d list --layer-groups
        """
    )
    
    parser.add_argument('--database', '-d', default='c2db.db',
                       help='Path to c2db database (default: c2db.db)')
    
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Search command
    search_parser = subparsers.add_parser('search', 
                                         help='Search for polar stackings')
    search_parser.add_argument('--uid', '-u', type=str,
                              help='Material UID (e.g., MoS2-165798ab5e18)')
    search_parser.add_argument('--formula', '-f', type=str,
                              help='Chemical formula (e.g., MoS2)')
    search_parser.add_argument('--grid', '-g', type=int, default=50,
                              help='Grid size for scanning (default: 50)')
    search_parser.add_argument('--output', '-o', type=str,
                              help='Output path (directory for CIF, file for JSON)')
    search_parser.add_argument('--output-format', '-of', type=str,
                              choices=['cif', 'json'], default='cif',
                              help='Output format (default: cif)')
    search_parser.add_argument('--auto-select', action='store_true',
                              help='Auto-select first match when multiple found')
    search_parser.add_argument('--polar-direction', '-pd', type=str,
                              choices=['x', 'y', 'z', 'xy', 'general', 'all'],
                              default='all',
                              help='Filter by polarization direction (default: all)')
    search_parser.add_argument('--interlayer-distance', '-d', type=float, default=None,
                              help='Interlayer distance in Angstroms (default: auto-estimate)')
    
    # List command
    list_parser = subparsers.add_parser('list',
                                       help='List materials in database')
    list_parser.add_argument('--formula', '-f', type=str,
                            help='Filter by formula')
    list_parser.add_argument('--layer-group', '-lg', type=str,
                            help='Filter by layer group')
    list_parser.add_argument('--layer-groups', action='store_true',
                            help='List all layer groups')
    list_parser.add_argument('--limit', '-n', type=int, default=20,
                            help='Maximum number of results (default: 20)')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    if args.command == 'search':
        return search_material(args)
    elif args.command == 'list':
        return list_materials(args)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())