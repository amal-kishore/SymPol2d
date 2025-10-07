"""
Interface to c2db database for material extraction
"""

import sqlite3
import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import json


@dataclass
class Material2D:
    """Represents a 2D material from c2db"""
    uid: str
    formula: str
    layer_group: str
    lattice: np.ndarray  # 3x3 matrix
    positions: np.ndarray  # Nx3 atomic positions
    numbers: np.ndarray  # Atomic numbers
    natoms: int
    
    def __repr__(self):
        return f"Material2D({self.formula}, {self.layer_group})"
    
    def get_chemical_symbols(self):
        """Get chemical symbols from atomic numbers"""
        # Simple atomic number to symbol mapping
        ATOMIC_SYMBOLS = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
            11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
            19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe',
            27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se',
            35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
            43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
            51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce',
            59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy',
            67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W',
            75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
            83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
            91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf'
        }
        return [ATOMIC_SYMBOLS.get(int(num), f'X{int(num)}') for num in self.numbers]


class C2DBInterface:
    """Interface to c2db SQLite database"""
    
    def __init__(self, db_path: str = 'c2db.db'):
        """Initialize connection to c2db database"""
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()
    
    def __del__(self):
        """Close database connection"""
        if hasattr(self, 'conn'):
            self.conn.close()
    
    def get_material_by_uid(self, uid: str) -> Optional[Material2D]:
        """
        Retrieve a material by its unique ID.
        
        Args:
            uid: Unique identifier (e.g., '1MoS2-3')
            
        Returns:
            Material2D object or None if not found
        """
        # First get the system data
        query = """
        SELECT s.id, s.cell, s.positions, s.numbers, s.natoms
        FROM systems s
        JOIN text_key_values t ON s.id = t.id
        WHERE t.key = 'uid' AND t.value = ?
        """
        
        self.cursor.execute(query, (uid,))
        result = self.cursor.fetchone()
        
        if not result:
            print(f"Material with uid '{uid}' not found")
            return None
        
        system_id, cell_blob, pos_blob, num_blob, natoms_raw = result
        
        # Get additional properties
        self.cursor.execute("""
        SELECT key, value FROM text_key_values 
        WHERE id = ? AND key IN ('layergroup', 'uid')
        """, (system_id,))
        
        properties = dict(self.cursor.fetchall())
        
        # Parse binary data
        if isinstance(natoms_raw, bytes):
            # Handle packed natoms if needed
            natoms = int.from_bytes(natoms_raw[:8], 'little')
        else:
            natoms = natoms_raw
            
        lattice = np.frombuffer(cell_blob, dtype=float).reshape(3, 3)
        positions = np.frombuffer(pos_blob, dtype=float).reshape(-1, 3)
        
        # Handle atomic numbers - they might be stored as int32 or int64
        try:
            numbers = np.frombuffer(num_blob, dtype=np.int32)
        except:
            try:
                numbers = np.frombuffer(num_blob, dtype=np.int64)
            except:
                # Fallback - assume one atom type
                numbers = np.ones(natoms, dtype=int)
        
        # Extract formula from uid (e.g., '1MoS2-3' -> 'MoS2')
        formula = uid
        if '-' in uid:
            parts = uid.split('-')
            if len(parts[0]) > 0 and parts[0][0].isdigit():
                # Remove leading digit
                formula = parts[0][1:]
            else:
                formula = parts[0]
        
        return Material2D(
            uid=uid,
            formula=formula,
            layer_group=properties.get('layergroup', 'p1'),
            lattice=lattice,
            positions=positions,
            numbers=numbers,
            natoms=natoms
        )
    
    def search_materials(self, formula: Optional[str] = None, 
                        layer_group: Optional[str] = None,
                        limit: int = 10) -> List[Dict[str, str]]:
        """
        Search for materials by formula or layer group.
        
        Args:
            formula: Chemical formula (e.g., 'MoS2')
            layer_group: Layer group symbol (e.g., 'p-4m2')
            limit: Maximum number of results
            
        Returns:
            List of material info dictionaries
        """
        conditions = []
        params = []
        
        if formula:
            conditions.append("t1.key = 'uid' AND t1.value LIKE ?")
            params.append(f'%{formula}%')
        
        if layer_group:
            conditions.append("t2.key = 'layergroup' AND t2.value = ?")
            params.append(layer_group)
        
        if not conditions:
            # Get some samples if no criteria specified
            query = """
            SELECT DISTINCT t1.value as uid, t2.value as layergroup
            FROM text_key_values t1
            JOIN text_key_values t2 ON t1.id = t2.id AND t2.key = 'layergroup'
            WHERE t1.key = 'uid'
            LIMIT ?
            """
            params = [limit]
        else:
            if formula and layer_group:
                query = f"""
                SELECT DISTINCT t1.value as uid, t2.value as layergroup
                FROM text_key_values t1
                JOIN text_key_values t2 ON t1.id = t2.id
                WHERE {' AND '.join(conditions)}
                LIMIT ?
                """
            elif formula:
                query = """
                SELECT DISTINCT t1.value as uid, 
                       (SELECT value FROM text_key_values WHERE id = t1.id AND key = 'layergroup') as layergroup
                FROM text_key_values t1
                WHERE t1.key = 'uid' AND t1.value LIKE ?
                LIMIT ?
                """
            else:  # layer_group only
                query = """
                SELECT DISTINCT t1.value as uid, t2.value as layergroup
                FROM text_key_values t1
                JOIN text_key_values t2 ON t1.id = t2.id
                WHERE t1.key = 'uid' AND t2.key = 'layergroup' AND t2.value = ?
                LIMIT ?
                """
            params.append(limit)
        
        self.cursor.execute(query, params)
        results = []
        
        for row in self.cursor.fetchall():
            uid = row[0]
            # Extract formula from uid
            formula_extracted = uid
            if '-' in uid:
                parts = uid.split('-')
                if len(parts[0]) > 0 and parts[0][0].isdigit():
                    formula_extracted = parts[0][1:]
                else:
                    formula_extracted = parts[0]
            
            results.append({
                'uid': uid,
                'formula': formula_extracted,
                'layer_group': row[1] or 'Unknown'
            })
        
        return results
    
    def get_all_layer_groups(self) -> List[Tuple[str, int]]:
        """Get all unique layer groups and their counts"""
        query = """
        SELECT value, COUNT(*) as count
        FROM text_key_values
        WHERE key = 'layergroup'
        GROUP BY value
        ORDER BY count DESC
        """

        self.cursor.execute(query)
        return self.cursor.fetchall()


def fetch_by_uid(uid: str, dbpath: str = "raw/c2db.db") -> Optional[Dict]:
    """
    Fetch material info by UID for the new CLI interface.

    Args:
        uid: C2DB unique identifier (e.g., '4P-1')
        dbpath: Path to c2db.db file

    Returns:
        Dictionary with 'uid', 'formula', 'layer_group' keys, or None if not found
    """
    import os
    import sqlite3

    if not os.path.exists(dbpath):
        return None

    try:
        conn = sqlite3.connect(dbpath)
        cursor = conn.cursor()

        # Get system ID from UID
        cursor.execute("""
            SELECT id FROM text_key_values
            WHERE key = 'uid' AND value = ?
        """, (uid,))

        result = cursor.fetchone()
        if not result:
            conn.close()
            return None

        system_id = result[0]

        # Get layer group
        cursor.execute("""
            SELECT value FROM text_key_values
            WHERE id = ? AND key = 'layergroup'
        """, (system_id,))

        lg_result = cursor.fetchone()
        layer_group = lg_result[0] if lg_result else "p1"

        # Extract formula from uid (e.g., '4P-1' -> 'P4')
        formula = uid.split('-')[0] if '-' in uid else uid
        # Handle leading digit: '4P' -> 'P4'
        if formula and formula[0].isdigit():
            i = 0
            while i < len(formula) and formula[i].isdigit():
                i += 1
            if i < len(formula):
                count = formula[:i]
                element = formula[i:]
                formula = element + count

        conn.close()

        return {
            "uid": uid,
            "formula": formula,
            "layer_group": layer_group
        }

    except Exception as e:
        print(f"Error fetching UID {uid}: {e}")
        return None