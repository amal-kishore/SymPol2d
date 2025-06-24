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