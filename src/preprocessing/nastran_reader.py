import os
from pathlib import Path
from typing import Union

# Importamos las clases definidas en tu data_model
from .data_model import (
    FEM_MODEL, 
    NODE, 
    ISOTROPIC_MATERIAL, 
    ROD_PROPERTY, 
    ROD_ELEMENT, 
    SINGLE_POINT_LOAD, 
    SINGLE_MOMENT_LOAD, 
    SPC_CONSTRAINT
)

"""
@file nastran_reader.py
@brief Module for reading, parsing, and interpreting NASTRAN files (.dat, .bdf).
@details Supports:
            - Relative path handling
            - Recursion for INCLUDE cards.
            - Fixed Field Format (Small Field, 8 characters) '$1234567'
            - Free Field Format (Free Field, separated by commas).
            - Comment lines
"""


def _parse_nastran_number(s: str) -> Union[float, None]:
    """
    @brief Converts a NASTRAN string to a Python float.
    @details Handles compressed scientific notation (e.g., '1.-4' -> 1.0e-4) and empty fields.
    @param s Raw text string from the field (8 characters). "$1234567"
    @return Float value or None if it is an empty field.
    """
    s = s.strip()
    if not s:
        return None
    
    try:
        # Detection of scientific format without "E" (ej: 2.1-1 -> 2.1E-1)
        # Search for a "-" or "+" that is not at the beginning and does not have an "E" or "D" in front of it.
        if ('-' in s[1:] or '+' in s[1:]) and 'E' not in s and 'D' not in s:
            # We scan to find the sign that acts as an exponent.
            for i, char in enumerate(s):
                if i > 0 and (char == '-' or char == '+'):
                    # If the previous character is not "E" or "D", insert "E".
                    if s[i-1] not in ['E', 'D', 'e', 'd']:
                        s = s[:i] + 'E' + s[i:]
                        break
        
        # NASTRAN sometimes uses "D" for double precision (e.g., 1.0D+4), Python uses "E".
        s = s.replace('D', 'E').replace('d', 'e')
        return float(s)
    
    except ValueError:
        return ValueError(f"The value '{s}' can not be converted into a float number.")

def _get_entry_fields(line: str) -> list[str]:
    """
    @brief Convert a text line into a list of fields.
    @details Automatically detects whether it is free format (with commas) 
            or fixed format (8 characters). 
            Delete Comment lines ($).
    """
    
    if '$' in line: # 1º Comment lines are errased (everything following an '$')
        line = line.split('$')[0]
    
    line = line.rstrip() # Remove trailing spaces
    
    if ',' in line:
        # --- Free Field, separated by commas (CSV) ---
        return [f.strip() for f in line.split(',')]
    else:
        # --- Small Field, 8 characters ---
        return [line[i:i+8].strip() for i in range(0, len(line), 8)] # Pieces of 8 spaces range(0, len, 8)



def _safe_get(fields: list[str], index: int, default: str = "") -> str:
    """
    @brief Recovers a raw field by index or returns an empty string ("").
    @details If the line is short (e.g., “GRID, 1”), requesting field 7 would return an error without this.
    """
    return fields[index] if index < len(fields) else ""

# region --- SEMANTIC EXTRACTION HELPERS (VALIDATION) ---
def _req_float(fields: list[str], idx: int, name: str) -> float:
    """Extracts a REQUIRED float. Throws an error if missing."""
    val = _parse_nastran_number(_safe_get(fields, idx))
    if val is None:
        raise ValueError(f"Required a float field, '{name}' is missing or empty.")
    return val

def _opt_float(fields: list[str], idx: int, default: float = 0.0) -> float:
    """Extracts an OPTIONAL float. Returns default if missing."""
    val = _parse_nastran_number(_safe_get(fields, idx))
    return val if val is not None else default

def _req_int(fields: list[str], idx: int, name: str) -> int:
    """Returns an integer REQUIRED. Throws an error if missing."""
    s = _safe_get(fields, idx)
    if not s:
        raise ValueError(f"Required a float field, '{name}' is missing or empty.")
    try:
        return int(s)
    except ValueError:
        raise ValueError(f"Required a integer fieldfor '{name}', recived value '{s}' ")

def _opt_int(fields: list[str], idx: int, default: int = 0) -> int:
    """Extracts an OPTIONAL int. Returns default if missing."""
    s = _safe_get(fields, idx)
    if not s: return default
    try:
        return int(s)
    except ValueError:
        return default # Or raise a warning

# endregion --- SEMANTIC EXTRACTION HELPERS (VALIDATION) ---

# region --- CARD READING FUNCTIONS ---

def _read_grid(fields: list[str], model: FEM_MODEL):
    # GRID(0) ID(1) CP(2) X1(3) X2(4) X3(5) CD(6) PS(7) SEID(8)

    nid = _req_int(fields, 1, "GRID ID")
    cp = _opt_int(fields, 2, 0)
    
    x1 = _req_float(fields, 3, f"X1 (Nodo {nid})")
    x2 = _req_float(fields, 4, f"X2 (Nodo {nid})")
    x3 = _req_float(fields, 5, f"X3 (Nodo {nid})")
    
    cd = _opt_int(fields, 6, 0)
    ps = _safe_get(fields, 7)
    seid = _opt_int(fields, 8, 0)
    
    new_node = NODE(nid=nid, coord_id=cp, x1=x1, x2=x2, x3=x3, cd_id=cd, ps=ps, seid=seid)
    
    # Quick pre-processing (assuming CP=0)
    if cp == 0: # no user coord system defined
        new_node.x, new_node.y, new_node.z = x1, x2, x3 
        
    model.nodes[nid] = new_node

def _read_crod(fields: list[str], model: FEM_MODEL):
    # CROD(0) EID(1) PID(2) G1(3) G2(4)
    eid = _req_int(fields, 1, "CROD EID")
    pid = _req_int(fields, 2, "CROD PID")
    n1 = _req_int(fields, 3, "CROD Node 1")
    n2 = _req_int(fields, 4, "CROD Node 2")
    
    model.elements[eid] = ROD_ELEMENT(eid=eid, prop_id=pid, nodes=[n1, n2])

def _read_prod(fields: list[str], model: FEM_MODEL):
    # PROD(0) PID(1) MID(2) A(3) J(4) C(5) NSM(6)
    pid = _req_int(fields, 1, "PROD PID")
    mid = _req_int(fields, 2, "PROD MID")
    A = _req_float(fields, 3, f"Area (A) en PROD {pid}")
    J = _opt_float(fields, 4, 0.0)
    C = _opt_float(fields, 5, 0.0)
    NSM = _opt_float(fields, 6, 0.0)
    
    model.properties[pid] = ROD_PROPERTY(pid=pid, mat_id=mid, A=A, J=J, C=C, NSM=NSM)

def _read_mat1(fields: list[str], model: FEM_MODEL):
    # MAT1(0) MID(1) E(2) G(3) NU(4) RHO(5)
    mid = _req_int(fields, 1, "MAT1 MID")
    # At least two fields between E,G and NU are mandatory. ISOTROPIC_MATERIAL will check it in '__post_init__'.
    E = _parse_nastran_number(_safe_get(fields, 2)) # Puede ser None
    G = _parse_nastran_number(_safe_get(fields, 3))
    nu = _parse_nastran_number(_safe_get(fields, 4))
    rho = _opt_float(fields, 5, 0.0)
    
    model.materials[mid] = ISOTROPIC_MATERIAL(mid=mid, E=E, nu=nu, G=G, rho=rho)

def _read_spc1(fields: list[str], model: FEM_MODEL):
    # SPC1(0) SID(1) C(2) G1(2) G2(3) ...
    sid = _req_int(fields, 1, "SPC1 SID")
    components = _safe_get(fields, 2)
    if not components:
        raise ValueError(f"SPC1 {sid} without restricted components.")
    
    # Nodos variables desde índice 3
    for i in range(3, len(fields)):
        node_str = fields[i]
        if not node_str: continue
        try:
            nid = int(node_str)
            model.constraints.append(
                SPC_CONSTRAINT(sid=sid, node_id=nid, components=components, value=0.0)
            )
        except ValueError:
            pass # Ignore THRU for now

def _read_force(fields: list[str], model: FEM_MODEL):
    # FORCE(0) SID(1) G(2) CID(3) F(4) N1(5) N2(6) N3(7)
    sid = _req_int(fields, 1, "FORCE SID")
    nid = _req_int(fields, 2, "FORCE Node ID")
    cid = _opt_int(fields, 3, 0)
    module = _req_float(fields, 4, f"Factor (F) in FORCE {sid}")
    
    n1 = _req_float(fields, 5, "Vector component  N1")
    n2 = _req_float(fields, 6, "Vector component N2")
    n3 = _req_float(fields, 7, "Vector component N3")
    
    model.loads.append(
        SINGLE_POINT_LOAD(sid=sid, appnode_id=nid, module=module, vector=[n1, n2, n3], cid=cid)
    )

def _read_moment(fields: list[str], model: FEM_MODEL):
    # MOMENT(0) SID(1) G(2) CID(3) F(4) N1(5) N2(6) N3(7)
    sid = _req_int(fields, 1, "MOMENT SID")
    nid = _req_int(fields, 2, "MOMENT Node ID")
    cid = _opt_int(fields, 3, 0)
    module = _req_float(fields, 4, f"Magnitud (M) en MOMENT {sid}")
    
    n1 = _req_float(fields, 5, "Vector component N1")
    n2 = _req_float(fields, 6, "Vector component N2")
    n3 = _req_float(fields, 7, "Vector component N3")
    
    model.loads.append(
        SINGLE_MOMENT_LOAD(sid=sid, appnode_id=nid, module=module, vector=[n1, n2, n3], cid=cid)
    )

# endregion --- CARD READING FUNCTIONS ---


# region --- RELATIONSHIP DICTIONARIES ---
# Map the card name to the function that processes it.
_CARD_HANDLERS: dict[str, callable[[list[str], FEM_MODEL], None]] = {
    'GRID': _read_grid,
    'CROD': _read_crod,
    'PROD': _read_prod,
    'MAT1': _read_mat1,
    'SPC1': _read_spc1,
    'FORCE': _read_force,
    'MOMENT': _read_moment
}

# endregion --- RELATIONSHIP DICTIONARIES ---

# region --- MAIN PROCESS ---
def _process_file(filepath: Path, model: FEM_MODEL) -> None:
    print(f"[INFO] Reading: {filepath.name}")
    if not filepath.exists():
        raise FileNotFoundError(f"[ERROR] file not found: {filepath}")

    with open(filepath, 'r') as f:
        for line_idx, line in enumerate(f):
            line_raw = line.strip()
            if not line_raw or line_raw.startswith('$'):
                continue
            
            # INCLUDE management (special recursive case)
            if line_raw.upper().startswith('INCLUDE'):
                quote = "'" if "'" in line_raw else '"'
                parts = line_raw.split(quote)
                if len(parts) >= 2:
                    inc_path = filepath.parent / parts[1]
                    _process_file(inc_path, model)
                continue

            # Get fields and card name
            fields = _get_entry_fields(line)
            if not fields: continue
            
            card_name = fields[0].upper()

            
            if card_name in _CARD_HANDLERS: # Search for known  NASTRAN cards.
                try:
                    _CARD_HANDLERS[card_name](fields, model)
                except Exception as e:
                    print(f"[ERROR] Error in line {line_idx+1} ({filepath.name}): {e}")
                    print(f"        Line content: {line_raw}")
                    # Aquí decides si haces 'raise' para abortar o 'continue' para intentar seguir
                    # raise e 
            else:
                print(f"[DEBUG] Card not supported or ignored: {card_name}")
                print(f"        Line {line_idx+1} ({filepath.name})")
                pass

def parse_nastran_input(entry_file_path: str) -> FEM_MODEL:
    """
    @brief Main entry point of nastran_reader.py file
    """
    model = FEM_MODEL()
    path_obj = Path(entry_file_path)
    
    try:
        _process_file(path_obj, model)
        print(f"\n[SUCCESS] Model successfully loaded.")
        print(f"  - Nodes: {len(model.nodes)}")
        print(f"  - Elements: {len(model.elements)}")
    except Exception as e:
        print(f"\n[FAIL] Error while trying to read nastran input file ({entry_file_path}): {e}")
        
    return model

# endregion --- MAIN PROCESS ---












