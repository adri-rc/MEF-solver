from typing import Optional
import numpy as np

# Intentamos importar pyNastran. Si no está, avisamos al usuario.
try:
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.errors import MissingDeckSections
    from pyNastran.bdf.cards.materials import MAT1
    from pyNastran.bdf.cards.properties.rods import PROD
    from pyNastran.bdf.cards.properties.bars import PBAR
    from pyNastran.bdf.cards.elements.rods import CROD
    from pyNastran.bdf.cards.elements.bars import CBAR
    from pyNastran.bdf.cards.properties.shell import PSHELL
    from pyNastran.bdf.cards.elements.shell import CTRIA3, CQUAD4
except ImportError:
    raise ImportError(
        "'pyNastran' is required to read Nastran files. Install it with: pip install pyNastran"
    )

from .data_model import (
    FEMModel, 
    Node, 
    IsotropicMaterial, 
    RodProperty,
    BarProperty,
    RodElement,
    BarElement,
    SinglePointLoad, 
    SinglePointMoment, 
    SPCConstraint,
    ShellProperty,
    CTRIA3Element,
    CQUAD4Element
)

"""
@file nastran_reader_v2.py
@brief Nastran reader module using the pyNastran library.
@details
    This version replaces manual text parsing from old nastran_reader.py with the pyNastran API (BDF class).
    It acts as a 'mapper' between pyNastran internal objects 
    and the user's data_model.py.
"""

def _map_materials(bdf_model: BDF, model: FEMModel):
    """
    @brief Translates pyNastran materials to IsotropicMaterial project class.
    """
    for mid, mat in bdf_model.materials.items():
        if isinstance(mat, MAT1):
            # pyNastran returns 0.0 or None if the field is empty.
            # The data_model requires 'None' to trigger automatic calculation logic.
            # We convert 0.0 to None for E, G, nu if necessary.
            
            E_val = mat.e if mat.e != 0.0 else None
            G_val = mat.g if mat.g != 0.0 else None
            # Poisson's ratio (nu) CAN be 0.0 (ej. corcho).
            # pyNastran returns None if the field is empty in the file.
            # Accept the value as it is (None or float).
            nu_val = mat.nu
            
            # Note: pyNastran handles rho as mass density directly.
            rho_val = mat.rho

            try:
                model.materials[mid] = IsotropicMaterial(
                    mid=mid, E=E_val, G=G_val, nu=nu_val, rho=rho_val
                )
            except ValueError as e:
                print(f"\n[ERROR] Material ID {mid} is invalid/inconsistent.")
                print(f"              Reason: {e}")
                print("              Execution stopped to prevent incorrect results.\n")
                raise e

def _map_properties(bdf_model: BDF, model: FEMModel):
    """
    @brief Translates properties (PROD, PBAR) to data_model classes project.
    """
    for pid, prop in bdf_model.properties.items():
        
        # --- PROD Mapping ---
        if isinstance(prop, PROD):
            model.properties[pid] = RodProperty(
                pid=pid,
                mat_id=prop.Mid(),
                A=prop.Area(),
                J=prop.J(),
                C=prop.c,
                NSM=prop.nsm
            )
        
        # --- PBAR Mapping ---
        elif isinstance(prop, PBAR):
            # EXTRACTION OF STRESS RECOVERY POINTS (C, D, E, F)
            # pyNastran stores them in attributes .c, .d, .e, .f 
            # c1, c2 -> Y, Z coordinates  of point C
            # d1, d2 -> Y, Z coordinates  of point D, etc.
            
            extracted_points = []
            
            # Helper para extraer tupla (y, z) de forma segura
            def get_stress_point(prefix: str):
                # Intentamos leer prop.c1 y prop.c2 (por ejemplo)
                # Usamos getattr por seguridad, aunque suelen existir con default 0.0
                y = getattr(prop, f"{prefix}1", 0.0) or 0.0
                z = getattr(prop, f"{prefix}2", 0.0) or 0.0
                return (float(y), float(z))

            extracted_points = []
            # Iteramos ordenadamente C, D, E, F
            for char in ['c', 'd', 'e', 'f']:
                extracted_points.append(get_stress_point(char))
            
            # build the property class
            model.properties[pid] = BarProperty(
                pid=pid,
                mat_id=prop.Mid(),
                A=prop.Area(),
                I1=prop.I11(),
                I2=prop.I22(),
                J=prop.J(),
                I12=prop.I12(),
                NSM=prop.nsm,
                # K1 y K2 no siempre están en todas las versiones de tarjetas PBAR, 
                # pero pyNastran suele tenerlos.
                # Nota: pyNastran a veces los llama k1, k2.
                # Validamos si existen atributos, si no 0.0
                K1=getattr(prop, 'k1', 0.0) or 0.0,
                K2=getattr(prop, 'k2', 0.0) or 0.0,
                stress_points=extracted_points
            )

        # --- PSHELL Mapping ---
        elif isinstance(prop, PSHELL):
            
            mid2_val = prop.Mid2()
            if mid2_val == 0:
                mid2_val = -1
            
            mid3_val = prop.Mid3()
            if mid3_val == 0 or mid3_val is None:
                mid3_val = -1
                
            model.properties[pid] = ShellProperty(
                pid=pid,
                mat_id=prop.Mid1(),
                t=prop.t,
                mat_id_bending=mid2_val,
                mat_id_shear=None if mid3_val == -1 else mid3_val,
                bending_inertia_ratio=prop.twelveIt3,
                shear_thickness_ratio=prop.tst,
            )



def _map_nodes(bdf_model: BDF, model: FEMModel):
    """
    @brief Translates pyNastran nodes to Node project class.
    """
    for nid, node in bdf_model.nodes.items():
        # node.get_position() returns coordinates in the Global Coordinate System,
        # resolving any local coordinate system (CP) logic automatically.
        
        # node.xyz returns coordinates in the global system if requested,
        # but “node.xyz” raw is in the local CD system.
        # Project Node class requests x1, x2, x3 raw and then calculates x,y,z.
        
        # pyNastran already solves the global position
        global_xyz = node.get_position()
        
        # Raw coords (as written in the file, relative to CP)
        raw_xyz = node.xyz
        
        new_node = Node(
            nid=nid,
            coord_id=node.cp,
            x1=raw_xyz[0],
            x2=raw_xyz[1],
            x3=raw_xyz[2],
            cd_id=node.cd,
            ps=node.ps if node.ps else "", # pyNastran may return None or ''
            seid=node.seid
        )
        
        # Assign the global coordinates pre-calculated by pyNastran.
        # Safer than assuming CP=0 as in manual nastran_reader.py
        new_node.x = global_xyz[0]
        new_node.y = global_xyz[1]
        new_node.z = global_xyz[2]

        if node.cd not in (None, 0):
            coord = node.cd_ref
            basis_global = np.column_stack([
                coord.transform_vector_to_global(np.eye(3)[axis]) for axis in range(3)
            ])
            new_node.dof_transform = basis_global.T
        
        model.nodes[nid] = new_node

def _map_elements(bdf_model: BDF, model: FEMModel):
    """
    @brief Translates elements (CROD, CBAR, CTRIA3, CQUAD4) to data_model project classes.
    @details Raises ValueError if unsupported elements are found.
    """
    # Set to avoid spamming the same warning 1000 times if there are 1000 identical elements
    unsupported_types = set()

    def material_angle_degrees(shell_element) -> float:
        _, _, element_x, element_y, _ = shell_element.element_coordinate_system()
        _, _, material_x, _, _ = shell_element.material_coordinate_system()
        return round(float(np.degrees(np.arctan2(
            np.dot(material_x, element_y), np.dot(material_x, element_x)
        ))), 12)

    for eid, elem in bdf_model.elements.items():
        
        # --- CROD ---
        if isinstance(elem, CROD):
            model.elements[eid] = RodElement(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids
            )
            
        # --- CBAR ---
        elif isinstance(elem, CBAR):
            # Orientation: pyNastran distinguishes between vector (x) or node (g0)
            orient_vec = None
            g0 = None
            
            if elem.g0 is not None:
                g0 = elem.g0
            else:
                # elem.x is the orientation vector
                if elem.x is not None:
                    orient_vec = list(elem.x)
            
            # Let pyNastran resolve G0/X, OFFT and the WA/WB coordinate flags.
            # The returned orientation and offsets are expressed in basic axes.
            axes_failed, axes = elem.get_axes(bdf_model)
            if axes_failed:
                raise ValueError(f"Unable to determine axes for CBAR {eid}.")
            _, _, yhat, _, wa_basic, wb_basic = axes
            orient_vec = list(yhat)
            wa = list(wa_basic)
            wb = list(wb_basic)
            
            model.elements[eid] = BarElement(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids,
                orientation_vector=orient_vec,
                g0=g0,
                offt=str(elem.offt), # usually a string or int in pyNastran
                pa=elem.pa,
                pb=elem.pb,
                wa=wa,
                wb=wb
            )

        # --- CTRIA3 ---
        elif isinstance(elem, CTRIA3):
            model.elements[eid] = CTRIA3Element(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids,
                theta=material_angle_degrees(elem)
            )
        
        # --- CQUAD4 ---
        elif isinstance(elem, CQUAD4):
            model.elements[eid] = CQUAD4Element(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids,
                theta=material_angle_degrees(elem)
            )

        # --- CATCH-ALL FOR UNSUPPORTED ELEMENTS ---
        else:
            # store the type of element
            unsupported_types.add(elem.type)
        
        # CRITICAL: Stop execution if the mesh contains elements we cannot solve.
    if unsupported_types:
        error_msg = (
            f"\n[FATAL ERROR] The model contains element types not supported by this solver:\n"
            f"              Found: {', '.join(unsupported_types)}\n"
            f"              Supported types are only: CROD, CBAR, CTRIA3, CQUAD4.\n"
            f"              Execution stopped to prevent incorrect physics."
        )
        print(error_msg)
        raise ValueError(error_msg)

def _map_loads(
    bdf_model: BDF,
    model: FEMModel,
    load_sid: Optional[int] = None,
    map_all_base_cards: bool = True,
):
    """
    @brief Translates loads (FORCE, MOMENT) to data_model objects.
    @details 
    Raises ValueError if unsupported load types (GRAV, PLOAD, etc.) are found.
    pyNastran groups loads by Set ID (SID) in bdf_model.loads.
    """
    unsupported_loads = set()
    if load_sid is not None and not (bdf_model.loads or bdf_model.load_combinations):
        selected_loads = []
    elif load_sid is not None:
        loads_list, scale_factors, _ = bdf_model.get_reduced_loads(
            load_sid, consider_load_combinations=True
        )
        selected_loads = zip(loads_list, scale_factors)
    elif map_all_base_cards:
        selected_loads = (
            (load, 1.0)
            for loads_list in bdf_model.loads.values()
            for load in loads_list
        )
    else:
        selected_loads = []

    for load, combination_scale in selected_loads:
            sid = load_sid if load_sid is not None else load.sid
            if load.type not in ('FORCE', 'MOMENT'):
                unsupported_loads.add(load.type)
                continue
            vector_global = np.asarray(load.xyz, dtype=float)
            if load.cid not in (None, 0):
                vector_global = load.cid_ref.transform_vector_to_global(vector_global)
            
            # pyNastran has a “type” attribute to identify the card
            if load.type == 'FORCE':
                # vector force: load.mag * load.xyz
                model.loads.append(SinglePointLoad(
                    sid=sid,
                    appnode_id=load.node_id,
                    module=combination_scale * load.mag,
                    vector=list(vector_global),
                    cid=0
                ))
            
            elif load.type == 'MOMENT':
                model.loads.append(SinglePointMoment(
                    sid=sid,
                    appnode_id=load.node_id,
                    module=combination_scale * load.mag,
                    vector=list(vector_global),
                    cid=0
                ))
            # --- CATCH-ALL FOR UNSUPPORTED LOADS ---
            else:
                unsupported_loads.add(load.type)

    if unsupported_loads:
        error_msg = (
            f"\n[FATAL ERROR] The model contains load types not supported by this solver:\n"
            f"              Found: {', '.join(unsupported_loads)}\n"
            f"              Supported types are only: FORCE, MOMENT.\n"
            f"              Execution stopped."
        )
        print(error_msg)
        raise ValueError(error_msg)

def _map_constraints(
    bdf_model: BDF,
    model: FEMModel,
    spc_sid: Optional[int] = None,
    map_all_base_cards: bool = True,
):
    """
    @brief Translates SPCs (Single Point Constraints) to data_model project classes.
    @details 
        Strict validation:
        1. Raises Error if unsupported SPC types are found (e.g., SPCAX).
        2. Raises Error if ANY MPC (Multi-Point Constraint) or Rigid Element (RBE2, RBE3) is found.
    """
    # pyNastran groups SPCs in bdf_model.spcs = {sid: [SPC1, SPC, ...]}
    # NOTE: Depending on the version and configuration, you can put them in bdf.spc_cards
    # Standard grouping interface used
    
    unsupported_constraints = set()

    # --- 1. CHECK FOR MULTI-POINT CONSTRAINTS (MPCs / RBEs) ---
    # pyNastran stores RBE2, RBE3, and MPC cards in the .mpcs attribute.
    # Since our solver (likely) currently lacks logic for Lagrange multipliers 
    # or master/slave elimination, we must forbid them.
    if bdf_model.mpcs or bdf_model.rigid_elements:
        unsupported_constraints.add("MPC / Rigid Elements (RBE2, RBE3)")

    for node in model.nodes.values():
        if node.ps:
            model.constraints.append(SPCConstraint(
                sid=0, node_id=node.nid, components=str(node.ps), value=0.0
            ))
    
    # --- 2. CHECK FOR SINGLE-POINT CONSTRAINTS (SPCs) ---
    if spc_sid is not None and not (bdf_model.spcs or bdf_model.spcadds):
        spc_groups = []
    elif spc_sid is not None:
        spc_groups = [(spc_sid, bdf_model.get_reduced_spcs(spc_sid, consider_spcadd=True))]
    elif map_all_base_cards:
        spc_groups = bdf_model.spcs.items()
    else:
        spc_groups = []

    for sid, spcs_list in spc_groups:
        for spc in spcs_list:
            if spc.type == 'SPC1':
                # SPC1 applies to a list of nodes with 0 enforced displacement
                comps = str(spc.components) # ej: "123"
                for nid in spc.node_ids:
                    model.constraints.append(SPCConstraint(
                        sid=sid,
                        node_id=nid,
                        components=comps,
                        value=0.0
                    ))
            
            elif spc.type == 'SPC':
                # SPC applies to individual nodes with specific enforced displacement
                # In pyNastran, an SPC object can contain multiple nodes.
                for i, nid in enumerate(spc.node_ids):
                    comps = str(spc.components[i])
                    val = spc.enforced[i]
                    model.constraints.append(SPCConstraint(
                        sid=sid,
                        node_id=nid,
                        components=comps,
                        value=val
                    ))
            # --- CATCH-ALL FOR UNSUPPORTED SPC TYPES ---
            else:
                unsupported_constraints.add(spc.type)
    
    # --- FINAL VALIDATION ---
    if unsupported_constraints:
        error_msg = (
            f"\n[FATAL ERROR] The model contains constraint types not supported by this solver:\n"
            f"              Found: {', '.join(unsupported_constraints)}\n"
            f"              Supported types are only: SPC, SPC1.\n"
            f"              (Note: Rigid elements/MPCs are not currently supported).\n"
            f"              Execution stopped."
        )
        print(error_msg)
        raise ValueError(error_msg)

def _select_case_control(bdf: BDF, requested_subcase_id: Optional[int]):
    """Select one analysis SUBCASE and return its LOAD/SPC set IDs."""
    analysis_subcases = sorted(sid for sid in bdf.subcases if sid != 0)
    if requested_subcase_id is not None:
        if requested_subcase_id not in analysis_subcases:
            raise ValueError(
                f"SUBCASE {requested_subcase_id} does not exist. Available: {analysis_subcases}"
            )
        selected_id = requested_subcase_id
    elif len(analysis_subcases) == 1:
        selected_id = analysis_subcases[0]
    elif len(analysis_subcases) > 1:
        raise ValueError(
            "The input contains multiple analysis subcases. Pass subcase_id explicitly "
            f"to parse_nastran_input(). Available: {analysis_subcases}"
        )
    else:
        global_subcase = bdf.subcases.get(0)
        if global_subcase is None or not ({'LOAD', 'SPC'} & set(global_subcase.params)):
            return None, None, None
        selected_id = 0

    subcase = bdf.subcases[selected_id]
    load_sid = subcase.get_parameter('LOAD')[0] if 'LOAD' in subcase.params else None
    spc_sid = subcase.get_parameter('SPC')[0] if 'SPC' in subcase.params else None
    return selected_id, load_sid, spc_sid


def parse_nastran_input(entry_file_path: str, subcase_id: Optional[int] = None) -> FEMModel:
    """
    @brief Main entry point using pyNastran.
    @param entry_file_path Path to the .bdf or .dat file.
    @return Populated FEMModel instance.
    """
    print(f"[INFO] Starting read with pyNastran: {entry_file_path}")
    
    # 1. Instantiate BDF and read file
    # debug=False reduces console noise
    bdf = BDF(debug=False)
    

    try:
        # ATTEMPT 1: Read as a Standard Complete Nastran File (Default)
        # This handles files with 'SOL 101', 'TIME', 'CEND', 'BEGIN BULK'.
        bdf.read_bdf(entry_file_path, validate=True, xref=True, punch=False)
        
    except MissingDeckSections:
        # ATTEMPT 2: Fallback for "Snippet" or "Bulk-Only" files
        # This handles unit test files that lack headers.
        # We must create a new instance to avoid internal state pollution.
        # validate=True checks card integrity
        # xref=True is CRITICAL: resolves relationships (Material -> Prop -> Element)
        # punch = True allow pynastran to analyce a non compleate bdf, without the Executive control in it.
        #         the same as giiving a BULK DATA block. But BULK DATA must not be written in the .bdf file.          
        #         Nastran would raise an error while reading this file, pyNastran also stops the execution 
        #         if this is not marked.
        print("[INFO] Missing control decks detected. Retrying as 'punch' file (Bulk Data only)...")
        bdf = BDF(debug=False)
        bdf.read_bdf(entry_file_path, validate=True, xref=True, punch=True)
        
    except Exception as e:
        # Any other critical error (syntax, missing file, etc.)
        print(f"[FAIL] Error while trying to read nastran input file ({entry_file_path}): {e}")
        raise e
    
    
    print(f"[INFO] File {entry_file_path} parsed successfully. Mapping to FEMModel...")
    
    # 2. Initialize our container
    model = FEMModel()
    selected_subcase, load_sid, spc_sid = _select_case_control(bdf, subcase_id)
    model.active_subcase_id = selected_subcase
    model.active_load_sid = load_sid
    model.active_spc_sid = spc_sid

    # --- EXTRAER PARÁMETROS DEL DECK DE NASTRAN ---
    if 'K6ROT' in bdf.params:
        # pyNastran guarda el valor del parámetro en la lista '.values'
        model.k6rot = float(bdf.params['K6ROT'].values[0])
        print(f"[INFO] PARAM, K6ROT detectado en el archivo de entrada. Valor: {model.k6rot}")
    
    # 3. Entity Mapping (Logical Order: Mat -> Prop -> Node -> Elem -> Loads)
    try:
        _map_materials(bdf, model)
        _map_properties(bdf, model)
        _map_nodes(bdf, model)
        _map_elements(bdf, model)
        bulk_only = selected_subcase is None
        _map_loads(bdf, model, load_sid=load_sid, map_all_base_cards=bulk_only)
        _map_constraints(bdf, model, spc_sid=spc_sid, map_all_base_cards=bulk_only)
        model.validate(require_elements=False)
        
    except Exception as e:
        print(f"[ERROR] Error during data mapping: {e}")
        raise e

    print(f"[SUCCESS] Model successfully loaded.")
    print(f"  - Nodes: {len(model.nodes)}")
    print(f"  - Elements: {len(model.elements)}")
    
    return model



