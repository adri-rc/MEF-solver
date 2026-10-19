"""Read supported Nastran decks into the format-independent FEM model."""

from typing import Optional

import numpy as np

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
except ImportError as exc:
    raise ImportError(
        "'pyNastran' is required to read Nastran files. Install it with: pip install pyNastran"
    ) from exc

from .data_model import (
    AnalysisConfiguration,
    AnalysisType,
    BarProperty,
    BarElement,
    CQUAD4Element,
    CTRIA3Element,
    FEMModel,
    IsotropicMaterial,
    Node,
    NonlinearStepControl,
    RodElement,
    RodProperty,
    ShellProperty,
    SinglePointLoad,
    SinglePointMoment,
    SPCConstraint,
)


def _get_source_entry_counts(bdf_model: BDF) -> dict[str, int]:
    """Return F06-like logical Bulk Data counts for report presentation.

    MSC consolidates SPC1 cards that share the same SID and component string
    in MODEL SUMMARY.  pyNastran retains the physical-card count, so SPC1 is
    normalized here without affecting the constraints used by the solver.
    """
    counts = {
        str(name).upper(): int(count)
        for name, count in bdf_model.card_count.items()
        if name != "ENDDATA" and count
    }
    spc1_entries = {
        (int(sid), str(spc.components))
        for sid, cards in bdf_model.spcs.items()
        for spc in cards
        if spc.type == "SPC1"
    }
    if spc1_entries:
        counts["SPC1"] = len(spc1_entries)
    return dict(sorted(counts.items()))

def _map_materials(bdf_model: BDF, model: FEMModel) -> None:
    """Map supported pyNastran materials to project data classes."""
    for mid, mat in bdf_model.materials.items():
        if isinstance(mat, MAT1):
            # E/G use None to trigger completion from the other constants.
            E_val = mat.e if mat.e != 0.0 else None
            G_val = mat.g if mat.g != 0.0 else None
            # Unlike E/G, a zero Poisson ratio is valid and must be preserved.
            nu_val = mat.nu
            rho_val = mat.rho

            try:
                model.materials[mid] = IsotropicMaterial(
                    mid=mid, E=E_val, G=G_val, nu=nu_val, rho=rho_val
                )
            except ValueError as exc:
                print(f"\n[ERROR] Material ID {mid} is invalid/inconsistent.")
                print(f"              Reason: {exc}")
                print("              Execution stopped to prevent incorrect results.\n")
                raise


def _map_properties(bdf_model: BDF, model: FEMModel) -> None:
    """Map supported PROD, PBAR and PSHELL properties."""
    for pid, prop in bdf_model.properties.items():
        if isinstance(prop, PROD):
            model.properties[pid] = RodProperty(
                pid=pid,
                mat_id=prop.Mid(),
                A=prop.Area(),
                J=prop.J(),
                C=prop.c,
                NSM=prop.nsm,
            )
        elif isinstance(prop, PBAR):
            def get_stress_point(prefix: str) -> tuple[float, float]:
                """Return a PBAR C--F recovery point as section (y, z)."""
                y = getattr(prop, f"{prefix}1", 0.0) or 0.0
                z = getattr(prop, f"{prefix}2", 0.0) or 0.0
                return float(y), float(z)

            stress_points = [get_stress_point(label) for label in "cdef"]
            model.properties[pid] = BarProperty(
                pid=pid,
                mat_id=prop.Mid(),
                A=prop.Area(),
                I1=prop.I11(),
                I2=prop.I22(),
                J=prop.J(),
                I12=prop.I12(),
                NSM=prop.nsm,
                K1=getattr(prop, "k1", 0.0) or 0.0,
                K2=getattr(prop, "k2", 0.0) or 0.0,
                stress_points=stress_points,
            )

        elif isinstance(prop, PSHELL):
            mid1_val = prop.Mid1() or None
            mid2_val = prop.Mid2() or None
            mid3_val = prop.Mid3() or None
            model.properties[pid] = ShellProperty(
                pid=pid,
                mat_id=mid1_val,
                t=prop.t,
                mat_id_bending=mid2_val,
                mat_id_shear=mid3_val,
                bending_inertia_ratio=prop.twelveIt3,
                shear_thickness_ratio=prop.tst,
            )


def _map_nodes(bdf_model: BDF, model: FEMModel) -> None:
    """Map GRID entries and their CP/CD transformations."""
    for nid, node in bdf_model.nodes.items():
        # get_position() resolves GRID.CP into the basic coordinate system;
        # xyz retains the coordinates as entered in CP.
        global_xyz = node.get_position()
        raw_xyz = node.xyz

        new_node = Node(
            nid=nid,
            coord_id=node.cp,
            x1=raw_xyz[0],
            x2=raw_xyz[1],
            x3=raw_xyz[2],
            cd_id=node.cd,
            ps=node.ps if node.ps else "",
            seid=node.seid,
        )

        new_node.x = global_xyz[0]
        new_node.y = global_xyz[1]
        new_node.z = global_xyz[2]

        if node.cd not in (None, 0):
            coord = node.cd_ref
            basis_global = np.column_stack(
                [
                    coord.transform_vector_to_global(np.eye(3)[axis])
                    for axis in range(3)
                ]
            )
            new_node.dof_transform = basis_global.T

        model.nodes[nid] = new_node


def _map_elements(bdf_model: BDF, model: FEMModel) -> None:
    """Map supported elements and reject every unsupported element family."""
    unsupported_types = set()

    def material_angle_degrees(shell_element) -> float:
        _, _, element_x, element_y, _ = shell_element.element_coordinate_system()
        _, _, material_x, _, _ = shell_element.material_coordinate_system()
        angle = np.arctan2(
            np.dot(material_x, element_y), np.dot(material_x, element_x)
        )
        return round(float(np.degrees(angle)), 12)

    for eid, elem in bdf_model.elements.items():
        if isinstance(elem, CROD):
            model.elements[eid] = RodElement(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids,
            )
        elif isinstance(elem, CBAR):
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
                g0=elem.g0,
                offt=str(elem.offt),
                pa=elem.pa,
                pb=elem.pb,
                wa=wa,
                wb=wb,
            )

        elif isinstance(elem, CTRIA3):
            model.elements[eid] = CTRIA3Element(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids,
                theta=material_angle_degrees(elem),
            )
        elif isinstance(elem, CQUAD4):
            model.elements[eid] = CQUAD4Element(
                eid=eid,
                prop_id=elem.Pid(),
                nodes=elem.node_ids,
                theta=material_angle_degrees(elem),
            )
        else:
            unsupported_types.add(elem.type)

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
) -> None:
    """Map the selected FORCE/MOMENT entries and reject other load types."""
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
        if load.type not in {"FORCE", "MOMENT"}:
            unsupported_loads.add(load.type)
            continue

        vector_global = np.asarray(load.xyz, dtype=float)
        if load.cid not in (None, 0):
            vector_global = load.cid_ref.transform_vector_to_global(vector_global)

        load_arguments = {
            "sid": sid,
            "appnode_id": load.node_id,
            "module": combination_scale * load.mag,
            "vector": list(vector_global),
            "cid": 0,
        }
        if load.type == "FORCE":
            model.loads.append(SinglePointLoad(**load_arguments))
        else:
            model.loads.append(SinglePointMoment(**load_arguments))

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
) -> None:
    """Map SPC/SPC1 entries and reject unsupported kinematic constraints."""
    unsupported_constraints = set()

    # The solver has neither generalized DOF-set reduction, Lagrange
    # multipliers, nor master/slave elimination. pyNastran stores these cards
    # outside ``spcs``, so each collection must be checked explicitly.
    unsupported_collections = {
        "MPC / MPCADD": ("mpcs", "mpcadds"),
        "Rigid Elements (RBE/RBAR families)": ("rigid_elements",),
        "SUPORT / SUPORT1": ("suport", "suport1"),
        "ASET / BSET / CSET / QSET / USET / OMIT": (
            "asets",
            "bsets",
            "csets",
            "qsets",
            "usets",
            "omits",
        ),
        "SPCOFF / SPCOFF1": ("spcoffs",),
        "Superelement DOF sets": (
            "se_bsets",
            "se_csets",
            "se_qsets",
            "se_sets",
            "se_suport",
            "se_usets",
            "senqset",
            "release",
        ),
    }
    for label, attribute_names in unsupported_collections.items():
        if any(getattr(bdf_model, name, None) for name in attribute_names):
            unsupported_constraints.add(label)

    for node in model.nodes.values():
        if node.ps:
            model.constraints.append(
                SPCConstraint(
                    sid=0,
                    node_id=node.nid,
                    components=str(node.ps),
                    value=0.0,
                )
            )

    if spc_sid is not None and not (bdf_model.spcs or bdf_model.spcadds):
        spc_groups = []
    elif spc_sid is not None:
        spc_groups = [
            (
                spc_sid,
                bdf_model.get_reduced_spcs(spc_sid, consider_spcadd=True),
            )
        ]
    elif map_all_base_cards:
        spc_groups = bdf_model.spcs.items()
    else:
        spc_groups = []

    for sid, spcs_list in spc_groups:
        for spc in spcs_list:
            if spc.type == "SPC1":
                components = str(spc.components)
                for nid in spc.node_ids:
                    model.constraints.append(
                        SPCConstraint(
                            sid=sid,
                            node_id=nid,
                            components=components,
                            value=0.0,
                        )
                    )
            elif spc.type == "SPC":
                for index, nid in enumerate(spc.node_ids):
                    model.constraints.append(
                        SPCConstraint(
                            sid=sid,
                            node_id=nid,
                            components=str(spc.components[index]),
                            value=spc.enforced[index],
                        )
                    )
            else:
                unsupported_constraints.add(spc.type)

    if unsupported_constraints:
        error_msg = (
            f"\n[FATAL ERROR] The model contains constraint types not supported by this solver:\n"
            f"              Found: {', '.join(unsupported_constraints)}\n"
            f"              Supported kinematic constraints are only: SPC, SPC1, "
            f"SPCADD, and GRID.PS.\n"
            f"              Execution stopped."
        )
        print(error_msg)
        raise ValueError(error_msg)

def _select_case_control(
    bdf: BDF, requested_subcase_id: Optional[int]
) -> tuple[Optional[int], Optional[int], Optional[int]]:
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
        if global_subcase is None or not (
            {"LOAD", "SPC"} & set(global_subcase.params)
        ):
            return None, None, None
        selected_id = 0

    subcase = bdf.subcases[selected_id]
    load_sid = (
        subcase.get_parameter("LOAD")[0] if "LOAD" in subcase.params else None
    )
    spc_sid = (
        subcase.get_parameter("SPC")[0] if "SPC" in subcase.params else None
    )
    return selected_id, load_sid, spc_sid


def _nastran_float(value: str) -> float:
    """Parse ordinary and compact Nastran floating-point notation."""
    normalized = value.strip().upper().replace("D", "E")
    if "E" not in normalized:
        for index in range(1, len(normalized)):
            if normalized[index] in "+-" and any(
                character.isdigit() for character in normalized[:index]
            ):
                normalized = f"{normalized[:index]}E{normalized[index:]}"
                break
    return float(normalized)


def _rejected_card_rows(card_lines) -> list[list[str]]:
    """Return one field row per line of an unsupported pyNastran card."""
    rows = []
    for raw_line in card_lines:
        line = raw_line.split("$", 1)[0].rstrip()
        if not line.strip():
            continue
        if "," in line:
            fields = [field.strip() for field in line.split(",")]
        else:
            fields = [
                line[index : index + 8].strip()
                for index in range(0, len(line), 8)
            ]
        if fields:
            fields = fields[1:]
        while fields and fields[-1] == "":
            fields.pop()
        rows.append(fields)
    return rows


def _field(row: list[str], index: int, default=None):
    if index >= len(row) or row[index] == "":
        return default
    return row[index]


def _parse_nlstep_cards(bdf: BDF) -> dict[int, NonlinearStepControl]:
    """Parse MSC NLSTEP cards rejected by pyNastran 1.4.1."""
    controls = {}
    for rejected_card in bdf.reject_lines:
        rows = _rejected_card_rows(rejected_card)
        if not rows:
            continue
        # pyNastran may attach comments immediately preceding an unsupported
        # card to the same reject_lines group.  Identify the first physical
        # card line after removing comments instead of treating '$' as the
        # rejected card name.
        first_raw = next(
            (
                content
                for raw_line in rejected_card
                if (content := raw_line.split("$", 1)[0].strip())
            ),
            "",
        )
        card_name = first_raw.split(",", 1)[0].strip().split()[0].upper()
        if card_name != "NLSTEP":
            continue

        header = rows[0]
        if not header or _field(header, 0) is None:
            raise ValueError("NLSTEP requires a positive ID.")
        sid = int(_field(header, 0))
        values = {
            "sid": sid,
            "source": "NLSTEP",
            "total_time": _nastran_float(_field(header, 1, "1.0")),
            "control_definition": _field(header, 2),
            "increment_method": "FIXED",
            "number_of_increments": 50,
            "output_intervals": 0,
            "maximum_iterations": 10,
            "minimum_iterations": 1,
            "maximum_bisections": 10,
            "convergence": "PV",
            "stiffness_update_method": "PFNT",
        }

        increment_keyword_seen = False
        for row in rows[1:]:
            if not row:
                continue
            keyword = row[0].upper()
            if keyword == "GENERAL":
                values["maximum_iterations"] = int(_field(row, 1, "10"))
                values["minimum_iterations"] = int(_field(row, 2, "1"))
                values["maximum_bisections"] = int(_field(row, 3, "10"))
                creep = int(_field(row, 4, "0"))
                if creep != 0:
                    raise ValueError(f"NLSTEP {sid} requests unsupported creep analysis.")
            elif keyword == "FIXED":
                if increment_keyword_seen:
                    raise ValueError(f"NLSTEP {sid} defines more than one stepping method.")
                increment_keyword_seen = True
                values["increment_method"] = "FIXED"
                values["number_of_increments"] = int(_field(row, 1, "50"))
            elif keyword == "ADAPT":
                if increment_keyword_seen:
                    raise ValueError(f"NLSTEP {sid} defines more than one stepping method.")
                increment_keyword_seen = True
                values["increment_method"] = "ADAPT"
                values["number_of_increments"] = None
                values["initial_step_fraction"] = _nastran_float(_field(row, 1, "0.01"))
                values["minimum_step_fraction"] = _nastran_float(_field(row, 2, "1.E-5"))
                values["maximum_step_fraction"] = _nastran_float(_field(row, 3, "0.5"))
                values["desired_iterations"] = int(_field(row, 4, "4"))
                values["step_scale_factor"] = _nastran_float(_field(row, 5, "1.2"))
                values["output_intervals"] = int(_field(row, 6, "0"))
                values["maximum_increments"] = int(_field(row, 7, "99999"))
            elif keyword == "MECH":
                values["convergence"] = _field(row, 1, "PV")
                eps_u = _field(row, 2)
                eps_p = _field(row, 3)
                eps_w = _field(row, 4)
                values["displacement_tolerance"] = (
                    None if eps_u is None else _nastran_float(eps_u)
                )
                values["residual_tolerance"] = (
                    None if eps_p is None else _nastran_float(eps_p)
                )
                values["work_tolerance"] = (
                    None if eps_w is None else _nastran_float(eps_w)
                )
                values["stiffness_update_method"] = _field(row, 5, "PFNT")
                kstep = _field(row, 6)
                values["stiffness_update_interval"] = (
                    None if kstep is None else int(kstep)
                )
            else:
                raise ValueError(
                    f"NLSTEP {sid} uses unsupported continuation keyword {keyword}. "
                    "Supported keywords are GENERAL, FIXED, ADAPT and MECH."
                )
        if sid in controls:
            raise ValueError(f"Duplicate NLSTEP ID {sid}.")
        controls[sid] = NonlinearStepControl(**values)
    return controls


def _map_pbarn1_extensions(bdf: BDF, model: FEMModel) -> None:
    """Attach the SOL 400 PBARN1 options rejected by pyNastran 1.4.1.

    PBARN1 supplements, rather than replaces, the PBAR with the same PID. Its
    fields are retained in the project model so that unsupported nonlinear
    section/integration variants cannot be mistaken for an ordinary PBAR.
    """
    for rejected_card in bdf.reject_lines:
        rows = _rejected_card_rows(rejected_card)
        if not rows:
            continue
        first_raw = next(
            (
                content
                for raw_line in rejected_card
                if (content := raw_line.split("$", 1)[0].strip())
            ),
            "",
        )
        card_name = first_raw.split(",", 1)[0].strip().split()[0].upper()
        if card_name != "PBARN1":
            continue

        header = rows[0]
        pid_value = _field(header, 0)
        if pid_value is None:
            raise ValueError("PBARN1 requires the PID of its associated PBAR.")
        pid = int(pid_value)
        prop = model.properties.get(pid)
        if not isinstance(prop, BarProperty):
            raise ValueError(f"PBARN1 {pid} does not reference a mapped PBAR {pid}.")

        prop.nonlinear_section = _field(header, 3)
        for row in rows[1:]:
            if not row:
                continue
            keyword = row[0].upper()
            if keyword != "C2":
                raise ValueError(
                    f"PBARN1 {pid} uses unsupported continuation keyword {keyword}; "
                    "only C2 is currently mapped."
                )
            prop.nonlinear_behavior = _field(row, 1)
            prop.nonlinear_integration = _field(row, 2)

        for attribute in (
            "nonlinear_section",
            "nonlinear_behavior",
            "nonlinear_integration",
        ):
            value = getattr(prop, attribute)
            if value is not None:
                setattr(prop, attribute, str(value).upper())


def _legacy_nlparm_control(nlparm) -> NonlinearStepControl:
    output_intervals = -1 if nlparm.int_out == "NO" else 0
    return NonlinearStepControl(
        sid=nlparm.nlparm_id,
        source="NLPARM",
        increment_method="FIXED",
        number_of_increments=nlparm.ninc or 10,
        output_intervals=output_intervals,
        maximum_iterations=nlparm.max_iter,
        minimum_iterations=1,
        maximum_bisections=abs(nlparm.max_bisect),
        convergence=nlparm.conv,
        displacement_tolerance=nlparm.eps_u,
        residual_tolerance=nlparm.eps_p,
        work_tolerance=nlparm.eps_w,
        stiffness_update_method=nlparm.kmethod,
        stiffness_update_interval=nlparm.kstep,
    )


def _map_analysis_configuration(
    bdf: BDF,
    model: FEMModel,
    selected_subcase_id: Optional[int],
) -> None:
    """Map a supported SOL 101 or SOL 400 analysis configuration.

    Bulk-only punch files have no solution sequence and retain the linear
    default. Any other sequence or inconsistent ``SOL``/``ANALYSIS`` pair is
    rejected before it can be routed to the wrong solver.
    """
    solution_sequence = int(bdf.sol) if bdf.sol is not None else None
    subcase = (
        bdf.subcases.get(selected_subcase_id)
        if selected_subcase_id is not None
        else None
    )
    analysis_keyword = None
    if subcase is not None and "ANALYSIS" in subcase.params:
        analysis_keyword = str(subcase.get_parameter("ANALYSIS")[0]).upper()

    if solution_sequence is None and not bdf.punch:
        raise ValueError(
            "A complete Nastran deck requires SOL 101 or SOL 400; only a "
            "Bulk Data-only punch file may omit the solution sequence."
        )

    if solution_sequence not in {None, 101, 400}:
        raise ValueError(
            f"Unsupported solution sequence SOL {solution_sequence}. "
            "Supported sequences are SOL 101 and SOL 400."
        )

    nonlinear_aliases = {"NLSTAT", "NLSTATIC", "NLSTATICS"}
    if solution_sequence == 400:
        if analysis_keyword not in nonlinear_aliases:
            raise ValueError(
                "SOL 400 requires ANALYSIS=NLSTAT for the supported nonlinear "
                "static subset."
            )
        is_nonlinear_static = True
    else:
        if analysis_keyword in nonlinear_aliases:
            raise ValueError(
                f"ANALYSIS={analysis_keyword} requires SOL 400; "
                f"found SOL {solution_sequence}."
            )
        if analysis_keyword not in {None, "STATICS"}:
            raise ValueError(
                f"Unsupported ANALYSIS={analysis_keyword} for SOL {solution_sequence}."
            )
        is_nonlinear_static = False

    analysis_type = (
        AnalysisType.NONLINEAR_STATIC if is_nonlinear_static else AnalysisType.LINEAR_STATIC
    )

    lgdisp_value = None
    large_displacements = False
    if "LGDISP" in bdf.params:
        lgdisp_value = int(bdf.params["LGDISP"].values[0])
        if lgdisp_value not in {-1, 0, 1, 2}:
            raise ValueError(f"Unsupported PARAM,LGDISP value: {lgdisp_value}")
        large_displacements = lgdisp_value in {1, 2}
    follower_force_stiffness = lgdisp_value == 1
    if "FOLLOWK" in bdf.params:
        followk_value = str(bdf.params["FOLLOWK"].values[0]).strip().upper()
        if followk_value in {"YES", "Y", "1"}:
            follower_force_stiffness = True
        elif followk_value in {"NO", "N", "0"}:
            follower_force_stiffness = False
        else:
            raise ValueError(f"Unsupported PARAM,FOLLOWK value: {followk_value}")

    nonlinear_step = None
    if is_nonlinear_static:
        nlstep_sid = None
        nlparm_sid = None
        if subcase is not None and "NLSTEP" in subcase.params:
            nlstep_sid = int(subcase.get_parameter("NLSTEP")[0])
        if subcase is not None and "NLPARM" in subcase.params:
            nlparm_sid = int(subcase.get_parameter("NLPARM")[0])
        if nlstep_sid is not None and nlparm_sid is not None:
            raise ValueError(
                "Select either NLSTEP or NLPARM, not both, for one nonlinear subcase."
            )
        if nlstep_sid is not None:
            nlsteps = _parse_nlstep_cards(bdf)
            if nlstep_sid not in nlsteps:
                raise ValueError(
                    f"Case Control selects NLSTEP {nlstep_sid}, but that Bulk Data entry was not found."
                )
            nonlinear_step = nlsteps[nlstep_sid]
        elif nlparm_sid is not None:
            if nlparm_sid not in bdf.nlparms:
                raise ValueError(
                    f"Case Control selects NLPARM {nlparm_sid}, but that Bulk Data entry was not found."
                )
            nonlinear_step = _legacy_nlparm_control(bdf.nlparms[nlparm_sid])
        else:
            raise ValueError(
                "Nonlinear static analysis requires an NLSTEP (preferred for SOL 400) "
                "or NLPARM Case Control selection."
            )

    model.analysis = AnalysisConfiguration(
        solution_sequence=solution_sequence,
        analysis_type=analysis_type,
        analysis_keyword=analysis_keyword,
        large_displacement_parameter=lgdisp_value,
        large_displacements=large_displacements,
        follower_force_stiffness=follower_force_stiffness,
        nonlinear_step=nonlinear_step,
    )


def parse_nastran_input(
    entry_file_path: str, subcase_id: Optional[int] = None
) -> FEMModel:
    """Parse a Nastran ``.dat``/``.bdf`` file into an :class:`FEMModel`."""
    print(f"[INFO] Starting read with pyNastran: {entry_file_path}")

    bdf = BDF(debug=False)

    try:
        bdf.read_bdf(entry_file_path, validate=True, xref=True, punch=False)
    except MissingDeckSections:
        # Recreate the BDF to avoid retaining state from the failed full-deck read.
        print(
            "[INFO] Missing control decks detected. Retrying as 'punch' file "
            "(Bulk Data only)..."
        )
        bdf = BDF(debug=False)
        bdf.read_bdf(entry_file_path, validate=True, xref=True, punch=True)
    except Exception as exc:
        print(
            "[FAIL] Error while trying to read Nastran input file "
            f"({entry_file_path}): {exc}"
        )
        raise

    print(f"[INFO] File {entry_file_path} parsed successfully. Mapping to FEMModel...")

    model = FEMModel()
    selected_subcase, load_sid, spc_sid = _select_case_control(bdf, subcase_id)
    model.active_subcase_id = selected_subcase
    model.active_load_sid = load_sid
    model.active_spc_sid = spc_sid
    model.source_entry_counts = _get_source_entry_counts(bdf)
    _map_analysis_configuration(bdf, model, selected_subcase)

    if "K6ROT" in bdf.params:
        model.k6rot = float(bdf.params["K6ROT"].values[0])
        print(f"[INFO] Input PARAM,K6ROT value: {model.k6rot}")

    try:
        _map_materials(bdf, model)
        _map_properties(bdf, model)
        _map_pbarn1_extensions(bdf, model)
        _map_nodes(bdf, model)
        _map_elements(bdf, model)
        bulk_only = selected_subcase is None
        _map_loads(bdf, model, load_sid=load_sid, map_all_base_cards=bulk_only)
        _map_constraints(bdf, model, spc_sid=spc_sid, map_all_base_cards=bulk_only)
        model.validate(require_elements=False)
    except Exception as exc:
        print(f"[ERROR] Error during data mapping: {exc}")
        raise

    print("[SUCCESS] Model successfully loaded.")
    print(f"  - Nodes: {len(model.nodes)}")
    print(f"  - Elements: {len(model.elements)}")
    return model
