"""Generate the updated final nonlinear mixed CROD/CBAR bridge deck.

The model is a spatial three-span through-arch bridge with explicit abutment and
pier foundations.  It stays inside the implemented SOL 400 subset: elastic
MAT1 materials, PBAR/PBARN1/CBAR members, PROD/CROD members, nodal conservative
loads and explicit SPC constraints.  Units are mm, N and MPa.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_OUTPUT = Path(__file__).with_name("tfm_final_nonlinear_bridge_updated.dat")
STATION_SPACING = 1_000.0
N_STATIONS = 31
DECK_Y = (-1_800.0, -600.0, 600.0, 1_800.0)
SUPPORT_INDICES = (0, 8, 20, 30)
SPAN_RISE = {(0, 8): 1_900.0, (8, 20): 3_200.0, (20, 30): 2_400.0}
SUPPORT_LOWER_Z = {0: -900.0, 8: -1_800.0, 20: -2_200.0, 30: -1_100.0}
SUPPORT_GROUND_Z = {0: -1_800.0, 8: -3_200.0, 20: -3_600.0, 30: -2_100.0}
SUPPORT_SETTLEMENT = {0: 0.0, 8: -2.0, 20: -6.0, 30: -10.0}


@dataclass
class BridgeModel:
    """Deterministic geometry, connectivity and named groups."""

    nodes: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    bars: list[tuple[int, int, int, int, tuple[float, float, float]]] = field(
        default_factory=list
    )
    rods: list[tuple[int, int, int, int]] = field(default_factory=list)
    deck_lines: list[list[int]] = field(default_factory=lambda: [[], [], [], []])
    arch_negative: list[int] = field(default_factory=list)
    arch_positive: list[int] = field(default_factory=list)
    support_lower: dict[int, tuple[int, int]] = field(default_factory=dict)
    support_ground: dict[int, list[int]] = field(default_factory=dict)
    inner_deck_nodes: list[int] = field(default_factory=list)

    def add_bar(
        self,
        element_id: int,
        property_id: int,
        node_a: int,
        node_b: int,
        orientation: tuple[float, float, float],
    ) -> None:
        self.bars.append((element_id, property_id, node_a, node_b, orientation))

    def add_rod(
        self, element_id: int, property_id: int, node_a: int, node_b: int
    ) -> None:
        self.rods.append((element_id, property_id, node_a, node_b))


def _arch_height(station_index: int) -> float:
    for (start, end), rise in SPAN_RISE.items():
        if start <= station_index <= end:
            xi = (station_index - start) / (end - start)
            return 4.0 * rise * xi * (1.0 - xi)
    raise ValueError(f"Station {station_index} does not belong to an arch span")


def _deck_grid_id(line_index: int, station_index: int) -> int:
    return 1001 + 100 * line_index + station_index


def build_model() -> BridgeModel:
    model = BridgeModel()

    # Four deck lines provide a recognizable rail/deck platform.  The two
    # outer tie girders are CBAR; the two inner stringers are pin-jointed CROD.
    for line_index, y in enumerate(DECK_Y):
        for station_index in range(N_STATIONS):
            grid_id = _deck_grid_id(line_index, station_index)
            x = station_index * STATION_SPACING
            model.nodes[grid_id] = (x, y, 0.0)
            model.deck_lines[line_index].append(grid_id)
            if line_index in (1, 2):
                model.inner_deck_nodes.append(grid_id)

    # Two arch planes share the outer deck GRID at the four springing points.
    for station_index in range(N_STATIONS):
        if station_index in SUPPORT_INDICES:
            model.arch_negative.append(model.deck_lines[0][station_index])
            model.arch_positive.append(model.deck_lines[3][station_index])
            continue
        x = station_index * STATION_SPACING
        z = _arch_height(station_index)
        negative_id = 2001 + station_index
        positive_id = 2101 + station_index
        model.nodes[negative_id] = (x, DECK_Y[0], z)
        model.nodes[positive_id] = (x, DECK_Y[3], z)
        model.arch_negative.append(negative_id)
        model.arch_positive.append(positive_id)

    # Abutments and piers: each deck side rests on a vertical CBAR column that
    # opens into two inclined CBAR legs and four explicit ground-footing GRID.
    for support_number, station_index in enumerate(SUPPORT_INDICES):
        x = station_index * STATION_SPACING
        lower_z = SUPPORT_LOWER_Z[station_index]
        negative_lower = 3001 + 10 * support_number
        positive_lower = negative_lower + 1
        model.nodes[negative_lower] = (x, DECK_Y[0], lower_z)
        model.nodes[positive_lower] = (x, DECK_Y[3], lower_z)
        model.support_lower[station_index] = (negative_lower, positive_lower)

        ground_ids: list[int] = []
        for ground_number, y in enumerate((-2_800.0, -1_000.0, 1_000.0, 2_800.0)):
            ground_id = 4001 + 10 * support_number + ground_number
            model.nodes[ground_id] = (x, y, SUPPORT_GROUND_Z[station_index])
            ground_ids.append(ground_id)
        model.support_ground[station_index] = ground_ids

    # Outer deck tie girders.
    for station_index in range(N_STATIONS - 1):
        model.add_bar(
            10001 + station_index,
            101,
            model.deck_lines[0][station_index],
            model.deck_lines[0][station_index + 1],
            (0.0, 0.0, 1.0),
        )
        model.add_bar(
            10101 + station_index,
            101,
            model.deck_lines[3][station_index],
            model.deck_lines[3][station_index + 1],
            (0.0, 0.0, 1.0),
        )

    # Three consecutive parabolic spans in both arch planes.
    for station_index in range(N_STATIONS - 1):
        model.add_bar(
            11001 + station_index,
            103,
            model.arch_negative[station_index],
            model.arch_negative[station_index + 1],
            (0.0, 1.0, 0.0),
        )
        model.add_bar(
            11101 + station_index,
            103,
            model.arch_positive[station_index],
            model.arch_positive[station_index + 1],
            (0.0, 1.0, 0.0),
        )

    # Pier/abutment columns and their splayed foundation legs.
    column_eid = 12001
    leg_eid = 12101
    for station_index in SUPPORT_INDICES:
        lower_negative, lower_positive = model.support_lower[station_index]
        model.add_bar(
            column_eid,
            105,
            model.deck_lines[0][station_index],
            lower_negative,
            (1.0, 0.0, 0.0),
        )
        column_eid += 1
        model.add_bar(
            column_eid,
            105,
            model.deck_lines[3][station_index],
            lower_positive,
            (1.0, 0.0, 0.0),
        )
        column_eid += 1

        ground_0, ground_1, ground_2, ground_3 = model.support_ground[station_index]
        for lower, ground in (
            (lower_negative, ground_0),
            (lower_negative, ground_1),
            (lower_positive, ground_2),
            (lower_positive, ground_3),
        ):
            model.add_bar(leg_eid, 106, lower, ground, (1.0, 0.0, 0.0))
            leg_eid += 1

    # Inner deck stringers.
    for station_index in range(N_STATIONS - 1):
        model.add_rod(
            20001 + station_index,
            202,
            model.deck_lines[1][station_index],
            model.deck_lines[1][station_index + 1],
        )
        model.add_rod(
            20101 + station_index,
            202,
            model.deck_lines[2][station_index],
            model.deck_lines[2][station_index + 1],
        )

    # Transverse deck members connect the four longitudinal lines.
    transverse_eid = 20201
    for station_index in range(N_STATIONS):
        for line_index in range(3):
            model.add_rod(
                transverse_eid,
                203,
                model.deck_lines[line_index][station_index],
                model.deck_lines[line_index + 1][station_index],
            )
            transverse_eid += 1

    # Bilateral vertical hangers; CROD is not treated as a tension-only cable.
    hanger_eid = 20401
    for station_index in range(N_STATIONS):
        if station_index in SUPPORT_INDICES:
            continue
        model.add_rod(
            hanger_eid,
            201,
            model.deck_lines[0][station_index],
            model.arch_negative[station_index],
        )
        hanger_eid += 1
        model.add_rod(
            hanger_eid,
            201,
            model.deck_lines[3][station_index],
            model.arch_positive[station_index],
        )
        hanger_eid += 1

    # Transverse ties and alternating wind diagonals between the arch planes.
    tie_eid = 20501
    for station_index in range(N_STATIONS):
        if station_index in SUPPORT_INDICES:
            continue
        model.add_rod(
            tie_eid,
            205,
            model.arch_negative[station_index],
            model.arch_positive[station_index],
        )
        tie_eid += 1

    for station_index in range(N_STATIONS - 1):
        if station_index % 2 == 0:
            node_a = model.arch_negative[station_index]
            node_b = model.arch_positive[station_index + 1]
        else:
            node_a = model.arch_positive[station_index]
            node_b = model.arch_negative[station_index + 1]
        model.add_rod(20601 + station_index, 205, node_a, node_b)

    # Alternating lower plan bracing closes the deck rectangles in shear.
    for station_index in range(N_STATIONS - 1):
        if station_index % 2 == 0:
            node_a = model.deck_lines[0][station_index]
            node_b = model.deck_lines[3][station_index + 1]
        else:
            node_a = model.deck_lines[3][station_index]
            node_b = model.deck_lines[0][station_index + 1]
        model.add_rod(20701 + station_index, 204, node_a, node_b)

    # Each inner stringer is tied longitudinally to its adjacent exterior
    # girder.  These diagonals remove the otherwise free rigid sliding mode of
    # a CROD-only line while retaining a pin-jointed deck-grid idealization.
    inner_plan_eid = 21001
    for station_index in range(N_STATIONS - 1):
        model.add_rod(
            inner_plan_eid,
            204,
            model.deck_lines[0][station_index],
            model.deck_lines[1][station_index + 1],
        )
        inner_plan_eid += 1
        model.add_rod(
            inner_plan_eid,
            204,
            model.deck_lines[3][station_index],
            model.deck_lines[2][station_index + 1],
        )
        inner_plan_eid += 1

    # Network hangers give both inner deck lines a direct out-of-plane load
    # path.  At springing stations they connect to the lower pier capitals.
    network_eid = 20901
    for station_index in range(N_STATIONS):
        if station_index in SUPPORT_INDICES:
            lower_negative, lower_positive = model.support_lower[station_index]
            negative_target = lower_negative
            positive_target = lower_positive
        else:
            negative_target = model.arch_negative[station_index]
            positive_target = model.arch_positive[station_index]
        model.add_rod(
            network_eid,
            207,
            model.deck_lines[1][station_index],
            negative_target,
        )
        network_eid += 1
        model.add_rod(
            network_eid,
            207,
            model.deck_lines[2][station_index],
            positive_target,
        )
        network_eid += 1

    # Ground beams make each idealized foundation base explicit in the model.
    footing_eid = 20801
    for station_index in SUPPORT_INDICES:
        ground_ids = model.support_ground[station_index]
        for node_a, node_b in zip(ground_ids[:-1], ground_ids[1:]):
            model.add_rod(footing_eid, 206, node_a, node_b)
            footing_eid += 1

    _validate_model(model)
    return model


def _validate_model(model: BridgeModel) -> None:
    if len(model.nodes) != 202:
        raise ValueError(f"Expected 202 GRID, found {len(model.nodes)}")
    if len(model.bars) != 144:
        raise ValueError(f"Expected 144 CBAR, found {len(model.bars)}")
    if len(model.rods) != 428:
        raise ValueError(f"Expected 428 CROD, found {len(model.rods)}")

    element_ids = [element[0] for element in model.bars + model.rods]
    if len(element_ids) != len(set(element_ids)):
        raise ValueError("Element IDs must be unique")
    for _eid, _pid, node_a, node_b, *_rest in model.bars + model.rods:
        if node_a == node_b or node_a not in model.nodes or node_b not in model.nodes:
            raise ValueError("Every element must connect two existing distinct GRID entries")


def _append_properties(lines: list[str]) -> None:
    lines.extend(
        [
            "$ Materials: structural steel CBAR and high-strength rod steel.",
            "MAT1,1,210000.,80769.2308,0.30",
            "MAT1,2,195000.,75000.,0.30",
            "$ PBAR: PID,MID,A,I1,I2,J,NSM followed by four recovery points.",
            "$ PID 101 - exterior longitudinal tie/deck girders.",
            "PBAR,101,1,18000.,1.20E9,4.50E8,9.00E7,0.",
            "+,300.,450.,-300.,450.,-300.,-450.,300.,-450.",
            "$ PID 103 - three-span arch ribs.",
            "PBAR,103,1,12000.,7.00E8,3.20E8,5.00E7,0.",
            "+,250.,360.,-250.,360.,-250.,-360.,250.,-360.",
            "$ PID 105 - vertical pier and abutment columns.",
            "PBAR,105,1,24000.,1.60E9,1.00E9,1.50E8,0.",
            "+,350.,500.,-350.,500.,-350.,-500.,350.,-500.",
            "$ PID 106 - inclined foundation legs.",
            "PBAR,106,1,16000.,8.00E8,5.00E8,9.00E7,0.",
            "+,280.,380.,-280.,380.,-280.,-380.,280.,-380.",
            "$ PBARN1 selects the validated nonlinear LC interpolation.",
        ]
    )
    for property_id in (101, 103, 105, 106):
        lines.extend((f"PBARN1,{property_id},,,S", ",C2,BAR,LC"))
    lines.extend(
        [
            "$ PROD: bilateral hangers and pin-jointed spatial truss members.",
            "PROD,201,2,800.,0.,0.,0.",
            "PROD,202,1,2500.,0.,0.,0.",
            "PROD,203,1,1400.,0.,0.,0.",
            "PROD,204,2,900.,0.,0.,0.",
            "PROD,205,2,650.,0.,0.,0.",
            "PROD,206,1,3000.,0.,0.,0.",
            "PROD,207,2,700.,0.,0.,0.",
        ]
    )


def _append_constraints(lines: list[str], model: BridgeModel) -> None:
    lines.append("$ Ground bases: fixed T1/T2 and all rotations; T3 is explicit below.")
    for ground_ids in model.support_ground.values():
        for grid_id in ground_ids:
            lines.append(f"SPC1,10,12456,{grid_id}")

    lines.append("$ Left base has zero T3; the other foundations settle differentially.")
    for grid_id in model.support_ground[SUPPORT_INDICES[0]]:
        lines.append(f"SPC1,11,3,{grid_id}")
    for station_index in SUPPORT_INDICES[1:]:
        settlement = SUPPORT_SETTLEMENT[station_index]
        for grid_id in model.support_ground[station_index]:
            lines.append(f"SPC,20,{grid_id},3,{settlement:.6f}")

    lines.append("$ CROD-only inner deck GRID: inactive rotations constrained explicitly.")
    for grid_id in model.inner_deck_nodes:
        lines.append(f"SPC1,30,456,{grid_id}")
    lines.append("SPCADD,100,10,11,20,30")


def _append_loads(lines: list[str], model: BridgeModel) -> None:
    traffic = (
        (2, 85_000.0),
        (4, 105_000.0),
        (6, 90_000.0),
        (10, 115_000.0),
        (12, 130_000.0),
        (14, 145_000.0),
        (16, 135_000.0),
        (18, 110_000.0),
        (22, 100_000.0),
        (24, 125_000.0),
        (26, 115_000.0),
        (28, 95_000.0),
    )
    lines.append("$ Load set 501 - asymmetric vertical traffic on the two inner rails.")
    for station_index, magnitude in traffic:
        negative_grid = model.deck_lines[1][station_index]
        positive_grid = model.deck_lines[2][station_index]
        lines.append(f"FORCE,501,{negative_grid},0,{magnitude:.1f},0.,0.,-1.")
        lines.append(
            f"FORCE,501,{positive_grid},0,{0.92 * magnitude:.1f},0.,0.,-1."
        )

    wind_stations = (3, 5, 11, 14, 17, 23, 26, 28)
    lines.append("$ Load set 502 - lateral wind on selected arch nodes.")
    for load_number, station_index in enumerate(wind_stations):
        magnitude = 24_000.0 + 2_500.0 * load_number
        lines.append(
            f"FORCE,502,{model.arch_negative[station_index]},0,{magnitude:.1f},0.,1.,0."
        )
        lines.append(
            f"FORCE,502,{model.arch_positive[station_index]},0,{0.85 * magnitude:.1f},0.,1.,0."
        )

    lines.extend(
        [
            "$ Load set 503 - longitudinal braking and nodal couples.",
            f"FORCE,503,{model.deck_lines[0][13]},0,80000.,1.,0.,0.",
            f"FORCE,503,{model.deck_lines[3][13]},0,80000.,1.,0.,0.",
            f"FORCE,503,{model.deck_lines[0][16]},0,95000.,1.,0.,0.",
            f"FORCE,503,{model.deck_lines[3][16]},0,95000.,1.,0.,0.",
            f"MOMENT,503,{model.deck_lines[0][14]},0,1.10E8,1.,0.,0.",
            f"MOMENT,503,{model.deck_lines[3][14]},0,-8.00E7,1.,0.,0.",
            f"MOMENT,503,{model.deck_lines[0][22]},0,1.40E8,0.,-1.,0.",
            "$ Wind enters at 75 percent; traffic and braking enter at full value.",
            "LOAD,500,1.0,1.0,501,0.75,502,1.0,503",
        ]
    )


def render_deck(model: BridgeModel) -> str:
    lines = [
        "$ UPDATED FINAL NONLINEAR VALIDATION - THREE-SPAN THROUGH-ARCH BRIDGE",
        "$ Units: mm, N, MPa. Geometric nonlinearity only; elastic MAT1.",
        "$ 202 GRID; 144 CBAR; 428 CROD; explicit abutment/pier foundations.",
        "$ AUTOSPC is disabled: every inactive or supported DOF is constrained here.",
        "TIME 9999",
        "SOL 400",
        "CEND",
        "TITLE = TFM UPDATED FINAL NONLINEAR THREE SPAN ARCH BRIDGE",
        "SUBTITLE = TRAFFIC WIND BRAKING MOMENTS AND DIFFERENTIAL FOUNDATION SETTLEMENT",
        "ECHO = NONE",
        "NLOPRM OUTCTRL = (SOLUTION,INTERM) NLDBG = ADVDBG",
        "SUBCASE 1",
        "  LABEL = UPDATED THREE SPAN BRIDGE COMBINED ACTIONS",
        "  ANALYSIS = NLSTAT",
        "  SPC = 100",
        "  LOAD = 500",
        "  NLSTEP = 500",
        "  DISPLACEMENT(PRINT,PUNCH) = ALL",
        "  SPCFORCES(PRINT,PUNCH) = ALL",
        "  OLOAD(PRINT,PUNCH) = ALL",
        "  GPFORCE(PRINT,PUNCH) = ALL",
        "  STRESS(SORT1,REAL,PRINT,PUNCH) = ALL",
        "  STRAIN(SORT1,REAL,PRINT,PUNCH) = ALL",
        "  ELFORCE(SORT1,REAL,PRINT,PUNCH) = ALL",
        "  NLSTRESS(SORT1) = ALL",
        "  BOUTPUT(SORT1,REAL) = ALL",
        "  ESE(PRINT,PUNCH) = ALL",
        "BEGIN BULK",
        "PARAM,LGDISP,1",
        "PARAM,FOLLOWK,NO",
        "PARAM,AUTOSPC,NO",
        "PARAM,POST,-1",
    ]
    _append_properties(lines)

    lines.append("$ GRID: ID,CP,X,Y,Z,CD,PS,SEID")
    for node_id, (x, y, z) in sorted(model.nodes.items()):
        lines.append(f"GRID,{node_id},,{x:.6f},{y:.6f},{z:.6f}")

    lines.append("$ CBAR groups use separate EID ranges for Patran grouping.")
    for element_id, property_id, node_a, node_b, orientation in model.bars:
        x1, x2, x3 = orientation
        lines.append(
            f"CBAR,{element_id},{property_id},{node_a},{node_b},"
            f"{x1:.1f},{x2:.1f},{x3:.1f}"
        )

    lines.append("$ CROD: EID,PID,G1,G2")
    for element_id, property_id, node_a, node_b in model.rods:
        lines.append(f"CROD,{element_id},{property_id},{node_a},{node_b}")

    _append_constraints(lines, model)
    _append_loads(lines, model)
    lines.extend(
        [
            "$ Fixed stepping makes the final MSC correlation reproducible.",
            "NLSTEP,500,1.0,SEVERELY",
            ",GENERAL,80,1,24",
            ",FIXED,20",
            ",MECH,PV,3.E-7,3.E-7,3.E-7,PFNT",
            "ENDDATA",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    model = build_model()
    args.output.write_text(render_deck(model), encoding="ascii", newline="\n")
    print(
        f"Wrote {args.output} with {len(model.nodes)} GRID, "
        f"{len(model.bars)} CBAR and {len(model.rods)} CROD."
    )


if __name__ == "__main__":
    main()
