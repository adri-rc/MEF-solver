"""Generate the final mixed CROD/CBAR SOL 400 bridge validation deck.

The bridge is a spatial tied-arch idealization.  CBAR elements form the deck,
cross-girders, arch ribs and upper struts; CROD elements form hangers and plan
bracing.  Units are mm, N and MPa.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_OUTPUT = Path(__file__).with_name("tfm_final_nonlinear_bridge.dat")
SPAN = 12_000.0
WIDTH = 2_400.0
ARCH_RISE = 2_400.0
N_STATIONS = 9


@dataclass
class BridgeModel:
    """Deterministic geometry and connectivity for the validation bridge."""

    nodes: dict[int, tuple[float, float, float]] = field(default_factory=dict)
    bars: list[tuple[int, int, int, int, tuple[float, float, float]]] = field(
        default_factory=list
    )
    rods: list[tuple[int, int, int, int]] = field(default_factory=list)
    deck_a: list[int] = field(default_factory=list)
    deck_b: list[int] = field(default_factory=list)
    arch_a: list[int] = field(default_factory=list)
    arch_b: list[int] = field(default_factory=list)

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
    xi = station_index / (N_STATIONS - 1)
    return 4.0 * ARCH_RISE * xi * (1.0 - xi)


def build_model() -> BridgeModel:
    model = BridgeModel()
    spacing = SPAN / (N_STATIONS - 1)
    half_width = WIDTH / 2.0

    for index in range(N_STATIONS):
        x = index * spacing
        deck_a = 101 + index
        deck_b = 201 + index
        model.nodes[deck_a] = (x, -half_width, 0.0)
        model.nodes[deck_b] = (x, half_width, 0.0)
        model.deck_a.append(deck_a)
        model.deck_b.append(deck_b)

        if index == 0 or index == N_STATIONS - 1:
            model.arch_a.append(deck_a)
            model.arch_b.append(deck_b)
        else:
            arch_a = 300 + index
            arch_b = 400 + index
            z = _arch_height(index)
            model.nodes[arch_a] = (x, -half_width, z)
            model.nodes[arch_b] = (x, half_width, z)
            model.arch_a.append(arch_a)
            model.arch_b.append(arch_b)

    # Longitudinal deck girders.
    for index in range(N_STATIONS - 1):
        model.add_bar(1001 + index, 101, model.deck_a[index], model.deck_a[index + 1], (0.0, 0.0, 1.0))
        model.add_bar(1101 + index, 101, model.deck_b[index], model.deck_b[index + 1], (0.0, 0.0, 1.0))

    # Transverse deck cross-girders.
    for index in range(N_STATIONS):
        model.add_bar(1201 + index, 102, model.deck_a[index], model.deck_b[index], (0.0, 0.0, 1.0))

    # Two parabolic arch ribs and their transverse upper struts.
    for index in range(N_STATIONS - 1):
        model.add_bar(2001 + index, 103, model.arch_a[index], model.arch_a[index + 1], (0.0, 1.0, 0.0))
        model.add_bar(2101 + index, 103, model.arch_b[index], model.arch_b[index + 1], (0.0, 1.0, 0.0))
    for index in range(1, N_STATIONS - 1):
        model.add_bar(2200 + index, 104, model.arch_a[index], model.arch_b[index], (0.0, 0.0, 1.0))

    # Vertical hangers. CROD is deliberately bilateral: this is not a
    # tension-only cable or a slackening/contact model.
    rod_id = 3001
    for index in range(1, N_STATIONS - 1):
        model.add_rod(rod_id, 201, model.deck_a[index], model.arch_a[index])
        rod_id += 1
        model.add_rod(rod_id, 201, model.deck_b[index], model.arch_b[index])
        rod_id += 1

    # Alternating plan braces under the deck.
    for index in range(N_STATIONS - 1):
        if index % 2 == 0:
            node_a, node_b = model.deck_a[index], model.deck_b[index + 1]
        else:
            node_a, node_b = model.deck_b[index], model.deck_a[index + 1]
        model.add_rod(rod_id, 202, node_a, node_b)
        rod_id += 1

    # Alternating wind braces between the arch ribs.
    for index in range(N_STATIONS - 1):
        if index % 2 == 0:
            node_a, node_b = model.arch_a[index], model.arch_b[index + 1]
        else:
            node_a, node_b = model.arch_b[index], model.arch_a[index + 1]
        model.add_rod(rod_id, 203, node_a, node_b)
        rod_id += 1

    _validate_model(model)
    return model


def _validate_model(model: BridgeModel) -> None:
    if len(model.nodes) != 32 or len(model.bars) != 48 or len(model.rods) != 30:
        raise ValueError("Unexpected bridge model size")
    element_ids = [element[0] for element in model.bars + model.rods]
    if len(element_ids) != len(set(element_ids)):
        raise ValueError("Element IDs must be unique")
    for _eid, _pid, node_a, node_b, *_rest in model.bars + model.rods:
        if node_a == node_b or node_a not in model.nodes or node_b not in model.nodes:
            raise ValueError("Every element must connect two existing distinct GRID entries")


def render_deck(model: BridgeModel) -> str:
    lines = [
        "$ FINAL NONLINEAR VALIDATION - SPATIAL TIED-ARCH BRIDGE",
        "$ Units: mm, N, MPa. Geometric nonlinearity only; elastic MAT1.",
        "$ CBAR: deck, cross-girders, arch ribs and upper struts.",
        "$ CROD: bilateral axial hangers and alternating spatial braces.",
        "TIME 9999",
        "SOL 400",
        "CEND",
        "TITLE = TFM FINAL NONLINEAR MIXED CBAR CROD BRIDGE",
        "SUBTITLE = TRAFFIC WIND BRAKING MOMENTS AND DIFFERENTIAL SETTLEMENT",
        "ECHO = NONE",
        "NLOPRM OUTCTRL = (SOLUTION,INTERM) NLDBG = ADVDBG",
        "SUBCASE 1",
        "  LABEL = COMBINED SERVICE ACTIONS AND PRESCRIBED BEARING SETTLEMENT",
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
        "PARAM,POST,-1",
        "$ Materials: structural steel CBAR and high-strength steel CROD.",
        "MAT1,1,210000.,80769.2308,0.30",
        "MAT1,2,195000.,75000.,0.30",
        "$ PBAR: PID,MID,A,I1,I2,J,NSM",
        "$ PID 101 - twin longitudinal deck box girders.",
        "PBAR,101,1,12000.,7.5E8,2.4E8,5.0E7,0.",
        "+,220.,320.,-220.,320.,-220.,-320.,220.,-320.",
        "$ PID 102 - transverse deck cross-girders.",
        "PBAR,102,1,6000.,1.8E8,8.0E7,2.0E7,0.",
        "+,160.,220.,-160.,220.,-160.,-220.,160.,-220.",
        "$ PID 103 - parabolic arch ribs.",
        "PBAR,103,1,9000.,4.5E8,2.2E8,4.0E7,0.",
        "+,190.,270.,-190.,270.,-190.,-270.,190.,-270.",
        "$ PID 104 - transverse upper struts.",
        "PBAR,104,1,3500.,6.0E7,3.0E7,8.0E6,0.",
        "+,110.,150.,-110.,150.,-110.,-150.,110.,-150.",
        "$ PBARN1 selects the validated nonlinear LC interpolation.",
    ]
    for property_id in (101, 102, 103, 104):
        lines.extend(
            [
                f"PBARN1,{property_id},,,S",
                ",C2,BAR,LC",
            ]
        )
    lines.extend(
        [
            "$ PROD: hangers, deck plan braces and arch wind braces.",
            "PROD,201,2,900.,0.,0.,0.",
            "PROD,202,2,650.,0.,0.,0.",
            "PROD,203,2,500.,0.,0.,0.",
            "$ GRID: ID,CP,X,Y,Z,CD,PS,SEID",
        ]
    )
    for node_id, (x, y, z) in sorted(model.nodes.items()):
        lines.append(f"GRID,{node_id},,{x:.6f},{y:.6f},{z:.6f}")

    lines.append("$ CBAR groups are separated by EID range for Patran grouping.")
    for element_id, property_id, node_a, node_b, orientation in model.bars:
        x1, x2, x3 = orientation
        lines.append(
            f"CBAR,{element_id},{property_id},{node_a},{node_b},"
            f"{x1:.1f},{x2:.1f},{x3:.1f}"
        )

    lines.append("$ CROD: EID,PID,G1,G2")
    for element_id, property_id, node_a, node_b in model.rods:
        lines.append(f"CROD,{element_id},{property_id},{node_a},{node_b}")

    lines.extend(
        [
            "$ Left abutment: two full CBAR end-section clamps.",
            "SPC1,10,123456,101,201",
            "$ Right bearings: transverse translation fixed, axial motion free.",
            "SPC1,20,2,109,209",
            "$ Differential vertical settlement, ramped with the load factor.",
            "SPC,30,109,3,-12.0",
            "SPC,30,209,3,-6.0",
            "SPCADD,100,10,20,30",
            "$ Load set 501 - asymmetric vertical traffic actions.",
            "FORCE,501,104,0,2.20E5,0.,0.,-1.",
            "FORCE,501,204,0,2.60E5,0.,0.,-1.",
            "FORCE,501,105,0,3.20E5,0.,0.,-1.",
            "FORCE,501,205,0,2.80E5,0.,0.,-1.",
            "FORCE,501,106,0,2.50E5,0.,0.,-1.",
            "FORCE,501,206,0,2.10E5,0.,0.,-1.",
            "FORCE,501,107,0,1.80E5,0.,0.,-1.",
            "$ Load set 502 - lateral wind actions on both arch ribs.",
            "FORCE,502,303,0,3.50E4,0.,1.,0.",
            "FORCE,502,304,0,5.00E4,0.,1.,0.",
            "FORCE,502,305,0,4.00E4,0.,1.,0.",
            "FORCE,502,403,0,3.00E4,0.,1.,0.",
            "FORCE,502,404,0,4.50E4,0.,1.,0.",
            "FORCE,502,405,0,3.50E4,0.,1.,0.",
            "$ Load set 503 - braking forces and nodal couples.",
            "FORCE,503,105,0,9.00E4,1.,0.,0.",
            "FORCE,503,205,0,9.00E4,1.,0.,0.",
            "MOMENT,503,105,0,1.20E8,1.,0.,0.",
            "MOMENT,503,205,0,-8.00E7,1.,0.,0.",
            "MOMENT,503,106,0,1.50E8,0.,-1.,0.",
            "$ Wind is combined at 80 percent; all other sets at full value.",
            "LOAD,500,1.0,1.0,501,0.8,502,1.0,503",
            "$ Fixed stepping makes the final MSC correlation reproducible.",
            "NLSTEP,500,1.0,SEVERELY",
            ",GENERAL,80,1,24",
            ",FIXED,50",
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
