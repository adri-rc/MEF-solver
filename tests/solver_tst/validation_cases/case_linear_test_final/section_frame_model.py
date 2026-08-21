"""Generate the final SOL 101 airframe-frame validation deck.

The model is deliberately generated instead of hand-written so that the two
O-grid meshes, element IDs and load-introduction nodes remain reproducible.
Units are mm, N and MPa.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


DEFAULT_OUTPUT = Path(__file__).with_name("tfm_final_airframe_frame.dat")
N_PERIMETER = 32
N_RADIAL_LAYERS = 5
BLOCK_HALF_SIZE = 300.0
HOLE_RADIUS = 150.0


def _chunks(values: list[int], size: int = 8):
    for start in range(0, len(values), size):
        yield values[start : start + size]


class DeckBuilder:
    """Small deterministic mesh/deck builder for this validation model."""

    def __init__(self) -> None:
        self.nodes: dict[int, tuple[float, float, float]] = {}
        self._node_by_coordinates: dict[tuple[float, float, float], int] = {}
        self.quads: list[tuple[int, int, int, int, int]] = []
        self.trias: list[tuple[int, int, int, int]] = []
        self.bars: list[tuple[int, int, int, int]] = []
        self.rods: list[tuple[int, int, int]] = []
        self._next_quad = 1001
        self._next_tria = 5001
        self._next_bar = 7001
        self._next_rod = 8001

    def add_node(self, x: float, y: float, z: float = 0.0) -> int:
        key = (round(x, 8), round(y, 8), round(z, 8))
        existing = self._node_by_coordinates.get(key)
        if existing is not None:
            return existing
        nid = len(self.nodes) + 1
        self.nodes[nid] = key
        self._node_by_coordinates[key] = nid
        return nid

    def add_quad(self, nodes: tuple[int, int, int, int], pid: int = 10) -> None:
        self.quads.append((self._next_quad, pid, *nodes))
        self._next_quad += 1

    def add_tria(self, nodes: tuple[int, int, int], pid: int = 11) -> None:
        self.trias.append((self._next_tria, pid, *nodes))
        self._next_tria += 1

    def add_bar(self, node_a: int, node_b: int, pid: int) -> None:
        if node_a == node_b:
            raise ValueError("A CBAR cannot have zero length")
        self.bars.append((self._next_bar, pid, node_a, node_b))
        self._next_bar += 1

    def add_rod(self, node_a: int, node_b: int, pid: int = 30) -> None:
        self.rods.append((self._next_rod, pid, node_a, node_b))
        self._next_rod += 1

    def node_at(self, x: float, y: float, tolerance: float = 1.0e-6) -> int:
        matches = [
            nid
            for nid, (nx, ny, _nz) in self.nodes.items()
            if abs(nx - x) <= tolerance and abs(ny - y) <= tolerance
        ]
        if len(matches) != 1:
            raise ValueError(f"Expected one GRID at ({x}, {y}); found {matches}")
        return matches[0]

    def boundary_nodes(self, axis: int, value: float) -> list[int]:
        matches = [
            nid
            for nid, coordinates in self.nodes.items()
            if math.isclose(coordinates[axis], value, abs_tol=1.0e-6)
        ]
        other_axis = 1 - axis
        return sorted(matches, key=lambda nid: self.nodes[nid][other_axis])


def _square_boundary_point(
    center_x: float, center_y: float, theta: float
) -> tuple[float, float]:
    direction_x = math.cos(theta)
    direction_y = math.sin(theta)
    scale = BLOCK_HALF_SIZE / max(abs(direction_x), abs(direction_y))
    return center_x + scale * direction_x, center_y + scale * direction_y


def _build_web_block(
    builder: DeckBuilder, center_x: float, center_y: float
) -> list[list[int]]:
    rings: list[list[int]] = []
    for radial_index in range(N_RADIAL_LAYERS + 1):
        blend = radial_index / N_RADIAL_LAYERS
        ring: list[int] = []
        for perimeter_index in range(N_PERIMETER):
            theta = 2.0 * math.pi * perimeter_index / N_PERIMETER
            inner_x = center_x + HOLE_RADIUS * math.cos(theta)
            inner_y = center_y + HOLE_RADIUS * math.sin(theta)
            outer_x, outer_y = _square_boundary_point(center_x, center_y, theta)
            x = inner_x + blend * (outer_x - inner_x)
            y = inner_y + blend * (outer_y - inner_y)
            ring.append(builder.add_node(x, y))
        rings.append(ring)

    for radial_index in range(N_RADIAL_LAYERS):
        for perimeter_index in range(N_PERIMETER):
            next_index = (perimeter_index + 1) % N_PERIMETER
            n0 = rings[radial_index][perimeter_index]
            n1 = rings[radial_index + 1][perimeter_index]
            n2 = rings[radial_index + 1][next_index]
            n3 = rings[radial_index][next_index]
            if radial_index == 0:
                # CTRIA3s provide the transition immediately around each
                # circular mouse hole; all remaining web cells are CQUAD4.
                builder.add_tria((n0, n1, n2))
                builder.add_tria((n0, n2, n3))
            else:
                builder.add_quad((n0, n1, n2, n3))
    return rings


def _add_bar_polyline(builder: DeckBuilder, nodes: list[int], pid: int) -> None:
    for node_a, node_b in zip(nodes, nodes[1:]):
        builder.add_bar(node_a, node_b, pid)


def build_model() -> tuple[DeckBuilder, dict[str, object]]:
    builder = DeckBuilder()
    left_rings = _build_web_block(builder, center_x=300.0, center_y=300.0)
    right_rings = _build_web_block(builder, center_x=900.0, center_y=300.0)

    # T-cap idealization around the coupon and on its central frame member.
    bottom_nodes = builder.boundary_nodes(axis=1, value=0.0)
    top_nodes = builder.boundary_nodes(axis=1, value=600.0)
    left_nodes = builder.boundary_nodes(axis=0, value=0.0)
    right_nodes = builder.boundary_nodes(axis=0, value=1200.0)
    center_nodes = builder.boundary_nodes(axis=0, value=600.0)
    for polyline in (bottom_nodes, top_nodes, left_nodes, right_nodes, center_nodes):
        _add_bar_polyline(builder, polyline, pid=20)

    # Four radial secondary ribs bound the upper and lower ligaments of the
    # holes. They use existing shell nodes, without MPC/RBE connections.
    for rings in (left_rings, right_rings):
        for perimeter_index in (8, 24):
            rib_nodes = [ring[perimeter_index] for ring in rings]
            _add_bar_polyline(builder, rib_nodes, pid=21)

    left_clevis = builder.add_node(300.0, 800.0)
    right_clevis = builder.add_node(900.0, 800.0)
    for target_x in (200.0, 400.0):
        attachment = min(
            top_nodes, key=lambda nid: abs(builder.nodes[nid][0] - target_x)
        )
        builder.add_rod(left_clevis, attachment)
    for target_x in (800.0, 1000.0):
        attachment = min(
            top_nodes, key=lambda nid: abs(builder.nodes[nid][0] - target_x)
        )
        builder.add_rod(right_clevis, attachment)

    top_center = builder.node_at(600.0, 600.0)
    rivet_indices = (0, 3, 6, 10, 13, 16)
    rivet_nodes = [bottom_nodes[index] for index in rivet_indices]
    metadata: dict[str, object] = {
        "left_clevis": left_clevis,
        "right_clevis": right_clevis,
        "top_center": top_center,
        "rivet_nodes": rivet_nodes,
        "mesh_nodes": [nid for nid in builder.nodes if nid not in (left_clevis, right_clevis)],
    }
    return builder, metadata


def _format_spc1(sid: int, components: str, node_ids: list[int]) -> list[str]:
    return [
        ",".join(["SPC1", str(sid), components, *(str(nid) for nid in chunk)])
        for chunk in _chunks(node_ids)
    ]


def render_deck(builder: DeckBuilder, metadata: dict[str, object]) -> str:
    left_clevis = int(metadata["left_clevis"])
    right_clevis = int(metadata["right_clevis"])
    top_center = int(metadata["top_center"])
    rivet_nodes = list(metadata["rivet_nodes"])
    mesh_nodes = list(metadata["mesh_nodes"])

    lines = [
        "$ FINAL LINEAR VALIDATION - AIRFRAME FRAME/WEB COUPON",
        "$ Units: mm, N, MPa. Small-displacement SOL 101.",
        "$ Two circular mouse holes use CTRIA3 transition rings; the remaining",
        "$ web uses CQUAD4. CBARs idealize T-caps/ribs and CRODs introduce load.",
        "TIME 9999",
        "SOL 101",
        "CEND",
        "TITLE = TFM FINAL LINEAR AIRFRAME FRAME VALIDATION",
        "SUBTITLE = CQUAD4 CTRIA3 CBAR CROD FORCE MOMENT ENFORCED DISPLACEMENT",
        "ECHO = NONE",
        "SUBCASE 1",
        "  LABEL = COMBINED PULL MOMENT AND PRESCRIBED CLEVIS MOTION",
        "  SPC = 130",
        "  LOAD = 200",
        "  DISPLACEMENT(PRINT,PLOT) = ALL",
        "  SPCFORCES(PRINT,PLOT) = ALL",
        "  OLOAD(PRINT,PLOT) = ALL",
        "  STRESS(PRINT,PLOT) = ALL",
        "  FORCE(PRINT,PLOT) = ALL",
        "BEGIN BULK",
        "PARAM,K6ROT,100.0",
        "PARAM,AUTOSPC,NO",
        "PARAM,POST,-1",
        "$ Elements and Element Properties for region : WEB_CQUAD4",
        "$ PID 10: 2 mm 7075-T6-like aluminum web, full shell behavior.",
        "PSHELL,10,1,2.0,1,,1,0.8333333333",
        "$ Elements and Element Properties for region : MOUSEHOLE_CTRIA3",
        "$ PID 11: same web; MID3=99 makes CTRIA3 transverse shear negligible.",
        "$ The test fixture constrains T3/R1/R2 globally, so validation is in-plane.",
        "PSHELL,11,1,2.0,1,,99,0.8333333333",
        "$ Elements and Element Properties for region : T_CAP_CBAR",
        "$ Equivalent T: flange 40x4 mm plus stem 25x3 mm (A=235 mm2).",
        "PBAR,20,1,235.0,21389.0,14854.0,1078.0",
        "+,20.0,-6.63,-20.0,-6.63,1.5,22.37,-1.5,22.37",
        "$ Elements and Element Properties for region : SECONDARY_RIB_CBAR",
        "PBAR,21,1,120.0,6000.0,3000.0,350.0",
        "+,12.0,-8.0,-12.0,-8.0,12.0,8.0,-12.0,8.0",
        "$ Elements and Element Properties for region : CLEVIS_CROD",
        "$ A=50 mm2, J=20 mm4; clevis inactive DOFs are fixture-guided.",
        "PROD,30,1,50.0,20.0,2.5",
        "$ Materials: structural aluminum and near-zero CTRIA3 MID3 surrogate.",
        "MAT1,1,71700.0,,0.33",
        "MAT1,99,1.0E-20,,0.33",
        "$ GRID: ID,CP,X,Y,Z,CD,PS,SEID",
    ]
    for nid, (x, y, z) in builder.nodes.items():
        lines.append(f"GRID,{nid},,{x:.8f},{y:.8f},{z:.8f}")

    lines.append("$ CQUAD4: EID,PID,G1,G2,G3,G4,THETA")
    for eid, pid, n1, n2, n3, n4 in builder.quads:
        lines.append(f"CQUAD4,{eid},{pid},{n1},{n2},{n3},{n4},0.0")
    lines.append("$ CTRIA3: EID,PID,G1,G2,G3,THETA")
    for eid, pid, n1, n2, n3 in builder.trias:
        lines.append(f"CTRIA3,{eid},{pid},{n1},{n2},{n3},0.0")
    lines.append("$ CBAR: EID,PID,GA,GB,X1,X2,X3 (orientation +Z)")
    for eid, pid, node_a, node_b in builder.bars:
        lines.append(f"CBAR,{eid},{pid},{node_a},{node_b},0.0,0.0,1.0")
    lines.append("$ CROD: EID,PID,G1,G2")
    for eid, pid, node_a, node_b in builder.rods:
        lines.append(f"CROD,{eid},{pid},{node_a},{node_b}")

    lines.extend(
        [
            "$ SPC set 100: planar test-fixture guidance.",
            *_format_spc1(100, "345", mesh_nodes),
            *_format_spc1(100, "3456", [left_clevis, right_clevis]),
            "$ SPC set 110: six idealized fully fixed rivet stations on lower cap.",
            *_format_spc1(110, "126", rivet_nodes),
            "$ SPC set 120: upward actuator stroke of 0.8 mm at right clevis.",
            f"SPC,120,{right_clevis},2,0.8",
            "SPCADD,130,100,110,120",
            "$ Load set 201: upward pull at left clevis plus lateral frame force.",
            f"FORCE,201,{left_clevis},0,12000.0,0.0,1.0,0.0",
            f"FORCE,201,{top_center},0,2500.0,1.0,0.0,0.0",
            "$ Load set 202: in-plane frame moment about global Z.",
            f"MOMENT,202,{top_center},0,250000.0,0.0,0.0,1.0",
            "LOAD,200,1.0,1.0,201,1.0,202",
            "ENDDATA",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    builder, metadata = build_model()
    args.output.write_text(render_deck(builder, metadata), encoding="utf-8")
    print(
        f"Wrote {args.output} with {len(builder.nodes)} GRID, "
        f"{len(builder.quads)} CQUAD4, {len(builder.trias)} CTRIA3, "
        f"{len(builder.bars)} CBAR and {len(builder.rods)} CROD."
    )


if __name__ == "__main__":
    main()
