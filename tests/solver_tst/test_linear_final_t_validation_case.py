"""Structural and MSC checks for the final three-dimensional T-frame model.

The external values below are transcribed from
``tfm_final_airframe_frame_T.f06``, generated with MSC Nastran Student
Edition 2026.1 on 2026-08-24.  This test deliberately does not read F06/OP2
at runtime, so the reference remains explicit and auditable in source.

Shell acceptance covers membrane resultants and the positive-fibre stress
row.  CTRIA3 centroidal bending moments and transverse shears are excluded:
their recovery does not correlate in this model and is documented in the
case README instead of being hidden behind a loose tolerance.
"""

from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.data_model import (
    BarElement,
    CQUAD4Element,
    CTRIA3Element,
    RodElement,
    SinglePointLoad,
    SinglePointMoment,
)
from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


CASE_FILE = (
    Path(__file__).parent
    / "validation_cases"
    / "case_linear_test_final_T"
    / "tfm_final_airframe_frame_T.dat"
)


MSC_DISPLACEMENTS = {
    9: (-4.409887e-1, 3.083221, -3.609286e-1, -7.644225e-4, -1.990196e-6, 7.442750e-4),
    27: (-1.468139e-1, 3.048432, -1.806144e-1, -1.058733e-3, -4.252649e-4, 1.064546e-3),
    51: (-1.752734e-1, 2.893064, -1.975081e-1, -7.641289e-4, 2.320034e-4, 8.427706e-4),
    165: (-5.612553e-1, 3.304615, -4.757499e-1, -7.943589e-4, 1.957877e-8, 3.458345e-3),
    376: (3.019772, 0.0, 0.0, 0.0, 0.0, 0.0),
    377: (2.911066, 0.0, 0.0, 0.0, 0.0, 0.0),
    405: (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    407: (-3.298431e-4, 1.216538, 1.807086e-5, -4.005495e-3, 9.942546e-6, 1.524455e-2),
    413: (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    415: (6.035173e-4, 1.366869, 2.850431e-4, -6.982817e-3, 7.441017e-6, 1.678660e-2),
    418: (-1.440211e-3, 2.884130, -7.410030e-4, -6.031379e-3, 6.230641e-3, -3.451217e-3),
    422: (0.0, 2.0e-1, 0.0, 0.0, 0.0, 0.0),
    430: (0.0, 2.0e-1, 0.0, 0.0, 0.0, 0.0),
    434: (1.757970e-3, 3.985038, -7.849077e-4, -7.022823e-3, -7.633911e-3, 6.283141e-3),
    440: (-1.525886e-1, 3.789752, -3.213940e-1, 3.888141e-4, -5.191459e-3, 4.346272e-3),
}

MSC_SPC_REACTIONS = {
    376: (0.0, 0.0, 0.0, -1.371455, 0.0, 0.0),
    377: (0.0, -1.665335e-12, 0.0, 1.371431, -3.806487e-16, 0.0),
    405: (8.905963e1, -2.097370e3, -1.753117e2, 1.495664e5, 3.684117e-2, -8.422031e2),
    413: (-8.905981e1, -2.509189e3, -2.242664e2, 1.788843e5, 6.267507e-2, -6.119595e2),
    422: (8.900810e1, -2.089024e3, 1.753116e2, -1.483155e5, -3.695332e-2, -8.442631e2),
    430: (-8.900793e1, -2.500844e3, 2.242665e2, -1.776335e5, -6.256934e-2, -6.099108e2),
}

# axial force, torque, axial stress and torsional stress
MSC_CROD_RESULTS = {
    8001: (6.000000e3, -1.371455, 1.200000e2, -1.714319e-1),
    8002: (6.000000e3, 1.371431, 1.200000e2, 1.714289e-1),
}

# F06 order: MA1, MA2, MB1, MB2, shear 1, shear 2, axial and torque.
MSC_CBAR_FORCES = {
    7001: (4.383186, 2.984130e2, -2.466016, -6.833267e2, 6.880411e-2, 9.862131, -5.995633, 8.074378),
    7025: (2.594579, 8.152888e4, 3.115252, -1.769260e4, -5.230459e-3, 9.967358e2, 8.462374e2, 2.709984),
    7026: (3.227453, -1.767426e4, 4.194040, -2.277374e3, -1.268661e-2, -2.020867e2, 8.162088e2, 3.046995),
    7030: (5.878735, -8.907751e3, 5.024330, 2.869711e4, 1.322805e-2, -5.822054e2, 3.456869e3, -2.116683),
    7056: (5.221014, 1.856010e4, 6.102160, -8.445209e4, -8.851618e-3, 1.034816e3, 5.163453, 4.528889e-3),
    7066: (1.598638e2, 3.374461e1, 2.326017e2, -2.799158e2, -2.424597, 1.045535e1, 1.526256e3, -1.779132e-1),
    7072: (-5.696044, 1.263707, 1.129040e1, -1.744475e2, -5.662148e-1, 5.857040, 3.472674e2, -3.478488e-1),
    7076: (1.598630e2, -4.251286e1, 2.325989e2, 3.034604e2, -2.424531, -1.153244e1, 1.825175e3, -4.262938e-2),
    7101: (3.301574e2, 1.175620e4, -1.048279e3, -8.751497e3, 9.731415, 1.447792e2, 4.183925e2, 3.363920e2),
    7102: (-1.048279e3, -8.751497e3, -2.426715e3, -2.925919e4, 9.731415, 1.447792e2, 4.183925e2, 3.363920e2),
    7105: (8.007095e1, 1.945693e4, -1.289018e2, -1.448147e4, 1.475295, 2.395966e2, 7.677404e2, 2.669413e1),
    7106: (-1.289018e2, -1.448147e4, -3.378746e2, -4.841988e4, 1.475295, 2.395966e2, 7.677404e2, 2.669413e1),
    7108: (-1.289542e2, 1.442385e4, -3.378672e2, 4.827597e4, 1.474874, -2.389875e2, 7.661993e2, -2.669237e1),
    7109: (-4.830369e2, 1.397104e4, 1.245382e3, -1.039888e4, -1.220220e1, 1.720455e2, 4.983369e2, -4.117372e2),
    7112: (1.267301e3, 1.035689e4, 2.982883e3, 3.466293e4, -1.211158e1, -1.715946e2, 4.970979e2, 4.126094e2),
}

# Four bending fibres and axial/max/min at A, then fibres and max/min at B.
MSC_CBAR_STRESSES = {
    7001: (1.290964e-1, 1.372935e-1, -4.497149e-1, -4.491001e-1, -2.551333e-2, 1.117802e-1, -4.752282e-1, -3.026932e-1, -3.073049e-1, 1.029257, 1.028911, 1.003744, -3.328183e-1),
    7025: (3.638754e1, 3.639239e1, -1.227820e2, -1.227816e2, 3.601010, 3.999340e1, -1.191810e2, -7.899904, -7.894078, 2.664468e1, 2.664512e1, 3.024613e1, -4.298894),
    7026: (-7.891825, -7.885789, 2.661706e1, 2.661751e1, 3.473229, 3.009074e1, -4.418596, -1.020415, -1.012572, 3.429412, 3.430001, 6.903230, 2.452814),
    7030: (-3.981422, -3.970428, 1.341459e1, 1.341541e1, 1.471008e1, 2.812549e1, 1.072866e1, 1.280410e1, 1.281349e1, -4.321796e1, -4.321725e1, 2.752357e1, -2.850787e1),
    7056: (8.279314, 8.289078, -2.795172e1, -2.795098e1, 2.197214e-2, 8.311049, -2.792974e1, -3.770042e1, -3.768901e1, 1.271837e2, 1.271846e2, 1.272065e2, -3.767845e1),
    7066: (-2.297419e-1, 4.097132e-1, -4.097132e-1, 2.297419e-1, 1.271880e1, 1.312852e1, 1.230909e1, -1.211645, -2.812388e-1, 2.812388e-1, 1.211645, 1.393045e1, 1.150716e1),
    7072: (1.476197e-2, -8.022203e-3, 8.022203e-3, -1.476197e-2, 2.893895, 2.908657, 2.879133, -4.877741e-1, -4.426125e-1, 4.426125e-1, 4.877741e-1, 3.381669, 2.406121),
    7076: (-4.330936e-1, 2.063584e-1, -2.063584e-1, 4.330936e-1, 1.520979e1, 1.564289e1, 1.477670e1, 3.440300e-1, 1.274426, -1.274426, -3.440300e-1, 1.648422e1, 1.393537e1),
    7101: (3.068954e1, 3.201017e1, -3.201017e1, -3.068954e1, 3.486604, 3.549677e1, -2.852357e1, -2.124077e1, -2.543388e1, 2.543388e1, 2.124077e1, 2.892049e1, -2.194728e1),
    7102: (-1.062038e1, -1.271694e1, 1.271694e1, 1.062038e1, 1.743302, 1.446024e1, -1.097364e1, -3.658554e1, -4.143897e1, 4.143897e1, 3.658554e1, 4.318227e1, -3.969567e1),
    7105: (5.172500e1, 5.204528e1, -5.204528e1, -5.172500e1, 6.397837, 5.844312e1, -4.564745e1, -3.835946e1, -3.887506e1, 3.887506e1, 3.835946e1, 4.527290e1, -3.247723e1),
    7106: (-1.917973e1, -1.943753e1, 1.943753e1, 1.917973e1, 3.198918, 2.263645e1, -1.623862e1, -6.422195e1, -6.489771e1, 6.489771e1, 6.422195e1, 6.809663e1, -6.169879e1),
    7108: (1.936075e1, 1.910284e1, -1.910284e1, -1.936075e1, 3.192497, 2.255325e1, -1.616825e1, 6.470583e1, 6.403009e1, -6.403009e1, -6.470583e1, 6.789832e1, -6.151333e1),
    7109: (3.822218e1, 3.629003e1, -3.629003e1, -3.822218e1, 4.152807, 4.237498e1, -3.406937e1, -3.022111e1, -2.523958e1, 2.523958e1, 3.022111e1, 3.437391e1, -2.606830e1),
    7112: (1.254189e1, 1.507649e1, -1.507649e1, -1.254189e1, 2.071241, 1.714773e1, -1.300525e1, 4.323436e1, 4.920013e1, -4.920013e1, -4.323436e1, 5.127137e1, -4.712889e1),
}

# Nx, Ny, Nxy and sx, sy, txy, von Mises at the positive shell fibre.
MSC_CQUAD4_RESULTS = {
    1001: (6.715000, 1.210035e1, -2.278013, 3.403141, 6.063821, -1.088615, 5.592004),
    1064: (8.012225, 9.254911, 3.336241e-1, 4.093295, 4.647607, 2.060890e-1, 4.411202),
    1096: (9.520939, 7.786807, -2.454638e-1, 4.873770, 3.926027, -9.968547e-2, 4.479127),
    1106: (3.194841e1, 2.569418e1, 2.292160e1, 1.598233e1, 1.284843e1, 1.146436e1, 2.468732e1),
    1128: (1.011816e1, 8.811956, -1.721577, 5.208381, 4.455737, -8.509609e-1, 5.093727),
    1151: (7.399487e-1, -3.881808, 1.653563e1, 2.783689e-1, -2.144658, 8.178474, 1.435048e1),
    1152: (2.259203, -1.201187e1, 1.141972e1, 9.403380e-1, -6.528524, 5.684102, 1.210669e1),
    1155: (2.328659, 1.614440, -1.507198e1, 1.148011, 7.818509e-1, -7.405272, 1.286646e1),
    1192: (4.801901, 1.351162e1, 1.926801, 2.480541, 6.772149, 1.001731, 6.182789),
    1230: (3.780788e1, 1.273233e1, 1.016028e1, 1.891001e1, 6.365088, 5.074481, 1.884118e1),
    1231: (3.188634e1, 2.656965e1, -2.112242e1, 1.595124e1, 1.328614e1, -1.056477e1, 2.353458e1),
    1256: (1.461237, -3.905959e-1, -2.182878e-1, 8.769323e-1, -1.452785e-1, -1.040949e-1, 9.746911e-1),
    2001: (3.281703e-1, -6.177590e-1, -6.844019e-1, -5.124774, 1.051409, 3.882257, 8.830219),
    2013: (4.381458, 8.655326e-1, 1.049993, 1.743851e1, 2.503325, 2.100358, 1.673166e1),
    2032: (4.788122e-1, -6.544561e-1, -8.999107e-1, -6.398695, 9.845112e-1, 4.840159, 1.088548e1),
}

MSC_CTRIA3_RESULTS = {
    5001: (2.282894, 1.238230e1, -1.490849, 1.164108, 6.215184, -6.790882e-1, 5.842260),
    5008: (9.204722, 3.337836e1, 1.148728e1, 4.633539, 1.665979e1, 5.752853, 1.791955e1),
    5016: (1.710870e1, -9.917123, 4.082125, 8.566326, -4.951472, 2.046102, 1.236424e1),
    5024: (-1.504183e1, -6.787716, -1.250633e1, -7.531513, -3.367955, -6.276157, 1.268351e1),
    5032: (1.855988e1, 3.179518e1, 1.770983e1, 9.235324, 1.601316e1, 8.848541, 2.070546e1),
    5040: (1.593967e1, 1.875634e1, 1.496562e1, 7.873923, 9.508310, 7.370407, 1.550831e1),
    5048: (-9.686528, -2.949691e1, -1.216357e1, -5.372273, -1.530767e1, -6.355631, 1.738192e1),
    5050: (-2.278585e1, -1.136446e1, -1.340770e1, -1.185329e1, -6.258149, -6.985918, 1.587118e1),
    5051: (3.608114, -4.354444, -1.069074e1, 1.763290, -2.407699, -5.265838, 9.815234),
    5064: (1.096055e1, 1.077981e1, 6.888840, 5.571722, 5.304486, 3.445131, 8.076715),
    5065: (-1.223495, 2.899547e1, 5.805656, -5.912412e-1, 1.451627e1, 2.969716, 1.568795e1),
    5096: (2.038524e1, 2.038499e1, 1.412040e1, 1.015398e1, 1.031397e1, 7.056507, 1.594165e1),
    5111: (-1.821839, -1.590698e1, 6.109672, -1.041148, -8.533697, 3.031708, 9.622727),
    5112: (-1.385217e1, -3.992323e1, -1.649518e1, -7.460320, -2.052220e1, -8.523336, 2.327352e1),
    5128: (3.245955e1, 3.464449e1, 2.505899e1, 1.631543e1, 1.723217e1, 1.252775e1, 2.743764e1),
}


def _relative_l2(actual, reference):
    actual_array = np.asarray(actual, dtype=float)
    reference_array = np.asarray(reference, dtype=float)
    denominator = np.linalg.norm(reference_array)
    assert denominator > 0.0
    return np.linalg.norm(actual_array - reference_array) / denominator


def _cbar_force_tuple(data):
    end_a = data["forces_n1"]
    end_b = data["forces_n2"]
    return (
        -end_a["Mz"],
        end_a["My"],
        end_b["Mz"],
        -end_b["My"],
        -end_a["Fy"],
        -end_a["Fz"],
        data["axial_force"],
        data["torque"],
    )


def _cbar_stress_tuple(data):
    end_a = data["stress_points_end_a"]
    end_b = data["stress_points_end_b"]
    values = tuple(end_a[f"P{i}"]["sigma_bending"] for i in range(4))
    values += (
        data["axial_stress"],
        max(end_a[f"P{i}"]["sigma_total"] for i in range(4)),
        min(end_a[f"P{i}"]["sigma_total"] for i in range(4)),
    )
    values += tuple(end_b[f"P{i}"]["sigma_bending"] for i in range(4))
    values += (
        max(end_b[f"P{i}"]["sigma_total"] for i in range(4)),
        min(end_b[f"P{i}"]["sigma_total"] for i in range(4)),
    )
    return values


def _shell_core_tuple(data):
    return (
        data["membrane_force_x"],
        data["membrane_force_y"],
        data["membrane_force_xy"],
        data["stress_top_fiber_x"],
        data["stress_top_fiber_y"],
        data["tau_top_fiber_xy"],
        data["von_mises_top"],
    )


@pytest.fixture(scope="module")
def solved_t_frame():
    model = parse_nastran_input(str(CASE_FILE))
    solver = LinearStaticSolver(model)
    solver.assemble()
    return model, solver, solver.solve()


def test_final_t_frame_contains_the_intended_supported_physics(solved_t_frame):
    model, _solver, _results = solved_t_frame
    type_counts = Counter(type(element) for element in model.elements.values())
    property_counts = Counter(element.prop_id for element in model.elements.values())

    assert len(model.nodes) == 417
    assert len(model.elements) == 506
    assert type_counts == {
        RodElement: 2,
        BarElement: 88,
        CTRIA3Element: 128,
        CQUAD4Element: 288,
    }
    assert property_counts == {
        10: 256,
        11: 128,
        12: 32,
        20: 56,
        21: 20,
        22: 6,
        23: 6,
        30: 2,
    }
    assert model.k6rot == pytest.approx(100.0)
    compact_deck = "".join(CASE_FILE.read_text(encoding="utf-8").split()).upper()
    assert "PARAM,AUTOSPC,NO" in compact_deck
    assert "PARAM,SNORM,0.0" in compact_deck


def test_final_t_frame_geometry_matches_requested_t_layout(solved_t_frame):
    model, _solver, _results = solved_t_frame

    web_nodes = {
        node_id
        for element in model.elements.values()
        if element.prop_id in {10, 11}
        for node_id in element.nodes
    }
    base_nodes = {
        node_id
        for element in model.elements.values()
        if element.prop_id == 12
        for node_id in element.nodes
    }
    web_coords = np.array([model.nodes[node_id].get_coords() for node_id in web_nodes])
    base_coords = np.array([model.nodes[node_id].get_coords() for node_id in base_nodes])

    assert web_coords[:, 0].min() == pytest.approx(0.0)
    assert web_coords[:, 0].max() == pytest.approx(1200.0)
    assert web_coords[:, 1].min() == pytest.approx(0.0)
    assert web_coords[:, 1].max() == pytest.approx(600.0)
    assert web_coords[:, 2] == pytest.approx(np.zeros(len(web_coords)))
    assert base_coords[:, 0].min() == pytest.approx(0.0)
    assert base_coords[:, 0].max() == pytest.approx(1200.0)
    assert base_coords[:, 1] == pytest.approx(np.zeros(len(base_coords)))
    assert base_coords[:, 2].min() == pytest.approx(-150.0)
    assert base_coords[:, 2].max() == pytest.approx(150.0)

    assert model.properties[12].t == pytest.approx(4.0)
    web_brace = model.properties[22]
    base_brace = model.properties[23]
    assert (base_brace.A, base_brace.I1, base_brace.I2, base_brace.J) == pytest.approx(
        tuple(2.0 * value for value in (web_brace.A, web_brace.I1, web_brace.I2, web_brace.J))
    )

    for web_eid, base_eid, mirror_web_eid, mirror_base_eid in (
        (7101, 7102, 7103, 7104),
        (7105, 7106, 7107, 7108),
        (7109, 7110, 7111, 7112),
    ):
        web_segment = model.elements[web_eid]
        base_segment = model.elements[base_eid]
        mirror_web = model.elements[mirror_web_eid]
        mirror_base = model.elements[mirror_base_eid]
        assert (web_segment.prop_id, base_segment.prop_id) == (22, 23)
        assert (mirror_web.prop_id, mirror_base.prop_id) == (22, 23)
        assert web_segment.nodes[0] == mirror_web.nodes[0]

        web_point = model.nodes[web_segment.nodes[0]].get_coords()
        split_point = model.nodes[web_segment.nodes[1]].get_coords()
        base_point = model.nodes[base_segment.nodes[1]].get_coords()
        mirror_split = model.nodes[mirror_web.nodes[1]].get_coords()
        mirror_base_point = model.nodes[mirror_base.nodes[1]].get_coords()
        assert split_point == pytest.approx((web_point + base_point) / 2.0)
        assert mirror_split == pytest.approx((web_point + mirror_base_point) / 2.0)
        assert mirror_split == pytest.approx(split_point * np.array([1.0, 1.0, -1.0]))
        assert mirror_base_point == pytest.approx(base_point * np.array([1.0, 1.0, -1.0]))
        assert np.linalg.norm(split_point - web_point) == pytest.approx(
            np.linalg.norm(base_point - split_point)
        )

    assert set(model.elements[8001].nodes) == {376, 171}
    assert set(model.elements[8002].nodes) == {377, 359}


def test_final_t_frame_loads_and_constraints_match_requested_layout(solved_t_frame):
    model, _solver, _results = solved_t_frame
    forces = [load for load in model.loads if isinstance(load, SinglePointLoad)]
    moments = [load for load in model.loads if isinstance(load, SinglePointMoment)]
    assert len(forces) == 2
    assert len(moments) == 1

    for load in forces:
        assert load.appnode_id in {376, 377}
        assert load.module == pytest.approx(6000.0)
        rod_id = 8001 if load.appnode_id == 376 else 8002
        rod = model.elements[rod_id]
        attached_node = next(node_id for node_id in rod.nodes if node_id != load.appnode_id)
        outward = (
            model.nodes[load.appnode_id].get_coords()
            - model.nodes[attached_node].get_coords()
        )
        outward /= np.linalg.norm(outward)
        assert np.asarray(load.vector) == pytest.approx(outward)

    assert moments[0].appnode_id == 165
    assert moments[0].module == pytest.approx(250000.0)
    assert moments[0].vector == pytest.approx((0.0, 0.0, 1.0))

    constraints = {
        (constraint.node_id, str(constraint.components), constraint.value)
        for constraint in model.constraints
    }
    assert constraints == {
        (376, "23456", 0.0),
        (377, "23456", 0.0),
        (405, "123456", 0.0),
        (413, "123456", 0.0),
        (422, "13456", 0.0),
        (430, "13456", 0.0),
        (422, "2", 0.2),
        (430, "2", 0.2),
    }
    assert model.nodes[405].get_coords() == pytest.approx((300.0, 0.0, -150.0))
    assert model.nodes[413].get_coords() == pytest.approx((900.0, 0.0, -150.0))
    assert model.nodes[422].get_coords() == pytest.approx((300.0, 0.0, 150.0))
    assert model.nodes[430].get_coords() == pytest.approx((900.0, 0.0, 150.0))


def test_final_t_frame_solves_with_prescribed_motion_and_global_equilibrium(
    solved_t_frame,
):
    model, solver, results = solved_t_frame
    displacements = results["displacements"]
    assert np.all(np.isfinite(displacements))
    assert np.max(np.abs(displacements)) < 5.0
    for node_id in (422, 430):
        t2 = solver.dof_mapper.get_dof_indices(node_id)[1]
        assert displacements[t2] == pytest.approx(0.2, abs=1.0e-12)

    applied_force = np.zeros(3)
    applied_moment = np.zeros(3)
    for load in model.loads:
        vector = load.module * np.asarray(load.vector, dtype=float)
        position = model.nodes[load.appnode_id].get_coords()
        if isinstance(load, SinglePointLoad):
            applied_force += vector
            applied_moment += np.cross(position, vector)
        else:
            applied_moment += vector

    reaction_force = np.zeros(3)
    reaction_moment = np.zeros(3)
    for node_id, local_reaction in results["reactions"].items():
        node = model.nodes[node_id]
        basic_force = node.dof_transform.T @ local_reaction[:3]
        basic_moment = node.dof_transform.T @ local_reaction[3:]
        position = node.get_coords()
        reaction_force += basic_force
        reaction_moment += basic_moment + np.cross(position, basic_force)

    assert reaction_force + applied_force == pytest.approx(
        np.zeros(3), rel=1.0e-9, abs=1.0e-6
    )
    assert reaction_moment + applied_moment == pytest.approx(
        np.zeros(3), rel=1.0e-9, abs=2.0e-3
    )


def test_final_t_frame_matches_transcribed_msc_nastran_reference(solved_t_frame):
    """Correlate selected GRID/EID values without reading external results."""
    _model, solver, results = solved_t_frame
    displacements = results["displacements"]

    actual_u = np.array(
        [displacements[solver.dof_mapper.get_dof_indices(node_id)] for node_id in MSC_DISPLACEMENTS]
    )
    reference_u = np.array(list(MSC_DISPLACEMENTS.values()))
    assert _relative_l2(actual_u[:, :3], reference_u[:, :3]) < 1.1e-2
    assert _relative_l2(actual_u[:, 3:], reference_u[:, 3:]) < 3.5e-2
    assert np.max(np.abs(actual_u[:, :3] - reference_u[:, :3])) < 5.5e-2
    assert np.max(np.abs(actual_u[:, 3:] - reference_u[:, 3:])) < 2.5e-4

    actual_reactions = np.array(
        [results["reactions"][node_id] for node_id in MSC_SPC_REACTIONS]
    )
    reference_reactions = np.array(list(MSC_SPC_REACTIONS.values()))
    assert _relative_l2(actual_reactions[:, :3], reference_reactions[:, :3]) < 3.0e-4
    assert _relative_l2(actual_reactions[:, 3:], reference_reactions[:, 3:]) < 6.0e-4
    assert np.max(np.abs(actual_reactions[:, :3] - reference_reactions[:, :3])) < 5.0e-1
    assert np.max(np.abs(actual_reactions[:, 3:] - reference_reactions[:, 3:])) < 6.5e1

    actual_crod = []
    for element_id in MSC_CROD_RESULTS:
        data = results["stresses"][element_id]
        actual_crod.append(
            (data["axial_force"], data["torque"], data["axial_stress"], data["torsional_stress"])
        )
    actual_crod = np.asarray(actual_crod)
    reference_crod = np.asarray(list(MSC_CROD_RESULTS.values()))
    assert _relative_l2(actual_crod[:, (0, 2)], reference_crod[:, (0, 2)]) < 2.0e-6
    assert _relative_l2(actual_crod[:, (1, 3)], reference_crod[:, (1, 3)]) < 1.2e-2

    actual_cbar_forces = np.array(
        [_cbar_force_tuple(results["stresses"][element_id]) for element_id in MSC_CBAR_FORCES]
    )
    reference_cbar_forces = np.asarray(list(MSC_CBAR_FORCES.values()))
    # MA2/MB2/V2/axial are the dominant in-plane actions.  The orthogonal
    # MA1/MB1/V1 and torque are retained as looser 3-D diagnostics.
    assert _relative_l2(
        actual_cbar_forces[:, (1, 3, 5, 6)], reference_cbar_forces[:, (1, 3, 5, 6)]
    ) < 5.0e-4
    assert _relative_l2(actual_cbar_forces[:, 0], reference_cbar_forces[:, 0]) < 1.8e-1
    assert _relative_l2(actual_cbar_forces[:, 2], reference_cbar_forces[:, 2]) < 1.0e-1
    assert _relative_l2(actual_cbar_forces[:, 4], reference_cbar_forces[:, 4]) < 2.0e-1
    assert _relative_l2(actual_cbar_forces[:, 7], reference_cbar_forces[:, 7]) < 7.0e-2

    actual_cbar_stress = np.array(
        [_cbar_stress_tuple(results["stresses"][element_id]) for element_id in MSC_CBAR_STRESSES]
    )
    reference_cbar_stress = np.asarray(list(MSC_CBAR_STRESSES.values()))
    assert _relative_l2(actual_cbar_stress, reference_cbar_stress) < 7.0e-3
    assert np.max(np.abs(actual_cbar_stress - reference_cbar_stress)) < 3.5e-1

    actual_cquad = np.array(
        [_shell_core_tuple(results["stresses"][element_id]) for element_id in MSC_CQUAD4_RESULTS]
    )
    reference_cquad = np.asarray(list(MSC_CQUAD4_RESULTS.values()))
    assert _relative_l2(actual_cquad[:, :3], reference_cquad[:, :3]) < 1.0e-4
    assert _relative_l2(actual_cquad, reference_cquad) < 1.2e-2
    assert np.max(np.abs(actual_cquad[:, 3:] - reference_cquad[:, 3:])) < 3.5e-1

    actual_ctria = np.array(
        [_shell_core_tuple(results["stresses"][element_id]) for element_id in MSC_CTRIA3_RESULTS]
    )
    reference_ctria = np.asarray(list(MSC_CTRIA3_RESULTS.values()))
    assert _relative_l2(actual_ctria[:, :3], reference_ctria[:, :3]) < 3.0e-5
    assert _relative_l2(actual_ctria, reference_ctria) < 2.0e-2
    assert np.max(np.abs(actual_ctria[:, 3:] - reference_ctria[:, 3:])) < 4.5e-1
