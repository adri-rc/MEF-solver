"""Tests for linear CQUAD4 stiffness and result recovery."""

import numpy as np
import pytest

from src.solver.elements.cquad4_nastran_element import CQUAD4NastranElementSolver
from src.solver.elements.cquad4_reissner_mindlin_element import (
    CQUAD4ReissnerMindlinElementSolver,
)


class MockNode:
    def __init__(self, nid, coords):
        self.nid = nid
        self._coords = np.array(coords, dtype=float)
    def get_coords(self):
        return self._coords

class MockMat:
    def __init__(self, E, nu, name="Mat"):
        self.E = E
        self.nu = nu

class MockProp:
    def __init__(self, t, mat_bending=None):
        self.pid = 1
        self.t = t
        self.mat_bending_ref = mat_bending
        self.mat_id_bending = 999 if mat_bending else None

class MockDofMapper:
    """Provide the contiguous six-DOF layout used by this four-node model."""

    def get_dof_indices(self, nid):
        start = (nid - 1) * 6
        return list(range(start, start + 6))


CASES = [
    {
        "id": "1_CANONICAL_SQUARE_XY",
        "nodes": [[0.,0.,0.], [2.,0.,0.], [2.,2.,0.], [0.,2.,0.]],
        "mat_mem": (210e9, 0.3),
        "mat_bend": None,
        "thick": 0.1,
        "desc": "Baseline square in the XY plane"
    },
    {
        "id": "2_ROTATED_3D",
        "nodes": [[1.,1.,1.], [3.,1.,1.], [3.,3.,1.], [1.,3.,1.]],
        "mat_mem": (70e9, 0.33),
        "mat_bend": None,
        "thick": 0.05,
        "desc": "Square translated in 3D"
    },
    {
        "id": "3_TRAPEZOID_DISTORTED",
        "nodes": [[0.,0.,0.], [3.,0.,0.], [2.5, 2., 0.], [0.5, 1.5, 0.]],
        "mat_mem": (110e9, 0.3),
        "mat_bend": None,
        "thick": 0.01,
        "desc": "Distorted quadrilateral with a variable Jacobian"
    },
    {
        "id": "4_COMPOSITE_DECOUPLED",
        "nodes": [[0.,0.,0.], [1.,0.,0.], [1.,1.,0.], [0.,1.,0.]],
        "mat_mem": (200e9, 0.25), # Carbon-like stiff membrane
        "mat_bend": (5e9, 0.4),   # Soft foam-like bending material
        "thick": 0.2,
        "desc": "Independent membrane and bending materials"
    },
    {
        "id": "5_FAR_OFFSET",
        "nodes": [[10000., 10000., 10000.],
                  [10002., 10000., 10000.],
                  [10002., 10002., 10000.],
                  [10000., 10002., 10000.]],
        "mat_mem": (210e9, 0.3),
        "mat_bend": None,
        "thick": 0.1,
        "desc": "Floating-point precision with large coordinates"
    },
    {
        "id": "6_INCLINED_PLANAR",
        # Plane z=0.2*x-0.1*y: 3D orientation without warping.
        "nodes": [[0.,0.,0.], [2.,0.,0.4], [2.,2.,0.2], [0.,2.,-0.2]],
        "mat_mem": (210e9, 0.3),
        "mat_bend": None,
        "thick": 0.15,
        "desc": "Inclined coplanar element in 3D"
    }
]

@pytest.fixture(params=CASES, ids=[c["id"] for c in CASES])
def element_setup(request):
    """Build the solver for each fixture case."""
    data = request.param

    n1 = MockNode(1, data["nodes"][0])
    n2 = MockNode(2, data["nodes"][1])
    n3 = MockNode(3, data["nodes"][2])
    n4 = MockNode(4, data["nodes"][3])

    mat_mem = MockMat(*data["mat_mem"], name="Membrane")

    mat_bend = mat_mem
    if data["mat_bend"]:
        mat_bend = MockMat(*data["mat_bend"], name="Bending")

    prop = MockProp(t=data["thick"], mat_bending=mat_bend)

    solver = CQUAD4NastranElementSolver(
        n1, n2, n3, n4, prop, mat_mem, mat_bending=mat_bend, k6rot=100.0
    )

    return {
        "solver": solver,
        "data": data,
        "mat_mem": mat_mem,
        "mat_bend": mat_bend
    }

def test_stiffness_matrix_shape(element_setup):
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    assert K.shape == (24, 24)

def test_stiffness_symmetry(element_setup):
    """Validate global stiffness symmetry (Maxwell-Betti reciprocity)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()

    np.testing.assert_allclose(K, K.T, rtol=1e-5, atol=1e-8,
                               err_msg="The stiffness matrix is not symmetric.")

def test_nastran_hourglass_eigenvalues_on_square():
    """The two bilinear membrane modes match the Nastran E/3 result."""
    nodes = [
        MockNode(1, [0.0, 0.0, 0.0]),
        MockNode(2, [2.0, 0.0, 0.0]),
        MockNode(3, [2.0, 2.0, 0.0]),
        MockNode(4, [0.0, 2.0, 0.0]),
    ]
    young_modulus = 210000.0
    element = CQUAD4NastranElementSolver(
        *nodes,
        MockProp(1.0),
        MockMat(young_modulus, 0.3),
        mat_bending=MockMat(0.0, 0.3),
        mat_shear=MockMat(0.0, 0.3),
    )
    assert element.formulation == "nastran"

    membrane_dofs = [0, 1, 6, 7, 12, 13, 18, 19]
    membrane = element._build_local_stiffness()[np.ix_(membrane_dofs, membrane_dofs)]
    eigenvalues = np.linalg.eigvalsh(membrane)
    np.testing.assert_allclose(
        eigenvalues[3:5], [young_modulus / 3.0] * 2, rtol=1e-12
    )


def test_oblique_torsion_metric_matches_dmig_control_points():
    """Lock the planar skew continuation identified from independent DMIGs."""
    metric = CQUAD4NastranElementSolver._oblique_torsion_metric

    assert metric(0.0) == pytest.approx(1.0)
    assert metric(0.0700468506) == pytest.approx(1.03068940806, rel=1e-9)
    assert metric(0.4028759652) == pytest.approx(1.26174688986, rel=1e-9)

    values = np.array([metric(value) for value in np.linspace(0.0, 1.0, 101)])
    assert np.all(np.isfinite(values))
    assert np.all(values > 0.0)
    assert np.all(np.diff(values) > 0.0)


def test_k6rot_uses_ctria_nodal_macro_rotation_penalty():
    """CQUAD4 extends CTRIA3's theta_z - omega penalty to four nodes."""
    nodes = [
        MockNode(1, [0.0, 0.0, 0.0]),
        MockNode(2, [2.0, 0.0, 0.0]),
        MockNode(3, [2.0, 2.0, 0.0]),
        MockNode(4, [0.0, 2.0, 0.0]),
    ]
    material = MockMat(210000.0, 0.3)
    element = CQUAD4NastranElementSolver(
        *nodes, MockProp(0.2), material, k6rot=100.0
    )

    stiffness = np.zeros((24, 24))
    element._add_drilling_stiffness(stiffness)

    corner_derivatives = (
        (
            np.array([-0.5, 0.5, 0.0, 0.0]),
            np.array([-0.5, 0.0, 0.0, 0.5]),
        ),
        (
            np.array([-0.5, 0.5, 0.0, 0.0]),
            np.array([0.0, -0.5, 0.5, 0.0]),
        ),
        (
            np.array([0.0, 0.0, 0.5, -0.5]),
            np.array([0.0, -0.5, 0.5, 0.0]),
        ),
        (
            np.array([0.0, 0.0, 0.5, -0.5]),
            np.array([-0.5, 0.0, 0.0, 0.5]),
        ),
    )
    expected_local = np.zeros((12, 12))
    shear_modulus = material.E / (2.0 * (1.0 + material.nu))
    nodal_area = 4.0
    k_theta = 100.0 * 1.0e-6 * shear_modulus * 0.2 * nodal_area
    for residual_node, (dN_dx, dN_dy) in enumerate(corner_derivatives):
        residual = np.zeros(12)
        for node_index in range(4):
            base = 3 * node_index
            residual[base] = 0.5 * dN_dy[node_index]
            residual[base + 1] = -0.5 * dN_dx[node_index]
        residual[3 * residual_node + 2] = 1.0
        expected_local += k_theta * np.outer(residual, residual)
    drill_dofs = [0, 1, 5, 6, 7, 11, 12, 13, 17, 18, 19, 23]

    np.testing.assert_allclose(
        stiffness[np.ix_(drill_dofs, drill_dofs)], expected_local
    )
    inactive_dofs = sorted(set(range(24)) - set(drill_dofs))
    assert not np.any(stiffness[inactive_dofs, :])


def test_reissner_mindlin_is_an_independent_full_integration_formulation():
    nodes = [
        MockNode(1, [0.0, 0.0, 0.0]),
        MockNode(2, [2.0, 0.0, 0.0]),
        MockNode(3, [2.2, 1.5, 0.0]),
        MockNode(4, [0.1, 1.8, 0.0]),
    ]
    prop = MockProp(0.2)
    material = MockMat(210000.0, 0.3)
    nastran = CQUAD4NastranElementSolver(
        *nodes, prop, material, mat_bending=material, k6rot=100.0
    )
    mindlin = CQUAD4ReissnerMindlinElementSolver(
        *nodes, prop, material, mat_bending=material, k6rot=100.0
    )

    assert nastran.formulation == "nastran"
    assert mindlin.formulation == "reissner-mindlin"
    np.testing.assert_allclose(
        mindlin.get_stiffness_matrix(), mindlin.get_stiffness_matrix().T
    )
    assert not np.allclose(
        mindlin.get_stiffness_matrix(), nastran.get_stiffness_matrix()
    )

def test_rigid_body_modes(element_setup):
    """Validate the six near-zero rigid-body modes."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()

    eigenvalues = np.sort(np.linalg.eigvalsh(K))

    # A shell with drilling stabilization retains six rigid-body modes.
    max_k = np.max(np.abs(K))
    tolerance = max_k * 1e-9

    zero_modes = np.count_nonzero(np.abs(eigenvalues) < tolerance)
    assert zero_modes == 6, \
        f"Expected 6 rigid-body modes, found {zero_modes}: {eigenvalues[:12]}"

def test_local_axes_orientation(element_setup):
    """Validate the Nastran element-axis convention."""
    solver = element_setup["solver"]
    data = element_setup["data"]

    vx = solver.local_x_axis
    vy = solver.local_y_axis
    vz = solver.local_z_axis

    assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)

    if data["id"] == "1_CANONICAL_SQUARE_XY":
        np.testing.assert_allclose(vx, [1,0,0], atol=1e-9)
        np.testing.assert_allclose(vy, [0,1,0], atol=1e-9)
        np.testing.assert_allclose(vz, [0,0,1], atol=1e-9)

def test_membrane_bending_decoupling(element_setup):
    """
    Verify that in-plane membrane and out-of-plane bending responses are
    uncoupled in element coordinates.
    """
    solver = element_setup["solver"]

    if solver.prop.t == 0.0:
        pytest.skip("Element has zero thickness.")

    K_local = solver._build_local_stiffness()

    # Node 1 Tx (membrane) versus node 1 Tz (bending).
    k_cross = K_local[0, 2]

    assert k_cross == pytest.approx(0.0, abs=1e-9), \
        "Membrane-bending coupling found in the 24x24 matrix."

def test_material_decoupling_check(element_setup):
    """
    Verify that bending stiffness uses the independent bending material.
    """
    setup = element_setup
    solver = setup["solver"]
    data = setup["data"]

    if data["id"] == "4_COMPOSITE_DECOUPLED":
        K_local = solver._build_local_stiffness()

        k_membrane_axial = K_local[0,0]
        k_bending_th = K_local[3,3]

        # The bending material is 40 times softer than the membrane material.
        assert k_bending_th < k_membrane_axial, \
            "Bending stiffness is too high; D_membrane may have replaced D_bending."

def test_stress_recovery_robustness(element_setup):
    """
    Validate stress recovery with a generic finite displacement vector.
    """
    solver = element_setup["solver"]
    np.random.seed(42)
    u_dummy_global = np.random.rand(100) * 0.001

    try:
        results = solver.calculate_stresses(u_dummy_global, MockDofMapper())
    except Exception as e:
        pytest.fail(f"Stress recovery raised an exception: {e}")

    required_keys = [
        'sigma_x', 'sigma_y', 'tau_xy', 'von_mises',
        'bending_moment_x', 'shear_force_x',
        'stress_top_fiber_x', 'stress_top_fiber_y', 'tau_top_fiber_xy',
        'stress_bot_fiber_x', 'stress_bot_fiber_y', 'tau_bot_fiber_xy'
    ]

    for key in required_keys:
        assert key in results, f"Missing result key: {key}."
        assert not np.isnan(results[key]), f"Result {key} is NaN."

    deprecated_keys = {
        'stress_fiber_top_x', 'stress_fiber_top_y', 'tau_fiber_top_y',
        'stress_fiber_bot_x', 'stress_fiber_bot_y', 'tau_fiber_bot_y',
        'stress_z1_fiber_x', 'stress_z1_fiber_y', 'tau_z1_fiber_xy',
        'stress_z2_fiber_x', 'stress_z2_fiber_y', 'tau_z2_fiber_xy'
    }
    assert deprecated_keys.isdisjoint(results)

def test_constant_strain_patch_x_tension():
    """Apply pure uniaxial strain to an isoparametric membrane patch."""
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [1,1,0])
    n4 = MockNode(4, [0,1,0])

    prop = MockProp(t=1.0)
    mat = MockMat(E=1000.0, nu=0.3)
    elem = CQUAD4NastranElementSolver(
        n1, n2, n3, n4, prop, mat, mat_bending=mat
    )

    eps_x = 0.001
    eps_y = -mat.nu * eps_x

    u_global = np.zeros(24)
    u_global[6] = eps_x
    u_global[12] = eps_x
    u_global[13] = eps_y
    u_global[19] = eps_y

    res = elem.calculate_stresses(u_global, MockDofMapper())

    expected_sig_x = mat.E * eps_x

    assert res['sigma_x'] == pytest.approx(expected_sig_x, rel=1e-5)
    assert res['sigma_y'] == pytest.approx(0.0, abs=1e-5)
    assert res['tau_xy'] == pytest.approx(0.0, abs=1e-5)

def test_degenerate_element_error():
    """Reject a zero-area quadrilateral with collinear nodes."""
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [2,0,0])
    n4 = MockNode(4, [3,0,0])

    prop = MockProp(t=1.0)
    mat = MockMat(E=1e5, nu=0.3)

    with pytest.raises(ValueError, match="Degenerate CQUAD4"):
        CQUAD4NastranElementSolver(n1, n2, n3, n4, prop, mat)


def test_non_coplanar_element_error():
    """Reject warped CQUAD4 geometry outside the solver scope."""
    nodes = [
        MockNode(1, [0.0, 0.0, 0.0]),
        MockNode(2, [2.0, 0.0, 0.0]),
        MockNode(3, [2.0, 2.0, 0.1]),
        MockNode(4, [0.0, 2.0, 0.0]),
    ]
    prop = MockProp(t=0.1)
    mat = MockMat(E=2.1e5, nu=0.3)

    with pytest.raises(ValueError, match="Non-coplanar CQUAD4"):
        CQUAD4NastranElementSolver(*nodes, prop, mat, mat_bending=mat)

def test_negative_jacobian_error():
    """Reject a misnumbered bow-tie quadrilateral."""
    n1 = MockNode(1, [0, 0, 0])
    n2 = MockNode(2, [1.5, 1, 0])
    n3 = MockNode(3, [1, 0, 0]) # Crossing edge order
    n4 = MockNode(4, [0, 1.5, 0])

    prop = MockProp(t=1.0)
    mat = MockMat(E=1e5, nu=0.3)

    # Geometry is validated during construction even when K6ROT=0 or the
    # PSHELL leaves plate blocks inactive.
    with pytest.raises(ValueError, match="Negative or zero Jacobian"):
        CQUAD4NastranElementSolver(
            n1, n2, n3, n4, prop, mat, mat_bending=mat
        )
