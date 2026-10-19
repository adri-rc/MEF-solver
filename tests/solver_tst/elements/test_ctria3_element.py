"""Tests for linear CTRIA3 stiffness and result recovery."""

import numpy as np
import pytest

from src.solver.elements.ctria3_element import CTRIA3ElementSolver


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
    """Provide the contiguous six-DOF layout used by this three-node model."""

    def get_dof_indices(self, nid):
        start = (nid - 1) * 6
        return list(range(start, start + 6))


CASES = [
    {
        "id": "1_CANONICAL_XY",
        "nodes": [[0.,0.,0.], [2.,0.,0.], [0.,2.,0.]],
        "mat_mem": (210e9, 0.3),
        "mat_bend": None,        # Explicit MID2/MID3 with MID1 properties
        "thick": 0.1,
        "area": 2.0,
        "desc": "Baseline case in the XY plane"
    },
    {
        "id": "2_ROTATED_3D",
        "nodes": [[1.,1.,1.], [2.,2.,1.], [1.,2.,2.]],
        "mat_mem": (70e9, 0.33),
        "mat_bend": None,
        "thick": 0.05,
        "area": 0.866025,
        "desc": "Triangle rotated in 3D"
    },
    {
        "id": "3_NEEDLE_SHAPE",
        "nodes": [[0.,0.,0.], [0.1, 0., 0.], [0.05, 10., 0.]],
        "mat_mem": (110e9, 0.3),
        "mat_bend": None,
        "thick": 0.01,
        "area": 0.5,
        "desc": "Numerical stability with a high aspect ratio"
    },
    {
        "id": "4_COMPOSITE_DECOUPLED",
        "nodes": [[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]],
        "mat_mem": (200e9, 0.25), # Carbon-like stiff membrane
        "mat_bend": (5e9, 0.4),   # Soft foam-like bending material
        "thick": 0.2,
        "area": 0.5,
        "desc": "Independent membrane and bending materials"
    },
    {
        "id": "5_FAR_OFFSET",
        "nodes": [[10000., 10000., 10000.],
                  [10002., 10000., 10000.],
                  [10000., 10002., 10000.]],
        "mat_mem": (210e9, 0.3),
        "mat_bend": None,
        "thick": 0.1,
        "area": 2.0,
        "desc": "Floating-point precision with large coordinates"
    },
    {
        "id": "6_SKEWED_SPACE_COMPOSITE",
        "nodes": [
            [1.2, 3.4, -0.5],  # N1: arbitrary point
            [5.6, 1.1, 2.3],   # N2: distant and rotated
            [2.2, 6.8, 4.1]    # N3: oblique scalene triangle
        ],
        # Stiff steel-like membrane and soft polymer-like bending response.
        "mat_mem": (210e9, 0.3),
        "mat_bend": (2e9, 0.45),
        "thick": 0.15,
        # Area from the cross product: 0.5 * ||Vec12 x Vec13||.
        # Vec12=[4.4, -2.3, 2.8], Vec13=[1.0, 3.4, 4.6]
        # Cross = [-20.1, -17.44, 17.26] -> Norm=31.718 -> Area=15.859...
        "area": 15.8593127216,
        "desc": "Oblique scalene triangle with independent materials (M>B)"
    },
    {
        "id": "7_MICRO_TILTED_HYBRID",
        "nodes": [
            [0.0, 0.0, 0.0],
            [1e-3, 2e-3, 3e-3],   # Millimetric coordinates
            [1e-3, -1e-3, 2e-3]
        ],
        # Soft aluminum membrane and artificially stiff ceramic bending response.
        "mat_mem": (70e9, 0.33),
        "mat_bend": (400e9, 0.2),
        "thick": 1e-4, # Very small thickness
        # Vec12=[1, 2, 3]e-3, Vec13=[1, -1, 2]e-3
        # Cross = [7, 1, -3]e-6 -> Norm=sqrt(59)*1e-6 -> Area=3.84057e-6
        "area": 3.840572874e-6,
        "desc": "Inclined micro-element with bending stiffness above membrane stiffness"
    }
]

@pytest.fixture(params=CASES, ids=[c["id"] for c in CASES])
def element_setup(request):
    """Build a solver and reference data for each fixture case."""
    data = request.param

    n1 = MockNode(1, data["nodes"][0])
    n2 = MockNode(2, data["nodes"][1])
    n3 = MockNode(3, data["nodes"][2])

    mat_mem = MockMat(*data["mat_mem"], name="Membrane")

    bending_data = data["mat_bend"] or data["mat_mem"]
    mat_bend = MockMat(*bending_data, name="Bending/Shear")

    prop = MockProp(t=data["thick"], mat_bending=mat_bend)

    # These cases exercise the complete formulation, so MID2 and MID3 are
    # explicit even when they reuse the MID1 properties.
    solver = CTRIA3ElementSolver(
        n1,
        n2,
        n3,
        prop,
        mat_mem,
        mat_bending=mat_bend,
        mat_shear=mat_bend,
    )

    return {
        "solver": solver,
        "data": data,
        "mat_mem": mat_mem,
        "mat_bend": mat_bend
    }


def _calculate_cst_analytical_local(solver, data):
    """
    Build the analytical CST membrane stiffness in element coordinates.

    Using the solver's local geometry also covers rotated and translated cases.
    """
    E, nu = data["mat_mem"]
    t = data["thick"]
    # Use the processed local geometry so the check is independent of the
    # triangle's position and orientation in basic coordinates.
    x2 = solver.x2
    x3 = solver.x3
    y3 = solver.y3
    Area = solver.A

    # Standard CST coefficients use b_i = y_j - y_k and g_i = x_k - x_j.
    y2 = 0.0
    b1 = y2 - y3
    g1 = x3 - x2

    x1, y1 = 0.0, 0.0
    b2 = y3 - y1
    g2 = x1 - x3

    b3 = y1 - y2
    g3 = x2 - x1

    factor_B = 1.0 / (2.0 * Area)

    # Rows are eps_x, eps_y, gamma_xy; columns are (u, v) per node.
    B = factor_B * np.array([
        [b1, 0,  b2, 0,  b3, 0],
        [0,  g1, 0,  g2, 0,  g3],
        [g1, b1, g2, b2, g3, b3]
    ])

    # Plane-stress constitutive matrix.
    factor_D = E / (1 - nu**2)
    D = factor_D * np.array([
        [1,  nu, 0],
        [nu, 1,  0],
        [0,  0, (1-nu)/2]
    ])

    # B is constant, so the area integral reduces to volume times B.T @ D @ B.
    K_mem = (Area * t) * (B.T @ D @ B)

    return K_mem

def test_stiffness_matrix_shape(element_setup):
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    assert K.shape == (18, 18)

def test_stiffness_symmetry(element_setup):
    """Validate global stiffness symmetry (Maxwell-Betti reciprocity)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()

    # Relative tolerance accommodates the different force and moment scales.
    np.testing.assert_allclose(K, K.T, rtol=1e-5, atol=1e-8,
                               err_msg="The stiffness matrix is not symmetric.")

def test_rigid_body_modes(element_setup):
    """Validate the near-zero rigid-body translation modes."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()

    eigenvalues = np.linalg.eigvalsh(K)

    # The first three eigenmodes are rigid translations. Rotational modes may
    # carry the artificial drilling penalty.
    translations = eigenvalues[:3]

    # Scale the zero threshold with the largest stiffness coefficient.
    max_k = np.max(np.abs(K))
    tolerance = max_k * 1e-9

    # The needle geometry is more poorly conditioned.
    if element_setup["data"]["id"] == "3_NEEDLE_SHAPE":
        tolerance = max_k * 1e-7

    assert np.all(np.abs(translations) < tolerance), \
        f"Rigid-body mode check failed. Eigenvalues: {translations}"

def test_area_calculation(element_setup):
    """Validate the internally computed area against the reference value."""
    solver = element_setup["solver"]
    expected_area = element_setup["data"]["area"]

    # Covers rounded reference areas for the rotated and needle cases.
    assert solver.A == pytest.approx(expected_area, rel=1e-4)

def test_local_axes_orientation(element_setup):
    """
    Validate the element-axis orientation.

    The canonical case has exact expected axes; the other cases check
    orthogonality.
    """
    solver = element_setup["solver"]
    data = element_setup["data"]

    vx = solver.local_x_axis
    vy = solver.local_y_axis
    vz = solver.local_z_axis

    assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)

    if data["id"] == "1_CANONICAL_XY":
        np.testing.assert_allclose(vx, [1,0,0], atol=1e-9)
        np.testing.assert_allclose(vy, [0,1,0], atol=1e-9)
        np.testing.assert_allclose(vz, [0,0,1], atol=1e-9)

def test_membrane_bending_decoupling(element_setup):
    """
    For an isotropic flat plate, in-plane forces should not cause
    out-of-plane deformations (Tz, Rx, Ry).
    """
    solver = element_setup["solver"]

    K_local = solver._build_local_stiffness()

    # Node 1 Tx and Tz sample the membrane-to-bending coupling block.
    k_cross = K_local[0, 2]

    assert k_cross == pytest.approx(0.0, abs=1e-9), \
        "Coupling detected between Membrane (x) and Bending (z). Check B-matrices."

def test_material_decoupling_check(element_setup):
    """
    Check that case 4 uses the bending material rather than the membrane one.
    """
    setup = element_setup
    solver = setup["solver"]
    data = setup["data"]

    # This targeted comparison applies only to the independent-material case.
    if data["id"] == "4_COMPOSITE_DECOUPLED":
        K_local = solver._build_local_stiffness()

        # Membrane stiffness scales with E_mem*t; bending stiffness scales with
        # E_bend*t^3.

        k_membrane_axial = K_local[0,0]
        k_bending_w = K_local[2,2]

        # A mistaken MID1 inheritance would make this coefficient much larger.
        assert k_bending_w < k_membrane_axial, \
            "Bending stiffness is too high; the membrane material may have leaked into MID2."

def test_bending_stiffness_generated(element_setup):
    """Check that the plate formulation produces nonzero bending stiffness."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    prop = element_setup["data"]["thick"]

    if element_setup["data"]["id"] == "3_NEEDLE_SHAPE":
        return

    assert abs(K[2,2]) > 0.0, "Bending stiffness K_zz should be non-zero"

def test_stress_recovery_robustness(element_setup):
    """Validate finite stress recovery and the public result keys."""
    solver = element_setup["solver"]

    # Small deterministic trial displacement vector.
    np.random.seed(42)
    u_dummy_global = np.random.rand(100) * 0.001 # The mapper selects the relevant entries.

    class MockMapper:
        def get_dof_indices(self, nid):
            base = (nid - 1) * 6
            return list(range(base, base+6))

    try:
        results = solver.calculate_stresses(u_dummy_global, MockMapper())
    except Exception as e:
        pytest.fail(f"Stress recovery raised an exception: {e}")

    # Public result contract.
    required_keys = [
        'sigma_x', 'sigma_y', 'tau_xy', 'von_mises',
        'bending_moment_x',
        'stress_top_fiber_x', 'stress_top_fiber_y', 'tau_top_fiber_xy',
        'stress_bot_fiber_x', 'stress_bot_fiber_y', 'tau_bot_fiber_xy'
    ]

    for key in required_keys:
        assert key in results, f"Missing stress-result key: {key}."
        assert not np.isnan(results[key]), f"Stress result {key} is NaN."

    deprecated_keys = {
        'stress_fiber_top_x', 'stress_fiber_top_y', 'tau_fiber_top_y',
        'stress_fiber_bot_x', 'stress_fiber_bot_y', 'tau_fiber_bot_y',
        'stress_z1_fiber_x', 'stress_z1_fiber_y', 'tau_z1_fiber_xy',
        'stress_z2_fiber_x', 'stress_z2_fiber_y', 'tau_z2_fiber_xy'
    }
    assert deprecated_keys.isdisjoint(results)

def test_bending_rigid_body_z(element_setup):
    """
    A constant transverse translation produces no plate stress. The tolerance
    scales with material stiffness to accommodate floating-point
    noise in high-stiffness or micro-scale elements.
    """
    solver = element_setup["solver"]
    data = element_setup["data"]

    u_global = np.zeros(100)
    u_global[2] = 0.1
    u_global[8] = 0.1
    u_global[14] = 0.1

    class MockMapper:
        def get_dof_indices(self, nid):
            base = (nid - 1) * 6
            return list(range(base, base+6))

    res = solver.calculate_stresses(u_global, MockMapper())

    # Floating-point stress noise scales with Young's modulus. E*1e-12 is a
    # conservative engineering-zero threshold.
    if element_setup["mat_bend"]:
        E_ref = element_setup["mat_bend"].E
    else:
        E_ref = element_setup["mat_mem"].E

    # At least 1e-5 Pa, otherwise relative to the material stiffness.
    tol = max(1e-5, E_ref * 1e-12)

    assert res['bending_moment_x'] == pytest.approx(0.0, abs=tol)
    assert res['stress_top_fiber_x'] == pytest.approx(0.0, abs=tol)

def test_constant_strain_patch_x_tension():
    """Apply uniform uniaxial strain to a right-triangle membrane patch."""
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [0,1,0])

    prop = MockProp(t=1.0)
    mat = MockMat(E=1000.0, nu=0.3)

    elem = CTRIA3ElementSolver(n1, n2, n3, prop, mat)

    eps_x = 0.001
    eps_y = -mat.nu * eps_x

    u_global = np.zeros(18)

    u_global[6] = eps_x
    u_global[13] = eps_y

    mapper = MockDofMapper()
    res = elem.calculate_stresses(u_global, mapper)

    # Compatible Poisson contraction gives sigma_x = E*eps_x and sigma_y = 0.
    expected_sig_x = mat.E * eps_x

    assert res['sigma_x'] == pytest.approx(expected_sig_x, rel=1e-5)
    assert res['sigma_y'] == pytest.approx(0.0, abs=1e-5)
    assert res['tau_xy'] == pytest.approx(0.0, abs=1e-5)
    assert res['von_mises'] == pytest.approx(expected_sig_x, rel=1e-5)

def test_zero_area_error():
    """Reject collinear triangle nodes."""
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [2,0,0]) # Collinear
    prop = MockProp(t=1.0)
    mat = MockMat(E=1e5, nu=0.3)

    with pytest.raises(ValueError):
        CTRIA3ElementSolver(n1, n2, n3, prop, mat)


def test_membrane_stiffness_analytical_check(element_setup):
    """
    Compare the membrane block term by term with the analytical CST matrix.
    """
    data = element_setup["data"]
    solver = element_setup["solver"]

    K_mem_ref = _calculate_cst_analytical_local(solver, data)

    K_solver_full = solver._build_local_stiffness()

    mem_indices = [0, 1, 6, 7, 12, 13]
    K_solver_mem = K_solver_full[np.ix_(mem_indices, mem_indices)]

    # A 1e-5 relative tolerance covers the extreme geometries in double precision.
    np.testing.assert_allclose(
        K_solver_mem, K_mem_ref, rtol=1e-5, atol=1e-8,
        err_msg=f"Membrane stiffness mismatch for case {data['id']}"
    )

def test_bending_stiffness_integration_check(element_setup):
    """
    Reconstruct the Reissner-Mindlin plate stiffness as bending plus shear and
    compare it with the solver matrix.
    """
    data = element_setup["data"]
    solver = element_setup["solver"]

    # Skip a membrane-only property.
    if not (data["thick"] > 0.0):
        pytest.skip("Element is pure membrane (t=0)")

    # Constant pure-bending contribution over the element area.
    B_b = solver._get_B_bending_Nastran()
    D_b = solver._get_D_bending()
    K_bend_ref = (B_b.T @ D_b @ B_b) * solver.A

    # One-point transverse-shear contribution.
    N_cent = 1.0 / 3.0
    B_s = solver._get_B_shear_Nastran(N_cent, N_cent, N_cent)
    D_s = solver._get_D_shear()
    K_shear_ref = (B_s.T @ D_s @ B_s) * solver.A

    K_plate_theoretical = K_bend_ref + K_shear_ref

    # Extract the solver's plate block for a like-for-like comparison.
    K_local_full = solver._build_local_stiffness()

    # Plate DOFs are w(2), Rx(3), and Ry(4) at each six-DOF node.
    # N1: 2, 3, 4 | N2: 8, 9, 10 | N3: 14, 15, 16
    bend_indices_map = [2, 3, 4, 8, 9, 10, 14, 15, 16]

    K_plate_solver = np.zeros((9, 9))
    for i, row in enumerate(bend_indices_map):
        for j, col in enumerate(bend_indices_map):
            K_plate_solver[i, j] = K_local_full[row, col]

    np.testing.assert_allclose(K_plate_solver, K_plate_theoretical,
                                rtol=1e-5, atol=1e-6,
                                err_msg="MIN3 integration mismatch against the reconstructed matrix.")
