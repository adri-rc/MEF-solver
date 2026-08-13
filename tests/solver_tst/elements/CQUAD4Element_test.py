"""
@file CQUAD4ElementSolver_test.py
@brief Unit tests for the CQUAD4 Element Solver.
@details
Validates the physics of the 4-node Quadrilateral Shell element:
1. Stiffness Matrix Assembly (Membrane + RM Bending + Shear).
2. Geometric Transformations (Nastran standard local axes).
3. Material Decoupling & Membrane-Bending decoupling.
4. Stress Recovery.
5. Constant Strain Patch Test.
6. Error Handling (Degenerate elements).
"""

import pytest
import numpy as np
from src.solver.elements.cquad4_nastran_element import CQUAD4NastranElementSolver
from src.solver.elements.cquad4_reissner_mindlin_element import (
    CQUAD4ReissnerMindlinElementSolver,
)


# --- MOCKS (Data Objects Simulation) ---
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
    """Simulates the global-to-local index mapping (24 DOFs for CQUAD4)."""
    def get_dof_indices(self, nid):
        # Maps Node ID to a block of 6 indices
        start = (nid - 1) * 6
        return list(range(start, start + 6))

# --- DEFINICIÓN DE CASOS (FIXTURE DATA) ---

CASES = [
    {
        "id": "1_CANONICAL_SQUARE_XY",
        "nodes": [[0.,0.,0.], [2.,0.,0.], [2.,2.,0.], [0.,2.,0.]], # Cuadrado perfecto
        "mat_mem": (210e9, 0.3), # Acero
        "mat_bend": None,
        "thick": 0.1,
        "desc": "Cuadrado base en el plano XY"
    },
    {
        "id": "2_ROTATED_3D",
        "nodes": [[1.,1.,1.], [3.,1.,1.], [3.,3.,1.], [1.,3.,1.]], # Cuadrado trasladado
        "mat_mem": (70e9, 0.33), # Aluminio
        "mat_bend": None,
        "thick": 0.05,
        "desc": "Cuadrado trasladado en 3D"
    },
    {
        "id": "3_TRAPEZOID_DISTORTED",
        "nodes": [[0.,0.,0.], [3.,0.,0.], [2.5, 2., 0.], [0.5, 1.5, 0.]], # Trapecio
        "mat_mem": (110e9, 0.3), # Titanio
        "mat_bend": None,
        "thick": 0.01,
        "desc": "Elemento cuadrangular distorsionado (Jacobiano variable)"
    },
    {
        "id": "4_COMPOSITE_DECOUPLED",
        "nodes": [[0.,0.,0.], [1.,0.,0.], [1.,1.,0.], [0.,1.,0.]],
        "mat_mem": (200e9, 0.25), # Carbono (Membrana rígida)
        "mat_bend": (5e9, 0.4),   # Foam Core (Flexión blanda)
        "thick": 0.2,
        "desc": "Materiales distintos para membrana y flexión"
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
        "desc": "Test de precisión flotante con coordenadas grandes"
    },
    {
        "id": "6_SLIGHTLY_WARPED",
        "nodes": [[0.,0.,0.], [2.,0.,0.], [2.,2.,0.1], [0.,2.,0.]], # Un nodo fuera del plano
        "mat_mem": (210e9, 0.3), 
        "mat_bend": None, 
        "thick": 0.15,
        "desc": "Elemento ligeramente alabeado (warped)"
    }
]

# --- FIXTURE PARAMETRIZADO ---

@pytest.fixture(params=CASES, ids=[c["id"] for c in CASES])
def element_setup(request):
    """
    Crea el solver para cada uno de los casos definidos.
    """
    data = request.param
    
    # 1. Crear Nodos
    n1 = MockNode(1, data["nodes"][0])
    n2 = MockNode(2, data["nodes"][1])
    n3 = MockNode(3, data["nodes"][2])
    n4 = MockNode(4, data["nodes"][3])
    
    # 2. Crear Materiales
    mat_mem = MockMat(*data["mat_mem"], name="Membrane")
    
    mat_bend = mat_mem
    if data["mat_bend"]:
        mat_bend = MockMat(*data["mat_bend"], name="Bending")
        
    # 3. Propiedad
    prop = MockProp(t=data["thick"], mat_bending=mat_bend)
    
    # 4. Instanciar Solver
    solver = CQUAD4NastranElementSolver(
        n1, n2, n3, n4, prop, mat_mem, mat_bending=mat_bend, k6rot=100.0
    )
    
    return {
        "solver": solver,
        "data": data,
        "mat_mem": mat_mem,
        "mat_bend": mat_bend
    }

# --- TESTS ---

def test_stiffness_matrix_shape(element_setup):
    """@test Validate Stiffness Matrix Shape (24x24 for CQUAD4)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    assert K.shape == (24, 24)

def test_stiffness_symmetry(element_setup):
    """Valida la simetría de la matriz K global (Teorema de Maxwell-Betti)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    
    np.testing.assert_allclose(K, K.T, rtol=1e-5, atol=1e-8, 
                               err_msg="La matriz de rigidez no es simétrica.")

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
    """Valida que existen 6 modos de cuerpo rígido (Autovalores cercanos a 0)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    
    eigenvalues = np.sort(np.linalg.eigvalsh(K))
    
    # Para elementos de placa (con drilling DOF artificial), debería haber 6 modos ~ 0.
    max_k = np.max(np.abs(K))
    tolerance = max_k * 1e-9

    zero_modes = np.count_nonzero(np.abs(eigenvalues) < tolerance)
    assert zero_modes == 6, \
        f"Expected 6 rigid-body modes, found {zero_modes}: {eigenvalues[:12]}"

def test_local_axes_orientation(element_setup):
    """
    @test Valida la orientación de los ejes locales (Normativa Nastran).
    """
    solver = element_setup["solver"]
    data = element_setup["data"]
    
    vx = solver.local_x_axis
    vy = solver.local_y_axis
    vz = solver.local_z_axis
    
    # 1. Ortogonalidad pura
    assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)
    
    # 2. Direcciones específicas para el cuadrado canónico
    if data["id"] == "1_CANONICAL_SQUARE_XY":
        np.testing.assert_allclose(vx, [1,0,0], atol=1e-9)
        np.testing.assert_allclose(vy, [0,1,0], atol=1e-9)
        np.testing.assert_allclose(vz, [0,0,1], atol=1e-9)

def test_membrane_bending_decoupling(element_setup):
    """
    @test Verifica que Membrana (in-plane) y Flexión (out-of-plane) 
    están desacopladas en el sistema local.
    """
    solver = element_setup["solver"]
    
    if solver.prop.t == 0.0:
        pytest.skip("Elemento sin espesor.")
        
    K_local = solver._build_local_stiffness()
    
    # Nodo 1 Tx (índice 0, membrana) cruzado con Nodo 1 Tz (índice 2, flexión)
    k_cross = K_local[0, 2]
    
    assert k_cross == pytest.approx(0.0, abs=1e-9), \
        "Acoplamiento detectado entre Membrana y Flexión. Revisa el ensamblaje de 24x24."

def test_material_decoupling_check(element_setup):
    """
    Verifica que la rigidez de flexión usa el módulo de flexión si es distinto.
    """
    setup = element_setup
    solver = setup["solver"]
    data = setup["data"]
    
    if data["id"] == "4_COMPOSITE_DECOUPLED":
        K_local = solver._build_local_stiffness()
        
        k_membrane_axial = K_local[0,0] # Membrana (E=200GPa)
        k_bending_th = K_local[3,3]     # Rotacional flexión (E=5GPa)
        
        # Como el carbono es 40x más rígido, la componente de flexión debe
        # reflejar la blandura relativa del núcleo.
        assert k_bending_th < k_membrane_axial, \
            "Rigidez a flexión excesiva. ¿Se usó D_membrane en vez de D_bending?"

def test_stress_recovery_robustness(element_setup):
    """
    @test Valida que el método de posproceso no falle con un vector genérico.
    """
    solver = element_setup["solver"]
    np.random.seed(42)
    u_dummy_global = np.random.rand(100) * 0.001
    
    try:
        results = solver.calculate_stresses(u_dummy_global, MockDofMapper())
    except Exception as e:
        pytest.fail(f"Stress Recovery falló con excepción: {e}")
    
    required_keys = ['sigma_x', 'sigma_y', 'tau_xy', 'von_mises', 
                     'bending_moment_x', 'shear_force_x']
    
    for key in required_keys:
        assert key in results, f"Falta la clave {key}."
        assert not np.isnan(results[key]), f"El resultado {key} es NaN."

# --- STANDALONE SPECIFIC TESTS ---

def test_constant_strain_patch_x_tension():
    """
    @test Patch Test de Membrana (CQUAD4 Isoparamétrico).
    @brief Aplica tensión uniaxial pura en X y verifica Tensiones = E * deformación.
    """
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
    # n1(0,0): u=0, v=0
    # n2(1,0): u=eps_x, v=0
    u_global[6] = eps_x
    # n3(1,1): u=eps_x, v=eps_y
    u_global[12] = eps_x
    u_global[13] = eps_y
    # n4(0,1): u=0, v=eps_y
    u_global[19] = eps_y
    
    res = elem.calculate_stresses(u_global, MockDofMapper())
    
    expected_sig_x = mat.E * eps_x
    
    assert res['sigma_x'] == pytest.approx(expected_sig_x, rel=1e-5)
    assert res['sigma_y'] == pytest.approx(0.0, abs=1e-5)
    assert res['tau_xy'] == pytest.approx(0.0, abs=1e-5)

def test_degenerate_element_error():
    """@test Asegura que un cuadrilátero con área cero (nodos colineales) lance error."""
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [2,0,0])
    n4 = MockNode(4, [3,0,0])
    
    prop = MockProp(t=1.0)
    mat = MockMat(E=1e5, nu=0.3)
    
    with pytest.raises(ValueError, match="Degenerate CQUAD4"):
        CQUAD4NastranElementSolver(n1, n2, n3, n4, prop, mat)

def test_negative_jacobian_error():
    """@test Detecta un cuadrilátero tipo 'pajarita' (bow-tie) mal numerado."""
    n1 = MockNode(1, [0, 0, 0])
    n2 = MockNode(2, [1.5, 1, 0])
    n3 = MockNode(3, [1, 0, 0]) # Nodos cruzados
    n4 = MockNode(4, [0, 1.5, 0])
    
    prop = MockProp(t=1.0)
    mat = MockMat(E=1e5, nu=0.3)
    
    # El eje Z y X se crearán bien (diagonales no paralelas), 
    # pero el _get_jacobian fallará en la cuadratura.
    elem = CQUAD4NastranElementSolver(
        n1, n2, n3, n4, prop, mat, mat_bending=mat
    )
    
    with pytest.raises(ValueError, match="Negative or zero Jacobian"):
        elem.get_stiffness_matrix()
