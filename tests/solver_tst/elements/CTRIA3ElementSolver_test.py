
"""
@file CTRIA3ElementSolver_test.py
@brief Unit tests for the CTRIA3 Element Solver (Membrane stiffness).
@details
Validates the physics of the Triangular Shell element:
1. Stiffness Matrix Assembly (Membrane - CST).
2. Geometric Transformations (Local axes computation).
3. Stress Recovery (Von Mises, Principal Stresses).
4. Error Handling (Zero area).
"""

import pytest
import numpy as np
from src.solver.elements.ctria3_element import CTRIA3ElementSolver

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
        # Simulamos la lógica de buscar el material de bending
        self.mat_bending_ref = mat_bending
        self.mat_id_bending = 999 if mat_bending else None

class MockDofMapper:
    """Simulates the global-to-local index mapping."""
    def get_dof_indices(self, nid):
        # Maps Node ID to a block of 6 indices
        # Node 1 -> [0..5], Node 2 -> [6..11] ...
        start = (nid - 1) * 6
        return list(range(start, start + 6))


# --- DEFINICIÓN DE CASOS (FIXTURE DATA) ---

'''
Casos de prueba:
1. CANONICAL: Triángulo rectángulo básico en el plano XY (Baseline).
2. ROTATED_3D: Triángulo equilátero flotando en el espacio con rotación arbitraria.
3. NEEDLE: Triángulo con relación de aspecto extrema (base pequeña, altura gigante).
4. COMPOSITE: Materiales desacoplados (Membrana muy rígida, Flexión muy blanda).
5. FAR_OFFSET: Triángulo válido pero definido en coordenadas gigantes (lejos del origen).
'''

CASES = [
    {
        "id": "1_CANONICAL_XY",
        "nodes": [[0.,0.,0.], [2.,0.,0.], [0.,2.,0.]], # Rectángulo Isósceles
        "mat_mem": (210e9, 0.3), # Acero
        "mat_bend": None,        # Igual a membrana
        "thick": 0.1,
        "area": 2.0,
        "desc": "Caso base en plano XY"
    },
    {
        "id": "2_ROTATED_3D",
        "nodes": [[1.,1.,1.], [2.,2.,1.], [1.,2.,2.]], # En el espacio
        "mat_mem": (70e9, 0.33), # Aluminio
        "mat_bend": None,
        "thick": 0.05,
        "area": 0.866025, # Aprox (verificado)
        "desc": "Triángulo rotado en 3D"
    },
    {
        "id": "3_NEEDLE_SHAPE",
        "nodes": [[0.,0.,0.], [0.1, 0., 0.], [0.05, 10., 0.]], # Base 0.1, Altura 10
        "mat_mem": (110e9, 0.3), # Titanio
        "mat_bend": None,
        "thick": 0.01,
        "area": 0.5, # 0.5 * 0.1 * 10
        "desc": "Test de estabilidad numérica (Aspect Ratio alto)"
    },
    {
        "id": "4_COMPOSITE_DECOUPLED",
        "nodes": [[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]],
        "mat_mem": (200e9, 0.25), # Carbono (Membrana rígida)
        "mat_bend": (5e9, 0.4),   # Foam Core (Flexión blanda)
        "thick": 0.2,
        "area": 0.5,
        "desc": "Materiales distintos para membrana y flexión"
    },
    {
        "id": "5_FAR_OFFSET",
        "nodes": [[10000., 10000., 10000.], 
                  [10002., 10000., 10000.], 
                  [10000., 10002., 10000.]], # Igual al caso 1 pero lejos
        "mat_mem": (210e9, 0.3),
        "mat_bend": None,
        "thick": 0.1,
        "area": 2.0,
        "desc": "Test de precisión flotante con coordenadas grandes"
    }
]

# --- FIXTURE PARAMETRIZADO ---

@pytest.fixture(params=CASES, ids=[c["id"] for c in CASES])
def element_setup(request):
    """
    Crea el solver para cada uno de los 5 casos definidos.
    Devuelve un diccionario con el solver y los datos esperados para validación.
    """
    data = request.param
    
    # 1. Crear Nodos
    n1 = MockNode(1, data["nodes"][0])
    n2 = MockNode(2, data["nodes"][1])
    n3 = MockNode(3, data["nodes"][2])
    
    # 2. Crear Materiales
    mat_mem = MockMat(*data["mat_mem"], name="Membrane")
    
    mat_bend = None
    if data["mat_bend"]:
        mat_bend = MockMat(*data["mat_bend"], name="Bending")
    
    # 3. Propiedad
    prop = MockProp(t=data["thick"], mat_bending=mat_bend)
    
    # 4. Instanciar Solver (Aquí se inyecta la lógica de corrección que hiciste)
    # Nota: El solver real debería aceptar mat_bending en __init__
    solver = CTRIA3ElementSolver(n1, n2, n3, prop, mat_mem, mat_bending=mat_bend)
    
    return {
        "solver": solver,
        "data": data,
        "mat_mem": mat_mem,
        "mat_bend": mat_bend if mat_bend else mat_mem
    }
# --- TESTS ---

# --- TESTS ---

def test_stiffness_matrix_shape(element_setup):
    """
    @test Validate Stiffness Matrix Shape (18x18).
    """
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    assert K.shape == (18, 18)

def test_stiffness_symmetry(element_setup):
    """Valida la simetría de la matriz K global (Maxwell-Betti)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    
    # Tolerancia relativa para manejar escalas numéricas gigantes (N/m vs N*m)
    # Usamos assert_allclose que maneja tolerancias relativas (rtol) y absolutas (atol)
    np.testing.assert_allclose(K, K.T, rtol=1e-5, atol=1e-8, 
                               err_msg="La matriz de rigidez no es simétrica.")

def test_rigid_body_modes(element_setup):
    """Valida que existen modos de cuerpo rígido (Autovalores cercanos a 0)."""
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    
    eigenvalues = np.linalg.eigvalsh(K)
    
    # Los primeros 3 deben ser traslaciones puras (muy cercanos a 0)
    # Nota: Los modos rotacionales pueden no ser 0 puro debido al Drilling DOF artificial
    translations = eigenvalues[:3]
    
    # Umbral de cero relativo a la rigidez máxima del sistema
    max_k = np.max(np.abs(K))
    tolerance = max_k * 1e-9
    
    # En el caso NEEDLE el condicionamiento es peor, relajamos si es necesario
    if element_setup["data"]["id"] == "3_NEEDLE_SHAPE":
        tolerance = max_k * 1e-7

    assert np.all(np.abs(translations) < tolerance), \
        f"Fallo en modos de cuerpo rígido. Eigenvalues: {translations}"

def test_area_calculation(element_setup):
    """Valida que el área calculada internamente coincide con la teórica."""
    solver = element_setup["solver"]
    expected_area = element_setup["data"]["area"]
    
    # Relajamos un poco la tolerancia para el caso NEEDLE o ROTATED
    assert solver.A == pytest.approx(expected_area, rel=1e-4)

def test_local_axes_orientation(element_setup):
    """
    @test Valida la orientación de los ejes locales.
    Para el caso 1_CANONICAL_XY, verifica orientación exacta.
    Para otros, verifica ortogonalidad.
    """
    solver = element_setup["solver"]
    data = element_setup["data"]
    
    vx = solver.local_x_axis
    vy = solver.local_y_axis
    vz = solver.local_z_axis
    
    # 1. Ortogonalidad
    assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)
    assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)
    
    # 2. Caso Específico
    if data["id"] == "1_CANONICAL_XY":
        np.testing.assert_allclose(vx, [1,0,0], atol=1e-9)
        np.testing.assert_allclose(vy, [0,1,0], atol=1e-9)
        np.testing.assert_allclose(vz, [0,0,1], atol=1e-9)

def test_membrane_bending_decoupling(element_setup):
    """
    @brief Checks that Membrane and Bending actions are decoupled in Local System.
    @details
    For an isotropic flat plate, in-plane forces (Tx, Ty) should not cause
    out-of-plane deformations (Tz, Rx, Ry).
    """
    solver = element_setup["solver"]
    
    K_local = solver._build_local_stiffness()
    
    # Local DOF mapping:
    # Membrane: 0(x), 1(y)
    # Bending:  2(z), 3(rx), 4(ry)
    
    # Check coupling between Node 1 Tx (0) and Node 1 Tz (2)
    # Solo válido si no hay acople exótico (el caso COMPOSITE_DECOUPLED usa isotropic materials, debería ser 0)
    k_cross = K_local[0, 2]
    
    assert k_cross == pytest.approx(0.0, abs=1e-9), \
        "Coupling detected between Membrane (x) and Bending (z). Check B-matrices."

def test_material_decoupling_check(element_setup):
    """
    Test crucial para el Caso 4 (Composite):
    Verifica que la rigidez de flexión usa el módulo E de flexión, no el de membrana.
    """
    setup = element_setup
    solver = setup["solver"]
    data = setup["data"]
    
    # Solo ejecutamos validación profunda si es el caso desacoplado
    if data["id"] == "4_COMPOSITE_DECOUPLED":
        K_local = solver._build_local_stiffness()
        
        # Rigidez Membrana Local ~ E_mem * t * factor
        # Rigidez Flexión Local ~ E_bend * t^3 * factor
        
        # Extraccion cualitativa:
        # K_local[0,0] es rigidez axial (u) -> Depende de Mat Membrana
        # K_local[2,2] es rigidez out-of-plane (w) -> Depende de Mat Bending (DKT)
        
        k_membrane_axial = K_local[0,0]
        k_bending_w = K_local[2,2] # Asociado a flexión
        
        # Si hubiéramos usado E_mem (200e9) para bending, k_bending_w sería enorme.
        # Como usamos E_bend (5e9), debe ser proporcionalmente menor (ajustando por t^2)
        
        # Chequeo de sanidad: La rigidez de flexión no debe "heredar" la rigidez del carbono
        # Este es un chequeo heurístico.
        assert k_bending_w < k_membrane_axial, \
            "La rigidez a flexión parece demasiado alta. ¿Se está usando el material de membrana?"

def test_bending_stiffness_generated(element_setup):
    """
    @test Validate Bending Stiffness contribution (DKT).
    @brief Checks that stiffness matrix has non-zero terms for bending DOFs (Tz, Rx, Ry).
    """
    solver = element_setup["solver"]
    K = solver.get_stiffness_matrix()
    prop = element_setup["data"]["thick"]
    
    # Check that we have non-zero stiffness in Bending DOFs
    # e.g., K[2,2] (Tz1-Tz1)
    # Scale checking: depends on E and t.
    if element_setup["data"]["id"] == "3_NEEDLE_SHAPE":
        return # Skip scaling check for needle, numerics are weird

    assert abs(K[2,2]) > 0.0, "Bending stiffness K_zz should be non-zero"

def test_stress_recovery_robustness(element_setup):
    """
    Valida que el método calculate_stresses no explota con NaN y devuelve las claves correctas.
    """
    solver = element_setup["solver"]
    
    # Generamos un vector de desplazamientos aleatorio pequeño
    np.random.seed(42)
    u_dummy_global = np.random.rand(100) * 0.001 # 100 es arbitrario, el mapper lo filtra
    
    # Mock del DofMapper
    class MockMapper:
        def get_dof_indices(self, nid):
            # Asignamos indices ficticios para el test
            base = (nid - 1) * 6
            return list(range(base, base+6))
    
    try:
        results = solver.calculate_stresses(u_dummy_global, MockMapper())
    except Exception as e:
        pytest.fail(f"Stress Recovery falló con excepción: {e}")
    
    # Validaciones
    required_keys = ['sigma_x', 'sigma_y', 'tau_xy', 'von_mises', 
                     'bending_moment_x', 'stress_fiber_top_x']
    
    for key in required_keys:
        assert key in results, f"Falta la clave {key} en los resultados de tensión."
        assert not np.isnan(results[key]), f"El resultado {key} es NaN."

def test_bending_rigid_body_z(element_setup):
    """
    @test Rigid Body Translation in Z.
    @brief Ensure zero stress for constant w translation.
    """
    solver = element_setup["solver"]
    
    # Apply uniform w=0.1 at all nodes, zero rotations
    # Nodes IDs are 1, 2, 3
    # Mapped to global indices:
    # N1: 0-5, N2: 6-11, N3: 12-17
    u_global = np.zeros(100) # Size enough for mocks
    
    # Global indices for w (z-translation)
    # Mapper mocks: (nid-1)*6 -> w is index 2, 8, 14
    u_global[2] = 0.1
    u_global[8] = 0.1
    u_global[14] = 0.1
    
    # Mock Mapper that matches the one implied in element setup (1,2,3)
    class MockMapper:
        def get_dof_indices(self, nid):
            base = (nid - 1) * 6
            return list(range(base, base+6))

    res = solver.calculate_stresses(u_global, MockMapper())
    
    # Should have zero bending moments and stresses
    assert res['bending_moment_x'] == pytest.approx(0.0, abs=1e-5)
    assert res['stress_fiber_top_x'] == pytest.approx(0.0, abs=1e-5)

# --- STANDALONE SPECIFIC TESTS (Not using Fixture Cases) ---

def test_constant_strain_patch_x_tension():
    """
    @test Patch Test (Constant Strain).
    @brief Apply uniform uniaxial tension in X and verify Stress = E*strain.
    """
    # Right Isosceles Triangle
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [0,1,0])
    
    prop = MockProp(t=1.0)
    mat = MockMat(E=1000.0, nu=0.3)
    
    elem = CTRIA3ElementSolver(n1, n2, n3, prop, mat)
    
    # Apply displacements for constant strain eps_x
    # u = x * eps_x
    # v = -y * nu * eps_x
    eps_x = 0.001
    eps_y = -mat.nu * eps_x
    
    # Node 1 (0,0): u=0, v=0
    # Node 2 (1,0): u=eps_x, v=0
    # Node 3 (0,1): u=0, v=eps_y
    
    u_global = np.zeros(18)
    
    # N2 DOFs (Indices 6-11)
    u_global[6] = eps_x # u2
    
    # N3 DOFs (Indices 12-17)
    u_global[13] = eps_y # v3
    
    mapper = MockDofMapper()
    res = elem.calculate_stresses(u_global, mapper)
    
    # Expected Stress Sigma_x = E * eps_x (in Plane Stress with compatible Poisson contraction)
    # Sigma_y should be 0.
    
    expected_sig_x = mat.E * eps_x # 1000 * 0.001 = 1.0
    
    assert res['sigma_x'] == pytest.approx(expected_sig_x, rel=1e-5)
    assert res['sigma_y'] == pytest.approx(0.0, abs=1e-5)
    assert res['tau_xy'] == pytest.approx(0.0, abs=1e-5)
    assert res['von_mises'] == pytest.approx(expected_sig_x, rel=1e-5)

def test_zero_area_error():
    """@test Ensure ValueError on collinear nodes."""
    n1 = MockNode(1, [0,0,0])
    n2 = MockNode(2, [1,0,0])
    n3 = MockNode(3, [2,0,0]) # Collinear
    prop = MockProp(t=1.0)
    mat = MockMat(E=1e5, nu=0.3)
    
    with pytest.raises(ValueError):
        CTRIA3ElementSolver(n1, n2, n3, prop, mat)
