"""Validación final CROD SOL 400 contra valores explícitos de MSC Nastran.

El test no lee el F06/PCH: las referencias se transcriben para que una prueba
sea reproducible sin disponer de Nastran. Para añadir un caso, buscar primero
el último bloque ``LOAD STEP = 1.00000E+00`` de su F06 y usar estas tablas:

* ``DISPLACEMENT VECTOR`` -> ``expected_displacements``;
* ``FORCES OF SINGLE-POINT CONSTRAINT`` -> ``expected_reactions``;
* ``NONLINEAR STRESSES IN ROD ELEMENTS (CROD)`` ->
  ``expected_element_results``;
* ``GLOBAL WORKS & ENERGIES / TOT STRAIN ENERGY`` ->
  ``expected_strain_energy``.

Las claves de cada referencia son:

``case_name``
    Identificador único y legible utilizado también como ID de pytest.
``folder`` y ``filename``
    Subcarpeta de ``validation_cases`` y nombre exacto del DAT de entrada.
``description`` y ``note``
    Objetivo físico y, sólo cuando sea necesario, una advertencia importante.
``rel_tol`` y ``abs_tol``
    Única tolerancia relativa y absoluta aplicadas a todas las magnitudes.
``expected_displacements`` y ``expected_reactions``
    ``{grid_id: {dof_index: value}}``. El índice es cero basado:
    0=T1, 1=T2, 2=T3, 3=R1, 4=R2 y 5=R3.
``expected_element_results``
    ``{element_id: {result_name: value}}``. Los nombres disponibles son
    ``axial_stress`` (AXIAL STRESS), ``axial_strain`` (TOTAL STRAIN),
    ``torsional_stress`` (LIN. TORSIONAL STRESS), ``axial_force`` y
    ``torque``. Si el F06 no imprime las dos últimas directamente, usar
    ``axial_force = area * axial_stress`` y
    ``torque = J/C * torsional_stress`` con la PROD del elemento.
``expected_strain_energy``
    Valor ``TOT STRAIN ENERGY``. Sólo se compara en trayectorias donde ese
    acumulado representa la energía elástica del estado final.
``expected_increments``
    Número de incrementos convergidos, contando los mensajes
    ``JOB CONVERGES FOR THE CURRENT STEP``; puede diferir del FIXED solicitado
    cuando se usa ADAPT o hay cutbacks.
``expected_load_factors``
    Secuencia opcional de factores convergidos para auditar un control ADAPT.

Las entradas históricas 92--122 conservan internamente los alias compactos
``case/u/r/s/energy/increments``. Para cualquier caso nuevo debe usarse el
formato descriptivo documentado arriba y ejemplificado en 123--127.
"""

from pathlib import Path

import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"

# Ejemplo del caso 112, siempre tomando el bloque final LOAD STEP = 1.0:
#
# * DISPLACEMENT VECTOR, GRID 3, T2 = -1.538383E+01
#     -> "u": {3: {1: -15.38383}}
# * FORCES OF SINGLE-POINT CONSTRAINT, GRID 1, T1/T2 = 1.181807E+06/1.E+06
#     -> "r": {1: {0: 1.181807e6, 1: 1.0e6}}
# * NONLINEAR STRESSES..., CROD 1: AXIAL STRESS = -1.548118E+04 y
#   TOTAL STRAIN = -7.371989E-02; con A=100, N=A*sigma=-1.548118E+06
#     -> "s": {1: {"axial_force": -1.548118e6,
#                    "axial_stress": -15481.18,
#                    "axial_strain": -0.07371989}}
#
# Formato recomendado para introducir un nuevo caso:
# {
#     "case_name": "case_NNN_short_unique_name",
#     "folder": "case_NNN",
#     "filename": "ejNNN_model.dat",
#     "description": "Qué fenómeno físico y resultados valida.",
#     "note": "Sólo para una salvedad importante; puede omitirse.",
#     "rel_tol": 1.0e-5,
#     "abs_tol": 1.0e-5,
#     "expected_displacements": {grid_id: {dof_index: value}},
#     "expected_reactions": {grid_id: {dof_index: value}},
#     "expected_element_results": {
#         element_id: {"axial_stress": value, "axial_strain": value}
#     },
#     "expected_strain_energy": value,  # Opcional.
#     "expected_increments": integer,   # Opcional.
#     "expected_load_factors": [value, ...],  # Opcional para ADAPT.
# }
NONLINEAR_VALIDATION_DATA = [
    {"case": 92, "u": {2: {0: 1.0}}, "r": {1: {0: -2.1e5}}},
    {"case": 93, "u": {2: {0: 1.0}}, "s": {1: {"axial_stress": 2100.0, "axial_strain": 0.01}}},
    {"case": 94, "u": {2: {0: 20.0}}, "r": {1: {0: -4.2e6}}, "s": {1: {"axial_force": 4.2e6, "axial_stress": 42000.0, "axial_strain": 0.2}}, "energy": 4.2e7, "increments": 20},
    {"case": 95, "u": {2: {0: 20.0}}, "s": {1: {"axial_stress": 42000.0}}, "increments": 9},
    {"case": 96, "u": {2: {0: -10.0}}, "r": {1: {0: 2.1e6}}, "s": {1: {"axial_stress": -21000.0, "axial_strain": -0.1}}},
    {"case": 97, "u": {2: {0: 12.24745}}, "r": {1: {0: -2.1e6}}, "s": {1: {"axial_stress": 21000.0, "axial_strain": 0.1}}},
    {"case": 98, "u": {2: {0: 10.0}, 3: {0: 25.0}}, "s": {1: {"axial_stress": 21000.0}, 2: {"axial_stress": 21000.0}}},
    {"case": 99, "u": {2: {0: 2.666667}, 3: {0: 20.66667}}, "r": {1: {0: -8.4e5}}, "s": {1: {"axial_stress": 7000.0}, 2: {"axial_stress": 10500.0}}},
    {"case": 100, "u": {2: {0: 10.0}}, "r": {1: {0: -2.45e6}}, "s": {1: {"axial_stress": 21000.0}, 2: {"axial_stress": 7000.0}}},
    {"case": 101, "u": {2: {0: 8.095238}}, "r": {1: {0: -1.7e6}}, "s": {1: {"axial_stress": 17000.0, "axial_strain": 0.08095238}}},
    {"case": 102, "u": {3: {1: -223.3105}}, "r": {1: {0: -1.621922e6, 1: 2.0e6}, 2: {0: 1.621922e6, 1: 2.0e6}}, "s": {1: {"axial_stress": 25750.01, "axial_strain": 0.1226191}, 2: {"axial_stress": 25750.01}}},
    {"case": 103, "u": {3: {0: 9.454499, 1: -39.27897}}, "r": {1: {0: 1.608725e6, 1: 7.545773e5}, 2: {0: -2.608725e6, 1: 2.245423e6}}, "s": {1: {"axial_stress": -17769.03}, 2: {"axial_stress": -57366.68}}},
    {"case": 104, "u": {3: {1: -6.221705}}, "r": {1: {0: 266264.9, 1: 50000.0}}, "s": {1: {"axial_stress": -2709.188, "axial_strain": -0.01290089}}, "increments": 14},
    {"case": 105, "u": {2: {0: -13.27239}, 3: {0: -6.636195, 1: -233.9893}}, "r": {1: {1: 2.0e6}, 2: {1: 2.0e6}}, "s": {1: {"axial_stress": -13936.01}, 2: {"axial_stress": 32501.96}, 3: {"axial_stress": 32501.96}}},
    {"case": 106, "u": {3: {0: 22.52493, 1: -12.27022}, 4: {0: 18.92084, 1: -4.190009}}, "r": {1: {0: -576573.9}, 2: {0: -423426.1, 1: 2.998984e6}}, "s": {3: {"axial_stress": 8229.271}, 6: {"axial_stress": -23623.47}}},
    {"case": 107, "u": {4: {0: 10.33237, 1: -5.384856, 2: -9.4136}}, "r": {2: {0: -797209.0, 1: 475349.0, 2: 1.453991e6}}, "s": {2: {"axial_stress": -24642.71, "axial_strain": -0.1173462}}},
    {"case": 108, "u": {5: {0: 10.92065, 1: 1.751269, 2: -21.62411}}, "r": {3: {0: -1.318034e6, 1: -1.492986e6, 2: 1.877012e6}}, "s": {3: {"axial_stress": -27366.76, "axial_strain": -0.1303179}}},
    {"case": 109, "u": {3: {1: -6.221705}}, "s": {1: {"axial_stress": -2709.188}}, "increments": 5},
    {"case": 110, "u": {3: {1: -6.221705}}, "s": {1: {"axial_stress": -2709.188}}, "increments": 40},
    {"case": 111, "u": {5: {0: 2.536503, 1: 0.03771908, 2: -5.326263}, 6: {0: -7.533348, 1: 44.10162, 2: -17.23535}}, "r": {3: {0: -870512.1, 1: -459667.5, 2: 1.442809e6}}, "s": {5: {"axial_stress": 120.2893}, 9: {"axial_stress": -22619.5}}, "increments": 15},
    {"case": 112, "u": {3: {1: -15.38383}}, "r": {1: {0: 1.181807e6, 1: 1.0e6}}, "s": {1: {"axial_force": -1.548118e6, "axial_stress": -15481.18, "axial_strain": -0.07371989}}},
    {"case": 113, "u": {3: {0: 0.0, 1: -15.38383, 2: 0.0}}, "r": {1: {0: 585445.7, 1: 1.189084e6, 2: 8.0e5}}, "s": {1: {"axial_stress": -15481.18}}},
    {"case": 114, "u": {3: {0: 0.006127186, 1: -15.38383}}, "s": {1: {"axial_stress": -15474.23}, 2: {"axial_stress": -15488.12}}},
    {"case": 115, "u": {3: {0: -0.006127186, 1: -15.38383}}, "s": {1: {"axial_stress": -15488.12}, 2: {"axial_stress": -15474.23}}},
    {"case": 116, "u": {2: {3: 0.9904762}}, "r": {1: {3: -40000.0}}, "s": {1: {"torque": 40000.0, "torsional_stress": 2000.0}}, "energy": 19809.52},
    {"case": 117, "u": {2: {0: 20.0, 3: 0.9904762}}, "r": {1: {0: -4.2e6, 3: -40000.0}}, "s": {1: {"axial_stress": 42000.0, "torque": 40000.0, "torsional_stress": 2000.0}}, "energy": 4.201981e7},
    {"case": 118, "u": {3: {1: -15.38383}}, "r": {2: {0: -1.181807e6, 1: 1.0e6}}, "s": {2: {"axial_stress": -15481.18}}},
    {"case": 119, "u": {3: {1: -40.32496}}, "r": {1: {0: 3.183911e6, 1: 1.9e6}}, "s": {1: {"axial_stress": -37077.34, "axial_strain": -0.1765587}}},
    {"case": 121, "u": {3: {1: -0.269374}}, "r": {1: {0: 20000.0, 1: 20000.0}}},
    {"case": 122, "u": {3: {1: -0.2699212}}, "r": {1: {0: 20054.13, 1: 20000.0}}, "s": {1: {"axial_stress": -283.2257, "axial_strain": -0.001348694}}},
    {
        "case_name": "case_123_planar_arch_truss_10_crods",
        "folder": "case_123",
        "filename": "ej123_planar_arch_truss_10_crods.dat",
        "description": (
            "Planar roof arch with 10 CRODs, three PROD definitions and "
            "multiple horizontal and vertical nodal loads."
        ),
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            3: {0: 21.43302, 1: -28.98271},
            4: {0: 2.384703, 1: -9.052676},
            5: {0: 2.718968, 1: -11.89651},
            6: {0: 14.18247, 1: 8.534641},
        },
        "expected_reactions": {
            1: {0: 1.592394e6, 1: 1.035651e6},
            2: {0: -1.892394e6, 1: 1.264349e6},
        },
        "expected_element_results": {
            1: {"axial_stress": -2313.415, "axial_strain": -0.01101626},
            6: {"axial_stress": -5344.242, "axial_strain": -0.02544877},
            7: {"axial_stress": -984.4266, "axial_strain": -0.01406324},
            10: {"axial_stress": -5916.356, "axial_strain": -0.02817313},
        },
        "expected_increments": 25,
    },
    {
        "case_name": "case_124_planar_bridge_enforced_support",
        "folder": "case_124",
        "filename": "ej124_planar_bridge_20_crods_enforced_support.dat",
        "description": (
            "Planar 20-CROD bridge with five PROD definitions, distributed "
            "loads and an imposed T1=-4 support displacement."
        ),
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            2: {0: -1.395999, 1: -8.651058},
            3: {0: -1.740947, 1: -12.40334},
            5: {0: -4.0},
            7: {0: 1.837685, 1: -8.662676},
            8: {0: -0.2736143, 1: -13.53867},
            9: {0: -3.180125, 1: -5.145882},
            10: {0: -3.884803, 1: -1.100777},
        },
        "expected_reactions": {
            1: {0: 1.803893e6, 1: 1.327292e6},
            5: {0: -1.903893e6, 1: 1.322708e6},
        },
        "expected_element_results": {
            1: {"axial_stress": -2136.171, "axial_strain": -0.01017224},
            5: {"axial_stress": -2027.556, "axial_strain": -0.009655028},
            9: {"axial_stress": -2996.685, "axial_strain": -0.01426993},
            14: {"axial_stress": -2288.936, "axial_strain": -0.03269909},
            18: {"axial_stress": 3202.249, "axial_strain": 0.0152488},
        },
        "expected_increments": 30,
    },
    {
        "case_name": "case_125_space_dome_snap_through",
        "folder": "case_125",
        "filename": "ej125_space_dome_30_crods.dat",
        "description": (
            "Shallow 30-CROD space dome with mixed materials, four supports "
            "and a geometrically nonlinear snap-through response."
        ),
        "note": (
            "MSC and the project match the final mechanical state, but MSC "
            "TotStrEn accumulates work through the displacement jump between "
            "load factors 0.533333 and 0.566667; it is not asserted as final "
            "elastic strain energy."
        ),
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            2: {0: -0.07741382, 1: 1.757908, 2: -2.591646},
            4: {0: 3.487501, 1: -0.4843325, 2: -1.679655},
            5: {0: 0.9462261, 1: -0.4191779, 2: -73.46688},
            6: {0: -1.614412, 1: -0.2643717, 2: -1.642472},
            8: {0: 1.907263, 1: -2.523066, 2: -2.781328},
            9: {0: 3.230258, 1: 0.8163792, 2: 0.2546659},
        },
        "expected_reactions": {
            10: {0: 265753.3, 1: 199545.8, 2: 342243.1},
            11: {0: -248069.9, 1: 180473.6, 2: 581075.1},
            12: {0: -406233.7, 1: -151868.7, 2: 474126.8},
            13: {0: 138550.3, 1: -128150.7, 2: 352554.9},
        },
        "expected_element_results": {
            1: {"axial_stress": 1685.539, "axial_strain": 0.008026378},
            13: {"axial_stress": -1942.757, "axial_strain": -0.02775367},
            19: {"axial_stress": 530.2653, "axial_strain": 0.002525073},
            23: {"axial_stress": -870.1577, "axial_strain": -0.01243082},
            30: {"axial_stress": -450.2757, "axial_strain": -0.00643251},
        },
        "expected_increments": 30,
    },
    {
        "case_name": "case_126_space_canopy_enforced_sway",
        "folder": "case_126",
        "filename": "ej126_space_canopy_24_crods_enforced_sway.dat",
        "description": (
            "Three-dimensional 24-CROD canopy with mixed properties, "
            "multiple loads and an imposed T1=+5 support displacement."
        ),
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            3: {0: 5.0},
            5: {0: 3.811938, 1: -0.1556497, 2: -1.091279},
            6: {0: 3.555248, 1: 2.676464, 2: -1.686319},
            7: {0: 2.930649, 1: 2.957228, 2: -0.6619032},
            8: {0: 2.57219, 1: 0.2532403, 2: -1.820616},
            9: {0: 2.358079, 1: 1.624651, 2: -6.567131},
        },
        "expected_reactions": {
            1: {0: -1.341084e6, 1: -822321.7, 2: 79257.18},
            2: {0: -68405.09, 1: -108987.0, 2: 617911.6},
            3: {0: 3.401829e6, 1: 843628.0, 2: 264624.0},
            4: {0: -2.062339e6, 1: -132319.3, 2: 738207.2},
        },
        "expected_element_results": {
            1: {"axial_stress": 0.0, "axial_strain": 0.0},
            7: {"axial_stress": -740.2198, "axial_strain": -0.003524856},
            11: {"axial_stress": 580.8184, "axial_strain": 0.008297405},
            15: {"axial_stress": -278.0962, "axial_strain": -0.001324268},
            21: {"axial_stress": -1526.071, "axial_strain": -0.02180101},
        },
        "expected_increments": 40,
    },
    {
        "case_name": "case_127_braced_space_tower_adaptive",
        "folder": "case_127",
        "filename": "ej127_braced_space_tower_30_crods_adaptive.dat",
        "description": (
            "Braced 30-CROD space tower with six PROD definitions, mixed "
            "materials, three-dimensional loading and adaptive stepping."
        ),
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            5: {0: 8.056498, 1: 1.093566, 2: -0.7330719},
            6: {0: 6.531606, 1: 7.155821, 2: -2.133779},
            7: {0: 1.863507, 1: 6.093887, 2: -1.432072},
            8: {0: 1.568026, 1: 1.020863, 2: -0.1321723},
            9: {0: 18.08114, 1: 1.461072, 2: -3.417495},
            10: {0: 17.40426, 1: 13.26354, 2: -4.560487},
            11: {0: 4.472591, 1: 11.79553, 2: -2.282022},
            12: {0: 5.401279, 1: 0.3733238, 2: -0.9814105},
        },
        "expected_reactions": {
            1: {0: -683320.4, 1: -171353.9, 2: -1.005912e6},
            2: {0: -76321.9, 1: -220847.0, 2: 899246.6},
            3: {0: -403673.1, 1: -198680.5, 2: 1.933434e6},
            4: {0: 13315.33, 1: -59118.61, 2: 73231.45},
        },
        "expected_element_results": {
            1: {"axial_stress": 0.0, "axial_strain": 0.0},
            5: {"axial_stress": -1841.171, "axial_strain": -0.008767482},
            13: {"axial_stress": 1167.11, "axial_strain": 0.005557666},
            21: {"axial_stress": 1414.152, "axial_strain": 0.02020217},
            25: {"axial_stress": 4841.163, "axial_strain": 0.02305316},
            29: {"axial_stress": -620.335, "axial_strain": -0.008861929},
        },
        "expected_increments": 14,
        "expected_load_factors": [
            0.025,
            0.055,
            0.091,
            0.1342,
            0.18604,
            0.248248,
            0.322898,
            0.412477,
            0.512477,
            0.612477,
            0.712477,
            0.812477,
            0.912477,
            1.0,
        ],
    },
]


def _validation_id(reference):
    return reference.get("case_name") or f"case_{reference['case']}"


@pytest.mark.parametrize(
    "reference", NONLINEAR_VALIDATION_DATA, ids=_validation_id
)
def test_nonlinear_crod_final_results(reference):
    case_name = _validation_id(reference)
    if "folder" in reference and "filename" in reference:
        deck = VALIDATION_ROOT / reference["folder"] / reference["filename"]
        assert deck.is_file(), f"Missing DAT for {case_name}: {deck}"
    else:
        deck = next((VALIDATION_ROOT / f"case_{reference['case']}").glob("*.dat"))
    model = parse_nastran_input(str(deck))
    solver = NonlinearStaticSolver(model)
    result = solver.solve()
    rel_tol = reference.get("rel_tol", 3.0e-5)
    abs_tol = reference.get("abs_tol", 1.0e-5)

    assert result["load_factor"] == pytest.approx(1.0)
    expected_increments = reference.get(
        "expected_increments", reference.get("increments")
    )
    if expected_increments is not None:
        assert len(result["increments"]) == expected_increments
    if "expected_load_factors" in reference:
        computed_factors = [
            increment.load_factor for increment in result["increments"]
        ]
        assert computed_factors == pytest.approx(
            reference["expected_load_factors"], rel=rel_tol, abs=abs_tol
        )

    expected_displacements = reference.get(
        "expected_displacements", reference.get("u", {})
    )
    for node_id, dofs in expected_displacements.items():
        nodal = result["displacements"][solver.dof_mapper.get_dof_indices(node_id)]
        for dof, expected in dofs.items():
            assert nodal[dof] == pytest.approx(expected, rel=rel_tol, abs=abs_tol)

    expected_reactions = reference.get("expected_reactions", reference.get("r", {}))
    for node_id, dofs in expected_reactions.items():
        nodal = result["reactions"][node_id]
        for dof, expected in dofs.items():
            assert nodal[dof] == pytest.approx(expected, rel=rel_tol, abs=abs_tol)

    expected_element_results = reference.get(
        "expected_element_results", reference.get("s", {})
    )
    for element_id, expected_results in expected_element_results.items():
        computed = result["stresses"][element_id]
        for key, expected in expected_results.items():
            assert computed[key] == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            )

    expected_energy = reference.get(
        "expected_strain_energy", reference.get("energy")
    )
    if expected_energy is not None:
        assert result["strain_energy"] == pytest.approx(
            expected_energy, rel=rel_tol, abs=abs_tol
        )
