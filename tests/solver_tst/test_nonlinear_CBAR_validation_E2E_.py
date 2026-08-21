"""Validación final CBAR SOL 400 contra referencias explícitas de MSC Nastran.

El test no lee los F06/PCH/OP2: los resultados se transcriben para que sea
reproducible sin una instalación de Nastran. Para añadir un caso se utiliza el
último bloque ``LOAD STEP = 1.00000E+00`` del F06:

* ``DISPLACEMENT VECTOR`` -> ``expected_displacements``;
* ``FORCES OF SINGLE-POINT CONSTRAINT`` -> ``expected_reactions``;
* ``NONLINEAR STRESSES IN BAR ELEMENTS (CBAR)`` -> resultados de elemento;
* ``GLOBAL WORKS & ENERGIES / TOT STRAIN ENERGY`` -> energía, sólo cuando el
  estado permite compararla directamente con la energía elástica recuperada.

Claves del diccionario de cada caso:

``case_name``
    ID único que también aparece como identificador de pytest.
``folder`` y ``filename``
    Subcarpeta dentro de ``validation_cases`` y DAT exacto que se resuelve.
``description`` y ``note``
    Fenómeno validado y, opcionalmente, una salvedad relevante.
``rel_tol`` y ``abs_tol``
    Una tolerancia relativa y otra absoluta para todas las magnitudes.
``expected_displacements`` y ``expected_reactions``
    ``{grid_id: {dof_index: value}}``; 0=T1, 1=T2, 2=T3, 3=R1, 4=R2,
    5=R3. Sólo es necesario transcribir los GDL relevantes para el caso.
``expected_element_results``
    ``{element_id: {result_name: value}}``. Entre las claves disponibles están
    ``axial_stress``, ``axial_strain``, ``axial_force``, ``torque``,
    ``twist_rate``, ``max_element_stress`` y ``min_element_stress``. Para una
    sección simétrica, ``axial_stress`` es la media de C/D/E/F de una estación;
    con puntos no simétricos debe obtenerse como la ordenada en el centroide
    del plano de tensiones, no como la media aritmética.
``expected_strain_energy`` y ``expected_increments``
    Referencias opcionales de energía y número de incrementos convergidos.

Ejemplo de alta rápida de una nueva referencia::

    {
        "case_name": "case_NNN_short_name",
        "folder": "case_NNN",
        "filename": "ejNNN_model.dat",
        "description": "Physical behavior validated by the model.",
        "rel_tol": 1.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {1: 12.345, 5: 0.0123}},
        "expected_reactions": {1: {1: -100000.0}},
        "expected_element_results": {1: {"axial_stress": 125.0}},
        "expected_increments": 20,
    }
"""

from pathlib import Path

import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"

NONLINEAR_CBAR_VALIDATION_DATA = [
    {
        "case_name": "case_128_small_axial_linear_limit",
        "folder": "case_128",
        "filename": "ej128_cbar_small_axial_sol400_baseline.dat",
        "description": "SOL 400 CBAR axial response with LGDISP=-1.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: 0.4761905}},
        "expected_reactions": {1: {0: -1.0e5}},
        "expected_element_results": {1: {"axial_stress": 100.0}},
        "expected_strain_energy": 23809.52,
    },
    {
        "case_name": "case_129_small_bending_linear_limit",
        "folder": "case_129",
        "filename": "ej129_cbar_small_bending_sol400_baseline.dat",
        "description": "SOL 400 bending signs and stiffness with LGDISP=-1.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {1: 1.984127, 5: 0.002976190}},
        "expected_reactions": {1: {1: -1.0e4, 5: -1.0e7}},
        "expected_element_results": {1: {"axial_stress": 0.0}},
        "expected_strain_energy": 9920.635,
    },
    {
        "case_name": "case_130_finite_axial_tension",
        "folder": "case_130",
        "filename": "ej130_cbar_large_axial_tension_fixed.dat",
        "description": "Finite axial stretch and logarithmic section strain.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: 20.20134}},
        "expected_reactions": {1: {0: -4.2e6}},
        "expected_element_results": {1: {"axial_stress": 4200.0}},
        "expected_increments": 20,
    },
    {
        "case_name": "case_131_finite_axial_compression",
        "folder": "case_131",
        "filename": "ej131_cbar_subcritical_axial_compression.dat",
        "description": "Finite subcritical compression and geometric tangent.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -1.902949}},
        "expected_reactions": {1: {0: 4.0e5}},
        "expected_element_results": {1: {"axial_stress": -400.0}},
        "expected_increments": 20,
    },
    {
        "case_name": "case_132_large_bending_local_y",
        "folder": "case_132",
        "filename": "ej132_cbar_large_bending_local_y_force.dat",
        "description": "Large I1-plane bending, rotation and P-delta coupling.",
        "rel_tol": 1.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -1.680474, 1: 59.39549, 5: 0.08913815}},
        "expected_reactions": {1: {1: -3.0e5, 5: -2.994959e8}},
        "expected_element_results": {1: {"axial_stress": 17.81715}},
        "expected_increments": 30,
    },
    {
        "case_name": "case_133_large_bending_local_z",
        "folder": "case_133",
        "filename": "ej133_cbar_large_bending_local_z_force.dat",
        "description": "Large I2-plane bending and three-dimensional signs.",
        "rel_tol": 1.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -4.437906, 2: 94.68011, 4: -0.1422255}},
        "expected_reactions": {1: {2: -1.2e5, 4: 1.194675e8}},
        "expected_element_results": {1: {"axial_stress": 11.361}},
        "expected_increments": 30,
    },
    {
        "case_name": "case_134_pure_end_moment_z",
        "folder": "case_134",
        "filename": "ej134_cbar_pure_end_moment_local_z.dat",
        "description": "Constant-curvature bending about local Z.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -1.771017, 1: 59.48864, 5: 0.1190476}},
        "expected_reactions": {1: {5: -2.0e8}},
        "expected_element_results": {1: {"axial_stress": 0.0}},
        "expected_strain_energy": 1.190476e7,
        "expected_increments": 30,
    },
    {
        "case_name": "case_135_pure_end_moment_y",
        "folder": "case_135",
        "filename": "ej135_cbar_pure_end_moment_local_y.dat",
        "description": "Constant-curvature bending about local Y.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -4.531709, 2: -95.09407, 4: 0.1904760}},
        "expected_reactions": {1: {4: -8.0e7}},
        "expected_element_results": {1: {"axial_stress": 0.0}},
        "expected_strain_energy": 7.619038e6,
        "expected_increments": 30,
    },
    {
        "case_name": "case_136_pure_torsion",
        "folder": "case_136",
        "filename": "ej136_cbar_pure_torsion.dat",
        "description": "Finite twist with the MAT1/PBAR GJ law.",
        "rel_tol": 1.0e-5,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {3: 1.114286}},
        "expected_reactions": {1: {3: -9.0e7}},
        "expected_element_results": {1: {"axial_stress": 0.0, "torque": 9.0e7}},
        "expected_strain_energy": 5.014286e7,
        "expected_increments": 24,
    },
    {
        "case_name": "case_137_combined_spatial_base",
        "folder": "case_137",
        "filename": "ej137_cbar_combined_load_base.dat",
        "description": "Combined finite axial, biaxial and torsional loading.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: 0.747466, 1: 53.21096, 2: -39.02246, 3: 0.2197431, 4: 0.05488243, 5: 0.08399875}},
        "expected_reactions": {1: {0: -6.0e5, 1: -2.0e5, 2: 8.0e4, 3: -3.354761e7, 4: -3.664632e7, 5: -2.282229e8}},
        "expected_element_results": {1: {"axial_stress": 612.4239}},
        "expected_increments": 40,
    },
    {
        "case_name": "case_138_positive_tangent_probe",
        "folder": "case_138",
        "filename": "ej138_cbar_combined_load_positive_probe.dat",
        "description": "Positive T2 tangent probe around the case 137 state.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: 0.7464024, 1: 53.22897, 2: -39.02622, 3: 0.2197381, 4: 0.05488959, 5: 0.08402530}},
        "expected_reactions": {1: {1: -2.001e5, 5: -2.283120e8}},
        "expected_element_results": {1: {"axial_stress": 612.4323}},
        "expected_increments": 40,
    },
    {
        "case_name": "case_139_negative_tangent_probe",
        "folder": "case_139",
        "filename": "ej139_cbar_combined_load_negative_probe.dat",
        "description": "Negative T2 tangent probe around the case 137 state.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: 0.7485292, 1: 53.19295, 2: -39.01869, 3: 0.2197482, 4: 0.05487528, 5: 0.08397220}},
        "expected_reactions": {1: {1: -1.999e5, 5: -2.281339e8}},
        "expected_element_results": {1: {"axial_stress": 612.4152}},
        "expected_increments": 40,
    },
    {
        "case_name": "case_140_inclined_g0_cid",
        "folder": "case_140",
        "filename": "ej140_cbar_inclined_3d_g0_cd_cid.dat",
        "description": "Inclined CBAR orientation by G0 with transformed loads.",
        "rel_tol": 6.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -229.0154, 1: 283.4722, 2: -65.63926, 3: 0.01922614, 4: 0.1559748, 5: 0.6531215}},
        "expected_reactions": {1: {0: 123743.7, 1: -229809.7, 2: 125000.0, 3: 1.578028e8, 4: -2.434119e7, 5: -2.681675e8}},
        "expected_element_results": {1: {"axial_stress": 54.9315}},
        "expected_increments": 40,
    },
    {
        "case_name": "case_141_asymmetric_i12",
        "folder": "case_141",
        "filename": "ej141_cbar_asymmetric_section_i12.dat",
        "description": "Coupled biaxial bending from nonzero PBAR I12.",
        "rel_tol": 4.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -3.918113, 1: 57.24796, 2: 64.31614, 4: -0.1073329, 5: 0.09553706}},
        "expected_reactions": {1: {1: -2.0e5, 2: -1.4e5, 4: 1.254515e8, 5: -1.792164e8}},
        "expected_increments": 36,
    },
    {
        "case_name": "case_142_pbarn1_ls_shear_flexible",
        "folder": "case_142",
        "filename": "ej142_cbar_short_shear_flexible.dat",
        "description": "Short beam with PBAR K1/K2 and PBARN1 INT2=LS.",
        "note": "Known mismatch: MSC's advanced LS section interpolation is not represented by the reference PBAR Timoshenko matrix.",
        "expected_failure": "PBARN1 INT2=LS advanced section interpolation remains unsupported.",
        "rel_tol": 2.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -0.07416349, 1: 5.127897, 2: 2.893774, 3: -3.142474e-5, 4: -0.01904017, 5: 0.03808233}},
        "expected_reactions": {1: {1: -8.0e5, 2: -4.0e5, 4: 7.997033e7, 5: -1.599407e8}},
        "expected_element_results": {1: {"axial_stress": 13.14873}},
        "expected_increments": 30,
    },
    {
        "case_name": "case_143_enforced_transverse_displacement",
        "folder": "case_143",
        "filename": "ej143_cbar_enforced_transverse_displacement.dat",
        "description": "Nonzero transverse SPC and its reactions without FORCE.",
        "rel_tol": 1.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -10.76035, 1: 150.0, 5: 0.2257281}},
        "expected_reactions": {1: {1: -766559.3, 5: -7.583108e8}, 2: {1: 766559.3}},
        "expected_element_results": {1: {"axial_stress": 114.9210}},
        "expected_increments": 30,
    },
    {
        "case_name": "case_144_multiple_enforced_dofs_and_force",
        "folder": "case_144",
        "filename": "ej144_cbar_enforced_axial_rotation_with_force.dat",
        "description": "Enforced axial displacement and R3 combined with FORCE.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: 50.0, 1: 53.89090, 5: 0.15}},
        "expected_reactions": {2: {0: 1.053087e7, 5: 4.770293e8}},
        "expected_element_results": {1: {"axial_stress": 10522.15}},
        "expected_increments": 36,
    },
    {
        "case_name": "case_145_two_element_multimaterial",
        "folder": "case_145",
        "filename": "ej145_cbar_two_element_multimaterial.dat",
        "description": "Two CBARs with different MAT1/PBAR definitions.",
        "rel_tol": 3.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -0.4754261, 1: 24.49795, 3: 0.05885688, 4: 0.01154607, 5: 0.07354454}, 3: {0: -4.603622, 1: 91.46203, 2: -24.18931, 3: 0.1462988, 4: 0.06837297, 5: 0.1349438}},
        "expected_reactions": {1: {1: -2.56e5, 2: -168937.9, 5: -3.083517e8}},
        "expected_element_results": {1: {"axial_stress": 8.7100}, 2: {"axial_stress": 12.8813}},
        "expected_increments": 42,
    },
    {
        "case_name": "case_146_refined_four_element_cantilever",
        "folder": "case_146",
        "filename": "ej146_cbar_four_element_refined_cantilever.dat",
        "description": "Four-element refinement of the case 132 cantilever.",
        "rel_tol": 1.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {3: {0: -0.3875478, 1: 18.5523, 5: 0.06679726}, 5: {0: -2.002073, 1: 59.32523, 5: 0.08904246}},
        "expected_reactions": {1: {1: -3.0e5, 5: -2.993994e8}},
        "expected_element_results": {1: {"axial_stress": 6.12475}, 4: {"axial_stress": 26.12360}},
        "expected_increments": 40,
    },
    {
        "case_name": "case_147_rectangular_frame_release",
        "folder": "case_147",
        "filename": "ej147_cbar_rectangular_frame_end_release.dat",
        "description": "Three-member spatial frame with PB=6 release.",
        "rel_tol": 7.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {2: {0: -0.03381342, 1: -1.588213, 2: 0.6387557, 5: -0.002464904}, 3: {0: 0.4606923, 1: -1.508575, 2: 1.538141, 5: -0.0005517536}},
        "expected_reactions": {1: {0: 10022.15, 1: 49476.3, 2: -8393.725, 4: 5.407157e6, 5: 1.916532e7}},
        "expected_element_results": {1: {"axial_stress": -10.98022}, 3: {"axial_stress": 162.5322}},
        "expected_increments": 22,
    },
    {
        "case_name": "case_148_spatial_frame_offsets_sway",
        "folder": "case_148",
        "filename": "ej148_cbar_spatial_frame_offsets_enforced_sway.dat",
        "description": "Spatial frame with offsets, mixed properties and sway SPC.",
        "rel_tol": 1.5e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {4: {0: 10.0, 1: -1.606804, 2: -0.3157204, 3: 0.003417011, 4: 0.02312770, 5: 0.02144647}, 6: {0: 49.92885, 1: 8.907005, 2: -114.3430, 3: 0.1342353, 4: 0.2601071, 5: 0.06726841}},
        "expected_reactions": {1: {0: -112076.0, 1: 50606.29, 2: 195677.2}, 2: {0: 64311.42, 1: -606.2875, 2: 99322.76}},
        "expected_element_results": {1: {"axial_stress": -107.0665}, 5: {"axial_stress": 40.86925}},
        "expected_increments": 100,
    },
    {
        "case_name": "case_149_two_story_planar_frame",
        "folder": "case_149",
        "filename": "ej149_cbar_two_story_planar_frame_15_elements.dat",
        "description": "Fifteen-member two-story frame with global P-delta.",
        "rel_tol": 3.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {10: {0: 14.22139, 1: -0.1779158, 5: -0.002963303}, 12: {0: 14.15392, 1: -0.7310958, 5: -0.006010935}},
        "expected_reactions": {1: {0: -198840.6, 1: -162189.6, 5: 1.075273e8}},
        "expected_element_results": {12: {"axial_stress": -38.75915}, 15: {"axial_stress": 254.4926}},
        "expected_increments": 40,
    },
    {
        "case_name": "case_150_two_level_space_frame_adaptive",
        "folder": "case_150",
        "filename": "ej150_cbar_two_level_space_frame_20_elements.dat",
        "description": "Twenty-member space frame with rotated loads and ADAPT.",
        "rel_tol": 6.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {10: {0: 6.369471, 1: 11.44745, 2: -0.7304186, 3: -0.004087561, 4: 0.002582058, 5: 0.008510091}, 12: {0: 0.6841936, 1: 3.127055, 2: -0.7045448, 3: -0.001979834, 4: 0.0003218793, 5: 0.006298014}},
        "expected_reactions": {1: {0: -218715.5, 1: -31390.03, 2: -305955.3, 3: 2.046953e7, 4: -1.140910e8}},
        "expected_element_results": {18: {"axial_stress": -19.14626}, 20: {"axial_stress": -13.96953}},
    },
    {
        "case_name": "case_151_semicircular_arch_settlement",
        "folder": "case_151",
        "filename": "ej151_cbar_semicircular_arch_20_elements_settlement.dat",
        "description": "Curved twenty-member arch with imposed support settlement.",
        "rel_tol": 5.0e-3,
        "abs_tol": 1.0e-5,
        "expected_displacements": {19: {0: 28.88042, 1: -16.93187, 5: -0.03804681}, 21: {1: -20.0, 5: -0.05100275}},
        "expected_reactions": {1: {0: 145642.7, 1: 399316.9, 5: 1.749313e7}, 21: {0: -215642.7, 1: 500683.1}},
        "expected_element_results": {18: {"axial_stress": -226.9025}, 20: {"axial_stress": -211.1142}},
        "expected_increments": 50,
    },
    {
        "case_name": "case_152_warehouse_offsets_adaptive",
        "folder": "case_152",
        "filename": "ej152_cbar_warehouse_space_frame_20_elements_offsets.dat",
        "description": "Braced space frame with offsets, combined load and ADAPT.",
        "rel_tol": 1.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {8: {0: 23.62955, 1: -4.125601, 2: -1.288910, 3: 0.005113720, 4: 0.02257335, 5: -0.008415081}, 10: {0: 25.0, 1: 2.242930, 2: -3.534115, 3: -0.0006486959, 4: -0.01213227, 5: -0.007062152}},
        "expected_reactions": {1: {0: -58660.88, 1: -96823.25, 2: 147718.0, 3: 8.182445e7, 4: -4.980801e7}},
        "expected_element_results": {20: {"axial_stress": 36.06776}},
    },
    {
        "case_name": "case_153_spatial_zigzag_boom",
        "folder": "case_153",
        "filename": "ej153_cbar_spatial_zigzag_boom_15_elements_mixed_sections.dat",
        "description": "Extreme 3D boom response with five sections and imposed R3.",
        "rel_tol": 4.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {16: {0: -1785.712, 1: 106.6657, 2: -7707.853, 3: 0.7276860, 4: 1.397984, 5: 0.1200001}},
        "expected_reactions": {1: {0: -274000.0, 1: 70000.0, 2: 810000.0, 4: -1.585066e9, 5: 1.297450e8}},
        "expected_element_results": {12: {"axial_stress": 329.9808}, 15: {"axial_stress": 260.9041}},
        "expected_increments": 60,
    },
]


def _validation_parameter(reference):
    mark = ()
    if reference.get("expected_failure"):
        mark = pytest.mark.xfail(
            strict=True,
            reason=reference["expected_failure"],
        )
    return pytest.param(reference, id=reference["case_name"], marks=mark)


@pytest.mark.parametrize(
    "reference",
    [_validation_parameter(reference) for reference in NONLINEAR_CBAR_VALIDATION_DATA],
)
def test_nonlinear_cbar_final_results(reference):
    deck = VALIDATION_ROOT / reference["folder"] / reference["filename"]
    assert deck.is_file(), f"Missing DAT for {reference['case_name']}: {deck}"

    solver = NonlinearStaticSolver(parse_nastran_input(str(deck)))
    result = solver.solve()
    rel_tol = reference["rel_tol"]
    abs_tol = reference["abs_tol"]

    assert result["load_factor"] == pytest.approx(1.0)
    if "expected_increments" in reference:
        assert len(result["increments"]) == reference["expected_increments"]

    for node_id, dofs in reference.get("expected_displacements", {}).items():
        computed = result["displacements"][solver.dof_mapper.get_dof_indices(node_id)]
        for dof, expected in dofs.items():
            assert computed[dof] == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            )

    for node_id, dofs in reference.get("expected_reactions", {}).items():
        computed = result["reactions"][node_id]
        for dof, expected in dofs.items():
            assert computed[dof] == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            )

    for element_id, expected_results in reference.get(
        "expected_element_results", {}
    ).items():
        computed = result["stresses"][element_id]
        for result_name, expected in expected_results.items():
            assert computed[result_name] == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            )

    if "expected_strain_energy" in reference:
        assert result["strain_energy"] == pytest.approx(
            reference["expected_strain_energy"], rel=rel_tol, abs=abs_tol
        )
