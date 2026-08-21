import pytest
import numpy as np
from pathlib import Path

import main as main_module
from main import run_analysis_job
from src.preprocessing.data_model import (
    AnalysisType,
    BarElement,
    BarProperty,
    RodElement,
)
from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


def _write_deck(tmp_path, case_control: str, bulk_control: str):
    deck = f"""SOL 400
CEND
SUBCASE 1
  ANALYSIS = NLSTAT
  SPC = 10
  LOAD = 20
{case_control}
BEGIN BULK
PARAM,LGDISP,1
PARAM,FOLLOWK,NO
GRID,1,,0.,0.,0.
GRID,2,,100.,0.,0.
MAT1,1,210000.,,0.3
PROD,1,1,100.
CROD,1,1,1,2
SPC1,10,123456,1
SPC1,10,23456,2
FORCE,20,2,,210000.,1.,0.,0.
{bulk_control}
ENDDATA
"""
    path = tmp_path / "nonlinear.dat"
    path.write_text(deck)
    return path


def test_sol400_fixed_nlstep_is_mapped(tmp_path):
    path = _write_deck(
        tmp_path,
        "  NLSTEP = 92",
        """NLSTEP,92,1.0,MILDLY
,GENERAL,30,1,12
,FIXED,20
,MECH,PV,1.E-8,2.E-8,3.E-8,PFNT""",
    )

    model = parse_nastran_input(str(path))

    assert model.analysis.solution_sequence == 400
    assert model.analysis.analysis_type is AnalysisType.NONLINEAR_STATIC
    assert model.analysis.analysis_keyword == "NLSTAT"
    assert model.analysis.large_displacement_parameter == 1
    assert model.analysis.large_displacements is True
    assert model.analysis.follower_force_stiffness is False
    control = model.analysis.nonlinear_step
    assert control.sid == 92
    assert control.source == "NLSTEP"
    assert control.control_definition == "MILDLY"
    assert control.increment_method == "FIXED"
    assert control.number_of_increments == 20
    assert control.maximum_iterations == 30
    assert control.maximum_bisections == 12
    assert control.convergence == "PV"
    assert control.displacement_tolerance == pytest.approx(1.0e-8)
    assert control.residual_tolerance == pytest.approx(2.0e-8)
    assert control.work_tolerance == pytest.approx(3.0e-8)
    assert control.stiffness_update_method == "PFNT"


def test_sol400_adaptive_nlstep_and_compact_exponents_are_mapped(tmp_path):
    path = _write_deck(
        tmp_path,
        "  NLSTEP = 108",
        """$ A comment immediately before rejected NLSTEP belongs to its group.
NLSTEP,108,2.0,SEVERELY
,GENERAL,40,2,15
,ADAPT,.05,1.-5,.2,4,1.25,20,300
,MECH,UPV,1.-7,2.-7,3.-7,ITER,1""",
    )

    control = parse_nastran_input(str(path)).analysis.nonlinear_step

    assert control.total_time == pytest.approx(2.0)
    assert control.increment_method == "ADAPT"
    assert control.number_of_increments is None
    assert control.initial_step_fraction == pytest.approx(0.05)
    assert control.minimum_step_fraction == pytest.approx(1.0e-5)
    assert control.maximum_step_fraction == pytest.approx(0.2)
    assert control.desired_iterations == 4
    assert control.step_scale_factor == pytest.approx(1.25)
    assert control.output_intervals == 20
    assert control.maximum_increments == 300
    assert control.stiffness_update_method == "ITER"
    assert control.stiffness_update_interval == 1


def test_sol400_legacy_nlparm_is_normalized(tmp_path):
    path = _write_deck(
        tmp_path,
        "  NLPARM = 7",
        "NLPARM,7,12,,ITER,1,25,PV,YES",
    )

    control = parse_nastran_input(str(path)).analysis.nonlinear_step

    assert control.sid == 7
    assert control.source == "NLPARM"
    assert control.increment_method == "FIXED"
    assert control.number_of_increments == 12
    assert control.maximum_iterations == 25
    assert control.output_intervals == 0
    assert control.stiffness_update_method == "ITER"
    assert control.stiffness_update_interval == 1


def test_sol400_requires_selected_nonlinear_control(tmp_path):
    path = _write_deck(tmp_path, "", "")

    with pytest.raises(ValueError, match="requires an NLSTEP"):
        parse_nastran_input(str(path))


def test_linear_punch_model_keeps_linear_default(tmp_path):
    path = tmp_path / "linear.bdf"
    path.write_text(
        """GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
MAT1,1,1000.,,0.3
PROD,1,1,1.
CROD,1,1,1,2
"""
    )

    model = parse_nastran_input(str(path))

    assert model.analysis.solution_sequence is None
    assert model.analysis.analysis_type is AnalysisType.LINEAR_STATIC
    assert model.analysis.nonlinear_step is None


EXTENDED_CAMPAIGN_COUNTS = {
    123: {"nodes": 6, "elements": 10, "properties": 3, "loads": 7},
    124: {"nodes": 10, "elements": 20, "properties": 5, "loads": 7},
    125: {"nodes": 13, "elements": 30, "properties": 4, "loads": 6},
    126: {"nodes": 9, "elements": 24, "properties": 5, "loads": 6},
    127: {"nodes": 12, "elements": 30, "properties": 6, "loads": 8},
}


@pytest.mark.parametrize("case_id", range(92, 128))
def test_nonlinear_crod_campaign_decks_are_parseable(case_id):
    validation_root = (
        Path(__file__).parents[1] / "solver_tst" / "validation_cases"
    )
    paths = list((validation_root / f"case_{case_id}").glob("*.dat"))
    assert len(paths) == 1

    model = parse_nastran_input(str(paths[0]))

    assert model.analysis.solution_sequence == 400
    assert model.analysis.analysis_type is AnalysisType.NONLINEAR_STATIC
    expected_lgdisp = -1 if case_id in {92, 121} else 2 if case_id == 118 else 1
    assert model.analysis.large_displacement_parameter == expected_lgdisp
    assert model.analysis.large_displacements is (case_id not in {92, 121})
    assert model.analysis.follower_force_stiffness is False
    assert model.analysis.nonlinear_step.sid == case_id
    assert model.active_load_sid == case_id
    assert model.active_spc_sid == case_id
    assert model.elements
    assert all(isinstance(element, RodElement) for element in model.elements.values())

    if case_id in EXTENDED_CAMPAIGN_COUNTS:
        expected = EXTENDED_CAMPAIGN_COUNTS[case_id]
        assert len(model.nodes) == expected["nodes"]
        assert len(model.elements) == expected["elements"]
        assert len(model.properties) == expected["properties"]
        assert len(model.loads) == expected["loads"]
        assert len(model.materials) == 2

    if case_id == 124:
        assert model.source_entry_counts["SPC"] == 1
        assert model.source_entry_counts["SPC1"] == 3
        assert model.source_entry_counts["FORCE"] == 7
    elif case_id == 127:
        assert model.source_entry_counts["SPC1"] == 2
        assert model.source_entry_counts["FORCE"] == 8

    prescribed_values = sorted(
        constraint.value for constraint in model.constraints if constraint.value != 0.0
    )
    if case_id == 124:
        assert prescribed_values == [-4.0]
    elif case_id == 126:
        assert prescribed_values == [5.0]


NONLINEAR_CBAR_CAMPAIGN_COUNTS = {
    **{
        case_id: {"nodes": 2, "elements": 1, "properties": 1}
        for case_id in range(128, 140)
    },
    140: {"nodes": 3, "elements": 1, "properties": 1},
    141: {"nodes": 2, "elements": 1, "properties": 1},
    142: {"nodes": 2, "elements": 1, "properties": 1},
    143: {"nodes": 2, "elements": 1, "properties": 1},
    144: {"nodes": 2, "elements": 1, "properties": 1},
    145: {"nodes": 3, "elements": 2, "properties": 2},
    146: {"nodes": 5, "elements": 4, "properties": 1},
    147: {"nodes": 4, "elements": 3, "properties": 1},
    148: {"nodes": 6, "elements": 5, "properties": 2},
    149: {"nodes": 12, "elements": 15, "properties": 3},
    150: {"nodes": 12, "elements": 20, "properties": 4},
    151: {"nodes": 21, "elements": 20, "properties": 4},
    152: {"nodes": 10, "elements": 20, "properties": 4},
    153: {"nodes": 16, "elements": 15, "properties": 5},
    **{
        case_id: {"nodes": 2, "elements": 1, "properties": 1}
        for case_id in range(154, 161)
    },
}


@pytest.mark.parametrize("case_id", range(128, 161))
def test_nonlinear_cbar_campaign_decks_are_parseable(case_id):
    validation_root = (
        Path(__file__).parents[1] / "solver_tst" / "validation_cases"
    )
    paths = list((validation_root / f"case_{case_id}").glob("*.dat"))
    assert len(paths) == 1

    model = parse_nastran_input(str(paths[0]))
    expected = NONLINEAR_CBAR_CAMPAIGN_COUNTS[case_id]

    assert model.analysis.solution_sequence == 400
    assert model.analysis.analysis_type is AnalysisType.NONLINEAR_STATIC
    assert model.analysis.large_displacement_parameter == (-1 if case_id < 130 else 1)
    assert model.analysis.large_displacements is (case_id >= 130)
    assert model.analysis.follower_force_stiffness is False
    assert model.analysis.nonlinear_step.sid == case_id
    assert model.active_spc_sid == case_id
    assert model.active_load_sid == (None if case_id == 143 else case_id)
    assert len(model.nodes) == expected["nodes"]
    assert len(model.elements) == expected["elements"]
    assert len(model.properties) == expected["properties"]
    assert all(isinstance(element, BarElement) for element in model.elements.values())
    assert all(isinstance(prop, BarProperty) for prop in model.properties.values())


@pytest.mark.parametrize("case_id", range(128, 161))
def test_nonlinear_cbar_campaign_activates_sol400_bar_output(case_id):
    validation_root = (
        Path(__file__).parents[1] / "solver_tst" / "validation_cases"
    )
    path = next((validation_root / f"case_{case_id}").glob("*.dat"))
    deck_text = path.read_text()

    # PBARN1 is the SOL 400 nonlinear extension associated with each PBAR.
    assert deck_text.count("\nPBARN1,") == deck_text.count("\nPBAR,")
    assert "NLSTRESS(SORT1) = ALL" in deck_text
    assert "BOUTPUT(SORT1,REAL) = ALL" in deck_text
    assert "ELFORCE(SORT1,REAL,PRINT,PUNCH) = ALL" in deck_text
    assert "STRAIN(SORT1,REAL,PRINT,PUNCH) = ALL" in deck_text


@pytest.mark.parametrize("case_id", range(149, 154))
def test_large_nonlinear_cbar_campaign_has_no_initial_linear_mechanism(case_id):
    validation_root = (
        Path(__file__).parents[1] / "solver_tst" / "validation_cases"
    )
    path = next((validation_root / f"case_{case_id}").glob("*.dat"))
    solver = LinearStaticSolver(parse_nastran_input(str(path)))
    solver.build_stiffness_matrix()
    free_dofs, _, _ = solver._get_dof_partition()
    free_stiffness = solver.K_global.tocsr()[free_dofs, :][:, free_dofs].toarray()

    assert np.linalg.matrix_rank(free_stiffness) == len(free_dofs)


def test_nonlinear_cbar_campaign_preserves_advanced_card_fields():
    validation_root = (
        Path(__file__).parents[1] / "solver_tst" / "validation_cases"
    )

    def load(case_id):
        path = next((validation_root / f"case_{case_id}").glob("*.dat"))
        return parse_nastran_input(str(path))

    inclined = load(140)
    assert inclined.elements[1].g0 == 99
    assert inclined.nodes[2].cd_id == 0
    assert inclined.loads[0].cid == 0
    assert inclined.loads[0].vector == pytest.approx(
        [-0.49497474683058323, 0.9192388155425117, -0.5]
    )
    assert inclined.loads[1].cid == 0
    assert inclined.loads[1].vector == pytest.approx(
        [0.4242640687119285, -0.1414213562373095, 1.0]
    )

    coupled = load(141)
    assert coupled.properties[1].I12 == pytest.approx(1.2e6)

    shear_flexible = load(142)
    assert shear_flexible.properties[1].K1 == pytest.approx(0.75)
    assert shear_flexible.properties[1].K2 == pytest.approx(0.50)
    assert shear_flexible.properties[1].nonlinear_section == "S"
    assert shear_flexible.properties[1].nonlinear_behavior == "BAR"
    assert shear_flexible.properties[1].nonlinear_integration == "LS"

    combined_loading = load(137)
    assert combined_loading.properties[1].nonlinear_section == "S"
    assert combined_loading.properties[1].nonlinear_behavior == "BAR"
    assert combined_loading.properties[1].nonlinear_integration == "LC"

    diagnostic_ls = load(155)
    diagnostic_lc = load(156)
    assert diagnostic_ls.properties[1].nonlinear_integration == "LS"
    assert diagnostic_lc.properties[1].nonlinear_integration == "LC"

    enforced = load(143)
    assert sorted(c.value for c in enforced.constraints if c.value) == [150.0]

    enforced_and_loaded = load(144)
    assert sorted(c.value for c in enforced_and_loaded.constraints if c.value) == [0.15, 50.0]

    released = load(147)
    assert released.elements[1].pb == 6

    offset_frame = load(148)
    for eid in (3, 4):
        assert offset_frame.elements[eid].wa == pytest.approx([0.0, 20.0, 0.0])
        assert offset_frame.elements[eid].wb == pytest.approx([0.0, 20.0, 0.0])


@pytest.mark.parametrize("case_id", range(128, 161))
def test_nonlinear_cbar_campaign_documents_each_bulk_card_header(case_id):
    validation_root = (
        Path(__file__).parents[1] / "solver_tst" / "validation_cases"
    )
    path = next((validation_root / f"case_{case_id}").glob("*.dat"))
    text = path.read_text()
    bulk_text = text.split("BEGIN BULK", maxsplit=1)[1].split("ENDDATA", maxsplit=1)[0]
    documented_cards = {
        line.split(",", maxsplit=1)[0]
        for line in (raw_line.strip() for raw_line in bulk_text.splitlines())
        if line and not line.startswith(("$", ",", "+"))
    }

    for card_name in documented_cards:
        assert f"$ {card_name}:" in text, f"Missing {card_name} field header in {path.name}"


def test_cli_routes_sol400_to_nonlinear_solver(tmp_path, monkeypatch):
    path = _write_deck(
        tmp_path,
        "  NLSTEP = 1",
        """NLSTEP,1,1.0,QLINEAR
,FIXED,1
,MECH,PV,1.E-8,1.E-8,1.E-8,PFNT""",
    )

    output_dir = tmp_path / "output"
    output_dir.mkdir()
    monkeypatch.setattr(main_module, "OUTPUT_DIR", output_dir)

    result = run_analysis_job(path)

    assert result["load_factor"] == pytest.approx(1.0)
    assert result["displacements"][6] == pytest.approx(1.0)
    report = (output_dir / "nonlinear.rpt").read_text()
    assert "NONLINEAR STATIC ANALYSIS REPORT" in report
    assert "ANALYSIS CONFIGURATION" in report
    assert "NONLINEAR STEP CONTROL" in report
    assert "NUMBER OF INCREMENTS" in report
    assert "CONVERGED INCREMENTS" in report
    assert "CONVERGED INCREMENT HISTORY" in report
    assert "NEWTON ITERATION HISTORY" in report
    assert "FINAL CONVERGED STATE RESULTS" in report
    assert report.index("FINAL CONVERGED STATE RESULTS") < report.index(
        "CONVERGED INCREMENT HISTORY"
    )
