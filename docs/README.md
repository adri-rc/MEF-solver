# MEF-solver documentation

MEF-solver is an academic finite element solver for static structural analysis
and comparison with MSC Nastran models. It includes linear statics for CROD,
CBAR, CTRIA3, and coplanar CQUAD4, plus geometric nonlinear `SOL 400` for
`MAT1` models made of CROD and/or CBAR. The project implements a deliberately
small subset of Nastran and is not a general replacement for it.

This documentation describes the repository state verified on
**19 September 2026**. The complete suite reports
`718 passed, 6 xfailed`.

## Contents

- [Usage](usage.md): environment, execution, supported input, and results.
- [Result report](result_report.md): report tables, F06 correspondence, and
  Newton history.
- [Architecture](architecture.md): modules, data model, and analysis flow.
- [Mathematical formulation](formulation.md): equations, conventions, and
  element formulations.
- [Validation](validation.md): actual coverage, test results, and Nastran
  references.
- [CQUAD4 development and validation](cquad4.md): detailed element status.
- [Complete CQUAD4 formulation](cquad4_formulation.md): mathematical
  derivation, parameters, code correspondence, evidence, and limitations.
- [Nonlinear CROD analysis](nonlinear_crod.md): formulation, `SOL 400`
  entries, incremental solver, and cases 92--127.
- [Nonlinear CBAR analysis](nonlinear_cbar.md): corotational formulation,
  limitations, cases 128--166, and mixed CROD/CBAR cases 167--176.
- [Roadmap](roadmap.md): consolidation work and possible future extensions.
- [Technical handover](HANDOVER.md): verified state and continuation context.

## Capability summary

| Capability | Read | Stiffness/tangent | Global solution | Recovery | Validation |
|---|---:|---:|---:|---:|---|
| CROD + PROD + MAT1 | Yes | Yes | Yes | Axial and torsional | Unit and Nastran E2E |
| CBAR + PBAR + MAT1 | Yes | Yes | Yes | End forces and fibre stresses | Unit and Nastran E2E |
| CTRIA3 + PSHELL + MAT1 | Yes | Yes | Yes | Membrane, bending, and shear | DMIG and Nastran E2E |
| CQUAD4 + PSHELL + MAT1 | Yes | Yes | Yes | Membrane, bending, shear, and fibres | DMIG and E2E; coplanar only |
| SOL 400 + NLSTEP/NLPARM | Yes | CROD/CBAR tangent | Yes | Increment history and CROD/CBAR results | MSC campaigns 92--176 and two integration bridges |
| FORCE / MOMENT / LOAD | Yes | N/A | Yes | Loads transformed from CID | Automated tests |
| SPC / SPC1 / SPCADD | Yes | N/A | Yes | Reactions and prescribed values | Automated tests |
| 3D solid elements | No | No | No | No | Outside scope |

## Current physical scope

- Small-deformation linear statics.
- Geometric nonlinear CROD/CBAR statics with `SOL 400`, `LGDISP`, full
  Newton, and `FIXED`/`ADAPT` controls.
- Isotropic linear elastic `MAT1`.
- Six DOFs per node: `[Tx, Ty, Tz, Rx, Ry, Rz]`.
- Sparse global assembly and direct SciPy solution.
- Nodal loads and single-point constraints.
- No dynamics, material nonlinearity, contact, follower or distributed loads,
  MPC/RBE, warped CQUAD4, or solid elements.

Start with [usage.md](usage.md) to run the program. Read
[HANDOVER.md](HANDOVER.md) before continuing development, and follow the
repository rules in [AGENTS.md](../AGENTS.md).
