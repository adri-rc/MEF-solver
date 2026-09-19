# Working guide for AI assistants

These instructions apply to the entire repository. Their purpose is to preserve
the numerical traceability of this academic solver against MSC Nastran.

## Project onboarding

Before proposing or modifying code:

1. run `git status --short --branch` and do not overwrite local changes;
2. read `README.md`, `docs/HANDOVER.md`, and `docs/roadmap.md`;
3. consult the technical document for the area being changed;
4. review the related tests and validation cases;
5. run `.\.venv\Scripts\python.exe -m pytest -q` as a baseline on Windows
   when the environment exists.

The verifiable state, recent history, and recommended reading order are in
`docs/HANDOVER.md`. Code and tests take precedence if a document is outdated;
update that document in the same change.

## Invariants that must be preserved

- The linear path covers small strains with isotropic `MAT1` and six DOFs per
  node for CROD, CBAR, CTRIA3, and coplanar CQUAD4.
- The `SOL 400` path covers geometric nonlinearity only with `MAT1`,
  `PROD/CROD`, and `PBAR/PBARN1/CBAR`, conservative nodal loads, full Newton,
  and `FIXED`/`ADAPT`. It does not include plasticity, contact, arc-length
  control, or follower loads.
- Do not present the project as a general replacement for Nastran.
- There is no general `AUTOSPC`: a real structural mechanism must fail rather
  than be stabilized silently.
- `GRID.CP` coordinates, DOFs in `GRID.CD`, and load vectors in their `CID`
  must be transformed at the stages documented in `docs/architecture.md`.
- For PSHELL, a blank `MID2` disables bending and a blank `MID3` disables the
  physical transverse-shear flexibility; neither field inherits another MID.
- CQUAD4 rejects warped geometries. The strict warped-geometry and `SNORM`
  `xfail` cases are explicit scope limits, not pending defects in the coplanar
  element.
- Do not recalibrate CQUAD4 stiffness solely to fit centroidal `Qx/Qy` without
  a new DMIG reference: the known discrepancy is in recovery.
- Result presentation belongs in `src/postprocessing/report_writer.py`; the
  solver retains numerical recovery and `main.py` only orchestrates the flow.

## Working practices

- Make small changes and support them with focused tests first. Run the full
  suite before delivery whenever the environment permits.
- Every new physical capability needs a minimal case, an auditable external
  reference, and an end-to-end test that fails before implementation.
- Do not relax tolerances or convert a failure to `xfail` without documenting
  the numerical cause and reference evidence.
- Keep code names, identifiers, comments, docstrings, and user documentation
  in English. Write delivery summaries in the user's language.
- `.dat`/`.bdf` models, `.f06` results, and other binaries use Git LFS. Confirm
  that they are not text pointers before diagnosing a read failure.
- `output/`, caches, and `draft_*` files are local artifacts; they are not a
  source of truth and must not be committed.
- `docs/TFM/` may contain local, untracked thesis material. Do not edit, delete,
  or add it to Git without an explicit user request.
- Do not create commits, push, or rewrite history unless explicitly requested.
  At completion, report tests, changed files, and risks.

## Documentation accompanying changes

- Usage, dependency, or CLI changes: `README.md` and `docs/usage.md`.
- Module or flow changes: `docs/architecture.md`.
- Formulation changes: `docs/formulation.md` and the element document.
- New references or tolerances: `docs/validation.md`.
- General status, priorities, or decisions: `docs/HANDOVER.md` and
  `docs/roadmap.md`.
