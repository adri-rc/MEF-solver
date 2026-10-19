# Result report reference

## Purpose and report order

`src/postprocessing/report_writer.py` writes `output/<job>.rpt`. The report
presents quantities calculated by MEF-solver and explains where to find their
counterparts in an MSC Nastran F06. It does not reproduce the F06 layout
verbatim. Every value uses the coherent unit system of the input deck.

The report order is:

1. linear or nonlinear heading;
2. `MODEL SUMMARY` and analysis configuration;
3. `SOL 400` control and convergence summary, when applicable;
4. final-state results;
5. converged-increment and Newton-correction histories.

## `MODEL SUMMARY`

This table counts logical entries in the source Bulk Data, like the F06
`M O D E L  S U M M A R Y`. It does not count the normalized objects used
internally by the solver.

Constraint counting deserves particular attention:

- `SPC` counts `SPC` entries with prescribed values;
- `SPC1` counts logical `SID`/component combinations; repeated `SPC1`
  lines with the same fields are consolidated like MSC Nastran;
- `SPCADD` remains a separate type when present;
- a permanent constraint in `GRID.PS` belongs to the `GRID` entry and does
  not create another `SPC` summary entry;
- `FORCE` and `MOMENT` are loads, remain separate types, and never
  contribute to the `SPC` count.

For example, case 127 expands to twelve constrained nodal components for the
solution, while both F06 and RPT show `SPC1 = 2`: the deck contains two
logical groups, `123456` and `456`, for one SID.

## Analysis configuration and `SOL 400` control

| RPT block | F06 location |
|---|---|
| `ANALYSIS CONFIGURATION` | `CASE CONTROL ECHO` for `ANALYSIS`, `SPC`, `LOAD`, and `NLSTEP`; Executive Control for `SOL 101` or `SOL 400`. |
| `PARAM,LGDISP` and `PARAM,FOLLOWK` | Bulk Data or deck echo; both also contribute to `PARAM` in `MODEL SUMMARY`. |
| `NONLINEAR STEP CONTROL` | `NON-LINEAR ITERATION SOLUTION CONTROL PARAMETERS`. MSC uses fields such as `TOTTIME`, `MAXITER`, `MAXBIS`, `DTINITF`, `DTMINF`, `DTMAXF`, `NDESIR`, `SFACT`, `EPSP`, and `KMETHOD`. |
| `NONLINEAR SOLUTION SUMMARY` | Derived from the last `LOAD STEP`, convergence messages, iteration table, and `GLOBAL WORKS & ENERGIES`; MSC does not collect the same five values in one table. |

`FINAL LOAD FACTOR = 1` means that the full reference load was reached. It
does not mean that forces equal one: the solver applies
`F_ext(lambda) = lambda F_ref`.

## Final results

`FINAL CONVERGED STATE RESULTS` identifies the last accepted equilibrium.
With intermediate output enabled, an F06 may repeat the same tables at every
load step. Compare the report with the final occurrence under the last
converged `LOAD STEP`.

| RPT result | Typical `SOL 400` F06 heading |
|---|---|
| `NODAL DISPLACEMENTS` | `D I S P L A C E M E N T  V E C T O R`. |
| `SPC REACTION FORCES AND MOMENTS` | `F O R C E S  O F  S I N G L E - P O I N T  C O N S T R A I N T`. |
| CROD stress/strain | `N O N L I N E A R  S T R E S S E S  I N  R O D  E L E M E N T S (CROD)`. |
| `SECTION RESPONSE IN CBAR ELEMENTS` | Centroidal axial stress and strain relate to `N O N L I N E A R  S T R E S S E S  I N  B A R  E L E M E N T S (CBAR)`. Under bending the F06 reports combined C--F values, not a separate centroidal value. The preserved F06 files have no direct columns for axial force, torque, twist rate, or torsional stress. |
| `NONLINEAR STATE IN CBAR ELEMENTS` | `CURRENT LENGTH` is reconstructed from displaced geometry and offsets; MSC does not print it directly in these listings. `ELEMENT STRAIN ENERGY` is the element's elastic contribution, while `GLOBAL WORKS & ENERGIES` reports only the model total. |
| `END FORCES AND MOMENTS IN CBAR ELEMENTS` | The reference `SOL 400` F06 files contain no textual OEF table. Checks use equilibrium with `OLOAD`, `SPCFORCES`, and PCH/OP2 when available. |
| `STRESSES AT RECOVERY POINTS IN CBAR ELEMENTS` | `N O N L I N E A R  S T R E S S E S  I N  B A R  E L E M E N T S (CBAR)`, using `GRID` for the end and C, D, E, or F for the point. |
| `STRAIN ENERGY` | `GLOBAL WORKS & ENERGIES`, field `TOT STRAIN ENERGY`; after a branch jump MSC may accumulate work rather than report the instantaneous elastic state energy used here. |
| CROD axial force and torque | A readable OEF table is not always present. Compare against stress, `OLOAD`, `SPCFORCES`, and PCH/OP2 when available. |

## Linear and nonlinear CBAR results

The report distributes CBAR response over four tables. Three are present in
both `SOL 101` and `SOL 400`. `NONLINEAR STATE IN CBAR ELEMENTS` appears
only when the element result contains converged geometric state.

### `SECTION RESPONSE IN CBAR ELEMENTS`

| Column | Meaning |
|---|---|
| `AXIAL FORCE` | Axial section resultant, positive in tension along the converged local axis. |
| `TORQUE` | Torsional resultant about the element local axis. |
| `AXIAL STRESS` | Uniform axial stress. In large-displacement `SOL 400`, it is $E\log(L/L_0)$. |
| `MAX TORSIONAL STRESS` | Largest absolute value of the torsional approximation at the configured recovery points; printed as nonnegative. |
| `AXIAL STRAIN` | Centroidal strain: linear in `SOL 101`, logarithmic in large-displacement `SOL 400`. |
| `TWIST RATE` | Relative section rotation per unit reference length. |

In a purely axial case, `AXIAL STRESS` and `AXIAL STRAIN` can be located
directly in the nonlinear F06 CBAR table because all C--F points receive the
same contribution. Under bending, MSC reports axial plus bending at each
point; derive the centroidal contribution from the section law or point
combinations rather than comparing it with one C--F row.

### `NONLINEAR STATE IN CBAR ELEMENTS`

| Column | Meaning |
|---|---|
| `CURRENT LENGTH` | Distance between current ends, including `WA/WB` offsets rotating with their nodes. |
| `ELEMENT STRAIN ENERGY` | Instantaneous elastic energy in the CBAR at the last converged state. Summing CROD and CBAR contributions gives the project global energy. |
| `MAX NORMAL STRESS` | Largest total normal stress over all C--F points and both ends. |
| `MIN NORMAL STRESS` | Smallest total normal stress over all C--F points and both ends. |

The table is omitted for linear CBAR results because their dictionaries do not
contain `current_length` or individual element energy. Its absence does not
mean that linear forces or stresses are missing; it avoids inventing zeros for
quantities recovered only by the nonlinear path.

### End actions and C--F points

`END FORCES AND MOMENTS IN CBAR ELEMENTS` gives local nodal actions
`Fx`, `Fy`, `Fz`, `Mx`, `My`, and `Mz` at ends A and B. They are
forces exerted by the element on its nodes, so their signs are not two section
cuts with the same outward normal. Offsets may add the moment
$\mathbf{r}\times\mathbf{f}$.

`STRESSES AT RECOVERY POINTS IN CBAR ELEMENTS` separates:

- `AXIAL`: uniform contribution;
- `BENDING`: contribution from `My`, `Mz`, `I1`, `I2`, and `I12`;
- `TOTAL`: `AXIAL + BENDING`;
- `TORSIONAL`: the $T r/J$ approximation, separate from normal stress.

The numerical dictionary uses internal keys `P0` through `P3`. The RPT
presents them as C through F in the order defined by `PBAR`. Cases 161--166
are the principal reference for locating these rows in the final F06 load
step.

## `CONVERGED INCREMENT HISTORY`

Each row is an increment that reached equilibrium and was accepted.

| Column | Meaning | F06 correspondence |
|---|---|---|
| `INCREMENT` | Consecutive accepted increment number in MEF-solver. | `INC` in `NON-LINEAR ITERATION MODULE OUTPUT`; numbering may differ if MSC adapts differently. |
| `LOAD FACTOR` | Total reached factor $\lambda$. | `LOAD STEP` in the iteration row and result block. |
| `STEP SIZE` | Difference between the current and previous factors. | `TIME INCREMENT HAS BEEN CHANGED` messages and adjacent load-step differences. |
| `ITERATIONS` | Newton corrections required by the accepted attempt. | Number of `ITR` rows in the attempt that converges. |
| `RESIDUAL NORM` | Final infinity norm over free DOFs. | Conceptually `MAXIMUM RESIDUAL FORCE`; exact identity is not expected with different internal systems/scaling. |
| `RELATIVE RESIDUAL` | Residual divided by the MEF-solver reference. | Analogous to `FORCE RESIDUAL CONVERGENCE RATIO`, but with a different unpublished MSC denominator. |
| `CUTBACKS` | Failed attempts immediately preceding the accepted reduced step. | Bisection/increment-change messages and iteration-control columns. |
| `STRAIN ENERGY` | Elastic energy of the accepted state. | `TOT STRAIN ENERGY`, subject to the branch-jump caveat above. |

## `NEWTON ITERATION HISTORY`

For one target load factor, full Newton-Raphson repeats:

```text
r_i = lambda F_ext - F_int(u_i)
K_T(u_i) Delta_u_i = r_i
u_(i+1) = u_i + Delta_u_i
```

Each history row is written after applying one correction $\Delta u_i$.

| Column | Meaning | F06 correspondence |
|---|---|---|
| `INCREMENT` | Next converged increment being attempted; the same number reappears after a cutback. | `INC`, with possible differences from MSC's adaptive strategy. |
| `LOAD FACTOR` | Attempt target; it may change when retrying after a cutback. | `LOAD STEP` in the row beginning with `%`. |
| `ITERATION` | Newton correction number within the attempt. | `ITR`. |
| `RESIDUAL NORM` | `max(abs(r))` over free DOFs after correction, including forces and moments according to DOF. | Conceptually `MAXIMUM RESIDUAL FORCE ... IS EQUAL TO ...`. |
| `RELATIVE RESIDUAL` | `RESIDUAL NORM / max(norm(lambda F_ext), norm(F_int), 1)`, using infinity norms on free DOFs. This is the project `P` criterion. | MSC's analogue is `FORCE RESIDUAL CONVERGENCE RATIO`, with a different internal reference. |
| `CORRECTION NORM` | `max(abs(Delta_u))` over free DOFs. | No direct column in the validation F06 files. MSC `DISP MAX` is total state displacement, not this correction. |

### Why the row count varies

- `FIXED` fixes the nominal number of increments, not the iterations within
  each increment.
- `ADAPT` also changes increment count and size according to
  `DESIRED ITERATIONS` and `STEP SCALE FACTOR`.
- Stronger nonlinearity generally needs more corrections.
- A cutback keeps the failed-attempt rows, reduces the step, and retries the
  same increment number at another target.
- A state already within tolerance may converge with zero correction rows.
- `MINIMUM ITERATIONS`, `MAXIMUM ITERATIONS`, and `MAXIMUM BISECTIONS`
  bound the process.

Consequently, `RECORDED NEWTON CORRECTIONS` may exceed the sum of
`CONVERGED NEWTON ITERATIONS`: the former includes attempts discarded by
cutbacks, while the latter includes only accepted attempts.

## Nonlinear CROD torsion

Nonlinear CROD accepts `MOMENT` loads and implements elastic torsion when
`MAT1.G` and `PROD.J > 0` are available:

```text
phi = (theta_2 - theta_1) dot n_0
twist_rate = phi / L_0
T = G J twist_rate
tau = T C / J
k_t = G J / L_0
```

The RPT exposes `TORQUE`, `TORSIONAL STRESS`, and `TWIST RATE`. In case
116, for example, the F06 reports `LIN. TORSIONAL STRESS`, rotation in
`R1`, applied moment under `OLOAD`, and torsional reaction under
`SPCFORCES`.

The deliberate limitation is that torsion remains elastic, linear, and
uncoupled from axial response. It is measured about reference axis $n_0$, and
extension does not change $GJ/L_0$. General finite section rotations and
bending-torsion coupling belong to nonlinear CBAR. If `J = 0`, CROD supplies
no torsional stiffness and its axial rotations must be constrained.

In the final `LOAD STEP = 1.00000E+00` of case 117, the F06 reports
`TOTAL STRAIN = 2.000000E-01` and
`LIN. TORSIONAL STRESS = 2.000000E+03`. It does not print `TWIST RATE`
directly. With `R1=0` at GRID 1, `R1=9.904762E-01` at GRID 2, and
`L_0=100`, the report value is:

```text
TWIST RATE = (R1_2 - R1_1) / L_0
           = 9.904762E-01 / 100
           = 9.904762E-03
```
