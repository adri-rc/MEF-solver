# Guía del informe de resultados `.rpt`

## Propósito y orden del informe

`src/postprocessing/report_writer.py` genera `output/<job>.rpt`. El informe
presenta magnitudes calculadas por MEF-solver y facilita su localización en un
F06 de MSC Nastran, pero no intenta reproducir literalmente su maquetación.
Todas las magnitudes usan el sistema coherente de unidades del deck.

El orden es:

1. cabecero lineal o no lineal;
2. `MODEL SUMMARY` y configuración;
3. para SOL 400, control y resumen de convergencia;
4. resultados del estado final;
5. al final, historial de incrementos y correcciones de Newton.

## `MODEL SUMMARY`

La tabla cuenta entradas lógicas del Bulk Data de origen, igual que la tabla
`M O D E L  S U M M A R Y` del F06. No cuenta los objetos normalizados que usa
internamente el solver.

Esto es especialmente importante para las restricciones:

- `SPC` cuenta entradas `SPC` con desplazamiento impuesto;
- `SPC1` cuenta combinaciones lógicas de `SID` y componentes; líneas `SPC1`
  repetidas con ambos campos iguales se consolidan como hace MSC;
- `SPCADD` se mantiene como tipo independiente si existe en el deck;
- una restricción permanente en `GRID.PS` pertenece a la entrada `GRID` y no
  crea una entrada `SPC` adicional en el resumen;
- `FORCE` y `MOMENT` son cargas, aparecen como tipos independientes y nunca se
  incluyen en el número de `SPC`.

Por ejemplo, el caso 127 tiene doce restricciones nodales una vez expandido
para resolver, pero el F06 y el `.rpt` muestran `SPC1 = 2`: existen dos grupos
lógicos, `123456` y `456`, para el mismo SID.

## Configuración y control SOL 400

| Bloque del `.rpt` | Dónde buscar en el F06 |
|---|---|
| `ANALYSIS CONFIGURATION` | `CASE CONTROL ECHO` para `ANALYSIS`, `SPC`, `LOAD` y `NLSTEP`; el Executive Control muestra `SOL 101` o `SOL 400`. |
| `PARAM,LGDISP` y `PARAM,FOLLOWK` | Bulk Data o eco del deck; también se contabilizan como `PARAM` en `MODEL SUMMARY`. |
| `NONLINEAR STEP CONTROL` | `NON-LINEAR ITERATION SOLUTION CONTROL PARAMETERS`. MSC usa nombres como `TOTTIME`, `MAXITER`, `MAXBIS`, `DTINITF`, `DTMINF`, `DTMAXF`, `NDESIR`, `SFACT`, `EPSP` y `KMETHOD`. |
| `NONLINEAR SOLUTION SUMMARY` | Se deduce del último `LOAD STEP`, de los mensajes de convergencia, de la tabla de iteraciones y de `GLOBAL WORKS & ENERGIES`; MSC no reúne exactamente estos cinco valores en una sola tabla. |

`FINAL LOAD FACTOR = 1` significa que se alcanzó toda la carga de referencia.
No significa que las fuerzas tengan valor uno: el solver aplica
`F_ext(lambda) = lambda F_ref`.

## Resultados finales

`FINAL CONVERGED STATE RESULTS` identifica el último equilibrio aceptado. En
un F06 con salida intermedia pueden aparecer las mismas tablas muchas veces;
para comparar el estado final hay que localizar su última aparición bajo el
último `LOAD STEP` alcanzado.

| Resultado del `.rpt` | Encabezado habitual del F06 SOL 400 |
|---|---|
| `NODAL DISPLACEMENTS` | `D I S P L A C E M E N T  V E C T O R`. |
| `SPC REACTION FORCES AND MOMENTS` | `F O R C E S  O F  S I N G L E - P O I N T  C O N S T R A I N T`. |
| tensión/deformación CROD | `N O N L I N E A R  S T R E S S E S  I N  R O D  E L E M E N T S (CROD)`. |
| `STRAIN ENERGY` | `GLOBAL WORKS & ENERGIES`, campo `TOT STRAIN ENERGY`; en saltos de rama MSC puede acumular trabajo y no representar la energía elástica instantánea usada por el proyecto. |
| fuerza axial y par CROD | No siempre aparece una tabla OEF legible en los F06 de referencia SOL 400. Se contrasta con tensión, cargas `OLOAD`, reacciones `SPCFORCES` y, cuando está disponible, PCH/OP2. |

## `CONVERGED INCREMENT HISTORY`

Una fila representa un incremento que terminó en equilibrio y fue aceptado.

| Columna | Significado | Correspondencia F06 |
|---|---|---|
| `INCREMENT` | Número consecutivo de incremento aceptado por MEF-solver. | Columna `INC` de `NON-LINEAR ITERATION MODULE OUTPUT`, con posibles diferencias si MSC realiza una adaptación distinta. |
| `LOAD FACTOR` | Factor total `lambda` alcanzado. | Primer campo `LOAD STEP` de la tabla iterativa y línea `LOAD STEP = ...` de cada bloque de resultados. |
| `STEP SIZE` | Diferencia entre el factor actual y el anterior. | Mensajes `TIME INCREMENT HAS BEEN CHANGED` y diferencia entre dos `LOAD STEP`. |
| `ITERATIONS` | Correcciones de Newton requeridas por el intento aceptado. | Número de filas `ITR` del incremento que finalmente converge. |
| `RESIDUAL NORM` | Norma infinito final del residuo en GDL libres. | Conceptualmente `MAXIMUM RESIDUAL FORCE`; no se espera identidad exacta si MSC usa otro sistema o escalado interno. |
| `RELATIVE RESIDUAL` | Residuo anterior dividido por la referencia definida por MEF-solver. | Análogo a `FORCE RESIDUAL CONVERGENCE RATIO`, pero no numéricamente idéntico porque MSC no publica aquí el mismo denominador. |
| `CUTBACKS` | Intentos fallidos inmediatamente anteriores que obligaron a reducir ese paso. | Mensajes de bisección/cambio de incremento y columnas de control de la tabla iterativa. |
| `STRAIN ENERGY` | Energía elástica del estado aceptado. | `TOT STRAIN ENERGY`, con la salvedad indicada para saltos de rama. |

## `NEWTON ITERATION HISTORY`

Para un factor de carga objetivo, Newton-Raphson repite:

```text
r_i = lambda F_ext - F_int(u_i)
K_T(u_i) Delta_u_i = r_i
u_(i+1) = u_i + Delta_u_i
```

Cada fila del historial se escribe después de aplicar una corrección
`Delta_u_i`.

| Columna | Significado | Correspondencia F06 |
|---|---|---|
| `INCREMENT` | Próximo incremento convergido que se está intentando. Si hay recorte, el mismo número reaparece. | `INC`, aunque la numeración puede diferir por la estrategia adaptativa de MSC. |
| `LOAD FACTOR` | Factor objetivo del intento. Puede cambiar al repetir un incremento tras un recorte. | Campo `LOAD STEP` de la fila que comienza por `%`. |
| `ITERATION` | Número de corrección de Newton dentro de ese intento. | Columna `ITR`. |
| `RESIDUAL NORM` | `max(abs(r))` sobre los GDL libres después de corregir. Incluye fuerzas y momentos según el GDL. | `MAXIMUM RESIDUAL FORCE ... IS EQUAL TO ...`, como magnitud conceptualmente equivalente. |
| `RELATIVE RESIDUAL` | `RESIDUAL NORM / max(norm(lambda F_ext), norm(F_int), 1)` usando norma infinito y GDL libres. Es el criterio `P` del proyecto. | `FORCE RESIDUAL CONVERGENCE RATIO` es el análogo MSC, pero emplea una referencia interna distinta y no debe compararse término a término. |
| `CORRECTION NORM` | `max(abs(Delta_u))` sobre los GDL libres. Mide cuánto cambió la solución en esa corrección. | No tiene una columna directamente equivalente en el F06 usado para validar. `DISP MAX` de MSC es el máximo desplazamiento total del estado, no esta corrección. |

### Por qué cambia el número de filas

No existe un número fijo de filas por caso:

- `FIXED` fija el número nominal de incrementos, no las iteraciones necesarias
  dentro de cada uno;
- `ADAPT` cambia además el número y tamaño de incrementos según
  `DESIRED ITERATIONS` y `STEP SCALE FACTOR`;
- una zona más no lineal suele requerir más correcciones;
- un recorte conserva las filas del intento fallido, reduce el paso y repite
  el mismo número de incremento con otro factor objetivo;
- un estado que ya satisface la tolerancia antes de corregir puede converger
  con cero filas de corrección;
- el proceso se limita mediante `MINIMUM ITERATIONS`, `MAXIMUM ITERATIONS` y
  `MAXIMUM BISECTIONS`.

Por eso `RECORDED NEWTON CORRECTIONS` puede superar la suma de
`CONVERGED NEWTON ITERATIONS`: la primera incluye intentos descartados por
recortes y la segunda sólo los intentos finalmente aceptados.

## Torsión del CROD no lineal

Sí está implementada y se aceptan cargas `MOMENT`. Requiere `G` en `MAT1` y
una constante torsional `J > 0` en `PROD`. La formulación actual usa:

```text
phi = (theta_2 - theta_1) dot n_0
twist_rate = phi / L_0
T = G J twist_rate
tau = T C / J
k_t = G J / L_0
```

El `.rpt` muestra `TORQUE`, `TORSIONAL STRESS` y `TWIST RATE`. En el F06 del
caso 116, por ejemplo, la tensión aparece como `LIN. TORSIONAL STRESS`, el giro
en `R1`, el momento aplicado en `OLOAD` y la reacción torsional en
`SPCFORCES`.

La capacidad tiene un límite deliberado: la torsión es elástica, lineal y
desacoplada del axil; se mide alrededor del eje de referencia `n_0`, y una
extensión axial no modifica `GJ/L_0`. No implementa cinemática general de
rotaciones finitas ni acoplamiento flexión-torsión; esa cinemática pertenece al
CBAR no lineal. Si `J = 0`, el CROD no aporta rigidez torsional y sus rotaciones
axiales deben quedar restringidas para evitar un mecanismo.

### Ejemplo de lectura del caso 117

En el último `LOAD STEP = 1.00000E+00`, la página 163 del F06 muestra
`TOTAL STRAIN = 2.000000E-01`; ese es el `AXIAL STRAIN` del `.rpt`. La misma
tabla muestra `LIN. TORSIONAL STRESS = 2.000000E+03`.

MSC no imprime directamente `TWIST RATE`. En la página 168, el
`DISPLACEMENT VECTOR` da `R1=0` para GRID 1 y `R1=9.904762E-01` para GRID 2.
Como el deck define `L_0=100`, el valor del `.rpt` se reconstruye mediante:

```text
TWIST RATE = (R1_2 - R1_1) / L_0
           = 9.904762E-01 / 100
           = 9.904762E-03
```
