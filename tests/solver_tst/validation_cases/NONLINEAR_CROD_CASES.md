# Campaña MSC Nastran SOL 400 para CROD

La campaña inicial contiene 20 casos, de `case_92` a `case_111`. Tras auditar
sus resultados se añadieron once controles focalizados, `case_112` a
`case_122`, para objetividad, tangente por diferencias centrales, torsión,
comportamiento próximo al punto límite y un límite lineal sensible a cambios
de geometría. Los casos 123--127 forman una campaña de escala estructural con
10--30 barras, propiedades mixtas, cargas múltiples y desplazamientos
impuestos; sus referencias MSC ya están auditadas e incorporadas al E2E.

## Objetivo de cada caso

| Caso | Modelo | Comprobación principal |
|---:|---|---|
| 92 | CROD axial, `LGDISP=-1` | Referencia SOL 400 sin no linealidad geométrica |
| 93 | Mismo CROD, carga pequeña, `LGDISP=1` | Límite de pequeñas deformaciones frente a 92 |
| 94 | CROD, tracción grande, 20 incrementos | Curva axial y estado final con `FIXED` |
| 95 | Réplica de 94 con `ADAPT` | Independencia del stepping, iteraciones y residuos |
| 96 | CROD, compresión finita | Signos, acortamiento, fuerza y tensión |
| 97 | CROD inclinado 3D con `GRID.CD`/`FORCE.CID` | Transformaciones antes de la iteración |
| 98 | Dos CROD iguales en serie | Ensamblaje y continuidad de fuerza interna |
| 99 | Serie con dos `MAT1` y dos `PROD` | Estado distinto por elemento y compatibilidad |
| 100 | Dos CROD coincidentes en paralelo | Reparto de fuerza no lineal |
| 101 | CROD axial con `LOAD` | Reducción de combinaciones antes de incrementar |
| 102 | V simétrica de dos CROD | Primera gran rotación y simetría |
| 103 | V asimétrica con carga combinada | Rotación, distinta sección y equilibrio 2D |
| 104 | V rebajada, `ADAPT`, diagnóstico avanzado | Convergencia cerca de pérdida de rigidez |
| 105 | Cercha triangular cerrada | Apoyo móvil, compatibilidad y reacción horizontal |
| 106 | Cuadro arriostrado en X | Ensamblaje 2D redundante y reparto por seis barras |
| 107 | Trípode tetraédrico | Equilibrio y rotación geométrica 3D |
| 108 | Pirámide de base cuadrada | Cuatro CROD, dos secciones y carga 3D |
| 109 | V rebajada con 5 incrementos | Sensibilidad de paso, pareja con 110 |
| 110 | V rebajada con 40 incrementos | Convergencia al refinar incrementos |
| 111 | Mástil espacial, 9 CROD, 2 materiales | Caso E2E final, `LOAD`, energía y stepping adaptativo |
| 112 | V simétrica estable a carga moderada | Referencia común para los casos 113--115 y 118 |
| 113 | Copia de 112 rotada rígidamente en 3D | Objetividad y transformaciones `GRID.CD`/`FORCE.CID` |
| 114 | Caso 112 con perturbación horizontal `+1000` | Rama positiva de una tangente secante central |
| 115 | Caso 112 con perturbación horizontal `-1000` | Rama negativa; simetría y tangente deformada |
| 116 | Un CROD a torsión pura con `J>0` | GDL rotacional, par interno y tensión torsional |
| 117 | Un CROD con tracción y torsión combinadas | Acoplamiento/independencia axial-torsional con longitud finita |
| 118 | Copia de 112 con `LGDISP=2` | Equivalencia con fuerzas no seguidoras |
| 119 | V simétrica a `3.8e6`, 60 incrementos | Estado estable refinado inmediatamente anterior al salto de 102 |
| 120 | V de 102 con `ADAPT` agresivo | Cambios de paso, cortes y salto; puede fallar de forma informativa |
| 121 | V simétrica con carga pequeña y `LGDISP=-1` | Límite lineal geométricamente sensible |
| 122 | Copia de 121 con `LGDISP=1` | Convergencia al límite lineal en una cercha que sí rota |
| 123 | Cubierta arqueada 2D, 10 CROD | Ensamblaje redundante, tres PROD y siete fuerzas nodales |
| 124 | Puente de celosía 2D, 20 CROD | Cinco PROD, cargas repartidas y apoyo con `T1=-4` impuesto |
| 125 | Cubierta espacial, 30 CROD | Gran desplazamiento 3D, cuatro apoyos y dos materiales |
| 126 | Marquesina espacial, 24 CROD | Fuerzas 3D combinadas y desplazamiento impuesto de un apoyo |
| 127 | Torre arriostrada 3D, 30 CROD | Seis PROD, dos materiales y stepping `ADAPT` |

## Auditoría completa

Resultados observados con **MSC Nastran Student Edition 2026.1**, motor
`NLPERF`, el 20 de agosto de 2026:

- los 30 trabajos convergidos 92--119/121--122 terminaron con `exit number 0`,
  alcanzaron factor de carga `1.0` y no contienen fatales ni valores no finitos;
- 92/93 coinciden en `u=1.0` y reacción `-2.1e5`; con `LGDISP=-1`, el caso 92
  no imprime la tabla `RODNL`, por lo que su control es desplazamiento/reacción;
- 94/95 coinciden exactamente en el estado final (`u=20`,
  `sigma=42000`, reacción `-4.2e6`); el 95 relanzado usa nueve incrementos
  adaptativos no uniformes;
- 104 y 110 coinciden en toda la curva impresa. El caso 109, con cinco pasos,
  llega al mismo desplazamiento, tensión y reacción, pero su energía/trabajo
  acumulado difiere un `0.65 %`, que es precisamente la sensibilidad buscada;
- 102 salta en el incremento final de `uy=-42.830` a `uy=-223.31`. No es una
  terminación errónea: convergió en 14 iteraciones con razón de residuo
  `4.66e-16`, pero pertenece a la rama posterior al *snap-through* y no debe
  usarse como primer caso de depuración del solver;
- 105 presenta otro salto convergido en el incremento 17, con factor de carga
  `0.48571`: `uy` pasa de `-33.975` a `-218.50` en 14 iteraciones. También es
  un control de *snap-through*, no un fallo del trabajo;
- los relanzamientos 95/104/109/110/111 ya no contienen el warning 602;
  95/104/111 usan 9/14/15 incrementos adaptativos reales;
- 112/113/118 coinciden al expresar el modelo rotado en sus ejes de salida;
  114/115 dan una flexibilidad tangente central `6.127186e-6`;
- 116/117 recuperan `R1=0.9904762`, par `40000` y tensión torsional `2000`;
- 119 converge antes del límite con `uy=-40.32496`; 120 termina anormalmente
  en `lambda=0.98387` tras alcanzar el paso mínimo, resultado que se conserva
  como referencia negativa;
- 121/122 dan `uy=-0.269374/-0.2699212`, una diferencia geométrica del `0.2031 %`.
- 123--127 terminaron con `exit number 0`, alcanzaron factor de carga `1.0` y
  reprodujeron `25/30/30/40/14` incrementos. La comparación exhaustiva de las
  tablas impresas obtuvo errores relativos máximos de `3.77e-7` en
  desplazamientos, `3.76e-7` en reacciones y `3.72e-7` en tensión/deformación;
- la secuencia adaptativa completa de 127 coincide desde `0.025` hasta `1.0`,
  incluidos los pasos crecientes y los cinco últimos incrementos de `0.1`;
- 124 y 126 confirman exactamente los desplazamientos impuestos finales
  `T1=-4.0` y `T1=+5.0`, junto con sus reacciones conjugadas;
- 125 presenta un *snap-through* de GRID 5 entre `lambda=0.533333` y
  `0.566667`, donde `T3` salta de `-13.5215` a `-69.40416`. El estado mecánico
  final coincide, pero `TotStrEn=4.7083e7` de MSC acumula trabajo a través del
  salto y no se compara con la energía elástica de estado `1.3756e7`.

No queda ningún trabajo pendiente de relanzamiento entre 92--127.

## Ejecución y archivos que deben conservarse

Ejecutar cada `.dat` desde su propia carpeta para que las salidas queden junto
al modelo. Conservar en el repositorio, como mínimo, el `.dat` y el `.f06`. El
F06 aporta los resultados impresos, el resumen `NLSTEP`, las iteraciones y la
versión de MSC. Conservar además el `.pch` sólo si se va a automatizar una
extracción de estados. Los `.op2`, `.f04`, `.itr`, `.log` y `.sts` son
artefactos de trabajo; para un caso fallido conviene mantener temporalmente
`.log` y `.sts` hasta diagnosticarlo, pero no son referencias permanentes por
defecto.

No renombrar la raíz del trabajo: los tests futuros asociarán las referencias
por el prefijo `ejNNN_...`.

En cada F06 registrar si MSC utilizó el motor clásico `CNLS` o el nuevo
`NLPERF`, además de la versión exacta. Si un caso no converge, no modificar
tolerancias ni carga todavía: devolver también el `.f06`, `.sts` y `.log` del
fallo para decidir si el deck o el fenómeno físico necesita otra estrategia.

## Datos intermedios prioritarios

- Casos 92/93: desplazamiento/fuerza final y diferencia entre ambos; la tabla
  de tensión `RODNL` sólo aparece en 93.
- Casos 94/95: todos los pares carga--desplazamiento y fuerza axial--incremento.
- Casos 104/109/110: tabla de `%` del F06, peor residuo, bisecciones y curva
  carga--desplazamiento del nodo 3.
- Caso 111: `GPFORCE`, energía, fuerzas de los nueve CROD y desplazamientos de
  nodos 5/6 en cada salida.
- Casos 112/113/118: igualdad de resultados al expresarlos en los ejes locales
  del modelo base.
- Casos 114/115: calcular `(u_plus-u_minus)/2000` y comprobar la media frente
  a 112 para obtener una referencia de flexibilidad tangente deformada.
- Casos 116/117: rotación `R1`, par, tensiones y deformaciones torsionales.
- Casos 119/120: curva carga--desplazamiento, número de ciclos, cortes y último
  residuo antes/después del punto límite.
- Casos 121/122: diferencia relativa de desplazamiento y reacciones para el
  límite de pequeñas rotaciones de una cercha transversalmente cargada.
- Casos 123--127: desplazamientos de todos los nodos libres, reacciones de
  todos los apoyos y tensiones de al menos un elemento por cada PROD; en
  124/126 comprobar además que el desplazamiento impuesto final sea exactamente
  `-4.0/+5.0` y registrar su reacción conjugada; en 127 registrar todos los
  factores de carga adaptativos.

La campaña no solicita una DMIG tangente por iteración: MSC Nastran no ofrece
una tarjeta estándar equivalente a la exportación lineal usada para las
láminas. Si el balance anterior no identifica la formulación, se añadirá un
deck DMAP específico sin alterar estos veinte controles.
