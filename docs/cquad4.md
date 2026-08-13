# Desarrollo y validación de CQUAD4

## Estado actual

El elemento se mantiene en dos implementaciones independientes:

- `CQUAD4NastranElementSolver`: reconstrucción de la rigidez de Nastran y
  formulación activa en `LinearStaticSolver`;
- `CQUAD4ReissnerMindlinElementSolver`: formulación isoparamétrica clásica de
  Reissner-Mindlin con integración completa `2 x 2`.

`CQUAD4ElementBase` contiene únicamente la geometría, el sistema local, las
matrices constitutivas, `K6ROT` y la recuperación de tensiones comunes. Las dos
clases producen matrices `24 x 24`, pero no comparten el cálculo de rigidez.

## Formulación Reissner-Mindlin

La implementación clásica integra con una regla `2 x 2` los bloques de
membrana, flexión y cortante transversal. Es una referencia teórica útil y no
pretende reproducir las decisiones internas de Nastran. Como es habitual en un
elemento de bajo orden, puede presentar bloqueo de cortante en placas delgadas.

## Reconstrucción Nastran

Para elementos planos, la ruta Nastran combina:

1. deformación uniforme evaluada en el centro;
2. dos modos bilineales ortogonales a traslaciones y campos lineales;
3. estabilización *hourglass* basada en las dos líneas medias del
   cuadrilátero, con coeficiente isotrópico `E/12`;
4. el mismo operador para membrana y curvaturas, con el mapeo
   `p = [theta_y, -theta_x]`;
5. el bloque transversal de MacNeal construido en flexibilidad aumentada:
   cuatro deformaciones de cortante selectivas se combinan con la
   flexibilidad residual de flexión y, si existe `MID3`, con su flexibilidad
   física de cortante;
6. penalización *drilling* opcional cuando `K6ROT > 0`.

En elementos alabeados se forma primero la rigidez sobre el plano que pasa por
los puntos medios de las aristas. Después se aplica una transformación de
equilibrio que traslada los resultantes a los nodos reales:

- las fuerzas de membrana se resuelven sobre las dos aristas incidentes y
  generan fuerzas transversales autoequilibradas (*kick-off forces*);
- los momentos se hacen tangentes al plano nodal y su componente normal se
  equilibra mediante fuerzas en el plano.

Si `W` transforma los resultantes del plano medio, la operación implementada es:

```text
K_warped = W K_mean W^T
K_global = T^T K_warped T
```

## Validación frente a DMIG

La prueba `test_E2E_CQUAD_NAstran_K_solution.py` lee cada deck, ensambla la
matriz y la compara con la DMIG `KAAX` exportada a CSV. Los casos planos comparan sus GDL
físicamente activos con tolerancia del `1 %`; el caso 54 compara los 24 GDL con
tolerancia del `3 %`, porque incluye la reconstrucción pública de alabeo y no
el código propietario de la versión comercial usada para generar la DMIG.

Error relativo de Frobenius observado:

| Caso | Geometría y bloque | Nastran | Reissner-Mindlin |
|---|---|---:|---:|
| 51 | distorsionada, membrana | `0.152382 %` | `12.557935 %` |
| 52 | distorsionada, flexión | `0.152382 %` | `12.557935 %` |
| 53 | plana sesgada, membrana + flexión | `0.372560 %` | `25.279406 %` |
| 54 | inclinada y alabeada, matriz `24 x 24` | `2.113059 %` | `45.074180 %` |
| 55 | inclinada 3D, `K6ROT=100`, matriz `24 x 24` | `1.077026 %` | -- |
| 65 | cortante aislado con `MID3`, `K6ROT=0` | `5.74e-9 %` | -- |
| 67 | MID1/MID2/MID3, `K6ROT=100` | `0.000298 %` | -- |
| 69-71, 73 | pares MID3/G3/TST rectangulares | `< 7.3e-9 %` | -- |
| 72 | MID3 vacío, `K6ROT=100` | `0.000121 %` | -- |
| 74-75 | cuadrilátero plano sesgado, con/sin MID3 | `0.72145 %` | -- |

Sin la corrección de alabeo, el error Nastran del caso 54 era `38.551436 %`.
La corrección lo reduce a `2.113059 %`; por tanto, el caso queda cubierto por el
criterio de aceptación del `3 %`, aunque no existe igualdad término a término.

### Flexibilidad transversal validada por los casos 65-75

La expresión implementada para el subespacio `[T3,R1,R2]` es:

```text
K_s = D_s^T (Z_b + Z_s)^-1 D_s
```

`D_s` contiene cuatro deformaciones selectivas: `gamma_xz` en
`(0,+/-1/sqrt(3))` y `gamma_yz` en `(+/-1/sqrt(3),0)`. `Z_b` es la
flexibilidad residual de flexión QUAD4 con corrección de aspecto
`epsilon=0.025`, identificada con las DMIG de MSC. Cuando `MID3` existe se
añade diagonalmente:

```text
Z_s(i,i) = 1 / (2 detJ_i G3 t (TS/T))
```

Por tanto, `MID3` aporta **flexibilidad en serie**, no una rigidez
Reissner-Mindlin en paralelo: reducir `G3` a la mitad duplica esa corrección y
reducir `TS/T` de `5/6` a `1/2` la multiplica por `5/3`. Los casos 70, 71 y 73
verifican independientemente estas escalas.

Un `MID3` vacío no hereda `MID2`. En ese caso `Z_s=0`, pero permanece `Z_b`;
por eso el elemento conserva la rigidez residual de placa observada en los
casos 69 y 72. Esta semántica coincide con el código público NASTRAN-95, que
activa `SHRFLX` únicamente si `MID3 > 0`.

Los casos rectangulares 65, 69-71 y 73 coinciden término a término al orden de
`1e-10`; los casos completos 67 y 72 quedan por debajo de `3e-6` de error
relativo. Los F06 66 y 68 validan además desplazamientos, cargas, reacciones,
fuerzas, momentos y tensiones. En 74-75 queda un error global del `0.72 %`,
atribuido a la reconstrucción geométrica pública del cuadrilátero irregular y
no a la ley de `MID3`, porque el error es prácticamente igual con y sin MID3.

La comparación puede reproducirse con:

```powershell
.\.venv\Scripts\python.exe tests/solver_tst/elements/run_cquad4_nastran_comparison.py
```

La opción `--export-dir ruta` guarda las matrices calculadas sin modificar las
referencias de validación.

## Fuentes y trazabilidad

- R. H. MacNeal, [A simple quadrilateral shell element](https://www.sciencedirect.com/science/article/pii/0045794978900202), *Computers & Structures*, 8(2), 175-183, 1978, DOI `10.1016/0045-7949(78)90020-2`.
- R. H. MacNeal, [Derivation of element stiffness matrices by assumed strain distributions](https://www.sciencedirect.com/science/article/pii/002954938290262X), *Nuclear Engineering and Design*, 70(1), 3-12, 1982, DOI `10.1016/0029-5493(82)90262-X`.
- NASA, [nasa/NASTRAN-95](https://github.com/nasa/NASTRAN-95): código histórico público; para CQUAD4 son relevantes `mis/quad4s.f`, `mis/quad4d.f`, `mis/q4bmgs.f`, `mis/q4bmgd.f`, `mis/trplms.f` y `mis/trplmd.f`.
- NASA NTRS, [Benchmarking the QUAD4/TRIA3 element](https://ntrs.nasa.gov/citations/19940013364): documenta diferencias entre variantes Nastran en endurecimiento de cortante e integración numérica.
- NASA NTRS, [informe 19900011106](https://ntrs.nasa.gov/citations/19900011106): describe el cortante transversal del MSC/NASTRAN QUAD4 mediante aproximación *string-net* y flexibilidad aumentada.
- B. P. Naganarayana y G. Prathap, [Force and moment corrections for the warped four-node quadrilateral plane shell element](https://repository.ias.ac.in/39268/1/138_PUB.pdf), *Computers & Structures*, 33(4), 1107-1115, 1989.

Las dos publicaciones de MacNeal describen la base conceptual de QUAD4 y las
deformaciones asumidas. NASTRAN-95 permite rastrear una implementación pública
histórica, pero no es el código de MSC Nastran actual. La ruta `nastran` de este
proyecto es, por ello, una reconstrucción contrastada contra DMIG, no una copia
literal de una formulación comercial propietaria.

## Validación pendiente

- casos alabeados adicionales para separar error de forma y error de alabeo;
- reproducir con mayor fidelidad la geometría irregular de QUAD4 para reducir
  el `0.72 %` residual de los casos 74-75;
- materiales ortótropos, laminados, `MID4`, `ZOFFS` y espesores nodales;
- *patch tests* multielemento y estudios de convergencia de placa;
- estudio de sensibilidad a aspecto, *skew*, *taper* y espesor.
