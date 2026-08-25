# Análisis de resultados y correlación con MSC Nastran

## Referencia y método

La referencia independiente es el último bloque `LOAD STEP = 1.00000E+00` de
`tfm_final_nonlinear_bridge_updated.f06`, generado con MSC Nastran Student
Edition 2026.1 el 25 de agosto de 2026. Los desplazamientos y reacciones se
contrastaron además con el bloque final del PCH. El modelo resuelto por MSC y
por MEF-solver es el mismo DAT, sin pasar por una regeneración de Patran.

Para un vector o familia completa se usa el error relativo L2

```text
epsilon_L2 = ||x_local - x_MSC||_2 / ||x_MSC||_2
```

Los errores relativos componente a componente sólo se interpretan cuando la
referencia tiene magnitud suficiente. Cerca de cero se informa también la
diferencia absoluta, porque un cociente porcentual puede resultar engañoso.

Las tensiones axiales centroidales `CBAR` se reconstruyen, como en el resto de
la campaña mixta, mediante la media de las tensiones C/D/E/F de ambos extremos.
Las tensiones de fibra se comparan directamente punto por punto. La longitud
actual se calcula con las coordenadas originales y el `DISPLACEMENT VECTOR`
final de MSC. Las referencias del test son literales escritos en el código: el
E2E no abre el F06 ni el PCH en tiempo de ejecución.

## Convergencia

| Magnitud | MSC Nastran | MEF-solver |
|---|---:|---:|
| Factor de carga final | 1.0 | 1.0 |
| Incrementos convergidos | 20 | 20 |
| Iteraciones de Newton | 40 | 40 |
| Recortes/bisecciones | 0 | 0 |
| Razón residual final MSC | 6.369e-8 | — |
| Energía de deformación [Nmm] | 1.092467e7 | 1.0924822116e7 |

MSC imprime dos iteraciones en cada incremento, `JOB CONVERGES FOR THE CURRENT
STEP` y `END OF JOB`. No hay mensajes `USER FATAL`, `SYSTEM FATAL`,
`USER WARNING` ni `SYSTEM WARNING`. El residuo absoluto final de MSC es
`2.919e-2 N`; la razón normalizada queda por debajo del criterio `3e-7`.
La energía difiere un `1.392e-5` relativo. Los valores de residuo se conservan
con la precisión visible en el F06 final.

## Respuesta global

La comparación completa, no sólo los puntos incluidos en el E2E, da:

| Familia completa | Error relativo L2 |
|---|---:|
| Traslaciones, 202 GRID | 3.280e-5 |
| Rotaciones, 202 GRID | 2.416e-5 |
| Reacciones de fuerza | 2.898e-5 |
| Reacciones de momento | 4.250e-5 |
| Energía de deformación | 1.392e-5 |

La máxima norma de traslación aparece en `GRID 1127`:

| Código | T1 [mm] | T2 [mm] | T3 [mm] | Norma [mm] |
|---|---:|---:|---:|---:|
| MSC | 1.36855459 | 1.83992660 | -15.23815823 | 15.40972867 |
| MEF-solver | 1.36860658 | 1.84006768 | -15.23830206 | 15.40989236 |

El giro de mayor magnitud de MSC es `R1=0.04409650 rad` en `GRID 1015`.
Los asientos `-2/-6/-10 mm` se recuperan exactamente en los GRID de las tres
bases afectadas. No se detecta deriva de los valores impuestos.

## Resultados CROD

La correlación sobre las tensiones axiales de los 428 `CROD` da
`epsilon_L2=2.125e-5`. La máxima diferencia absoluta es `0.006625 MPa` en
`EID 20610`. El peor porcentaje bruto, `0.780 %`, aparece en `EID 20122`, cuya
referencia es casi nula (`-0.01918 MPa`); por ello no gobierna la valoración
física. Limitando la lectura a `|sigma_MSC|>1 MPa`, el máximo error relativo es
`0.0537 %`.

| Resultado CROD completo | Error relativo L2 |
|---|---:|
| Tensión axial | 2.125056e-5 |
| Deformación axial | 2.125684e-5 |
| Fuerza axial | 2.123120e-5 |
| Longitud actual | 8.367439e-9 |

En `EID 20610`, que gobierna las diferencias absolutas, la fuerza reconstruida
es `-67029.495 N` en MSC y `-67033.8012 N` localmente; las longitudes son
`3818.973969` y `3818.973840 mm`.

La máxima tensión axial en valor absoluto es `+220.4665 MPa` en la péndola de
red `EID 20929`, `PID 207`. El signo indica tracción. El material del modelo es
elástico sin límite de fluencia: el valor sirve para comparar códigos, no para
aceptar un acero o dimensionar la péndola. Otras péndolas pueden aparecer en
compresión porque `CROD` es bilateral y no implementa destensado.

## Resultados CBAR

La tensión axial centroidal de los 144 `CBAR` presenta
`epsilon_L2=1.379e-5`; la mayor diferencia absoluta es `0.000661 MPa` en
`EID 12111`. Por tanto, la parte axial y la transferencia global de esfuerzos
correlacionan al mismo nivel que los desplazamientos.

| Resultado CBAR completo | Error relativo L2 |
|---|---:|
| Tensión axial centroidal | 1.379175e-5 |
| Deformación axial | 1.379135e-5 |
| Fuerza axial | 1.627389e-5 |
| Longitud actual | 1.047300e-9 |

En `EID 12111`, la fuerza axial es `-193312.564 N` en MSC frente a
`-193323.145 N` localmente; la diferencia de longitud es sólo
`-4.885e-6 mm`.

La recuperación C/D/E/F es más sensible. Comparando las 1152 fibras de ambos
extremos:

| Grupo CBAR | Error relativo L2 de fibras |
|---|---:|
| Todos los CBAR | 1.673e-2 |
| Largueros de tablero, PID 101 | 4.212e-2 |
| Arcos, PID 103 | 8.779e-4 |
| Columnas, PID 105 | 2.294e-4 |
| Patas de cimentación, PID 106 | 2.171e-4 |

La mayor diferencia absoluta es `2.529 MPa` en `EID 10014`, extremo B,
fibra E/P2. Entre fibras cuya referencia supera `10 MPa`, el peor porcentaje
es `11.57 %` en `EID 10015`, extremo B, fibra D/P1: MSC imprime
`19.48201 MPa` y MEF-solver recupera `21.73605 MPa`.

La discrepancia no se observa en los axiles centroidales y se concentra en los
largueros continuos de tablero, donde los momentos de extremo y la convención
corrotacional de recuperación de fibra tienen mayor influencia. No se ha
recalibrado la rigidez ni se ha relajado la formulación para forzar ese ajuste.
El E2E usa una tolerancia separada para las rutas `stress_points_*`: de este
modo, la envolvente de fibra no debilita las comprobaciones de desplazamientos,
reacciones, energía, fuerzas axiales ni estados `CROD`.

La tensión normal de fibra máxima del modelo MSC es `121.1869 MPa` en
`EID 12110`, extremo asociado a `GRID 3021`, fibra F. Es inferior al máximo
axial de las péndolas, aunque ambos valores siguen perteneciendo a una
comparación elástica y no a una comprobación resistente.

## Señal de no linealidad geométrica

Se resolvió localmente un control idéntico cambiando sólo la evaluación a
`LGDISP=-1`, sin modificar el DAT conservado. En `GRID 1015/T2`:

| Evaluación | T2 [mm] |
|---|---:|
| `LGDISP=1` | 2.20017237795 |
| `LGDISP=-1` | 2.15731057757 |

El cambio es `+1.9868164 %`. La norma máxima pasa de `15.40989236 mm` a
`15.45020968 mm`, y el error L2 entre todas las traslaciones libres de ambos
controles es `0.3092466 %`. La geometría actualizada ejercita así una señal no
lineal medible sin acercarse a un punto límite ni necesitar recortes.

## Selección incorporada a la regresión

La entrada de `test_nonlinear_mixed_CROD_CBAR_validation_E2E.py` conserva una
muestra distribuida entre los tres vanos, ambas caras, pilas y cimentaciones:

- GRID de tablero, arco, línea interior y base, incluidos desplazamientos y
  reacciones;
- `CBAR` de largueros, arcos, columnas y patas de cimentación;
- `CROD` de largueros interiores, travesaños, péndolas, arriostramientos,
  vigas de base y péndolas de red;
- energía final y 20 incrementos;
- las fibras `EID 10014/10015` que gobiernan las diferencias de recuperación.

En conjunto quedan escritos 230 valores numéricos de referencia: 54
componentes de desplazamiento, 24 de reacción, 76 resultados de 15 CBAR
(incluidos 16 puntos de fibra), 75 resultados de 15 CROD y la energía. Se
comprueban además 30 tipos de elemento, el factor final y los 20 incrementos.
El test no obtiene ninguno de ellos leyendo archivos externos.

La tolerancia común se mantiene independiente de la tolerancia de fibra. Esto
evita ocultar una regresión global bajo la envolvente del `11.57 %` que sólo
corresponde a una magnitud recuperada concreta.

Las envolventes de la entrada quedan separadas por familia: `0.1 %` para
desplazamientos y energía, `0.6 %` para reacciones (gobierna una reacción
vertical de sólo `236.6 N`), `1 %` para estados centroidales de elemento
(gobierna el `CROD 20122` casi descargado) y `12 %` exclusivamente para las
fibras dominantes `stress_points_*`. `abs_tol=1e-5` no encubre ninguna de esas
componentes gobernantes.

## Evidencia visual de postproceso

Patran 2026.1 importó el DAT y los 20 estados del OP2. Se conservaron seis
capturas reales: geometría, apoyos/asientos, tráfico, viento,
frenado/momentos y configuración final deformada. La deformada superpone la
geometría original azul y la respuesta blanca con amplificación gráfica:

![Configuración final deformada en Patran](tfm_final_nonlinear_bridge_updated_patran_deformed_final.png)

El conjunto completo y la explicación de qué SPC se muestra en cada vista se
encuentran en la sección [Vistas de Patran](README.md#vistas-de-patran) del
README del caso.

No se obtuvo un fringe de tensiones porque el registro que Patran construye
desde el OP2 sólo ofrece cargas, desplazamientos y fuerzas de restricción. El
OP2 consolidado carece de tablas `CROD/CBAR` de fuerza, tensión o deformación,
aunque esas tablas sí están impresas en el F06 y sustentan la comparación
anterior. Se ha preferido documentar esa limitación de exportación antes que
representar un mapa que no proceda de los resultados reales de elemento.

## Observaciones y límites de interpretación

- Las cuatro bases parten de cotas diferentes; representan alturas de
  estribos/pilas y terreno escalonado, no un único plano horizontal de suelo.
- Arco, tablero y pila comparten nudo en los arranques; no hay aparatos de
  apoyo, liberaciones ni contacto unilateral.
- Las líneas interiores, péndolas y arriostramientos `CROD` no aportan flexión
  ni pandean. Sus rotaciones inactivas se restringen explícitamente.
- No se incluyen peso propio automático, cargas repartidas, plasticidad,
  fatiga, pandeo, cable flojo, contacto ni interacción suelo-estructura.
- Las capturas de tensiones de Patran deben separar `CBAR` y `CROD`; no existe
  un campo superficial continuo de von Mises en una malla exclusivamente 1D.

Con estas salvedades, el caso es adecuado para cerrar una validación mixta de
ensamblaje, cinemática geométricamente no lineal y recuperación axial a escala
de aproximadamente 200 nodos. La recuperación de fibra del tablero queda
cuantificada y vigilada explícitamente, no presentada como coincidencia exacta.
