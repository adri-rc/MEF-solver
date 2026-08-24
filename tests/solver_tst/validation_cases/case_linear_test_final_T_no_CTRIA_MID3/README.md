# Caso lineal T sin rigidez física MID3 en CTRIA3

Este caso es el control A/B del
[`case_linear_test_final_T`](../case_linear_test_final_T/README.md). Conserva
su geometría, conectividad, propiedades estructurales, cargas y condiciones de
contorno, pero sustituye **sólo** el material de cortante transversal `MID3` de
los 128 CTRIA3 por un material auxiliar con `E=1e-20 MPa`. Su propósito es
aislar la discrepancia observada en la recuperación centroidal de momentos y
cortantes del CTRIA3; no representa una variante físicamente utilizable del
frame.

El nombre `no_CTRIA_MID3` es una abreviatura del experimento numérico. El
campo `MID3` de PSHELL no está vacío: apunta a `MAT1 99`, cuya rigidez es
despreciable. Se usa este sustituto porque reproduce la limitación ya
documentada en el proyecto sin cambiar el resto del modelo y permite que
ambos códigos recuperen los términos de lámina bajo la misma hipótesis casi
nula.

## Diferencia exacta frente al modelo físico

En el modelo T original, la propiedad de los triángulos es:

```text
PSHELL,11,1,2.0,1,,1,0.8333333333
```

En este control pasa a ser:

```text
PSHELL,11,1,2.0,1,,99,0.8333333333
MAT1,99,1.0E-20,,0.33
```

Por tanto, `MID1=1` y `MID2=1` conservan la membrana y la flexión del aluminio
estructural; únicamente `MID3` cambia de `MAT1 1` a `MAT1 99`. Los CQUAD4 de
PID 10 y 12 continúan usando `MID3=1`. No se modifican GRID, EID, PID,
espesores, orientación, coordenadas, cargas, SPC ni parámetros de solución.

## Modelo que se conserva

El deck sigue representando una sección de frame aeronáutico en T, en el
sistema coherente `mm - N - MPa` y con pequeñas deformaciones:

| Familia | Cantidad | PID | Función |
|---|---:|---:|---|
| CQUAD4 | 256 | 10 | Paños del web perforado, `t=2 mm` |
| CTRIA3 | 128 | 11 | Transición alrededor de los dos aligeramientos, `t=2 mm` |
| CQUAD4 | 32 | 12 | Ala de la T, `t=4 mm` |
| CBAR | 56 | 20 | Cap perimetral y montante central |
| CBAR | 20 | 21 | Refuerzos secundarios |
| CBAR | 6 | 22 | Mitades de los tres nervios próximas al web |
| CBAR | 6 | 23 | Mitades próximas al ala, con sección doble |
| CROD | 2 | 30 | Tirantes de introducción de carga |

Son 417 GRID y 506 elementos: 288 CQUAD4, 128 CTRIA3, 88 CBAR y 2
CROD. Se mantienen las dos fuerzas de `6000 N` en A, el momento de
`250000 Nmm` en B y el movimiento diferencial `T2=+0.2 mm`. Dos puntos del
ala están empotrados y otros dos representan remaches guiados con el
desplazamiento impuesto. Las guías de los extremos libres de los CROD también
son idénticas a las del caso original.

`PARAM,AUTOSPC,NO` permanece escrito expresamente. En consecuencia, ni MSC
Nastran ni el solver local deben ocultar mediante estabilización automática la
pérdida de rigidez transversal introducida deliberadamente.

## Ejecución de referencia MSC Nastran

El
[`tfm_final_airframe_frame_T_no_CTRIA_MID3.dat`](tfm_final_airframe_frame_T_no_CTRIA_MID3.dat)
se ejecutó con MSC Nastran Student Edition 2026.1 el 24 de agosto de 2026. El
[`F06`](tfm_final_airframe_frame_t_no_ctria_mid3.f06) alcanza `END OF JOB`, no
contiene mensajes `FATAL` y el lanzador termina con `NSEXIT: EXIT(0)`. MSC
informa:

- número de condición de la matriz de rigidez: `5.660208e28`;
- epsilon energético: `1.2989111e-13`;
- máximo módulo de traslación: `3.982957 mm` en GRID 434;
- máximo `|T3|`: `0.466976792 mm` en GRID 364.

El número de condición, casi veinte órdenes de magnitud mayor que el del
modelo físico (`3.172218e9`), es parte esencial del resultado: la terminación
normal de MSC no convierte este deck en una estructura bien condicionada.

| Fichero | SHA-256 |
|---|---|
| `tfm_final_airframe_frame_T_no_CTRIA_MID3.dat` | `9088444b921ee252f53310b76aa2e2bc9f9bfeab951c4dd63e216dfc803a789e` |
| `tfm_final_airframe_frame_t_no_ctria_mid3.f06` | `5d82a82c9a2bf9951e8b4f1211674f9e503726e8dab244eb9d06eadf600baec0` |
| `tfm_final_airframe_frame_t_no_ctria_mid3.op2` | `f52d9e98fe108fc60882c82d1efdd3501dbdeeda042a27c1d439d2476e8f4236` |

## Qué demuestra la comparación

La tabla siguiente usa campos completos del OP2 y el error relativo
`||resultado_MEF - resultado_MSC||2 / ||resultado_MSC||2`. El `core 7` de
CTRIA3 agrupa `Nx, Ny, Nxy` y las cuatro magnitudes de tensión de la fibra
positiva usadas en la validación. Los resultados impresos seleccionados se
transcriben como literales en
[`test_linear_final_t_validation_case.py`](../../test_linear_final_t_validation_case.py);
el test no lee el F06 ni el OP2 durante pytest.

| Magnitud | T con MID3 físico | T con `E(MID3)=1e-20` | Lectura |
|---|---:|---:|---|
| Traslaciones, vector completo | `8.206e-3` | `1.745e9` | El segundo valor está gobernado por el casi-mecanismo `T3` local |
| Rotaciones | `2.546e-2` | `5.079e-6` | Correlación muy estrecha en el control A/B |
| Reacciones de fuerza | `1.705e-4` | `2.298e-7` | Cierre global estrecho |
| Reacciones de momento | `3.553e-4` | `4.897e-8` | Cierre global estrecho |
| CTRIA3, membrana `N` | `7.218e-6` | `7.218e-6` | El cambio MID3 no altera la membrana |
| CTRIA3, momentos `M` | `6.820e-1` | `7.106e-5` | Mejora aproximada de `9.6e3` veces |
| CTRIA3, cortantes `Q` relativos | `1.730e1` | `4.612e12` | El cociente de la segunda columna no es interpretable porque MSC es casi cero |
| CTRIA3, `core 7` | `5.820e-3` | `6.880e-6` | Membrana y fibra positiva correlacionan |

En el caso casi nulo, las componentes de desplazamiento que siguen
controladas también coinciden: los errores L2 son `1.795e-7` en `T1`,
`3.233e-8` en `T2` y `5.079e-6` en las rotaciones. El máximo error absoluto
es `2.806e-7 mm` en T1 (GRID 164), `3.831e-7 mm` en T2 (GRID 364) y menor de
`1.40e-7 rad` en las tres rotaciones. En las seis estaciones de reacción, las
máximas diferencias son `4.692e-4 N` en F3 (GRID 430) y `1.213e-2 Nmm` en M1
(GRID 422).

Los dos CROD dan un error agregado `1.248e-10`; su máxima diferencia es
`7.618e-7 Nmm` en torsión y `9.522e-8 MPa` en tensión torsional. En los 15
CBAR seleccionados, los errores L2 de fuerzas y tensiones son
`2.237e-6/2.140e-6`, con máximos absolutos de `0.20584 Nmm` en `MA2` del EID
7056 y `3.608e-4 MPa` en `MINA` del mismo elemento. Los 15 CQUAD4 dan
`2.160e-5` en resultantes completas y `5.269e-6` en la fibra positiva; sus
máximos son `1.579e-3 N/mm` en `Qy` y `1.418e-4 MPa` en von Mises, ambos en
el EID 1151.

En los 15 CTRIA3, la membrana difiere como máximo `2.247e-4 N/mm` en `Ny` del
EID 5065. Los errores L2 por componente de momento son
`4.11e-5/4.64e-5/7.36e-5` para `Mx/My/Mxy`; la máxima diferencia seleccionada
es `1.455e-5 Nmm/mm` en `My` del EID 5112. La fibra positiva da
`3.675e-6`, con máximo `1.122e-4 MPa` en `smin` del EID 5111.

Para `Q`, la norma local del campo completo es aproximadamente
`4.002e-11 N/mm` y la de MSC `8.677e-24 N/mm`; la máxima diferencia absoluta
es sólo `1.867e-11 N/mm`, en `Qy` del EID 5084. En la selección del test el
máximo local baja a `1.166e-11 N/mm`, frente a valores MSC del orden de
`1e-25`. Dividir por una referencia casi nula produce el aparente error
`4.612e12`, pero no expresa una diferencia física relevante. Normalizada con
la escala de `Q` del caso T original, la diferencia es `6.817e-10`. La
evidencia correcta es que ambos recuperan cortante de ingeniería cero, no que
exista una correlación relativa de doce órdenes.

## Casi-mecanismo transversal y alcance de la evidencia

El solver local obtiene un máximo `|T3|=7.3586443229e10 mm` en GRID 203,
mientras que MSC conserva allí `T3=-0.342126399 mm` y mantiene todo su campo
transversal finito. El error relativo L2 de `T3` es `2.163e10`. Esta
divergencia no debe ocultarse ni aceptarse mediante una tolerancia grande. Al
retirar en la práctica la contribución MID3, la formulación local deja casi
sin control el desplazamiento transversal de los triángulos. Que MSC conserve
un campo finito es compatible con una estabilización o flexibilidad residual
distinta, pero el DAT y el F06 no permiten identificar de forma única ese
tratamiento interno. El enorme número de condición confirma que se está
interrogando un límite numérico, no una configuración estructural robusta.

Por ello, este caso aporta evidencia causal sobre la **recuperación de
resultados de elemento**: los momentos CTRIA3 pasan de una discrepancia del
orden de `68 %` a `7.1e-5`, y las membranas, fibras y restantes resultados
seleccionados permanecen correlacionados. No valida el desplazamiento `T3`, no
reemplaza el modelo con MID3 físico y no amplía el alcance declarado del
solver. La validación lineal sólo puede considerarse cerrada dentro de esas
limitaciones explícitas.

## Ejecución local

Desde la raíz del repositorio:

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_linear_test_final_T_no_CTRIA_MID3\tfm_final_airframe_frame_T_no_CTRIA_MID3.dat
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\test_linear_final_t_validation_case.py
```

Las capturas de geometría, apoyos, cargas, deformada y tensiones del modelo T
físico se conservan en el caso original. No se presentan como imágenes de
validación de esta variante porque la deformada `T3` local está dominada por
el casi-mecanismo deliberadamente introducido.
