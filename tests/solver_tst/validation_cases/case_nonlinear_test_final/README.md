# Caso no lineal final: puente arco atirantado espacial

Este caso es el modelo mixto `SOL 400` para la defensa del TFM. Representa un
puente arco atirantado de `12 000 mm` de luz, `2 400 mm` de ancho y
`2 400 mm` de flecha. No pretende dimensionar un puente real: es un modelo de
validación geométrica no lineal, elástico y reproducible dentro del alcance
actual de MEF-solver.

El mismo DAT se ejecutó sin modificar con MSC Nastran Student Edition
2026.1/NLPERF el 21 de agosto de 2026. Alcanzó los 50 incrementos, declaró
convergencia en `LOAD STEP=1.0`, terminó con `END OF JOB`/`EXIT(0)` y no
contiene `FATAL`. El [F06 auditado](tfm_final_nonlinear_bridge.f06) y el PCH se
conservan junto al deck; las referencias seleccionadas están transcritas en
el E2E mixto, que no lee resultados externos durante su ejecución.

## Modelo estructural

El sistema coherente de unidades es `mm - N - MPa`. Contiene 32 GRID y 78
elementos:

| Parte | Elementos | Propiedad |
|---|---:|---|
| Dos vigas longitudinales del tablero | 16 CBAR | PBAR 101, acero, `A=12 000 mm2` |
| Vigas transversales del tablero | 9 CBAR | PBAR 102, `A=6 000 mm2` |
| Dos arcos parabólicos | 16 CBAR | PBAR 103, `A=9 000 mm2` |
| Riostras transversales superiores | 7 CBAR | PBAR 104, `A=3 500 mm2` |
| Péndolas verticales | 14 CROD | PROD 201, `A=900 mm2` |
| Arriostramiento inferior alterno | 8 CROD | PROD 202, `A=650 mm2` |
| Arriostramiento entre arcos | 8 CROD | PROD 203, `A=500 mm2` |

Los CBAR usan `PBARN1,SECT=S,BEH2=BAR,INT2=LC`, que selecciona la rama
corrotacional validada. Los CROD son barras axiales bilaterales: las péndolas
no son cables *tension-only*, no pierden tensión y no incluyen pandeo.

## Acciones y apoyos

- El estribo izquierdo, GRID `101/201`, queda empotrado mediante `SPC1 10` en
  `123456`.
- Los apoyos derechos, GRID `109/209`, conservan libre el movimiento
  longitudinal y fijan `T2`.
- `SPC 30` prescribe asientos diferenciales `T3=-12 mm` y `T3=-6 mm` en esos
  apoyos. MEF-solver los incrementa con el factor de carga.
- `LOAD 501` contiene siete fuerzas verticales asimétricas sobre el tablero.
- `LOAD 502` contiene seis fuerzas laterales sobre los arcos y entra en la
  combinación con factor `0.8`.
- `LOAD 503` contiene dos fuerzas longitudinales de frenado y tres momentos
  nodales alrededor de los ejes globales X/Y.
- `LOAD 500` combina las 15 FORCE y 3 MOMENT en el único subcaso.

`PARAM,LGDISP,1` activa grandes desplazamientos y rotaciones;
`PARAM,FOLLOWK,NO` mantiene todas las cargas conservativas. `NLSTEP 500`
solicita 50 incrementos fijos y Newton completo `PFNT`, lo que facilita
comparar exactamente los mismos estados con MSC.

Quedan fuera del modelo la plasticidad, pandeo, contacto, pretensado,
gravedad automática, carga distribuida, cables, cargas seguidoras y arco de
longitud. Por ello debe presentarse como validación del solver mixto
CBAR/CROD, no como certificación estructural de un puente.

## Generación y ejecución local

Desde la raíz del repositorio:

```powershell
.\.venv\Scripts\python.exe tests\solver_tst\validation_cases\case_nonlinear_test_final\bridge_model.py
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_nonlinear_test_final\tfm_final_nonlinear_bridge.dat
```

El generador es la fuente editable del modelo. La solución local alcanza
`LOAD STEP=1.0` en 50 incrementos, sin recortes, con desplazamiento máximo
absoluto `21.512689 mm` y energía elástica final
`1.167100324e7 Nmm`.

### Correlación con MSC Nastran 2026.1

| Magnitud | MEF-solver | MSC Nastran |
|---|---:|---:|
| GRID 105, `T1/T2/T3` (mm) | `1.698907 / 0.598820 / -20.030602` | `1.699041 / 0.598616 / -20.031761` |
| GRID 205, `T1/T2/T3` (mm) | `1.840572 / 0.580943 / -15.778262` | `1.840642 / 0.580739 / -15.779271` |
| GRID 304, `T1/T2/T3` (mm) | `2.354343 / 2.784631 / -17.267393` | `2.354375 / 2.782789 / -17.268139` |
| Reacción vertical GRID 101 (N) | `442379.190` | `442384.531` |
| CBAR 1004, tensión axial (MPa) | `59.947538` | `59.952720` |
| CBAR 2004, tensión mínima de fibra (MPa) | `-198.581021` | `-198.588500` |
| CROD 3007, tensión axial (MPa) | `224.598734` | `224.632200` |
| CROD 3029, tensión axial (MPa) | `-152.441700` | `-152.414100` |
| Energía final (Nmm) | `1.167100324e7` | `1.167160e7` |

La comparación del estado completo disponible en el OP2 y de las referencias
F06 da:

| Campo | Error observado |
|---|---:|
| Desplazamientos, error relativo L2 | `8.238e-5` |
| Reacciones SPC, error relativo L2 | `1.619e-4` |
| Energía final, error relativo | `5.11e-5` |
| Tensiones axiales CROD seleccionadas, máximo relativo | `2.39e-4` |
| Estados centroidales CBAR seleccionados, máximo relativo | `< 3e-4` |
| Fibras CBAR seleccionadas, máximo relativo | `1.36e-2` |

El E2E usa una envolvente relativa del `2 %` para conservar explícitamente la
fibra C de CBAR 1004, cuya diferencia es `1.36 %`; no se ha eliminado ni se ha
ocultado con una tolerancia global amplia. Desplazamientos, reacciones,
energía y resultados centroidales quedan varios órdenes de magnitud por
debajo de ese límite.

Como control, al resolver la misma carga con `LGDISP=-1`, `T1` del GRID 105
es `1.733100 mm`; la solución con grandes desplazamientos es un `1.97 %`
menor. Esta diferencia comprueba que el caso no se limita a repetir la ruta
lineal, aunque permanece en una rama estable y apropiada para Newton por
control de carga.

## Auditoría MSC y apertura en Patran

La ejecución auditada corresponde exactamente a
`tfm_final_nonlinear_bridge.dat`, cuyo SHA-256 es
`221618ef61e366f86c0f01465bf968bdbccec288f8f09c6c6ba4aae89ac71acd`.
El SHA-256 del F06 es
`f7fb9b9c48762ae9b8bcd430fb3168e4e222487ef79aebd5e776e9d6863044fc`
y el del OP2 consolidado es
`8e17c666382208a627c2e60d9482333e0e7064800550949824e5271575be764a`.

La auditoría realizada fue:

1. Confirmar en LOG/F06 que no existe `FATAL`, que aparece
   `JOB CONVERGES FOR THE CURRENT STEP` y que el último estado es
   `LOAD STEP = 1.00000E+00`.
2. Conservar F06, OP2 y PCH. El F06/PCH son la fuente auditable y los valores
   seleccionados ya están transcritos en
   `test_nonlinear_mixed_CROD_CBAR_validation_E2E.py`.
3. Comparar los historiales de los GRID `105, 205, 109, 209, 304, 404`, las
   reacciones de `101, 201, 109, 209`, los CBAR `1004, 1104, 1205, 2004,
   2104, 2204` y los CROD `3007, 3008, 3015, 3029`.
4. Comparar también el número de incrementos, la energía final y los arrays
   completos disponibles del OP2 antes de elegir máximos para la presentación.

Para visualizarlo en Patran: crear una base nueva con preferencia MSC
Nastran, mantener unidades `mm-N-MPa` y usar `Analysis > Read Input File`.
Tras revisar `Reject Cards`, ejecutar el DAT original y acceder al OP2 mediante
`Analysis > Access Results > Read Output2`. Los rangos de EID y PID permiten
crear grupos separados para tablero, arcos, vigas transversales, péndolas y
arriostramientos.
