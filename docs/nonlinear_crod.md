# Análisis no lineal geométrico de CROD

## Estado y alcance

Esta fase incorpora el contrato de entrada, el elemento corrotacional y un
solver incremental para estática no lineal geométrica con CROD.
`LinearStaticSolver` permanece sin cambios físicos; `main.py` dirige SOL 400 a
`src/solver/nonlinear_analysis.py` y SOL 101 a su ruta lineal.

El alcance implementado es:

- `MAT1` elástico isótropo, `PROD` y `CROD`;
- grandes desplazamientos y rotaciones con deformación axial finita;
- cargas nodales conservativas `FORCE`/`MOMENT` y combinaciones `LOAD`;
- restricciones `SPC`, `SPC1`, `SPCADD` y `GRID.PS` ya normalizadas,
  incluidos valores no nulos de `SPC` aplicados de forma incremental;
- incremento `FIXED`/`ADAPT`, Newton-Raphson completo, residuo, recortes y
  tangente consistente;
- sin plasticidad, contacto, arc-length, cargas seguidoras ni CBAR.

La ampliación a CBAR será posterior y no debe reutilizar el CROD como si ya
incluyera flexión o una cinemática de sección.

La guía de SOL 400 indica que sus elementos 1D admiten grandes desplazamientos
y rotaciones. La gran deformación material es una opción distinta y no se ha
activado en esta primera fase. Los resultados 94--96 observados con `LGDISP=1`
dan `epsilon=(L-L0)/L0`, `sigma=E*epsilon` y `N=A*sigma` usando el área de
referencia; no muestran una actualización isócora de `A`. Esa es la ley que
reproduce el solver, limitado a no linealidad geométrica.

## Tarjetas principales de MSC Nastran

### Selección de análisis

```text
SOL 400
CEND
SUBCASE 1
  ANALYSIS = NLSTAT
  SPC = 10
  LOAD = 20
  NLSTEP = 30
```

- `SOL 400` selecciona el solver implícito no lineal de MSC Nastran.
- `ANALYSIS=NLSTAT` identifica una etapa estática no lineal. El lector
  normaliza también los alias `NLSTATIC` y `NLSTATICS`.
- `SPC` y `LOAD` seleccionan las condiciones de contorno y el valor **total**
  de carga de la etapa. `NLSTEP` define cómo se alcanza ese valor.
- El entero de Case Control `NLSTEP=30` debe coincidir con el `SID` de una
  entrada `NLSTEP` de Bulk Data.

### No linealidad geométrica

```text
PARAM,LGDISP,1
PARAM,FOLLOWK,NO
```

`PARAM,LGDISP,1` activa grandes desplazamientos. Esto obliga a actualizar la
geometría, las fuerzas internas y la rigidez tangente; no equivale a resolver
repetidamente la matriz elástica lineal sin término geométrico.

El preproceso conserva el valor completo: `-1` desactiva no linealidad
geométrica, `0` conserva sólo rigidez diferencial, `1` activa grandes
desplazamientos y cargas seguidoras, y `2` activa grandes desplazamientos sin
cargas seguidoras. Los casos usan `-1` en los controles 92/121 y `1` en los
restantes salvo el control 118, que usa `2`.

MSC activa también tratamiento de cargas seguidoras al usar grandes
desplazamientos. La campaña inicial fija `PARAM,FOLLOWK,NO` para excluir la
rigidez de carga seguidora y aislar la tangente estructural. Las `FORCE` se
definen en el sistema básico o en un `CID` fijo. La primera implementación del
proyecto se limitará a cargas nodales conservativas; cualquier ampliación a
cargas seguidoras necesitará referencias independientes.

### Desplazamientos impuestos

El preprocesador y `NonlinearStaticSolver` admiten desplazamientos impuestos
mediante un valor no nulo en una entrada `SPC`. Por ejemplo:

```text
SPC1,124,23456,5
SPC,124,5,1,-4.0
```

El GRID 5 queda restringido en los seis GDL y alcanza `T1=-4.0` cuando el
factor de carga llega a uno. Durante la solución se prescribe
`T1(lambda)=-4.0*lambda`, de manera coherente con el resto de acciones de la
etapa. Los casos 124 y 126 ejercitan conjuntamente fuerzas nodales y este tipo
de movimiento. La tarjeta de carga `SPCD` no forma parte todavía del
subconjunto leído; usarla no debe confundirse con el valor no nulo de `SPC`.

### Por qué `K6ROT` no aparece en los casos CROD

`PARAM,K6ROT` es una rigidez numérica para el giro de perforación, es decir,
la rotación alrededor de la normal de elementos de lámina como CTRIA3 y
CQUAD4. No añade rigidez física ni estabilización a un CROD. En una cercha de
CROD con `J=0`, los GDL rotacionales nodales no intervienen en la energía y se
restringen explícitamente con `SPC`; dejarlos libres sería un mecanismo real
del modelo discreto. Si `PROD.J>0`, sólo la rotación relativa proyectada sobre
el eje de la barra produce torsión, como comprueban los casos 116/117.

Por esa razón añadir `K6ROT` a un deck compuesto únicamente por CROD sería
inerte y podría ocultar la razón de las restricciones `456`. Los cinco decks
123--127 incluyen además un comentario explícito recordándolo.

### Control incremental e iterativo

Ejemplo con incrementos fijos:

```text
NLSTEP,30,1.0,MILDLY
,GENERAL,30,1,10
,FIXED,20
,MECH,PV,1.E-8,1.E-8,1.E-8,PFNT
```

Ejemplo adaptativo:

```text
NLSTEP,30,1.0,SEVERELY
,GENERAL,50,1,20
,ADAPT,.025,1.E-6,.10,4,1.20,0,500
,MECH,PV,1.E-9,1.E-9,1.E-9,PFNT
```

El subconjunto leído por el preproceso es:

- cabecera: `ID`, tiempo total y `CTRLDEF` (`QLINEAR`, `MILDLY` o
  `SEVERELY`);
- `GENERAL`: `MAXITER`, `MINITER` y `MAXBIS`;
- `FIXED`: número `NINC` de incrementos iguales;
- `ADAPT`: fracciones inicial/mínima/máxima, iteraciones deseadas, factor de
  crecimiento, intervalos de salida y máximo de incrementos;
- `MECH`: criterios `U/P/W/V/N/A`, tolerancias y política de actualización de
  rigidez.

`PFNT` solicita Newton-Raphson completo y es la opción usada en la campaña.
En `ADAPT`, `INTOUT=0` solicita salida en cada incremento convergido y `-1`
sólo al final; los casos de diagnóstico usan `0` para alimentar los OP2
intermedios de `NLOPRM`.
`ITER,KSTEP=1` se normaliza también, pero no se usa como referencia principal.
Los métodos de arco `ARCLN`, control térmico y acoplamiento se rechazan porque
están fuera del alcance acordado.

El solver implementado exige `PFNT` y que `MECH` incluya el criterio de residuo
de fuerzas `P`. Las banderas adicionales se conservan en el modelo, pero la
aceptación del incremento se decide por ese residuo; una entrada sin `P` se
rechaza en lugar de ignorar silenciosamente su criterio.

El lector admite además `NLPARM` como alternativa heredada y lo convierte al
mismo modelo de datos. Los nuevos decks utilizan `NLSTEP`, que es el control
preferido para SOL 400.

### Salida y diagnóstico

```text
NLOPRM OUTCTRL = (SOLUTION,INTERM) NLDBG = ADVDBG
DISPLACEMENT(PRINT,PUNCH) = ALL
SPCFORCES(PRINT,PUNCH) = ALL
GPFORCE(PRINT,PUNCH) = ALL
STRESS(PRINT,PUNCH) = ALL
FORCES(PRINT,PUNCH) = ALL
```

`NLOPRM` debe aparecer una sola vez por caso de control; sus opciones se
yuxtaponen en la misma sentencia. `NLOPRM` y las solicitudes anteriores no
cambian la física. Permiten conservar
los estados convergidos, la tabla de iteraciones/residuos y los balances en
puntos de malla. `NLOUT` se reserva normalmente para seleccionar resultados no
lineales avanzados; no es necesario para las tablas convencionales de CROD de
esta campaña.

#### Cómo localizar los resultados finales en el F06

Un F06 de SOL 400 repite las tablas para cada estado convergido. Antes de
transcribir una referencia E2E hay que localizar el último bloque
`LOAD STEP = 1.00000E+00` y confirmar el mensaje
`JOB CONVERGES FOR THE CURRENT STEP`. Dentro de ese bloque:

- `DISPLACEMENT VECTOR` contiene `T1 T2 T3 R1 R2 R3` por `POINT ID`;
- `FORCES OF SINGLE-POINT CONSTRAINT` contiene las reacciones con las mismas
  seis columnas;
- `NONLINEAR STRESSES IN ROD ELEMENTS (CROD)` contiene por `ELEMENT ID`
  `AXIAL STRESS`, `TOTAL STRAIN` y `LIN. TORSIONAL STRESS`;
- `GLOBAL WORKS & ENERGIES` contiene `TOT STRAIN ENERGY`;
- las líneas `FORCE RESIDUAL CONVERGENCE RATIO` y los mensajes de convergencia
  permiten auditar iteraciones, residuos y número de incrementos.

Por ejemplo, en el bloque final del caso 112, GRID 3 imprime
`T2=-1.538383E+01`; GRID 1 imprime reacciones
`T1=1.181807E+06`, `T2=1.000000E+06`; y CROD 1 imprime
`AXIAL STRESS=-1.548118E+04` y `TOTAL STRAIN=-7.371989E-02`.
Esos valores se convierten respectivamente en las claves `u`, `r` y `s` de
`test_nonlinear_validation_end_to_end.py`. La fuerza axial se obtiene como
`N=A*sigma` si no aparece directamente, y el par como
`T=(J/C)*tau_torsion`. El propio módulo de test contiene el ejemplo completo
y la correspondencia de todas sus claves.

No existe una salida estándar equivalente a una DMIG de **cada** tangente de
Newton. Por ello la primera auditoría de rigidez se apoya en:

1. el par 92/93 para el límite de carga pequeña;
2. equilibrio `F_ext - F_int` mediante `GPFORCE`, reacciones y fuerzas CROD;
3. factores de convergencia y peor residuo de `.f06`/`.sts`;
4. sensibilidad a discretización incremental en 94/95 y 104/109/110.

Si estas capas no bastan, se preparará una extracción DMAP específica para la
versión de MSC Nastran usada; no debe asumirse que una DMIG lineal representa
la tangente convergida.

La auditoría de MSC Nastran 2026.1 confirmó además que 114/115 pueden usarse
como sonda numérica: la diferencia central de desplazamientos ante cargas
horizontales `+/-1000` aproxima la flexibilidad tangente alrededor del estado
deformado de 112. Los casos 116/117 separan torsión pura y tracción+torsión;
119 fija una referencia estable anterior al salto observado en 102, 120
ejercita de forma deliberada el control adaptativo a través de ese salto y el
par 121/122 comprueba el límite lineal en una estructura que cambia de
orientación.

## Referencias primarias

- [MSC Nastran 2021.1 Nonlinear User's Guide SOL 400](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2021.1_Nonlinear_SOL_400_User_Guide/raw/resource/enus/MSC_Nastran_2021.1_Nonlinear_SOL_400_User_Guide.pdf): tipos de análisis, `NLSTEP`, convergencia y `LGDISP`.
- [MSC Nastran 2021.3 SOL 400 Getting Started Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2021.3_SOL_400_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2021.3_SOL_400_Getting_Started_Guide.pdf): elementos 1D, diagnóstico e intermedios `NLOPRM`.
- [MSC Nastran 2020 Reference Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2020_Reference_Guide/raw/resource/enus/MSC_Nastran_2020_Reference_Guide.pdf): ejemplos de análisis geométrico `SOL 400` con `PARAM,LGDISP`.
- [MSC Nastran 2023.4 Reference Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2023.4_Reference_Guide/raw/resource/enus/MSC_Nastran_2023.4_Reference_Guide.pdf): alcance de `K6ROT` como rigidez numérica del giro de perforación de láminas.

## Contrato del preproceso

`FEMModel.analysis` contiene un `AnalysisConfiguration` con:

- secuencia de solución y tipo de análisis;
- palabra `ANALYSIS` activa;
- valor `LGDISP`, activación derivada de grandes desplazamientos y rigidez de
  carga seguidora;
- `NonlinearStepControl` normalizado desde `NLSTEP` o `NLPARM`.

pyNastran 1.4.1 no implementa la entrada Bulk Data `NLSTEP`. La conserva en
`reject_lines`; `nastran_reader_v2.py` interpreta allí sus filas libres o de
campo pequeño y aplica validación estricta. Una continuación desconocida no se
ignora: detiene el preproceso.

## Formulación implementada

En cada incremento de factor de carga `lambda`, el solver busca:

```text
r(u, lambda) = lambda F_ext - F_int(u) = 0
K_T(u) Delta_u = r(u, lambda)
u_(i+1) = u_i + Delta_u
```

Para una cuerda actual `d`, longitud `l` y dirección `n`, usa:

```text
epsilon = (l-L0)/L0
N = E*A*epsilon
f_int = [-N*n, +N*n]
k = E*A/L0*(n outer n) + N/l*(I - n outer n)
```

El primer término es material y el segundo geométrico. El test elemental
compara la matriz analítica con la derivada central de `f_int`; el test global
compara la flexibilidad deformada `K_T^-1[ux,ux]` con la sonda Nastran 114/115.
La torsión usa `GJ/L0` sobre el eje de referencia, reproduce `R1`, par y
`tau=T*C/J` de 116/117 y permanece desacoplada de la extensión axial.

## Estado de validación

Los 35 casos que alcanzan carga completa (92--119 y 121--127) se resuelven de
extremo a extremo y se comparan contra valores Nastran escritos explícitamente
en `test_nonlinear_validation_end_to_end.py`. Se cubren desplazamientos,
reacciones, tensión/deformación axial, torsión, sistemas CD/CID, LOAD,
ensamblajes 2D/3D, desplazamientos impuestos y controles fijos/adaptativos.
Las entradas nuevas usan claves descriptivas (`case_name`, `folder`,
`filename`, `description`, tolerancias y resultados esperados) documentadas en
el propio módulo de test.

El historial de 112 coincide en 25, 50, 75 y 100 % de carga. La flexibilidad
tangente horizontal obtenida en el estado final es `6.127186e-6`, igual a la
diferencia central de 114/115. Los controles adaptativos 95, 104 y 111
reproducen 9, 14 y 15 incrementos, respectivamente.

El informe `.rpt` del proyecto conserva ahora ese tipo de trazabilidad durante
una ejecución ordinaria. Distingue SOL 400 en el cabecero, imprime la
configuración normalizada y separa el historial de incrementos convergidos del
historial de correcciones de Newton. Los desplazamientos, reacciones y estados
CROD bajo `FINAL CONVERGED STATE RESULTS` pertenecen exclusivamente al último
equilibrio convergido; los dos historiales se imprimen después. La guía
[`result_report.md`](result_report.md) explica cada columna y su localización
en F06.

Los casos estructurales 123--127 reproducen `25/30/30/40/14` incrementos; la
secuencia completa de factores adaptativos de 127 también coincide. La
auditoría completa de sus F06 obtuvo errores relativos máximos de `3.77e-7` en
desplazamientos, `3.76e-7` en reacciones y `3.72e-7` en resultados CROD. En
125, MSC salta de `T3=-13.5215` a `-69.40416` entre factores `0.533333` y
`0.566667`; por ello su `TotStrEn` acumulado no se usa como energía elástica
final, aunque todo el estado mecánico final coincide.

El caso 120 no alcanza carga completa bajo control de fuerza: tanto MSC como
el solver se detienen alrededor de `lambda=0.98387`, al aproximarse al punto
límite. El número exacto de recortes no se iguala porque la heurística interna
de adaptación de MSC no está publicada; el fenómeno y el último estado estable
sí coinciden. No se implementa arco de longitud.

La energía devuelta por el proyecto es la energía elástica del estado actual.
En ramas continuas coincide con Nastran salvo el error de integración
incremental. En 102/105, el campo global `TotStrEn` de Nastran acumula trabajo
a través del salto y no representa esa misma magnitud; por ello no se usa como
referencia final en esos dos casos.
