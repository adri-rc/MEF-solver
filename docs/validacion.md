# Validación

## Estado observado

Comando:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Resultado verificado el 25 de agosto de 2026 después de incorporar el puente
no lineal actualizado de tres vanos:

```text
688 passed, 6 xfailed, 2 warnings
```

Las dos advertencias proceden de que pytest no puede escribir su caché.

## Capas de comprobación

1. Modelo de datos, lector, referencias y sistemas de coordenadas.
2. Funciones de forma, jacobianos, transformaciones y cuadratura.
3. Propiedades matemáticas de cada elemento: simetría, modos rígidos,
   orientación, desacoplamientos y geometrías inválidas.
4. Matrices DMIG frente a Nastran.
5. Validación E2E de desplazamientos, cargas, reacciones y resultados de
   elemento con referencias transcritas en los tests.

El contrato del informe se comprueba aparte en
`tests/postprocessing_tst/test_report_writer.py`. Para CBAR existe un estado
SOL 400 explícito que exige las tablas de respuesta axial/torsional, longitud y
energía no lineales, acciones A/B y tensiones en puntos de recuperación. La
misma prueba comprueba que `P0` se presenta como C y que un CBAR lineal no
recibe una tabla `NONLINEAR STATE` ficticia. Como comprobación funcional, el
caso 128 genera las cuatro tablas a partir del resultado real del solver.

El caso 80 añade una comprobación integrada con CBAR, CTRIA3 y CQUAD4 en el
mismo modelo. Contrasta desplazamientos, cargas, reacciones, resultados de las
dos familias de lámina y esfuerzos de CBAR. Es la referencia más directa para
detectar regresiones de ensamblaje o recuperación entre tipos de elemento.

## Modelo lineal final para la defensa

`case_linear_test_final` es un cupón plano de frame/web de `1200 x 600 mm`
con dos huecos circulares. Su malla reproducible contiene 256 CQUAD4, 128
CTRIA3 alrededor de los huecos, 76 CBAR de cap/nervio y 4 CROD de introducción
de carga. En un único subcaso combina dos FORCE, un MOMENT, seis estaciones
empotradas y un desplazamiento impuesto de `0.8 mm`.

`test_linear_final_validation_case.py` comprueba el recuento y subconjunto de
tarjetas, la solución finita, el valor impuesto, el equilibrio global de
fuerzas/momentos y referencias MSC explícitas de desplazamientos, reacciones,
CROD, CBAR, CTRIA3 y CQUAD4. El CTRIA3 se usa en membrana: el útil fija
`T3/R1/R2` y su MID3 casi nulo declara la limitación de cortante transversal
sin estabilización silenciosa.

El DAT sin modificar terminó normalmente con MSC Nastran Student Edition
2026.1 el 21 de agosto de 2026. El F06 se conserva con el caso y sus valores
seleccionados están transcritos en el test, que no lee F06/OP2 en tiempo de
ejecución. La comparación de los arrays completos del OP2 dio errores relativos
L2 de `1.052e-6` en desplazamientos, `1.450e-5` en reacciones, inferiores a
`5.4e-6` en CROD/CBAR y a `1.1e-5` en resultados CTRIA3/CQUAD4. Es una
validación cruzada numérica del deck, no una correlación experimental ni una
validación general de Nastran. El dimensionado, tablas y flujo Patran están en el
[`README` del caso](../tests/solver_tst/validation_cases/case_linear_test_final/README.md).

## Modelo lineal final tridimensional en T

`case_linear_test_final_T` conserva el web perforado de `1200 x 600 mm` y
añade un ala CQUAD4 de `1200 x 300 mm` en el plano perpendicular. Seis
diagonales, partidas en doce CBAR, forman tres nervios simétricos; la mitad
próxima al ala duplica `A`, `I1`, `I2` y `J`. El inventario completo es 417
GRID, 288 CQUAD4, 128 CTRIA3, 88 CBAR y 2 CROD. Cuatro nodos del ala idealizan
uniones remachadas discretas: dos están empotrados y dos reciben `T2=0.2 mm`.
Dos fuerzas axiales de `6000 N` actúan en los CROD y un momento de
`250000 Nmm` en el punto medio superior. El deck declara `AUTOSPC=NO`.

`test_linear_final_t_validation_case.py` incluye controles de inventario,
geometría, acciones/contorno, equilibrio y referencia MSC. Los literales del
F06 cubren 15 GRID, seis reacciones, ambos CROD y 15 CBAR, 15 CTRIA3 y 15
CQUAD4; el test no lee resultados externos durante pytest. Frente al OP2
completo, los errores L2 son `8.206e-3` en traslaciones, `2.547e-2` en
rotaciones, `1.705e-4` en reacciones de fuerza, `3.553e-4` en reacciones de
momento y `3.930e-3` en tensiones normales CBAR. La membrana CQUAD4/CTRIA3
queda en `3.268e-5/7.218e-6` y las tensiones de fibra en
`6.582e-3/1.020e-2`.

La recuperación CTRIA3 centroidal de momentos y cortantes no se acepta como
correlacionada: sus errores L2 son `0.682` y `17.30`. Se conserva como límite
explícito y no se cubre con una tolerancia amplia ni motiva una recalibración
de rigidez. El F06 de MSC Nastran Student Edition 2026.1 acaba en `END OF JOB`
con `AUTOSPC=NO`. Diseño, propiedades, selección de EID, hashes, discrepancias
y capturas Patran están en el
[`README` del caso](../tests/solver_tst/validation_cases/case_linear_test_final_T/README.md).

### Diagnóstico complementario sin MID3 efectivo

`case_linear_test_final_T_no_CTRIA_MID3` repite exactamente la geometría,
cargas y restricciones del frame T. La única variación constitutiva es que
`PSHELL 11.MID3` referencia el nuevo `MAT1 99`, cuyo módulo es `E=1e-20`. MSC Nastran
2026.1 alcanza `END OF JOB` con `AUTOSPC=NO`, pero la matriz queda
extremadamente mal condicionada (`cond(K)=5.660208e28`). No es un modelo
físico alternativo, sino una sonda para aislar la contribución de cortante
transversal del CTRIA3.

En esta sonda, el error L2 de los momentos CTRIA3 del campo completo baja de
`0.682` a `7.106e-5`. Los cortantes son cero ingenieril en ambos resultados:
la norma local es menor que `4.1e-11` y la de MSC menor que `8.7e-24`. Las
demás familias de elemento, reacciones y grados de libertad no asociados a
`T3` conservan la correlación. En cambio, la traslación transversal del
CTRIA3 se convierte en un casi-mecanismo: en GRID 203 el solver obtiene
`7.3586e10 mm`, mientras MSC imprime `-0.3421 mm`. Por ello el deck degenerado
no valida desplazamientos globales ni una respuesta estructural físicamente
utilizable.

La pareja de casos cierra la investigación causal dentro de la limitación
formulativa actual: el frame original valida respuesta global, membrana y
fibras, y la sonda demuestra que la discrepancia `M/Q` depende del tratamiento
de cortante. No demuestra equivalencia física de la formulación ni incorpora
la recuperación física de momentos/cortantes CTRIA3 de Nastran; esa capacidad
permanece explícitamente fuera de alcance. El propósito, referencias y
advertencias completas están en el
[`README` de la sonda](../tests/solver_tst/validation_cases/case_linear_test_final_T_no_CTRIA_MID3/README.md).

## Modelo no lineal final para la defensa

`case_nonlinear_test_final` es un puente arco atirantado espacial de
`12 000 x 2 400 mm`, con 32 GRID, 48 CBAR y 30 CROD. Los CBAR representan dos
vigas de tablero, vigas transversales, dos arcos parabólicos y sus riostras;
los CROD representan 14 péndolas y 16 diagonales espaciales. El único subcaso
`SOL 400` combina 15 FORCE, 3 MOMENT, dos empotramientos nodales y asientos
diferenciales impuestos de `-12/-6 mm` en el estribo derecho.

`test_nonlinear_final_validation_case.py` comprueba el inventario, el alcance
`MAT1/PBAR/PBARN1/PROD`, el control `NLSTEP,FIXED,50`, la convergencia a factor
uno sin recortes, los valores impuestos y el residuo de equilibrio en todos
los GDL libres. Un control con `LGDISP=-1` da en GRID 105 una respuesta
longitudinal un `1.97 %` distinta, por lo que el modelo ejercita de forma
medible la cinemática no lineal sin aproximarse a un punto límite. Las
referencias MSC pertenecen al módulo E2E mixto común.

El DAT sin modificar terminó sin `FATAL` con MSC Nastran Student Edition
2026.1/NLPERF y alcanzó los mismos 50 incrementos. Los errores relativos L2
son `8.238e-5` en desplazamientos y `1.619e-4` en reacciones; la energía
difiere `5.11e-5`, los estados axiales seleccionados menos de `3e-4` y la
mayor diferencia de fibra CBAR seleccionada es `1.36e-2`. El E2E conserva una
envolvente del `2 %` debido a esa fibra, sin ampliar la tolerancia de las
campañas puras. El diseño, hashes, puntos GRID/EID y flujo Patran están en el
[`README` del caso](../tests/solver_tst/validation_cases/case_nonlinear_test_final/README.md).

### Puente actualizado de tres vanos

`case_nonlinear_test_final_updated` eleva la integración a 202 GRID y 572
elementos: 144 CBAR para largueros exteriores, tres pares de arcos y
cimentaciones, y 428 CROD para el tablero articulado, péndolas y
arriostramientos. La geometría tiene luces `8/12/10 m`, flechas distintas y
cuatro bases de terreno escalonadas. El único subcaso combina 44 FORCE, 3
MOMENT y asientos `-2/-6/-10 mm`, con `AUTOSPC=NO` y 20 incrementos fijos.

MSC Nastran 2026.1 y MEF-solver alcanzan factor uno sin recortes. Sobre las
familias completas, los errores L2 son `3.280e-5` en traslaciones,
`2.416e-5` en rotaciones, `2.898e-5` en reacciones de fuerza, `2.125e-5` en
tensión CROD y `1.379e-5` en tensión axial CBAR. La energía difiere
`1.392e-5`. Las fibras CBAR tienen `1.673e-2` L2 y concentran su diferencia en
los largueros de tablero: el peor punto con `|sigma_MSC|>10 MPa` alcanza
`11.57 %`. El E2E mantiene envolventes separadas para que esa recuperación no
relaje la respuesta global o axial.

El inventario, resultados completos, selección literal F06/PCH y alcance de
las vistas Patran se conservan en el
[`README` actualizado](../tests/solver_tst/validation_cases/case_nonlinear_test_final_updated/README.md)
y en su
[`análisis de resultados`](../tests/solver_tst/validation_cases/case_nonlinear_test_final_updated/RESULTS_ANALYSIS.md).

Las validaciones E2E de CBAR transcriben las tablas F06, sin leerlas durante la
prueba. Los casos 13--17 cubren voladizo, orientación 3D, vigas en serie,
marco espacial y trípode híbrido; el caso 15 comprueba las tensiones PBAR C/D
en ambos extremos A/B. El caso 18 amplía la recuperación a C--F y el 19 a una
geometría espacial con varias propiedades. Los signos de `forces_n1` y
`forces_n2` son acciones nodales del elemento, mientras que las tablas
`FORCES IN BAR ELEMENTS` del F06 son acciones internas de sección; cada
assertion conserva explícitamente esa conversión de convención.

## Matrices CQUAD4

El test `test_E2E_CQUAD_Nastran_K_solution.py` aplica dos criterios
complementarios:

- error relativo de Frobenius para medir la discrepancia global;
- comparación término a término mediante `allclose`.

Además repite la comparación en ejes del elemento para membrana,
placa/cortante, drilling y sus acoplamientos. Esta descomposición evita que la
mayor magnitud de la membrana oculte un término incorrecto de cortante.

Los casos activos cubren:

- membrana, flexión y combinación de ambas;
- `K6ROT=0/100`;
- `MID3` vacío y activo;
- cambios de `G3` y `TS/T`;
- rectángulos, trapecios y cuadriláteros afines sesgados;
- orientación arbitraria de un plano en 3D;
- DMIG completas de oblicuidad moderada y alta.

El caso 90 reveló que la ley oblicua anterior era demasiado blanda en alta
oblicuidad. Se mantuvieron el operador por aristas, la flexibilidad residual y
la contribución MID3 y se ajustó únicamente la métrica torsional racional. El
resultado cumple:

- matriz completa: error de Frobenius menor que `1e-6`;
- placa/cortante: menor que `4e-6`;
- términos del bloque: `rtol=1e-5`, `atol=2e-3`.

## E2E CQUAD4

`test_validation_end_to_end_cquad4.py` no lee F06 durante su ejecución. Los
valores de Nastran están escritos directamente en los casos, lo que permite
auditar cada referencia.

El caso 91, pareja cargada del DMIG 90, confirma el ajuste con dos paneles y
propiedades distintas:

- error vectorial de desplazamientos: `2.07e-5`;
- error vectorial de reacciones: `6.30e-5`;
- `N`, `M` y tensiones de fibra: menos de 0.1 % por componente;
- mayor discrepancia `Qx/Qy`: 47.1 %, con 0.0216 absoluto.

La diferencia de `Q` pertenece a la recuperación centroidal, no a la matriz de
rigidez. No se observa una discrepancia asociada en desplazamientos, reacciones
o momentos.

### Casos 57–59

- **57:** `cond(K)=2.90e25`; errores vectoriales `U=1.72e-3`, `R=4.07e-3` y
  resultados `9.28e-4`. La mayor diferencia relativa de desplazamiento está
  en una rotación casi nula.
- **58:** `cond(K)=8.52e26`; traslaciones próximas al 0.3 %, pero rotaciones,
  momentos secundarios y `Qy` sensibles. Error vectorial de resultados
  `2.94e-2`.
- **59:** geometría alabeada, fuera de alcance. El test es `xfail` estricto y
  espera que el solver la rechace.

## Significado de los seis xfail

| Casos | Motivo | Interpretación |
|---|---|---|
| 54, 59, 62, 63, 64 | CQUAD4 alabeados | Fuera del alcance plano; rechazo esperado |
| 78 | suavizado `SNORM` por defecto | No implementado; usar la variante `SNORM,0` |

Son límites declarados, no errores pendientes del CQUAD4 plano. Los casos 82 y
88 confirman ensamblajes facetados con normales independientes. El antiguo
`xfail` 142 se cerró con la interpolación lineal-cortante validada por
154--160.

## Criterio de cierre

La implementación CQUAD4 se considera cerrada para elementos coplanares porque:

- sus invariantes y modos rígidos están cubiertos;
- las DMIG activas coinciden globalmente, por bloques y término a término;
- los modelos con varios elementos, materiales y orientaciones reproducen la
  respuesta global de Nastran;
- el dominio no soportado se rechaza en código y se documenta en tests.

La única mejora específica no bloqueante es reproducir con mayor precisión la
convención de impresión centroidal de `Qx/Qy`.

## Campaña no lineal geométrica de CROD

Los casos 92--111 forman la campaña inicial de 20 decks `SOL 400`. Cubren el
límite lineal, tracción y compresión finitas, control fijo/adaptativo,
`CD/CID`, series, paralelo, `LOAD`, cerchas 2D, estructuras 3D y sensibilidad
al tamaño de incremento. La auditoría añadió los casos 112--122 para
objetividad 3D, tangente por perturbación central, torsión, `LGDISP=2` y
respuesta próxima al punto límite, además de un par lineal/no lineal sensible
a la rotación de la estructura.

Los casos 123--127 amplían la escala a estructuras de 10--30 CROD: una
cubierta arqueada, un puente de celosía, una cubierta espacial, una marquesina
y una torre 3D. Mezclan PROD/MAT1, orientaciones, cargas y apoyos; 124/126
incluyen valores no nulos de SPC. Los cinco convergen con MSC y con el solver,
reproducen `25/30/30/40/14` incrementos y forman parte del diccionario E2E.

El test `test_nonlinear_analysis_control.py` verifica que todos los decks:

- seleccionan `SOL 400`, `ANALYSIS=NLSTAT` y el `NLSTEP` correcto;
- sólo contienen CROD dentro del alcance de esta campaña;
- mapean `LGDISP` y `FOLLOWK` sin contaminar la ruta lineal;
- se enrutan desde el CLI a `NonlinearStaticSolver` sin pasar por
  `LinearStaticSolver`.

Como control adicional, los modelos se ensamblan y resuelven con su rigidez
lineal inicial cuando ese control es físicamente aplicable; los decks de
torsión conservan expresamente el GDL `R1` al que aporta rigidez el `J` de
`PROD`.

La primera ejecución 92--111 con MSC Nastran 2026.1 (`NLPERF`) terminó
normalmente y alcanzó el 100 % de carga en todos los casos. No hubo cortes ni
separaciones. Los pares 92/93 y 94/95 coinciden en sus estados finales; los
casos 104/109/110 coinciden en desplazamiento, tensión y reacción, mientras el
trabajo/energía del mallado incremental grueso 109 difiere un `0.65 %`. El
caso 102 converge a la rama posterior al *snap-through* en su último
incremento (`uy: -42.830 -> -223.31`); 105 hace lo mismo en el incremento 17
(`uy: -33.975 -> -218.50`). Ambos se conservan como controles avanzados.

La campaña completa se relanzó después de corregir `NLOPRM`. Los 35 casos
convergidos se comparan con diccionarios explícitos en
`test_nonlinear_validation_end_to_end.py`; ningún test lee F06/PCH durante la
ejecución. Los casos 95/104/111 reproducen 9/14/15 incrementos adaptativos.

`test_nonlinear_intermediate_validation.py` transcribe cuatro estados de 112 y
compara además la flexibilidad tangente horizontal con la diferencia central
de 114/115 (`6.127186e-6`). `test_crod_non_linear.py` verifica por diferencias
finitas todos los términos de la tangente elemental.

Para 123--127 se contrastaron todas las filas finales impresas: los máximos
errores relativos fueron `3.77e-7` en desplazamientos, `3.76e-7` en reacciones
y `3.72e-7` en tensión/deformación. El caso 125 atraviesa un *snap-through*
entre `lambda=0.533333/0.566667`; su estado final coincide, pero el acumulado
`TotStrEn` de MSC no se usa como energía elástica de estado.

El caso 120 conserva su terminación no convergida en `lambda≈0.98387`: es una
validación del punto límite bajo control de fuerza, no un fallo pendiente. Los
casos 102/105 convergen mediante salto a una rama posterior; sus resultados
mecánicos coinciden, pero el acumulado `TotStrEn` de Nastran no es la energía
elástica del estado final y no se usa como aserción.
Véanse
[`nonlinear_crod.md`](nonlinear_crod.md) y la
[guía de campaña](../tests/solver_tst/validation_cases/NONLINEAR_CROD_CASES.md).

## Campaña CBAR no lineal

La campaña CBAR comprende 39 decks, casos 128--166. Cubre límite lineal, axial,
flexión biaxial, momentos, torsión, carga combinada, una sonda tangente
central, orientación `G0`, cargas en `CID`, `I12`, `K1/K2`, SPC no nulos, varios
materiales/elementos, refinamiento, release y offsets con `LROFF`. Los casos
149--153 añaden estructuras de 15 o 20 CBAR con geometrías planas, espaciales y
curvas, varias propiedades y materiales, fuerzas, momentos y desplazamientos o
giros impuestos. Los casos 154--160 aíslan `PBARN1,INT2=LC/LS`, los dos planos
de cortante y la semántica de `K1/K2`; 161--166 añaden recuperación C--F bajo
rotación 3D, `I12`, `G0`, `CID`, offsets y giro impuesto.

La ejecución original de 128--148 se descartó porque faltaba `PBARN1`. Tras
corregir los decks y relanzar 140, 148 y 150, los 26 F06/PCH alcanzan factor de
carga 1.0 y contienen resultados no lineales auditables.

Todos los decks solicitan resultados intermedios. El triplete 137--139 comparte
40 factores de carga y permite calcular por diferencia central la flexibilidad
tangente en cada estado convergido; los demás historiales aportan
desplazamientos, reacciones, energía y `NLSTRESS/BOUTPUT` por incremento.

`test_nonlinear_analysis_control.py` comprueba que cada carpeta contiene un
único DAT, que los 39 decks CBAR y los diez mixtos se leen como SOL 400, que
los campos avanzados
llegan al modelo interno y que toda tarjeta Bulk Data usada tiene un cabecero
de campos comentado. También comprueba que existe una extensión `PBARN1` por
cada `PBAR` y que las salidas específicas de SOL 400 están solicitadas. La
rigidez lineal inicial de los cinco modelos nuevos tiene rango completo.

El plan de implementación, el orden de cierre y la correspondencia entre
casos están en [`nonlinear_cbar.md`](nonlinear_cbar.md) y en la
[guía de campaña](../tests/solver_tst/validation_cases/NONLINEAR_CBAR_CASES.md).

La implementación se valida mediante:

- once pruebas unitarias de `CBARNonlinear`, incluidas las matrices `LC/LS`,
  objetividad y derivada central de la tangente;
- siete pruebas intermedias con las trayectorias 132/137, las sondas 138/139 y
  estados convergidos de 155--159;
- 39 entradas E2E CBAR con referencias explícitas de desplazamiento, reacción,
  energía y esfuerzos seleccionados. Todas pasan, incluido el caso 142;
- doce entradas E2E mixtas: los diez casos 167--176 y los dos puentes finales,
  con referencias explícitas CBAR/CROD y de respuesta global.

Los decks 154--160 terminan sin `FATAL`. Los máximos errores vectoriales de la
implementación frente a esos resultados son `8.04e-5` en desplazamientos,
`6.03e-5` en reacciones y `1.52e-4` en energía. `LC` coincide con
Euler--Bernoulli e ignora `K1/K2`; `LS` usa interpolación lineal con un punto
de cortante y toma `K=1` cuando el factor PBAR es cero.

Los casos axiales 130/131 identifican la deformación logarítmica de MSC. La
flexibilidad tangente de 137 difiere menos de 0.85 % a lo largo de cinco
estados. Los modelos 149--153 quedan dentro de 0.3--4 % en las magnitudes
seleccionadas. La tensión axial centroidal está validada; la convención de
tensión de fibra bajo rotación 3D combinada sigue siendo una mejora abierta.

Para estudiar esa mejora se ejecutaron los casos 161--166: comienzan por una
sección simétrica y añaden inclinación 3D, puntos C--F asimétricos con `I12`,
orientación `G0`, cargas en `CID`, varios elementos/materiales, offsets y giro
impuesto. Los casos 167--176 forman una segunda campaña de integración con 2 a
30 elementos y combinan CROD y CBAR, propiedades/materiales distintos,
`FORCE`, `MOMENT`, `LOAD`, empotramientos y SPC no nulos. Los 16 DAT pasan el
preprocesado, tienen rigidez inicial de rango completo y alcanzan carga completa
en MSC y con el solver local. Los F06 no contienen `FATAL`.

La solución mixta global queda cerrada por un módulo E2E independiente con
tolerancia relativa del 3 %. Los errores vectoriales observados antes de
seleccionar componentes son menores de 0.04 % en desplazamientos, 0.03 % en
reacciones, 2.2 % en fibras CBAR y 1.3 % en CROD. En 161--166 la tensión axial
centroidal ajustada al plano C--F coincide entre aproximadamente 0.01 % y 1.7 %.
Las fibras alejadas del eje neutro seleccionadas pasan entre 1 % y 7 %; las
fibras casi neutras y los tramos más severos de 165/166 se conservan como
diagnóstico abierto de reparto de momentos, sin `xfail` ni tolerancias laxas.
La matriz de casos y la auditoría requerida están en la
[guía de cierre](../tests/solver_tst/validation_cases/NONLINEAR_FINAL_CLOSURE_CASES.md).
