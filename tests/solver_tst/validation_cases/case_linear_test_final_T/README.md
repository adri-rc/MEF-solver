# Caso lineal final T: sección tridimensional de frame aeronáutico

Este caso `SOL 101` es el modelo lineal tridimensional de integración para el
TFM. Representa un tramo de frame aeronáutico mediante una sección en T: el
alma o *web* perforada es la geometría del caso plano histórico y una placa
transversal forma el ala inferior de la T. El objetivo no es reproducir una
pieza certificable ni sustituir una correlación experimental, sino comprobar
en un único deck el ensamblaje 3D de CROD, CBAR, CTRIA3 y CQUAD4, cargas en
sistemas locales, momento nodal y desplazamientos impuestos.

El modelo usa el sistema coherente `mm - N - MPa`, pequeñas deformaciones y
material elástico lineal isótropo. `PARAM,AUTOSPC,NO` está escrito de forma
explícita: los únicos grados de libertad restringidos son los declarados por
el modelo, sin estabilización automática silenciosa.

![Geometría mallada del frame en T](tfm_final_airframe_frame_T_patran.png)

## Geometría e interpretación estructural

El eje longitudinal del tramo es `X`. El web ocupa el plano global `XY`, con
`Z=0`, mide `1200 x 600 mm` y conserva los dos aligeramientos circulares de
radio `150 mm`, centrados en `(300, 300, 0)` y `(900, 300, 0)`. El ala inferior
ocupa el plano `XZ`, con `Y=0`, y mide `1200 x 300 mm`. Ambos planos se cortan
en la línea longitudinal `Y=Z=0`, por lo que su sección transversal es una T.

Los CBAR del contorno del web no representan una línea arbitrariamente más
rígida. Idealizan el cap o pestaña que acompaña al alma de un frame real. Esa
pestaña aporta área axial e inercia de flexión, evita que el borde de la chapa
sea la única vía resistente y conduce las cargas alrededor de los huecos. Los
CBAR secundarios completan los caminos de carga superior e inferior entre el
contorno de los aligeramientos y el cap. Es una idealización centrada, sin
offset cap-web; introducir la excentricidad requeriría una decisión adicional
de modelización.

Los doce CBAR nuevos forman seis barras diagonales, agrupadas en tres
estaciones longitudinales (`X=0, 600, 1200 mm`) y colocadas simétricamente en
las dos caras del web. Cada diagonal se parte en dos elementos en su nodo
medio: el tramo próximo al web usa PID 22 y el tramo próximo al ala usa PID
23. En PID 23 se duplican `A`, `I1`, `I2` y `J`, simulando el aumento de
espesor del nervio al acercarse a la base del frame.

Los dos CROD sustituyen el antiguo arriostramiento triangular superior. Están
contenidos en el plano del web y unen los puntos libres A, GRID 376 y 377, con
los nodos superiores 171 y 359. Los sistemas `CORD2R 101/102` alinean el GDL
1 y las fuerzas con cada tirante. Las guías `23456` en A eliminan únicamente
los movimientos que un CROD aislado no puede resistir; no son apoyos del
frame ni remaches.

## Inventario de malla

El modelo contiene 417 GRID y 506 elementos:

| Familia | Cantidad | PID | Función |
|---|---:|---:|---|
| CQUAD4 | 256 | 10 | Paños cuadriláteros del web, `t=2 mm` |
| CTRIA3 | 128 | 11 | Transición de malla alrededor de los dos huecos, `t=2 mm` |
| CQUAD4 | 32 | 12 | Ala o base de la T, `t=4 mm` |
| CBAR | 56 | 20 | Cap perimetral, bordes longitudinales y montante central |
| CBAR | 20 | 21 | Refuerzos secundarios alrededor de los aligeramientos |
| CBAR | 6 | 22 | Mitades de los tres nervios, lado próximo al web |
| CBAR | 6 | 23 | Mitades de los tres nervios, lado próximo al ala |
| CROD | 2 | 30 | Tirantes de introducción de carga en A |

Por familias, son 288 CQUAD4, 128 CTRIA3, 88 CBAR y 2 CROD. Los 12 CBAR de
PID 22/23 son los nervios tridimensionales añadidos; los otros 76 conservan
la rigidización del web del modelo plano de partida.

## Materiales y propiedades

Todos los elementos usan `MAT1 1`, con `E=71700 MPa`, `nu=0.33` y
`G=26954.887 MPa`.

| Tarjeta | Parámetros principales | Interpretación |
|---|---|---|
| PSHELL 10 | `t=2`, `MID1=MID2=MID3=1`, `TS/T=5/6` | Web CQUAD4 con membrana, flexión y cortante físicos |
| PSHELL 11 | `t=2`, `MID1=MID2=MID3=1`, `TS/T=5/6` | Anillos CTRIA3 con el mismo material |
| PSHELL 12 | `t=4`, `MID1=MID2=MID3=1`, `TS/T=5/6` | Ala CQUAD4 más gruesa |
| PBAR 20 | `A=235`, `I1=21389`, `I2=14854`, `J=1078` | T equivalente `40x4 + 25x3 mm` para el cap |
| PBAR 21 | `A=120`, `I1=6000`, `I2=3000`, `J=350` | Refuerzo secundario |
| PBAR 22 | `A=120`, `I1=6000`, `I2=3000`, `J=350` | Mitad del nervio junto al web |
| PBAR 23 | `A=240`, `I1=12000`, `I2=6000`, `J=700` | Mitad del nervio junto al ala, sección 2x |
| PROD 30 | `A=50`, `J=20`, `C=2.5` | CROD de introducción de carga |

Las propiedades de lámina usan flexión y cortante físicos porque este caso ya
no es coplanar. `PARAM,K6ROT,100` aporta la rigidez de drilling documentada y
`PARAM,SNORM,0` mantiene normales independientes sin suavizado entre planos.

## Cargas y condiciones de contorno

Las cuatro estaciones del ala representan uniones discretas mediante remaches
o puntos equivalentes de un útil de ensayo. No se empotra toda la arista:
concentrar la unión en cuatro GRID permite que el ala se deforme entre puntos,
como cabe esperar en una fijación remachada idealizada.

| Lado del ala | GRID | Coordenadas `(X,Y,Z)` mm | Condición |
|---|---|---|---|
| `Z=-150` | 405, 413 | `(300,0,-150)`, `(900,0,-150)` | Empotramiento `123456=0` |
| `Z=+150` | 422, 430 | `(300,0,150)`, `(900,0,150)` | `13456=0` y `T2=+0.2 mm` |

El desplazamiento igual en los dos remaches del segundo lado introduce un
pequeño movimiento diferencial del útil. No es un `AUTOSPC`; pertenece a SPC
120 y se combina con los empotramientos y guías mediante `SPCADD 130`.

El único subcaso aplica simultáneamente:

- `FORCE=6000 N` en cada punto A (GRID 376/377), axial y hacia fuera del CROD
  correspondiente, mediante CID 101/102;
- `MOMENT Mz=250000 Nmm` en el punto medio superior B, GRID 165;
- el desplazamiento impuesto `T2=+0.2 mm` en GRID 422/430.

![Apoyos, guías y desplazamientos impuestos en Patran](tfm_final_airframe_frame_T_patran_spc.png)

![Fuerzas en A y momento en B mostrados en Patran](tfm_final_airframe_frame_T_patran_force_moment.png)

La deformada siguiente superpone la configuración original en azul y la
deformada en blanco, con una amplificación automática aproximada de 9. La
escala es exclusivamente visual; el análisis sigue siendo lineal de pequeñas
deformaciones.

![Estructura original y deformada amplificada](tfm_final_airframe_frame_T_patran_deformed.png)

El mapa siguiente muestra `Stress Tensor / At Z2 / VONM` en las láminas. Es
un fringe nodal suavizado útil para interpretar la distribución; los valores
de aceptación se toman de las estaciones de elemento explícitas del F06, no
de la interpolación gráfica.

![Von Mises en la fibra positiva de las láminas](tfm_final_airframe_frame_T_patran_von_mises_top.png)

## Referencia MSC Nastran

El [DAT](tfm_final_airframe_frame_T.dat) se ejecutó sin modificar con MSC
Nastran Student Edition 2026.1 el 24 de agosto de 2026. El
[F06 auditado](tfm_final_airframe_frame_T.f06) alcanza `END OF JOB`, no
contiene mensajes `FATAL` ni avisos numéricos de usuario/sistema y fue
calculado con `AUTOSPC=NO`. La única aparición literal de `Warning` pertenece
al texto legal de la cabecera. MSC informa:

- número de condición de la matriz de rigidez: `3.172218e9`;
- epsilon energético: `-3.0259522e-13`;
- máximo módulo de traslación: `3.985038 mm` en GRID 434;
- máximo módulo de rotación: `0.0181810 rad` en GRID 415;
- máxima reacción de fuerza: `2520.765 N`, en GRID 413;
- máxima reacción de momento: `178885.39 Nmm`, en GRID 413.

Trazabilidad de los dos ficheros fuente:

| Fichero | SHA-256 |
|---|---|
| `tfm_final_airframe_frame_T.dat` | `1c3f4b0300815fc2106a28adfb5655494c459dc54424c2eb5ed6c804168b851a` |
| `tfm_final_airframe_frame_T.f06` | `edb202790970e5f9ec41a91d341112bdfcb9dfc55276cf36a466c46b699a99a9` |
| `tfm_final_airframe_frame_t.op2` (auditoría local) | `a3ccf11aa227ed9137dbdd2fc06e8e2222a5cf6e9f06757ae439e1950b944bac` |

## Incorporación a la suite de validación

[`test_linear_final_t_validation_case.py`](../../test_linear_final_t_validation_case.py)
añade cinco comprobaciones independientes: inventario físico, geometría en T,
cargas y restricciones, solución/equilibrio y correlación MSC. El test no lee
el F06 ni el OP2 durante pytest; todos los valores MSC seleccionados están
transcritos como literales en el módulo.

Se comparan 15 GRID de desplazamiento, las seis estaciones con reacción, los
dos CROD y 15 elementos de cada familia restante:

- CBAR: `7001, 7025, 7026, 7030, 7056, 7066, 7072, 7076, 7101, 7102, 7105,
  7106, 7108, 7109, 7112`;
- CTRIA3: `5001, 5008, 5016, 5024, 5032, 5040, 5048, 5050, 5051, 5064, 5065,
  5096, 5111, 5112, 5128`;
- CQUAD4: `1001, 1064, 1096, 1106, 1128, 1151, 1152, 1155, 1192, 1230,
  1231, 1256, 2001, 2013, 2032`.

En CBAR se transcriben las ocho acciones F06 y las tensiones normales de las
cuatro fibras en A/B. En las láminas se aceptan `Nx, Ny, Nxy` y la fila de
fibra positiva `sx, sy, txy, von Mises`. La selección cubre web, ala, contorno,
refuerzos secundarios y las dos mitades de los nervios nuevos.

## Resultados y discrepancias

El error relativo indicado es
`||resultado_MEF - resultado_MSC||2 / ||resultado_MSC||2`. La auditoría de
campos completos usa la precisión binaria del OP2; el test conserva los
valores impresos del F06 a seis cifras. El OP2 no es una dependencia de la
regresión.

| Campo completo | Error relativo L2 | Observación |
|---|---:|---|
| Traslaciones | `8.206e-3` | Gobernado por la respuesta transversal |
| Rotaciones | `2.547e-2` | Las rotaciones fuera del plano son las más sensibles |
| Reacciones de fuerza | `1.705e-4` | Cierre global correcto |
| Reacciones de momento | `3.553e-4` | Máxima diferencia absoluta `59.66 Nmm` |
| CROD, conjunto fuerza/tensión | `1.871e-6` | Axial coincide prácticamente a precisión de máquina |
| CROD, torsión | `8.125e-3` | Acción secundaria de aproximadamente `1.37 Nmm` |
| CBAR, vector completo de fuerzas | `2.744e-3` | Las acciones dominantes en el plano correlacionan mejor |
| CBAR, tensiones normales | `3.930e-3` | Máxima diferencia absoluta `0.303 MPa` |
| CQUAD4, membrana | `3.268e-5` | Web y ala |
| CQUAD4, momentos / cortantes | `6.177e-3 / 5.883e-3` | Recuperación centroidal 3D razonablemente próxima |
| CQUAD4, tensiones de fibra | `6.582e-3` | Core seleccionado `4.790e-3` |
| CTRIA3, membrana | `7.218e-6` | Coincidencia muy estrecha |
| CTRIA3, tensiones de fibra | `1.020e-2` | Core seleccionado `6.775e-3` |
| CTRIA3, momentos / cortantes | `6.820e-1 / 1.730e1` | No aceptados como correlación |

Las máximas diferencias absolutas de desplazamiento son `0.000179 mm` en T1,
`0.010416 mm` en T2 y `0.050705 mm` en T3; para rotaciones son menores de
`1.93e-4 rad`. Como puntos de lectura, GRID 434 da `T2=3.985038 mm` en MSC y
`3.993304 mm` localmente. En GRID 413, `Fy=-2509.189 N` y
`Mx=178884.3 Nmm` en MSC pasan a `-2509.618 N` y `178944.0 Nmm` en el solver.
Los dos CROD recuperan `6000 N` y `120 MPa` axiales en ambos códigos.

La discrepancia no uniforme es una observación importante del caso, no un
fallo que deba ocultarse ampliando una tolerancia global. En CBAR, `MA2`,
`MB2`, `V2` y axial presentan errores L2 inferiores a `5e-4`; las componentes
ortogonales pequeñas `MA1`, `MB1`, `V1` y torsión se mantienen como
diagnósticos con errores de hasta `0.163`. En CTRIA3, membrana y tensiones de
fibra sí forman parte de la aceptación, pero los momentos y cortantes
centroidales no quedan validados por este deck. Su divergencia debe
investigarse como recuperación de resultados antes de atribuirle precisión
Nastran; no justifica recalibrar la rigidez de los elementos.

El equilibrio del solver, transformando antes las reacciones desde cada
`GRID.CD` al sistema básico, cierra por debajo de `1e-6 N` y `2e-3 Nmm` en el
test automatizado.

## Ejecución local

Desde la raíz del repositorio:

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_linear_test_final_T\tfm_final_airframe_frame_T.dat
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\test_linear_final_t_validation_case.py
```

El procedimiento detallado para reconstruir en Patran las vistas de
geometría, apoyos, cargas, deformada y contornos está en la
[guía de capturas para el TFM](../../../../docs/TFM/guia_capturas_patran_case_linear_test_final_T.md).
