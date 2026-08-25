# Caso final no lineal actualizado: puente arco espacial de tres vanos

## Propósito y alcance

Este caso amplía `case_nonlinear_test_final` para validar el ensamblaje no
lineal mixto `CBAR/CROD` sobre una estructura sensiblemente mayor. El modelo
representa un fragmento idealizado de puente ferroviario de arco superior con
tres vanos, tablero espacial, péndolas, arriostramientos y bases explícitas de
estribos y pilas. Las fotografías proporcionadas se han usado únicamente como
referencia conceptual para la silueta de vanos sucesivos, las péndolas y las
cimentaciones; no aportan dimensiones, propiedades ni cargas al cálculo.

El modelo usa unidades `mm-N-MPa` y permanece dentro del subconjunto realmente
implementado por el proyecto:

- `SOL 400`, `ANALYSIS=NLSTAT` y no linealidad geométrica;
- materiales elásticos isótropos `MAT1`;
- `PBAR/PBARN1/CBAR` y `PROD/CROD`;
- cargas nodales conservativas `FORCE` y `MOMENT`;
- restricciones `SPC/SPC1/SPCADD` explícitas;
- Newton completo `PFNT` con incrementos fijos.

No se usan sólidos, láminas, contacto, plasticidad, gravedad automática,
cargas distribuidas, cargas seguidoras ni `AUTOSPC`.

## Geometría

El eje longitudinal es `X`, el ancho se desarrolla en `Y` y la vertical es
`Z`. El tablero tiene `30 000 mm` de longitud y `3 600 mm` de ancho. Hay 31
estaciones separadas `1 000 mm` y cuatro líneas longitudinales en
`Y=-1800,-600,+600,+1800 mm`.

Los apoyos se sitúan en `X=0, 8000, 20000, 30000 mm`, por lo que las tres
luces miden `8000`, `12000` y `10000 mm`. Cada plano exterior contiene tres
tramos parabólicos consecutivos con flechas respectivas de `1900`, `3200` y
`2400 mm`. Los nodos de arranque de los arcos coinciden con los largueros
exteriores del tablero.

Las cuatro líneas de apoyo incluyen dos columnas bajo los bordes del tablero,
cuatro patas inclinadas y cuatro nodos de terreno. Las cotas de terreno se
escalonan para que los estribos y las dos pilas tengan alturas distintas:

| Apoyo X [mm] | Cota inferior de columna Z [mm] | Cota de terreno Z [mm] |
|---:|---:|---:|
| 0 | -900 | -1800 |
| 8000 | -1800 | -3200 |
| 20000 | -2200 | -3600 |
| 30000 | -1100 | -2100 |

## Inventario de entidades

| Entidad | Cantidad |
|---|---:|
| GRID | 202 |
| CBAR | 144 |
| CROD | 428 |
| Elementos totales | 572 |
| MAT1 | 2 |
| PBAR + PBARN1 | 4 + 4 |
| PROD | 7 |
| FORCE | 44 |
| MOMENT | 3 |
| LOAD | 1 |

La numeración GRID se agrupa de la forma siguiente:

| Grupo | GRID | Cantidad |
|---|---|---:|
| Cuatro líneas del tablero | `1001-1031`, `1101-1131`, `1201-1231`, `1301-1331` | 124 |
| Arco exterior en `Y=-1800` | bloque `2002-2030`, salvo los arranques reutilizados | 27 |
| Arco exterior en `Y=+1800` | bloque `2102-2130`, salvo los arranques reutilizados | 27 |
| Capiteles inferiores | `3001/3002`, `3011/3012`, `3021/3022`, `3031/3032` | 8 |
| Bases de terreno | `4001-4004`, `4011-4014`, `4021-4024`, `4031-4034` | 16 |

## Elementos y propiedades asociadas

### CBAR

Los `CBAR` se reservan para los miembros cuya rigidez de flexión y torsión es
relevante. Todos los `PBAR` llevan su `PBARN1` asociado con
`SECT=S`, `BEH2=BAR` e `INT2=LC`.

| PID | EID principales | Cantidad | Representación | A [mm²] | I1 [mm⁴] | I2 [mm⁴] | J [mm⁴] |
|---:|---|---:|---|---:|---:|---:|---:|
| 101 | `10001-10030`, `10101-10130` | 60 | dos largueros/tirantes exteriores | 18000 | 1.20e9 | 4.50e8 | 9.00e7 |
| 103 | `11001-11030`, `11101-11130` | 60 | dos costillas de arco de tres vanos | 12000 | 7.00e8 | 3.20e8 | 5.00e7 |
| 105 | `12001-12008` | 8 | columnas de pilas y estribos | 24000 | 1.60e9 | 1.00e9 | 1.50e8 |
| 106 | `12101-12116` | 16 | patas inclinadas de cimentación | 16000 | 8.00e8 | 5.00e8 | 9.00e7 |

### CROD

Los `CROD` forman la parte articulada de la celosía. Son elementos axiales
bilaterales: una péndola puede entrar matemáticamente en compresión y no se
modela ni destensado ni comportamiento de cable `tension-only`.

| PID | EID principales | Cantidad | Representación | A [mm²] | MID |
|---:|---|---:|---|---:|---:|
| 201 | `20401-20454` | 54 | péndolas verticales exteriores | 800 | 2 |
| 202 | `20001-20030`, `20101-20130` | 60 | largueros interiores axiales | 2500 | 1 |
| 203 | `20201-20293` | 93 | enlaces transversales del tablero | 1400 | 1 |
| 204 | `20701-20730`, `21001-21060` | 90 | diagonales en planta exteriores e interiores | 900 | 2 |
| 205 | `20501-20527`, `20601-20630` | 57 | tirantes y diagonales entre planos de arco | 650 | 2 |
| 206 | `20801-20812` | 12 | vigas axiales de base sobre el terreno | 3000 | 1 |
| 207 | `20901-20962` | 62 | péndolas diagonales de red y enlaces a capiteles | 700 | 2 |

Los dos materiales son:

| MID | Uso | E [MPa] | G [MPa] | ν |
|---:|---|---:|---:|---:|
| 1 | miembros principales y celosía de tablero/base | 210000 | 80769.2308 | 0.30 |
| 2 | péndolas y arriostramientos | 195000 | 75000 | 0.30 |

## Idea estructural de la representación

Los largueros exteriores trabajan como tirantes del sistema de arco y dan
continuidad flexional al tablero. Las costillas superiores descargan mediante
dos familias de péndolas: una vertical en cada plano exterior y otra diagonal
que conecta las líneas interiores con el arco. Esta segunda familia crea un
camino resistente vertical directo y evita que las líneas interiores formen
un mecanismo fuera del plano.

Las diagonales en planta atan longitudinalmente los largueros interiores a los
exteriores; sin ellas, cada línea `CROD` interior podría deslizar rígidamente
en `X`. Entre los dos planos de arco se alternan tirantes transversales y
diagonales de viento. Las pilas y estribos se representan mediante columnas
`CBAR` y patas acampanadas ancladas en cuatro GRID de terreno. Es una
idealización 1D de la transferencia a la cimentación, no un modelo de hormigón
ni de interacción suelo-estructura.

## Condiciones de contorno

`PARAM,AUTOSPC,NO` queda escrito en el DAT. Todos los grados de libertad sin
rigidez física se tratan explícitamente:

- los 16 GRID de terreno tienen `T1`, `T2`, `R1`, `R2` y `R3` restringidos;
- los cuatro GRID del apoyo izquierdo tienen además `T3=0`;
- las tres bases restantes reciben asientos uniformes de `-2`, `-6` y
  `-10 mm` en `T3`, rampados con el factor de carga;
- los 62 GRID de los dos largueros interiores, conectados sólo mediante
  `CROD`, tienen restringidas sus rotaciones `R1-R3` porque esas rotaciones no
  pertenecen a la cinemática axial de la celosía.

Estas restricciones no estabilizan silenciosamente un mecanismo: forman
parte declarada de la idealización articulada y de los anclajes al terreno.
La primera versión de la malla, antes de añadir las péndolas de red y las
diagonales interiores, produjo correctamente una matriz singular con
`AUTOSPC=NO`; el DAT conservado es la topología estable corregida.

## Cargas

El conjunto `LOAD 500` combina tres familias:

| SID | Acciones | Cantidad | Factor en LOAD 500 |
|---:|---|---:|---:|
| 501 | tráfico vertical asimétrico sobre las dos líneas interiores | 24 FORCE | 1.00 |
| 502 | viento lateral en ocho estaciones de ambos arcos | 16 FORCE | 0.75 |
| 503 | frenado longitudinal y pares nodales | 4 FORCE + 3 MOMENT | 1.00 |

Las resultantes aplicadas que imprime MSC en el sistema básico son
`Fx=350000 N`, `Fy=363525 N` y `Fz=-2592000 N`. Los momentos nodales de SID
503 son dos pares alrededor de `X` (`+1.10e8` y `-8.00e7 Nmm`) y uno alrededor
de `Y` (`-1.40e8 Nmm`). Todas las direcciones son fijas en el sistema básico;
`PARAM,FOLLOWK,NO` excluye rigidez de cargas seguidoras.

## Control no lineal y ejecución

El paso usa `NLSTEP 500`, veinte incrementos fijos, hasta 80 iteraciones y 24
bisecciones permitidas por incremento, criterio residual `PV=3e-7` y Newton
completo `PFNT`:

```text
NLSTEP,500,1.0,SEVERELY
,GENERAL,80,1,24
,FIXED,20
,MECH,PV,3.E-7,3.E-7,3.E-7,PFNT
```

Para regenerar y resolver localmente:

```powershell
.\.venv\Scripts\python.exe tests\solver_tst\validation_cases\case_nonlinear_test_final_updated\bridge_model_updated.py
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_nonlinear_test_final_updated\tfm_final_nonlinear_bridge_updated.dat
```

Para MSC Nastran 2026.1 Student Edition, desde la carpeta del caso:

```powershell
& 'C:\Program Files\MSC.Software\NaPa_SE\20261\Nastran\bin\nastran.exe' `
  'tfm_final_nonlinear_bridge_updated.dat' 'scr=yes' 'old=no'
```

MSC alcanzó `LOAD STEP=1.0`, completó los 20 incrementos sin bisecciones,
imprimió `JOB CONVERGES FOR THE CURRENT STEP` y terminó con `END OF JOB`, sin
mensajes `USER FATAL` ni `SYSTEM FATAL`. La razón residual del último estado
es `6.369e-8`, por debajo del criterio solicitado. La comparación cuantitativa
y los puntos transcritos al E2E se detallan en
[`RESULTS_ANALYSIS.md`](RESULTS_ANALYSIS.md).

## Vistas de Patran

Las siguientes capturas se han obtenido realmente en Patran 2026.1 importando
el DAT y el OP2 de este caso. Se han ocultado los valores junto a los glifos
para evitar solapes; las magnitudes y resultantes exactas permanecen escritas
en las secciones anteriores y en el F06.

Vista general de la malla espacial y las bases escalonadas:

![Geometría del puente en Patran](tfm_final_nonlinear_bridge_updated_patran_geometry.png)

Apoyos físicos y asientos impuestos. Para hacer legible la imagen se muestran
los anclajes `SPC1 10/11` y los asientos `SPC 20`, pero se omite de esta vista
el `SPC1 30` que sólo bloquea las rotaciones cinemáticamente inactivas de los
62 GRID exclusivamente `CROD`:

![Apoyos y asientos en Patran](tfm_final_nonlinear_bridge_updated_patran_supports.png)

Cargas verticales de tráfico, en alzado:

![Cargas de tráfico en Patran](tfm_final_nonlinear_bridge_updated_patran_loads_traffic.png)

Cargas laterales de viento, en perspectiva para distinguir su dirección:

![Cargas de viento en Patran](tfm_final_nonlinear_bridge_updated_patran_loads_wind.png)

Frenado longitudinal y momentos nodales, también en perspectiva:

![Frenado y momentos en Patran](tfm_final_nonlinear_bridge_updated_patran_loads_braking_moments.png)

Configuración final deformada, con la geometría original en azul y la
deformada en blanco. Patran magnifica gráficamente los desplazamientos para
hacer visible la respuesta; la imagen no está a escala real de deformación:

![Deformada final en Patran](tfm_final_nonlinear_bridge_updated_patran_deformed_final.png)

No se incluye un mapa de tensiones: el registro de resultados que Patran
construyó desde el OP2 consolidado sólo expone cargas aplicadas,
desplazamientos y fuerzas de restricción en los 20 estados. Ese OP2 no contiene
tablas `force/stress/strain` para `CROD/CBAR`, aunque los resultados de elemento
sí están impresos en el F06 y son los usados en la correlación. Generar un
fringe desde esa fuente habría atribuido a Patran un campo que no está
disponible realmente.

## Archivos de evidencia

- `tfm_final_nonlinear_bridge_updated.dat`: modelo exacto resuelto por ambos
  códigos;
- `tfm_final_nonlinear_bridge_updated.f06`: listado de MSC Nastran con los 20
  estados incrementales;
- `tfm_final_nonlinear_bridge_updated.pch`: desplazamientos, cargas y
  reacciones en formato estructurado;
- `bridge_model_updated.py`: generador determinista del DAT;
- `RESULTS_ANALYSIS.md`: métricas de correlación y observaciones de alcance;
- los seis PNG `*_patran_*.png`: geometría, condiciones de contorno, tres
  familias de carga y deformada final verificadas en Patran.

Los tres artefactos numéricos de referencia quedan identificados por SHA-256:

| Archivo | SHA-256 |
|---|---|
| `tfm_final_nonlinear_bridge_updated.dat` | `7c5f374b77a0abf83aab1e4c31320425faac5bd16335d32c9c0e5f9beac05439` |
| `tfm_final_nonlinear_bridge_updated.f06` | `4180efc0678eeec374d8c14d9f7f712d57985c300564e1d1da25010a660a91be` |
| `tfm_final_nonlinear_bridge_updated.pch` | `232dce9c28494fbee56cadfbbed0d03afbf8a80172d00051b82109407311c75e` |

El OP2 se conserva localmente para construir las vistas de Patran, pero no se
considera fuente de verdad versionada. Las referencias del test quedan
escritas como literales extraídos del último bloque del F06/PCH y el test no
lee resultados externos durante su ejecución.

## Limitaciones de interpretación

Este caso valida el solver académico dentro de su alcance, no dimensiona un
puente real. El acero permanece elástico aunque alguna tensión recuperada
pueda superar un valor admisible; no se evalúan pandeo, fatiga, plastificación,
uniones, cable flojo, apoyo unilateral ni suelo. Las bases, travesaños y
péndolas son líneas sin volumen. Un mapa de Patran para esta malla debe mostrar
por separado magnitudes axiales o tensiones de fibra `CBAR/CROD`; no existe un
campo superficial continuo de von Mises equivalente al de un modelo de
láminas.
