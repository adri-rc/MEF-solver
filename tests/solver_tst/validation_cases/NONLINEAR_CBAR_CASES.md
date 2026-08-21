# Campaña MSC Nastran SOL 400 para CBAR

## Objetivo

Los decks 128--153 generan las referencias externas necesarias antes de
implementar el CBAR corrotacional. La campaña contiene 26 casos: los 21 casos
iniciales 128--148 y cinco modelos mayores 149--153.

Cada carpeta contiene un único `.dat`. Todas las tarjetas Bulk Data presentes
van precedidas por comentarios con sus campos para que el deck pueda auditarse
sin consultar otro archivo. Los modelos usan unidades coherentes N-mm-MPa.

## Matriz de casos

| Caso | Modelo o variación | Comprobación principal |
|---:|---|---|
| 128 | CBAR axial pequeño, `LGDISP=-1` | referencia SOL 400 en límite lineal |
| 129 | flexión pequeña, `LGDISP=-1` | referencia lineal de flexión y signos |
| 130 | tracción axial finita | extensión, fuerza axial y rigidez material |
| 131 | compresión subcrítica | rigidez geométrica de compresión |
| 132 | fuerza grande en Y local | flexión sobre `I1`, giro y efecto P-delta |
| 133 | fuerza grande en Z local | flexión sobre `I2` y convención de signos |
| 134 | momento puro sobre Z local | curvatura constante en el primer plano |
| 135 | momento puro sobre Y local | curvatura constante en el segundo plano |
| 136 | torsión pura | `GJ`, giro relativo y momento torsor |
| 137 | carga combinada base | estado finito axial/biaxial/torsional |
| 138 | base más sonda `+100 N` en Y | diferencia positiva de flexibilidad tangente |
| 139 | base más sonda `-100 N` en Y | diferencia central de flexibilidad tangente |
| 140 | CBAR 3D inclinado, `G0`, cargas en `CID` | objetividad, orientación y transformación de cargas |
| 141 | sección asimétrica, `I12 != 0` | acoplamiento de los dos planos de flexión |
| 142 | viga corta, `K1/K2` | flexibilidad de cortante Timoshenko |
| 143 | desplazamiento transversal impuesto | SPC no nulo, reacción sin carga externa |
| 144 | desplazamiento/giro impuestos y FORCE | varios SPC no nulos combinados con carga |
| 145 | dos elementos, MAT1/PBAR distintos | ensamblaje e interfaz entre propiedades |
| 146 | voladizo refinado a cuatro elementos | convergencia de malla frente al caso 132 |
| 147 | marco rectangular con `PB=6` | condensación estable de un release de extremo |
| 148 | marco espacial con offsets | integración de materiales, cargas, SPC y `LROFF` |
| 149 | pórtico plano de dos plantas, 15 CBAR | ensamblaje de tres propiedades y efecto P-delta global |
| 150 | pórtico espacial de dos niveles, 20 CBAR | cargas rotadas, arriostramiento 3D y paso adaptativo |
| 151 | arco semicircular, 20 CBAR | geometría curva, cuatro propiedades y asiento impuesto |
| 152 | nave espacial arriostrada, 20 CBAR | offsets `LROFF`, carga combinada y desplazamiento impuesto |
| 153 | brazo 3D en zigzag, 15 CBAR | cinco secciones, tres materiales, `I12`, cortante y giro impuesto |

## Diagnóstico pendiente de INT2=LS: casos 154--160

La implementación reveló que el caso 142 no puede reproducirse sustituyendo
la integración `LS` de MSC por la matriz Timoshenko de `PBAR`. Se añaden siete
decks mínimos, aún pendientes de ejecución con Nastran:

| Caso | INT2 | K1/K2 | Carga | Comparación |
|---:|---|---|---|---|
| 154 | LC | 0.75/0.50 | Y + 0.5 Z | pareja directa de 142 |
| 155 | LS | 0.75/0.50 | Y | cortante/plano K1 |
| 156 | LC | 0.75/0.50 | Y | pareja de 155 |
| 157 | LS | 0.75/0.50 | Z | cortante/plano K2 |
| 158 | LC | 0.75/0.50 | Z | pareja de 157 |
| 159 | LS | 0/0 | Y + 0.5 Z | efecto LS sin factores PBAR |
| 160 | LC | 0/0 | Y + 0.5 Z | pareja de 159 |

Todos conservan 30 incrementos comunes y `OUTCTRL=(SOLUTION,INTERM)`. Al
devolver resultados se necesitan al menos DAT, F06 y PCH; los OP2 sólo son
útiles si se quiere automatizar la extracción de toda la trayectoria.

## Auditoría de la primera ejecución (2026-08-20)

Los F06 originales de los casos 128--148 terminaron normalmente: no contienen
mensajes `USER/SYSTEM FATAL`, todos completan el paso 1 y alcanzan factor de
carga 1.0. Sin embargo, esa ejecución no constituye una referencia CBAR no
lineal válida. Los DAT originales sólo contenían `PBAR`; en SOL 400 es necesario
asociar `PBARN1` para activar la formulación avanzada de CBAR. La respuesta de
los casos perturbados 137--139 resultó exactamente lineal y los F06/OP2 no
incluyeron resultados de elemento CBAR, lo que permitió detectar el problema.

Los DAT 128--148 se corrigieron con `PBARN1` (`BEH2=BAR`, `INT2=LC`, salvo
`INT2=LS` en el caso 142) y con `NLSTRESS`/`BOUTPUT`, además de solicitudes
explícitas de `STRAIN` y `ELFORCE`. La ejecución corregida cerró 23 de los 26
casos a factor de carga 1.0. Los tres restantes se auditaron así:

- 140 y 150 resolvían el primer incremento, pero la exportación intermedia
  terminaba con `USER FATAL 2025`: los ficheros incrementales de SOL 400 no
  conservaban el `CORD2R` usado por `GRID.CD`. Sus DAT dejan ahora los
  resultados en el sistema básico y mantienen `FORCE/MOMENT.CID` para validar
  la transformación de entrada.
- 148 convergía cada incremento, pero el control adaptativo agotaba
  `MAXINC=400` en el factor 0.775568. Su DAT usa ahora 100 incrementos fijos.

Los casos 140, 148 y 150 se relanzaron posteriormente. Sus F06/PCH actuales
alcanzan factor 1.0 y sustituyen las versiones fallidas descritas arriba.

Los casos 137--139 deben ejecutarse sin modificar ningún otro dato. En el
estado final, una aproximación de la flexibilidad tangente en Y del extremo es:

```text
d(uy)/d(Fy) ~= (uy_138 - uy_139) / 200 N
```

El caso 137 permite comprobar que la media de las dos soluciones perturbadas
permanece próxima al estado base.

## Referencias intermedias para implementar la tangente

Todos los decks solicitan `NLOPRM OUTCTRL=(SOLUTION,INTERM) NLDBG=ADVDBG`.
Nastran escribe por cada incremento convergido desplazamientos, cargas,
reacciones, trabajo/energía, tensiones y deformaciones no lineales de CBAR,
además del historial de iteraciones y residuos. Los casos 132 y 133 aíslan los
dos planos de flexión; 136 aísla torsión; 145--146 aíslan ensamblaje y
refinamiento; 149--153 comprueban ensamblajes mayores.

Los casos 137--139 son especialmente útiles porque comparten exactamente 40
incrementos fijos. Como Nastran multiplica también la sonda por el factor de
carga `lambda`, en cualquier incremento común distinto de cero puede
calcularse:

```text
d(uy)/d(Fy) | lambda ~= (uy_138(lambda) - uy_139(lambda)) / (200*lambda N)
```

Así se dispone de una referencia externa de flexibilidad tangente a lo largo
de toda la trayectoria, no sólo al final. El F06 no imprime la matriz tangente
global completa ni, con `PBARN1`, la tabla clásica `FORCES IN BAR ELEMENTS`;
para depurar se combinarán esta sonda, las reacciones nodales y las tablas
`NLSTRESS/BOUTPUT`. La implementación usa esas sondas en el test intermedio.
Para investigar la discrepancia `INT2=LS`, conviene añadir una pareja `+/-`
mínima sobre cada plano de cortante antes de ampliar los modelos grandes.

## Qué conservar de la ejecución

Para la auditoría inicial devolver al menos `.f06` y `.pch` junto a cada
`.dat`. El `.pch` simplifica la transcripción de vectores; el `.f06` conserva
la convergencia, los avisos y las tablas de elemento. Antes de incorporar otros
artefactos grandes, se decidirá cuáles aportan información irreemplazable.

En cada F06 revisar:

- `FATAL MESSAGE`, `USER FATAL MESSAGE`, `WARNING MESSAGE` y terminación del
  trabajo;
- `NON-LINEAR ITERATION MODULE OUTPUT` y el porcentaje/factor de carga final;
- `D I S P L A C E M E N T   V E C T O R`;
- `F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T`;
- tablas no lineales CBAR con `STRESS`, `TOTAL STRAIN` y estaciones C--F;
- energía de deformación y cualquier tabla intermedia solicitada por
  `OUTCTRL=(SOLUTION,INTERM)`.

La salida de CBAR de SOL 400 usa ejes locales corrotados. Al transcribirla hay
que anotar si una columna es esfuerzo interno de sección o acción nodal del
elemento; sus signos no tienen por qué coincidir directamente con
`forces_n1/forces_n2` del código.

## Orden de lanzamiento recomendado

Relanzar primero 128--136 y comprobar que aparecen resultados de `BOUTPUT` o
`NLSTRESS` para CBAR. Si terminan normalmente, seguir con 137--146, después
147--148 y finalmente 149--153. Un fallo aislado en 141, 142, 147, 148, 151 o
152 apunta respectivamente a `I12`, integración `LS`, releases, offsets,
asiento impuesto u offsets con una malla mayor, sin invalidar los casos básicos.

No renombrar los jobs ni sus carpetas: el número del caso coincide con los SID
de `SPC`, `LOAD` y `NLSTEP`, salvo conjuntos auxiliares documentados en el
propio deck.
