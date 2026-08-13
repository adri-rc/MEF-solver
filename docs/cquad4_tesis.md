# CQUAD4: formulación, parámetros y alcance de la implementación

## Objetivo

Este capítulo describe la formulación cuadrilateral de cuatro nodos usada por
el proyecto y deja trazabilidad suficiente para reproducir o continuar su
validación. La ruta activa no es un elemento Reissner-Mindlin genérico: es una
reconstrucción numérica del CQUAD4/QUAD4 de Nastran contrastada con matrices
DMIG y resultados F06. La formulación clásica se conserva en otra clase como
referencia académica.

Cada elemento tiene 24 grados de libertad, seis por nodo:
`[Tx, Ty, Tz, Rx, Ry, Rz]`. En ejes locales se separan membrana `[Tx,Ty]`,
placa/cortante `[Tz,Rx,Ry]` y rotación *drilling* `[Rz]`.

## Geometría y ejes

El origen es el promedio de los cuatro nodos. La normal se obtiene del producto
vectorial de las diagonales y el eje local x une los puntos medios de los lados
1-4 y 2-3. Para un cuadrilátero coplanar toda la rigidez se forma en ese plano.
Para uno alabeado se usa primero el plano de puntos medios y después una
transformación de equilibrio que traslada fuerzas y momentos a los nodos
reales. Esta última ruta es aproximada y no está cerrada término a término.

## Membrana y flexión

La reconstrucción Nastran combina deformación uniforme, dos modos bilineales
ortogonales a los movimientos rígidos y una estabilización *hourglass* basada
en las líneas medias. Se usa el mismo operador geométrico para membrana y
curvaturas, con la convención de rotaciones de placa
`p=[theta_y,-theta_x]`. `MID1` gobierna membrana y `MID2` flexión; si `MID2`
está vacío el bloque completo `[Tz,Rx,Ry]` es nulo.

## Cortante transversal por flexibilidad aumentada

El bloque se expresa como:

```text
K_s = D_s^T (Z_b + Z_s)^-1 D_s
```

`D_s` contiene cuatro muestras selectivas: dos de `gamma_xz` en la línea
`xi=0` y dos de `gamma_yz` en `eta=0`. Las direcciones se reconstruyen con las
aristas incidentes (*edge-string shear*), equivalente conceptualmente a
`UEV/ANGLEI/EDGSHR` del `QUAD4D` público. `Z_b` es la flexibilidad residual de
flexión. Si PSHELL define `MID3`, se añade en serie:

```text
Z_s(i,i) = 1 / (2 det(J_i) G3 t (TS/T))
```

para un rectángulo. Por tanto, bajar `G3` o `TS/T` ablanda el elemento. Un
`MID3` vacío no activa `Z_s`, pero `MID2` sigue alimentando `Z_b`; esa es la
precisión necesaria al decir que Nastran “sustituye MID3 por MID2” internamente.

En un paralelogramo oblicuo las dos muestras de cada dirección están acopladas.
Las DMIG 81/83/89 identifican la métrica positiva:

```text
c = coseno entre las dos líneas medias
mu = 1 + 0.66204667784 c^2 / (1 + 7.29631879857 c^2)
P = 1/2 [[1+mu, 1-mu], [1-mu, 1+mu]]
```

`P` preserva el modo promedio y modifica el antisimétrico torsional. Las
constantes proceden de identificación contra DMIG de MSC, no de documentación
oficial; por eso los casos 90/91 comprueban una geometría de alta oblicuidad no
usada en el ajuste.

## K6ROT y SNORM

`K6ROT` estabiliza la rotación nodal normal mediante una penalización entre la
rotación *drilling* y la macrorrotación de membrana. No debe usarse para ocultar
mecanismos estructurales ni confundirse con rigidez física de placa.

`SNORM` es otra operación. Nastran promedia las normales de elementos shell
adyacentes a un GRID si su diferencia angular queda dentro de `beta` grados;
el valor predeterminado comercial es 20. La normal común define cómo se miden
las rotaciones compartidas y cambia el acoplamiento de una superficie
facetada. No modifica `MAT1`, `MID3` ni el módulo de cortante. El solver no
implementa todavía ese promedio; las comparaciones equivalentes usan
`PARAM,SNORM,0.0`. Los casos 78/82 y, con mayor precisión, 88 demuestran esta
diferencia.

## Recuperación de resultados

En el centro se calculan `Nx, Ny, Nxy`, `Mx, My, Mxy`, tensiones de fibras Z1
y Z2, principales y von Mises. `THETA` orienta la constitutiva, pero Nastran
informa los resultados de un `MAT1` isótropo en ejes del elemento.

`QX/QY` se recupera de las fuerzas modales del bloque aumentado. Es una etapa
de postproceso independiente de la rigidez: los casos 60, 61, 76 y 84 muestran
que componentes centroidales pequeñas pueden tener gran error relativo aunque
desplazamientos, reacciones, `N/M` y tensiones coincidan. La máxima diferencia
absoluta nueva es `0.00977` en el elemento 4 del caso 84. No debe alterarse la
rigidez ya validada para hacer coincidir artificialmente esta salida.

## Evidencia numérica

- Rectángulos 65, 69-71 y 73: coincidencia término a término en torno a
  `1e-10`; casos completos 67/72 por debajo de `3e-6`.
- Paralelogramos oblicuos 81/83/89: errores relativos completos `1.03e-7`,
  `3.14e-6` y `8.89e-7`, con comparación por bloques y por términos.
- E2E 84-88: errores vectoriales de desplazamiento entre `1.69e-8` y
  `1.22e-5`, incluidos materiales/PSHELL distintos, CBAR, CROD, CTRIA3,
  `K6ROT=50/100` y desplazamientos impuestos.
- Alabeados 54, 62-64: útiles pero aproximados; no deben presentarse como
  validación definitiva del CQUAD4 plano.

## Limitaciones y trabajo abierto

- `SNORM`, `MID4`, `ZOFFS`, `TFLAG`, espesores nodales, ortotropía y laminados
  no están implementados físicamente de forma completa.
- La corrección de alabeo necesita más DMIG independientes.
- La recuperación comercial exacta de `QX/QY` sigue abierta.
- Los casos 90/91 deben ejecutarse en Nastran antes de considerar predictiva la
  métrica angular fuera de los datos 75/89.

## Fuentes

- R. H. MacNeal, *A simple quadrilateral shell element*, 1978.
- R. H. MacNeal, *Derivation of element stiffness matrices by assumed strain
  distributions*, 1982.
- NASA, repositorio público NASTRAN-95, especialmente `mis/quad4d.f` y
  `mis/q4bmgd.f`.
- NASA NTRS 19900011106, descripción *string-net* y flexibilidad aumentada.
- MSC Nastran 2021.4 Reference Guide, secciones de shell normals y `SNORM`.

Las referencias enlazadas y el histórico detallado de casos se mantienen en
[cquad4.md](cquad4.md) y [validacion.md](validacion.md).
