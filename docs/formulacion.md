# Formulación matemática

## Hipótesis generales

El solver resuelve equilibrio estático lineal:

```text
K u = F
```

Se suponen pequeñas deformaciones y rotaciones, geometría sin actualización, elasticidad lineal isótropa y ausencia de efectos dependientes del tiempo. No se consideran rigidez geométrica, plasticidad, contacto ni grandes desplazamientos.

La ley isótropa completa las constantes mediante:

```text
E = 2 G (1 + nu)
```

Si se proporcionan las tres, el modelo rechaza diferencias superiores al 5 % respecto a esa relación.

## Integración y transformación

Las reglas de Gauss-Legendre disponibles son de 1, 2 y 3 puntos en `[-1, 1]`. Para los elementos 1D se integra:

```text
K_e = integral(Bᵀ D B detJ dxi)
```

Para elementos isoparamétricos 2D, la integración es tensorial en `(xi, eta)`. Las matrices locales se transforman con:

```text
u_local  = T u_global
K_global = Tᵀ K_local T
```

## CROD

Elemento recto de dos nodos, 12 GDL globales y cuatro GDL locales activos:

- desplazamiento axial `u` en ambos extremos;
- giro torsional `theta_x` en ambos extremos.

Las rigideces exactas resultantes son:

```text
k_axial   = (E A / L) [[ 1, -1], [-1, 1]]
k_torsion = (G J / L) [[ 1, -1], [-1, 1]]
```

Se integran numéricamente con un punto de Gauss y se dispersan en una matriz local `12 x 12`. La recuperación usa:

```text
epsilon_x = (u2 - u1) / L
sigma_x   = E epsilon_x
N         = A sigma_x
phi'      = (theta2 - theta1) / L
T         = G J phi'
tau       = T C / J
```

Los GDL transversales no reciben rigidez; deben quedar restringidos o estabilizados por otros elementos.

## CBAR

Viga recta de Euler-Bernoulli con dos nodos y 12 GDL. Incluye:

- axial `EA`;
- torsión de Saint-Venant `GJ`;
- flexión en los dos planos principales mediante funciones cúbicas de Hermite e inercias `I1`, `I2`.

La integración usa al menos dos puntos de Gauss. La orientación construye el triedro local a partir del eje de la barra y del vector auxiliar de CBAR.

La formulación **no** incluye todavía deformación por cortante de Timoshenko: los factores `K1` y `K2` se leen, pero no se usan. Tampoco se aplican en el solver global:

- orientación mediante `G0`;
- *offsets* `WA/WB` y semántica `OFFT`;
- liberaciones `PA/PB`;
- producto de inercia `I12` en la matriz constitutiva.

La recuperación calcula fuerzas nodales locales mediante `f_local = K_local u_local` y tensiones normales simplificadas en los puntos C/D/E/F de PBAR:

```text
sigma = sigma_axial - Mz y / I1 + My z / I2
```

La tensión torsional en esos puntos se aproxima con `T r / J`; no es una recuperación general para secciones arbitrarias.

## CTRIA3

Elemento triangular de lámina con tres nodos y 18 GDL:

- membrana CST: `[u, v]` por nodo;
- placa Reissner-Mindlin: `[w, theta_x, theta_y]` por nodo;
- rotación *drilling* `theta_z` estabilizada opcionalmente.

Para tensión plana:

```text
D_m = E/(1-nu²) [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
```

La membrana multiplica `D_m` por el espesor. La matriz de flexión incorpora `t³/12`; el cortante transversal usa el factor de corrección `5/6` y `G t`. Las matrices B son constantes para el triángulo y se evalúan en el baricentro.

`PSHELL` permite materiales separados:

- MID1: membrana;
- MID2: flexión;
- MID3: cortante transversal.

La implementación contiene una penalización *drilling* cuando `K6ROT > 0`, basada en la diferencia entre la rotación nodal y la macrorrotación de membrana. La compatibilidad numérica con Nastran se ha ajustado empíricamente y cuenta con casos específicos de comparación de `K_global`.

Advertencia: los propios comentarios del código indican que el cortante transversal no reproduce directamente Nastran. Los casos de referencia suelen usar un MID3 de rigidez prácticamente nula para aislar membrana y flexión. Esta es una condición de validación, no una equivalencia general de PSHELL.

## CQUAD4

Elemento cuadrilateral bilineal de cuatro nodos y 24 GDL. Su sistema local se define con:

- origen en el promedio de los cuatro nodos;
- `z` normal obtenida del producto vectorial de las diagonales;
- `x` entre los puntos medios de los lados 1-4 y 2-3;
- `y = z x x`.

La formulación Reissner-Mindlin suma tres bloques integrados con una regla
completa `2 x 2`:

```text
K_mem   = integral(B_mᵀ D_m B_m t dA)
K_bend  = integral(B_bᵀ D_b B_b t³/12 dA)
K_shear = integral(B_sᵀ D_s B_s (5/6)t dA)
```

- Membrana: integración completa `2 x 2`.
- Flexión: integración completa `2 x 2`.
- Cortante transversal: integración completa `2 x 2`.

La reconstrucción Nastran utiliza deformación uniforme y dos modos bilineales
estabilizados con una métrica de líneas medias. El bloque transversal emplea
la flexibilidad aumentada de MacNeal:

```text
K_s = D_s^T (Z_b + Z_s)^-1 D_s
Z_s(i,i) = 1 / (2 detJ_i G3 t (TS/T))
```

`D_s` muestrea cuatro deformaciones de cortante selectivas y `Z_b` aporta la
flexibilidad residual de flexión. `Z_s` solo existe cuando PSHELL define
`MID3`; un `MID3` vacío no hereda `MID2`. Para un elemento alabeado se
trasladan las fuerzas y momentos desde el plano medio a los nodos reales antes
de aplicar la transformación global.

Los bloques se dispersan en `[Tx, Ty]` y `[Tz, Rx, Ry]`. `PARAM,K6ROT` llega a
CQUAD4 y activa una penalización *drilling* cuando es positivo. Las clases
`CQUAD4NastranElementSolver` y `CQUAD4ReissnerMindlinElementSolver` mantienen
separadas las dos formulaciones; la geometría y la recuperación se comparten en
`CQUAD4ElementBase`.

Los casos DMIG 51-53 verifican los bloques planos con error inferior al `1 %`.
El caso 54 verifica la matriz global `24 x 24` de un elemento inclinado y
alabeado con error `2.113059 %` y criterio del `3 %`.

El caso 55 verifica una geometría 3D con `K6ROT=100` al `1.077026 %`. Los casos
65, 67 y 69-75 validan por separado `MID3`, `G3`, `TS/T`, `K6ROT` y la ausencia
de `MID3`: los rectángulos regulares quedan entre `1e-10` y `3e-6` de error
relativo y los cuadriláteros sesgados 74-75 en `0.72 %`. Los F06 66 y 68
validan también la respuesta completa.

El análisis detallado y el plan de cierre se encuentran en [cquad4.md](cquad4.md).

## Condiciones de contorno y reacciones

Cada componente de SPC se convierte en un índice restringido. Se admiten valores prescritos no nulos, aplicando el término `K_fs u_s` al lado derecho. Si dos SPC asignan valores diferentes al mismo GDL, el código actual conserva silenciosamente el último valor no nulo procesado; conviene añadir una validación de conflicto.

Las reacciones se evalúan con la matriz sin modificar:

```text
R = K u - F
```

Esto permite recuperar correctamente fuerzas y momentos en los GDL restringidos.
