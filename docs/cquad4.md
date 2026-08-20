# CQUAD4 plano: formulación y validación final

## Estado

La implementación activa es `CQUAD4NastranElementSolver`. Reconstruye la
rigidez del QUAD4 de Nastran para cuadriláteros **coplanares**, no degenerados
y con jacobiano positivo. El plano puede tener cualquier orientación 3D y los
elementos de una estructura pueden pertenecer a planos distintos.

Los CQUAD4 alabeados quedan deliberadamente fuera del proyecto. La clase base
comprueba la distancia de los cuatro nodos al plano local y genera
`ValueError` si su rango supera `max(1e-10, 1e-8*L)`, con `L` igual a la mayor
longitud de lado. No se proyecta la geometría ni se aplica una corrección de
alabeo.

La colección actual arroja `430 passed, 6 xfailed`. Los seis `xfail` son límites
conocidos: casos 54, 59 y 62–64 alabeados, y caso 78 con suavizado `SNORM`
predeterminado. Los casos facetados equivalentes con `SNORM,0` sí pasan.

## Clases

- `CQUAD4ElementBase`: geometría plana, funciones de forma, constitutivas,
  transformación, `K6ROT` y recuperación común.
- `CQUAD4NastranElementSolver`: formulación usada por el solver global.
- `CQUAD4ReissnerMindlinElementSolver`: formulación clásica independiente,
  conservada para comparación académica; no es el elemento activo.

## Geometría y sistema local

El origen es el promedio de los cuatro nodos. La normal es el producto de las
diagonales `G1-G3` y `G2-G4`; el eje local x es la bisectriz de las diagonales
unitarias `G1-G3` y `G4-G2`; y el eje y completa el triedro.

Las funciones bilineales en `[-1,1]²` son:

```text
N1 = 1/4 (1-xi)(1-eta)       N2 = 1/4 (1+xi)(1-eta)
N3 = 1/4 (1+xi)(1+eta)       N4 = 1/4 (1-xi)(1+eta)
```

El jacobiano se evalúa en todos los puntos y vértices requeridos. Un
determinante no positivo se rechaza.

## PSHELL

- `MID1`: material de membrana.
- `MID2`: material de flexión; vacío desactiva la flexión y **no** hereda
  `MID1`.
- `MID3`: material de flexibilidad física de cortante; vacío no hereda `MID2`.
  Permanece la flexibilidad residual de flexión propia de QUAD4.
- `12I/T^3`: multiplica `t^3/12`, con valor por defecto 1.
- `TS/T`: multiplica el espesor de cortante, con valor por defecto `5/6`.
- `THETA`: orienta el material; con `MAT1` isótropo no cambia la constitutiva y
  la salida se recupera en ejes del elemento.

## Membrana y flexión

La formulación Nastran no usa una cuadratura completa 2x2. Cada bloque combina:

1. una deformación uniforme calculada en el centro;
2. dos modos bilineales `hourglass` ortogonalizados frente a los campos
   lineales;
3. una normalización por distorsión;
4. una métrica construida con las líneas que unen puntos medios de lados
   opuestos.

En forma abreviada:

```text
K0  = A*s*B0.T*D*B0
h   = h0 - G0.T*X.T*h0,  h0 = [1,-1,1,-1]
chi = 1 / (1 - ||h-h0||²/12)
M   = A*(a⊗a/||a||⁴ + b⊗b/||b||⁴)
Khg = s*L*(E*chi*M/12)*L.T
K   = K0 + Khg
```

`s=t` en membrana y `s=(t³/12)*(12I/T³)` en flexión.

## Cortante transversal

El operador cinemático procede de las cuatro aristas, siguiendo la
construcción pública `UEV/ANGLEI/EDGSHR` de NASTRAN-95. Se toman dos muestras
en cada línea natural central:

```text
(0, -1/sqrt(3)), (0, +1/sqrt(3)),
(-1/sqrt(3), 0), (+1/sqrt(3), 0)
```

La rigidez se obtiene desde una flexibilidad aumentada:

```text
Z = Z_residual_bending + Z_MID3
Kshear = S.T * solve(Z, S)
```

`MID3` se suma en serie como flexibilidad, por lo que activarlo hace más blando
el elemento. La corrección de aspecto usa `epsilon=0.025`.

En cuadriláteros afines oblicuos, el modo antisimétrico de cada pareja se
escala con:

```text
m(c²) = 1 + 0.4100259495463487*c² / (1 - 0.9156556465380681*c²)
```

`c` es el coseno entre las líneas medias. La ley vale 1 para un rectángulo,
permanece positiva en todo el dominio y se identificó con DMIG de oblicuidad
moderada (casos 81/89) y alta (caso 90). El caso cargado 60 verifica un ángulo
intermedio. Es una identificación del proyecto, no una constante publicada
por MSC.

## K6ROT

La implementación es la extensión a cuatro nodos de la ya validada en CTRIA3.
En cada vértice penaliza:

```text
r_i = theta_z_i - 0.5*(dv/dx - du/dy)_i
K6 = sum_i K6ROT*1e-6*G*t*(4*detJ_i)*(r_i outer r_i)
```

No es una rigidez diagonal aislada. La forma conserva los movimientos rígidos
y acopla el drilling nodal con la macrorrotación de membrana.

## Recuperación

En el centro se calculan fuerzas de membrana `N`, momentos `M`, tensiones en
fibras superior/inferior, principales y von Mises. Estas magnitudes, junto con
rigidez, desplazamientos y reacciones, coinciden estrechamente con Nastran.

La interfaz de resultados de los elementos shell usa una sola nomenclatura
para cada componente de fibra: `stress_top_fiber_x`,
`stress_top_fiber_y`, `tau_top_fiber_xy`, `stress_bot_fiber_x`,
`stress_bot_fiber_y` y `tau_bot_fiber_xy`. `top` representa `z=+t/2` y
`bot`, `z=-t/2` en los ejes locales. No se exponen alias basados en `Z1/Z2`,
porque esa etiqueta depende de la convención de salida del solver externo.

`Qx/Qy` se recuperan desde las cuatro fuerzas generalizadas de cortante. La
convención centroidal propietaria de MSC no está reproducida exactamente: en
el caso 91 la mayor diferencia es 47.1 % en una componente pequeña, 0.0216 en
valor absoluto, mientras `N`, `M` y fibras quedan por debajo del 0.1 %. Es una
limitación de postproceso, no de rigidez ni equilibrio.

## Validación

`tests/solver_tst/elements/test_E2E_CQUAD_Nastran_K_solution.py` compara DMIG:

- error relativo de Frobenius global;
- comparación término a término con tolerancia relativa y absoluta;
- bloques en ejes de elemento: membrana, placa/cortante y drilling;
- acoplamiento drilling–membrana cuando procede;
- exportación opcional mediante `debug_export` y `output_k_file`.

`tests/solver_tst/test_validation_end_to_end_cquad4.py` contiene directamente
las referencias de desplazamientos, cargas, reacciones y resultados; no lee el
F06 durante el test.

El caso 90 obligó a ajustar sólo la continuación geométrica del modo torsional.
Con la ley final su DMIG completa queda por debajo de `1e-6` de Frobenius, el
bloque placa/cortante por debajo de `4e-6` y sus términos cumplen
`rtol=1e-5`, `atol=2e-3`. El caso 91 obtiene errores vectoriales de
`2.07e-5` en desplazamientos y `6.30e-5` en reacciones.

Los casos 57 y 58 no son controles limpios de formulación: contienen un MID3
con módulo del orden de `1e-20` y números de condición `2.90e25` y `8.52e26`.
Sus comentarios y tolerancias documentan ahora qué componentes pequeñas son
sensibles. El caso 59 es alabeado y se conserva sólo como `xfail` de alcance.

## Limitaciones

- Sólo cuadriláteros coplanares y no degenerados.
- `MAT1` isótropo; sin laminados, `MAT8` ni acoplamiento de composite.
- Sin `ZOFFS`, espesores nodales, `TFLAG`, `MID4` o cargas de cara.
- Sin suavizado `SNORM`; usar normales independientes, equivalentes a
  `SNORM,0`, en ensamblajes facetados.
- Recuperación de `Qx/Qy` aproximada respecto a la impresión de MSC Nastran.

La derivación matemática y la correspondencia detallada con el código se
encuentran en [`cquad4_tesis.md`](cquad4_tesis.md).

## Referencias

- R. H. MacNeal, *A Simple Quadrilateral Shell Element*, Computers &
  Structures 8, 175–183, 1978.
- R. H. MacNeal, [*Derivation of element stiffness matrices by assumed strain
  distributions*](https://doi.org/10.1016/0029-5493(82)90262-X), 1982.
- NASA, [NASTRAN-95 `quad4d.f`](https://github.com/nasa/NASTRAN-95/blob/master/mis/quad4d.f).
- Hexagon, *MSC Nastran Reference Guide*.
