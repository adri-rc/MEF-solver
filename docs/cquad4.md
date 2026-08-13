# Desarrollo y validación de CQUAD4

## Estado actual

El elemento se mantiene en dos implementaciones independientes:

- `CQUAD4NastranElementSolver`: reconstrucción de la rigidez de Nastran y
  formulación activa en `LinearStaticSolver`;
- `CQUAD4ReissnerMindlinElementSolver`: formulación isoparamétrica clásica de
  Reissner-Mindlin con integración completa `2 x 2`.

`CQUAD4ElementBase` contiene únicamente la geometría, el sistema local, las
matrices constitutivas, `K6ROT` y la recuperación de tensiones comunes. Las dos
clases producen matrices `24 x 24`, pero no comparten el cálculo de rigidez.

## Formulación Reissner-Mindlin

La implementación clásica integra con una regla `2 x 2` los bloques de
membrana, flexión y cortante transversal. Es una referencia teórica útil y no
pretende reproducir las decisiones internas de Nastran. Como es habitual en un
elemento de bajo orden, puede presentar bloqueo de cortante en placas delgadas.

## Reconstrucción Nastran

Para elementos planos, la ruta Nastran combina:

1. deformación uniforme evaluada en el centro;
2. dos modos bilineales ortogonales a traslaciones y campos lineales;
3. estabilización *hourglass* basada en las dos líneas medias del
   cuadrilátero, con coeficiente isotrópico `E/12` y normalización del modo
   bilineal cuando las diagonales no se bisecan;
4. el mismo operador para membrana y curvaturas, con el mapeo
   `p = [theta_y, -theta_x]`;
5. el bloque transversal de MacNeal construido en flexibilidad aumentada:
   cuatro deformaciones de cortante selectivas se reconstruyen con las
   direcciones duales de las aristas (*edge-string shear*) y se combinan con la
   flexibilidad residual de flexión y, si existe `MID3`, con su flexibilidad
   física de cortante;
6. penalización *drilling* opcional cuando `K6ROT > 0`.

En elementos alabeados se forma primero la rigidez sobre el plano que pasa por
los puntos medios de las aristas. Después se aplica una transformación de
equilibrio que traslada los resultantes a los nodos reales:

- las fuerzas de membrana se resuelven sobre las dos aristas incidentes y
  generan fuerzas transversales autoequilibradas (*kick-off forces*);
- los momentos se hacen tangentes al plano nodal y su componente normal se
  equilibra mediante fuerzas en el plano.

Si `W` transforma los resultantes del plano medio, la operación implementada es:

```text
K_warped = W K_mean W^T
K_global = T^T K_warped T
```

## Validación frente a DMIG

La prueba `test_E2E_CQUAD_NAstran_K_solution.py` lee cada deck, ensambla la
matriz y la compara con la DMIG `KAAX` exportada a CSV. Los casos planos comparan sus GDL
físicamente activos con tolerancia del `1 %`; el caso 54 compara los 24 GDL con
tolerancia del `3 %`, porque incluye la reconstrucción pública de alabeo y no
el código propietario de la versión comercial usada para generar la DMIG.

Error relativo de Frobenius observado:

| Caso | Geometría y bloque | Nastran | Reissner-Mindlin |
|---|---|---:|---:|
| 51 | distorsionada, membrana | `0.039670 %` | `12.557935 %` |
| 52 | distorsionada, flexión | `0.039670 %` | `12.557935 %` |
| 53 | plana sesgada, membrana + flexión | `0.023385 %` | `25.279406 %` |
| 54 | inclinada y alabeada, matriz `24 x 24` | `2.079993 %` | `45.074180 %` |
| 55 | inclinada 3D, `K6ROT=100`, matriz `24 x 24` | `0.028215 %` | -- |
| 65 | cortante aislado con `MID3`, `K6ROT=0` | `5.74e-9 %` | -- |
| 67 | MID1/MID2/MID3, `K6ROT=100` | `0.000298 %` | -- |
| 69-71, 73 | pares MID3/G3/TST rectangulares | `< 7.3e-9 %` | -- |
| 72 | MID3 vacío, `K6ROT=100` | `0.000121 %` | -- |
| 74-75 | cuadrilátero plano sesgado, con/sin MID3 | `0.05436 %` | -- |
| 81 | paralelogramo oblicuo 3D plano, propiedades completas | `0.0000103 %` | -- |
| 83 | mismo paralelogramo, MID3 aislado | `0.000314 %` | -- |
| 89 | mismo paralelogramo, MID3 vacío | `0.0000889 %` | -- |

Sin la corrección de alabeo, el error Nastran del caso 54 era `38.551436 %`.
La corrección lo reduce a `2.079993 %`; por tanto, el caso queda cubierto por el
criterio de aceptación del `3 %`, aunque no existe igualdad término a término.

### Flexibilidad transversal validada por los casos 65-75

La expresión implementada para el subespacio `[T3,R1,R2]` es:

```text
K_s = D_s^T (Z_b + Z_s)^-1 D_s
```

`D_s` contiene cuatro deformaciones selectivas: `gamma_xz` en
`(0,+/-1/sqrt(3))` y `gamma_yz` en `(+/-1/sqrt(3),0)`. `Z_b` es la
flexibilidad residual de flexión QUAD4 con corrección de aspecto
`epsilon=0.025`, identificada con las DMIG de MSC. Para geometría rectangular,
cuando `MID3` existe la flexibilidad física se añade diagonalmente:

```text
Z_s(i,i) = 1 / (2 detJ_i G3 t (TS/T))
```

Por tanto, `MID3` aporta **flexibilidad en serie**, no una rigidez
Reissner-Mindlin en paralelo: reducir `G3` a la mitad duplica esa corrección y
reducir `TS/T` de `5/6` a `1/2` la multiplica por `5/3`. Los casos 70, 71 y 73
verifican independientemente estas escalas.

Un `MID3` vacío no activa flexibilidad física de cortante: `SHRFLX` solo es
cierto si `MID3 > 0`, de modo que `Z_s=0`. Sin embargo, esto no significa que
`MID2` sea irrelevante. El `QUAD4D` público usa internamente `MID2` al recorrer
la tercera matriz material cuando `MID3` está vacío y existe flexión; sus
propiedades alimentan la corrección residual/torsional, no una flexibilidad
física `MID3`. Esta es la interpretación implementada: el bloque conserva
`Z_b` y no inventa `Z_s`.

Un `MID2` vacío tiene una consecuencia distinta: la propiedad es solo de
membrana y el bloque completo `[T3,R1,R2]` debe ser nulo. El caso 51 compara
ahora los 24 GDL, y una prueba unitaria protege expresamente esta semántica
para impedir que `MID1` se reutilice de forma implícita como material de
flexión.

Los casos rectangulares 65, 69-71 y 73 coinciden término a término al orden de
`1e-10`; los casos completos 67 y 72 quedan por debajo de `3e-6` de error
relativo. Los F06 66 y 68 validan además desplazamientos, cargas, reacciones,
fuerzas, momentos y tensiones. En 74-75 queda un error global del `0.054 %`,
prácticamente igual con y sin MID3. Esto muestra que el residual de esos dos
casos pertenece a la geometría irregular y no a la flexibilidad física de
`MID3`.

Los casos 81, 83 y 89 cierran el paralelogramo oblicuo 3D plano. Proyectar sus
DMIG sobre las cuatro deformaciones selectivas reveló términos fuera de la
diagonal entre las dos muestras de cada dirección, presentes tanto en la
flexibilidad residual como en la física. Se introdujo una métrica simétrica de
pareja: conserva el modo promedio y escala el modo antisimétrico relacionado
con torsión según el coseno entre las dos líneas medias. La continuación usada
es positiva y acotada:

```text
c = (a · b) / (|a| |b|)
mu(c) = 1 + 0.66204667784 c^2 / (1 + 7.29631879857 c^2)
P(mu) = 1/2 [[1+mu, 1-mu], [1-mu, 1+mu]]
```

`P(mu)` se aplica por congruencia a la flexibilidad física y `mu` multiplica
los parámetros residuales `alpha/beta`. Es identidad para un rectángulo y
mantiene positiva la flexibilidad para cualquier ángulo. Los errores completos
de 81/83/89 son ahora `1.03e-7`, `3.14e-6` y `8.89e-7`; sus bloques
`[T3,R1,R2]` quedan en el mismo orden. Los tests verifican además cada término
con tolerancia absoluta, no solo la norma de Frobenius. La estructura coincide
con los pares `SFCTX1/SFCTX2` y `SFCTY1/SFCTY2` del `QUAD4D` público, pero las
constantes son una identificación contra DMIG comercial, no una constante
documentada por MSC. Los casos 90/91 son una validación predictiva de
oblicuidad alta y no deben considerarse cerrados hasta incorporar sus PCH/F06.

La interpolación de cortadura no usa directamente las derivadas bilineales
cartesianas. Para cada arista se forma una dirección dual a las dos aristas
incidentes mediante los cosenos de sus ángulos internos, igual que las
variables `UEV/ANGLEI/EDGSHR` de `Q4BMGD` en NASTRAN-95. Este cambio reduce el
error de desplazamientos del panel plano sesgado 76 de `39.66 %` a `0.0116 %`
en norma de Frobenius.

La recuperación impresa de `QX/QY` es una operación distinta de la rigidez.
En los casos planos 60, 61 y 76, los desplazamientos tienen errores globales
de `0.055 %`, `0.031 %` y `0.0116 %`, mientras componentes pequeñas de cortante
centroidal alcanzan diferencias relativas de `48 %`, `8.4 %` y `34 %`. Por
ello los tests aplican tolerancias específicas solamente a `QX/QY`; membrana,
momentos y tensiones conservan criterios mucho más estrictos. El caso 84
refuerza esta conclusión: desplazamientos y reacciones difieren `0.00122 %` y
`0.00288 %`, y `N/M` y tensiones menos de `0.003 %`, mientras un `QX` pequeño
del elemento 4 difiere `0.00977` en valor absoluto (`72 %` relativo). No se ha
modificado la rigidez para forzar esa salida de postproceso.

## Parámetros relacionados

- `K6ROT`: multiplica la penalización de la rotación *drilling* `R3`. No es
  rigidez física de membrana/flexión ni corrige mecanismos estructurales.
- `SNORM=beta`: Nastran construye en cada GRID una normal común a las facetas
  adyacentes cuyo ángulo con ella sea inferior a `beta`; el valor comercial
  predeterminado es `20` grados. Cambia el sistema usado para las rotaciones y
  por ello el acoplamiento membrana-flexión de superficies facetadas. No cambia
  la constitutiva `MAT1`. Este solver aún usa una normal independiente por
  elemento; para comparación equivalente debe usarse `PARAM,SNORM,0.0`.
- `MID1/T`: material y espesor de membrana.
- `MID2` y `12I/T^3`: material y escala de inercia de flexión.
- `MID3` y `TS/T`: material y espesor efectivo de cortante físico. Vacío no
  activa `Z_s`; el valor vacío de `TS/T` se interpreta como `5/6`.
- `THETA/MCID`: orientación de la constitutiva. Con `MAT1` isótropo no debe
  rotar los resultantes impresos fuera de los ejes del elemento.

La comparación puede reproducirse con:

```powershell
.\.venv\Scripts\python.exe tests/solver_tst/elements/run_cquad4_nastran_comparison.py
```

La opción `--export-dir ruta` guarda las matrices calculadas sin modificar las
referencias de validación.

## Fuentes y trazabilidad

- R. H. MacNeal, [A simple quadrilateral shell element](https://www.sciencedirect.com/science/article/pii/0045794978900202), *Computers & Structures*, 8(2), 175-183, 1978, DOI `10.1016/0045-7949(78)90020-2`.
- R. H. MacNeal, [Derivation of element stiffness matrices by assumed strain distributions](https://www.sciencedirect.com/science/article/pii/002954938290262X), *Nuclear Engineering and Design*, 70(1), 3-12, 1982, DOI `10.1016/0029-5493(82)90262-X`.
- NASA, [nasa/NASTRAN-95](https://github.com/nasa/NASTRAN-95): código histórico público; para CQUAD4 son relevantes `mis/quad4s.f`, `mis/quad4d.f`, `mis/q4bmgs.f`, `mis/q4bmgd.f`, `mis/trplms.f` y `mis/trplmd.f`.
- Hexagon, [MSC Nastran 2021.4 Reference Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2021.4_Reference_Guide/raw/resource/enus/MSC_Nastran_2021.4_Reference_Guide.pdf): documenta que `SNORM` genera una normal única promediada en cada GRID compartido y que su valor predeterminado es `20.0`.
- NASA NTRS, [Benchmarking the QUAD4/TRIA3 element](https://ntrs.nasa.gov/citations/19940013364): documenta diferencias entre variantes Nastran en endurecimiento de cortante e integración numérica.
- NASA NTRS, [informe 19900011106](https://ntrs.nasa.gov/citations/19900011106): describe el cortante transversal del MSC/NASTRAN QUAD4 mediante aproximación *string-net* y flexibilidad aumentada.
- B. P. Naganarayana y G. Prathap, [Force and moment corrections for the warped four-node quadrilateral plane shell element](https://repository.ias.ac.in/39268/1/138_PUB.pdf), *Computers & Structures*, 33(4), 1107-1115, 1989.

Las dos publicaciones de MacNeal describen la base conceptual de QUAD4 y las
deformaciones asumidas. NASTRAN-95 permite rastrear una implementación pública
histórica, pero no es el código de MSC Nastran actual. La ruta `nastran` de este
proyecto es, por ello, una reconstrucción contrastada contra DMIG, no una copia
literal de una formulación comercial propietaria.

## Validación pendiente

- ejecutar los casos 90/91 para comprobar la extrapolación angular acotada en
  una geometría de oblicuidad alta no usada para identificar sus constantes;
- reproducir exactamente la recuperación centroidal comercial de `QX/QY` sin
  degradar la rigidez ya validada;
- casos alabeados adicionales para separar error de forma y error de alabeo;
- reproducir con mayor fidelidad la geometría irregular de QUAD4 para reducir
  el `0.054 %` residual de los casos 74-75;
- materiales ortótropos, laminados, `MID4`, `ZOFFS` y espesores nodales;
- *patch tests* multielemento y estudios de convergencia de placa;
- estudio de sensibilidad a aspecto, *skew*, *taper* y espesor.

El caso 82 ya cierra la hipótesis de normales del 78: con `SNORM=0` el solver
coincide con Nastran al `0.0227 %` en desplazamientos y `0.00695 %` en
reacciones. El comportamiento por defecto de `SNORM` no está implementado y
se mantiene documentado mediante el `xfail` del caso 78.
