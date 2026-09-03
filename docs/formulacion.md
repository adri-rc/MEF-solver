# Formulación matemática

## Hipótesis generales

El solver resuelve estática lineal con pequeñas deformaciones y materiales
isótropos:

```text
K u = F
```

Cada nodo reserva seis grados de libertad `[Tx, Ty, Tz, Rx, Ry, Rz]`. Cada
elemento calcula primero su rigidez local y la transforma mediante:

```text
u_local  = T u_global
K_global = T.T K_local T
```

La matriz global se ensambla de forma dispersa. No se aplica AUTOSPC general:
un mecanismo real o un GDL libre sin rigidez produce un error de singularidad.

## Material isótropo

`MAT1` admite dos de las constantes `E`, `G`, `nu` y completa la tercera con:

```text
G = E / (2*(1+nu))
```

Para membrana y flexión de láminas se usa tensión plana:

```text
D = E/(1-nu²) * [[1, nu, 0],
                  [nu, 1, 0],
                  [0, 0, (1-nu)/2]]
```

## CROD

CROD es una barra recta de dos nodos. En su eje local:

```text
K_axial   = E*A/L * [[1,-1],[-1,1]]
K_torsion = G*J/L * [[1,-1],[-1,1]]
```

Se recuperan deformación y tensión axial, fuerza normal, giro por unidad de
longitud, torsor y tensión torsional.

### CROD no lineal geométrico

`CRODNonlinear` usa las posiciones de referencia `X1`, `X2` y los
desplazamientos totales para construir la cuerda actual:

```text
d  = (X2 + u2) - (X1 + u1)
l  = ||d||
n  = d/l
L0 = ||X2-X1||
epsilon = (l-L0)/L0
sigma   = E*epsilon
N       = A*sigma
```

La fuerza interna traslacional en coordenadas básicas es
`[-N*n, +N*n]`. Su derivada consistente contiene las partes material y
geométrica:

```text
k = (E*A/L0) (n outer n) + (N/l) (I - n outer n)
K_e = [[ k, -k],
       [-k,  k]]
```

La segunda contribución puede ser negativa en compresión y es imprescindible
para reproducir pérdida de rigidez, grandes rotaciones y *snap-through*. El
área y el módulo permanecen referidos al estado inicial: no se implementa
gran deformación material ni actualización isócora de sección.

En cada incremento el solver resuelve:

```text
r(u, lambda) = lambda F_ext - F_int(u)
K_T(u) Delta_u = r(u, lambda)
```

La tangente se actualiza en cada iteración mediante Newton completo. La torsión
permanece desacoplada y usa `T=G*J*(theta2-theta1).n0/L0`, tal como confirman
116/117 incluso con extensión axial. CROD proyecta el giro sobre el eje de
referencia; la cinemática general de sección con rotaciones finitas se reserva
al CBAR no lineal descrito más abajo.

Con `LGDISP=-1` se conserva la dirección y rigidez de referencia; `LGDISP=1`
y `2` activan la cuerda actual. Como las cargas soportadas son conservativas y
`FOLLOWK=NO`, ambos valores producen la misma respuesta en 112/118. El diseño,
las tarjetas y límites están en [`nonlinear_crod.md`](nonlinear_crod.md).

## CBAR

CBAR usa una formulación de viga espacial con extensión axial, torsión y
flexión en los dos planos. El triedro de referencia se define con la cuerda y
el vector de orientación o `G0`. `PBAR` aporta `A`, `I1`, `I2`, `I12`, `J` y
`K1/K2`; se incluyen offsets y releases.

### CBAR no lineal geométrico

`CBARNonlinear` separa el movimiento rígido mediante una terna corrotada. Las
rotaciones nodales se calculan con el mapa exponencial. La terna sigue la
cuerda mediante la rotación más corta y añade el giro medio de las secciones
alrededor de ella. Los offsets giran con sus nodos.

La medida axial identificada con los casos 130/131 es logarítmica:

```text
epsilon_log = log(l/L0)
sigma       = E*epsilon_log
N           = A*sigma
```

Los giros relativos de extremo y `L0*epsilon_log` se introducen en la matriz
local PBAR ya validada. Las acciones locales se rotan a la base y los offsets
añaden `r x f`. La tangente se obtiene como derivada central de esta fuerza
espacial completa. Puede ser no simétrica al parametrizar Newton mediante
vectores de rotación totales.

`LGDISP=-1` conserva la cinemática de referencia; `LGDISP=1/2` activa la
cinemática anterior. `INT2=LC` usa interpolación axial lineal y transversal
cúbica Euler--Bernoulli, independientemente de `K1/K2`. `INT2=LS` usa
desplazamiento y giro lineales, curvatura constante y cortante en el punto
medio; sus factores de rigidez transversal son `K1*A*G` y `K2*A*G`. La
derivación, recuperación y validación están en
[`nonlinear_cbar.md`](nonlinear_cbar.md).

## CTRIA3

CTRIA3 combina una membrana de deformación constante con una placa
Reissner–Mindlin. Los grados de libertad físicos locales son `[u,v]` para
membrana y `[w,theta_x,theta_y]` para placa. La rotación normal se estabiliza,
si `K6ROT>0`, penalizando la diferencia entre `theta_z` y la macrorrotación de
membrana:

```text
omega = 0.5*(dv/dx - du/dy)
```

El espesor multiplica la constitutiva de membrana; `t³/12` multiplica la de
flexión; y `TS/T` interviene en el cortante transversal. PSHELL permite separar
MID1, MID2 y MID3.

Los tres campos se asignan sin herencias implícitas. `MID1` vacío anula el
bloque de membrana y `MID2` vacío anula el bloque completo de placa. `MID3` no
puede definirse sin `MID2`, de acuerdo con el `FATAL 9994` observado en MSC
Nastran. Cuando `MID2` está activo y `MID3` vacío, el solver conserva la matriz
de flexión pero deja a cero el término físico de cortante de su formulación
Reissner--Mindlin; no sustituye `MID3` por el material de flexión.

Esta corrección reproduce la selección de materiales y los bloques inactivos,
pero no convierte la formulación local en la placa nativa de MSC. En
particular, MSC conserva respuesta transversal con `MID2` y `MID3` vacío,
mientras el término de cortante explícito del modelo local es nulo. Las nuevas
referencias DMIG validan la semántica de activación, no una igualdad completa
del bloque de placa.

La recuperación centroidal tridimensional de momentos y cortantes no se
considera equivalente a MSC Nastran dentro del alcance actual. El frame T con
MID3 físico correlaciona en desplazamientos, membrana y tensiones de fibra,
pero no en `M/Q` CTRIA3. La sonda complementaria con `E(MID3)=1e-20` reduce el
error de momentos a `7.106e-5` y lleva ambos cortantes a cero ingenieril, a
costa de introducir un casi-mecanismo `T3`. Esta evidencia localiza el límite
en el tratamiento de cortante/recuperación, pero no valida físicamente el
deck degenerado ni justifica recalibrar la rigidez. Véase
[`docs/validacion.md`](validacion.md#diagnóstico-complementario-sin-mid3-efectivo).

## CQUAD4

CQUAD4 tiene cuatro nodos y 24 GDL. La clase activa reconstruye el QUAD4 de
Nastran para geometría **coplanar**. Su sistema local usa:

- origen en el promedio de los cuatro nodos;
- normal por producto de diagonales G1-G3 y G2-G4;
- eje x como bisectriz de las diagonales unitarias G1-G3 y G4-G2;
- eje y que completa el triedro.

Las funciones de forma son bilineales en `[-1,1]²`. La formulación de membrana
y flexión combina deformación uniforme central y dos modos bilineales
estabilizados:

```text
K0  = A*s*B0.T*D*B0
h   = h0 - G0.T*X.T*h0
chi = 1/(1 - ||h-h0||²/12)
M   = A*(a⊗a/||a||⁴ + b⊗b/||b||⁴)
Khg = s*L*(E*chi*M/12)*L.T
```

`s=t` en membrana y `s=(t³/12)*(12I/T³)` en flexión.

El cortante usa las aristas de MacNeal y cuatro muestras sobre las líneas
naturales centrales. Se expresa mediante una flexibilidad aumentada:

```text
Z = Z_residual_bending + Z_MID3
K_shear = S.T * solve(Z, S)
```

MID3 se añade en serie como flexibilidad. MID3 vacío no hereda MID2 y MID2
vacío no hereda MID1. La corrección de aspecto usa `epsilon=0.025`.

En elementos afines oblicuos, el modo torsional se escala con:

```text
m(c²) = 1 + 0.4100259495463487*c² /
             (1 - 0.9156556465380681*c²)
```

La función se identificó con DMIG de los casos 81/89/90 y se confirmó con el
modelo cargado 60. El caso 91 verifica la respuesta de dos paneles de alta
oblicuidad.

K6ROT extiende la misma penalización de CTRIA3 a los cuatro vértices:

```text
r_i = theta_z_i - omega_i
K6  = sum_i K6ROT*1e-6*G*t*(4*detJ_i)*(r_i outer r_i)
```

Los cuatro nodos deben cumplir una tolerancia de coplanaridad
`max(1e-10,1e-8*L)`. El alabeo se rechaza y no se corrige mediante proyección.
El suavizado SNORM tampoco está implementado; una estructura facetada debe
compararse usando normales independientes, equivalentes a `SNORM,0`.

En el centro se recuperan N, M, Q, tensiones de fibras, principales y von
Mises. La recuperación centroidal Qx/Qy puede diferir de la convención de
impresión de MSC Nastran; rigidez, desplazamientos, reacciones, N, M y fibras
están validados de forma independiente.

La derivación completa y autocontenida está en
[cquad4_tesis.md](cquad4_tesis.md).

## Condiciones de contorno y reacciones

Los GDL se dividen en libres `f` y prescritos `s`:

```text
K_ff u_f = F_f - K_fs u_s
```

Esto permite valores SPC distintos de cero. Dos valores incompatibles sobre el
mismo GDL generan un error explícito. Las reacciones se calculan con la matriz
original:

```text
R = K u - F
```
