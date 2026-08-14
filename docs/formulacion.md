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

## CBAR

CBAR usa una formulación de viga espacial con extensión axial, torsión y
flexión en los dos planos principales. El triedro local se define con el eje
entre nodos y el vector de orientación resuelto por pyNastran. La propiedad
PBAR aporta `A`, `I1`, `I2`, `I12` y `J`. Se incluyen los campos de
excentricidad y liberación implementados en el solver.

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
