# CQUAD4 plano: formulación, implementación y alcance

## 1. Objetivo del documento

Este documento describe de forma autocontenida la implementación del elemento
cuadrilateral de cuatro nodos del proyecto. La explicación procede del código
Python activo y de los casos de validación; no depende de ningún artefacto
externo a la documentación Markdown del repositorio.

El solver global utiliza `CQUAD4NastranElementSolver`, una reconstrucción de la
estructura del QUAD4 de Nastran para elementos **planos**. La formulación
clásica Reissner–Mindlin se conserva en una clase distinta únicamente como
referencia académica y no es la implementación que se ensambla por defecto.

Los archivos principales son:

- `src/solver/elements/cquad4_element.py`: geometría, transformación,
  constitutivas, `K6ROT` y recuperación compartida.
- `src/solver/elements/cquad4_nastran_element.py`: rigidez Nastran activa.
- `src/solver/elements/cquad4_reissner_mindlin_element.py`: formulación
  isoparamétrica clásica independiente.
- `tests/solver_tst/elements/CQUAD4Element_test.py`: propiedades unitarias y
  validación geométrica.
- `tests/solver_tst/elements/test_E2E_CQUAD_Nastran_K_solution.py`: comparación
  de matrices DMIG.
- `tests/solver_tst/test_validation_end_to_end_cquad4.py`: desplazamientos,
  cargas, reacciones y resultados de elemento frente a Nastran.

## 2. Hipótesis físicas y geométricas

La formulación adopta:

- pequeñas deformaciones y pequeños desplazamientos;
- material elástico lineal e isótropo `MAT1`;
- estado plano de tensiones para membrana y flexión;
- cuatro nodos ordenados consecutivamente alrededor del contorno;
- transformación isoparamétrica no degenerada y con jacobiano positivo;
- cuatro nodos coplanares dentro de la tolerancia numérica.

El elemento puede estar inclinado y orientado arbitrariamente en el espacio
tridimensional. Una estructura puede contener paneles pertenecientes a planos
distintos, siempre que cada CQUAD4 individual sea plano.

### 2.1. Restricción de coplanaridad

Si `x_i` son las coordenadas de los nodos, `x_bar` su promedio, `e_z` la normal
local y `L` la mayor longitud de lado, se comprueba:

$$
\max_i[(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_z]
-\min_i[(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_z]
\leq \max(10^{-10},10^{-8}L).
$$

Si no se cumple, el constructor genera `ValueError`. El solver no proyecta los
nodos sobre un plano medio ni aplica una corrección de alabeo. Esta limitación
es deliberada y forma parte del alcance final.

El jacobiano se valida además en los cuatro vértices naturales al construir el
elemento, con independencia de `K6ROT` y de qué bloques de `PSHELL` estén
activos. Esto permite rechazar cuadriláteros cóncavos, invertidos o con
numeración cruzada antes del ensamblaje.

## 3. Sistema local y grados de libertad

El origen local se sitúa en el promedio de las cuatro posiciones:

$$
\bar{\mathbf{x}}=\frac{1}{4}\sum_{i=1}^{4}\mathbf{x}_i.
$$

Se definen las diagonales:

$$
\mathbf{d}_{13}=\mathbf{x}_3-\mathbf{x}_1,
\qquad
\mathbf{d}_{24}=\mathbf{x}_4-\mathbf{x}_2,
\qquad
\mathbf{d}_{42}=\mathbf{x}_2-\mathbf{x}_4.
$$

La base ortonormal utilizada por el código es:

$$
\mathbf{e}_z=
\frac{\mathbf{d}_{13}\times\mathbf{d}_{24}}
     {\|\mathbf{d}_{13}\times\mathbf{d}_{24}\|},
$$

$$
\mathbf{e}_x=
\frac{\widehat{\mathbf{d}}_{13}+\widehat{\mathbf{d}}_{42}}
     {\|\widehat{\mathbf{d}}_{13}+\widehat{\mathbf{d}}_{42}\|},
\qquad
\mathbf{e}_y=\mathbf{e}_z\times\mathbf{e}_x.
$$

Es importante normalizar las dos diagonales antes de sumarlas. Si se sumaran
con su longitud original, la orientación local dependería indebidamente de la
diagonal más larga en un cuadrilátero irregular.

Las coordenadas planas del nodo `i` son:

$$
x_i=(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_x,
\qquad
y_i=(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_y.
$$

Cada nodo reserva seis grados de libertad. El vector local completo tiene 24
componentes:

$$
\boldsymbol{\delta}^{e}=
[u_1,v_1,w_1,\theta_{x1},\theta_{y1},\theta_{z1},\ldots,
u_4,v_4,w_4,\theta_{x4},\theta_{y4},\theta_{z4}]^T.
$$

- `u`, `v`: desplazamientos de membrana.
- `w`, `theta_x`, `theta_y`: desplazamiento y rotaciones de placa.
- `theta_z`: rotación normal o *drilling DOF*.

Membrana y placa quedan desacopladas en el sistema local. El acoplamiento con
`theta_z` procede únicamente de la estabilización `K6ROT`.

## 4. Interpolación isoparamétrica

El dominio natural es `(xi, eta) in [-1,1]²`, con nodos:

```text
G1 = (-1,-1)    G2 = ( 1,-1)
G4 = (-1, 1)    G3 = ( 1, 1)
```

Las funciones bilineales son:

$$
N_1=\frac14(1-\xi)(1-\eta),
\qquad
N_2=\frac14(1+\xi)(1-\eta),
$$

$$
N_3=\frac14(1+\xi)(1+\eta),
\qquad
N_4=\frac14(1-\xi)(1+\eta).
$$

La geometría y cualquier campo nodal bilineal se interpolan como:

$$
x(\xi,\eta)=\sum_{i=1}^{4}N_i x_i,
\qquad
y(\xi,\eta)=\sum_{i=1}^{4}N_i y_i.
$$

El jacobiano empleado en Python es:

$$
\mathbf{J}=
\begin{bmatrix}
\partial x/\partial\xi & \partial y/\partial\xi\\
\partial x/\partial\eta& \partial y/\partial\eta
\end{bmatrix}.
$$

Las derivadas físicas se calculan mediante:

$$
\begin{bmatrix}N_{i,x}\\N_{i,y}\end{bmatrix}
=\mathbf{J}^{-1}
\begin{bmatrix}N_{i,\xi}\\N_{i,\eta}\end{bmatrix}.
$$

Se exige `det(J) > 0` en todos los puntos usados por la formulación.

## 5. Material, constitutiva y semántica de PSHELL

Para un material isótropo se emplea la matriz de tensión plana:

$$
\mathbf{D}(E,\nu)=\frac{E}{1-\nu^2}
\begin{bmatrix}
1&\nu&0\\
\nu&1&0\\
0&0&(1-\nu)/2
\end{bmatrix}.
$$

`PSHELL` separa los mecanismos físicos:

- `MID1`: material de membrana.
- `MID2`: material de flexión.
- `MID3`: material de flexibilidad física de cortante transversal.
- `12I/T^3`: relación de inercia de flexión.
- `TS/T`: relación de espesor efectivo de cortante.

La semántica de los campos vacíos es relevante:

- `MID2` vacío desactiva la flexión; **no** hereda `MID1`.
- `MID3` vacío desactiva la flexibilidad física de cortante; **no** hereda
  `MID2`.
- Aunque `MID3` esté vacío, permanece la flexibilidad residual de flexión
  característica del QUAD4 cuando existe `MID2`.

Los factores de sección son:

$$
A_m=t,
\qquad
I_b=\frac{t^3}{12}r_I,
\qquad
t_s=t\,r_s,
$$

con:

$$
r_I=\frac{12I}{t^3},
\qquad
r_s=\frac{TS}{t}.
$$

Los valores por defecto son `r_I=1` y `r_s=5/6`.

El campo `THETA` orienta los ejes de material. Con `MAT1` isótropo la matriz
constitutiva es invariante frente a ese giro. Las fuerzas y tensiones CQUAD4 se
recuperan en los ejes propios del elemento, como en la salida usada de Nastran.

## 6. Matrices deformación–desplazamiento

### 6.1. Membrana

Con el vector:

$$
\mathbf{u}_m=[u_1,v_1,u_2,v_2,u_3,v_3,u_4,v_4]^T,
$$

la matriz de membrana es:

$$
\mathbf{B}_m=
\begin{bmatrix}
N_{1,x}&0&\cdots&N_{4,x}&0\\
0&N_{1,y}&\cdots&0&N_{4,y}\\
N_{1,y}&N_{1,x}&\cdots&N_{4,y}&N_{4,x}
\end{bmatrix}.
$$

Por tanto:

$$
\boldsymbol{\varepsilon}_m=\mathbf{B}_m\mathbf{u}_m,
\qquad
\boldsymbol{\sigma}_m=\mathbf{D}_m\boldsymbol{\varepsilon}_m.
$$

### 6.2. Curvatura de placa

Para:

$$
\mathbf{q}_p=[w_1,\theta_{x1},\theta_{y1},\ldots,
w_4,\theta_{x4},\theta_{y4}]^T,
$$

se usa la convención:

$$
\boldsymbol{\kappa}=
\begin{bmatrix}
\partial\theta_y/\partial x\\
-\partial\theta_x/\partial y\\
\partial\theta_y/\partial y-\partial\theta_x/\partial x
\end{bmatrix}
=\mathbf{B}_b\mathbf{q}_p.
$$

La rigidez de flexión se construye inicialmente como un bloque plano de ocho
grados de libertad auxiliares y después se transforma al orden
`[w, theta_x, theta_y]` mediante el mapeo de rotaciones usado en el código.

## 7. Rigidez de membrana y flexión tipo Nastran

La implementación activa **no** obtiene estos bloques con una integración de
Gauss completa `2x2`. Para membrana y curvatura se combina una deformación
uniforme central con dos modos bilineales estabilizados.

Sea:

- `B0 = B(0,0)`;
- `J0 = det(J(0,0))`;
- `A = 4 J0`;
- `s=t` para membrana;
- `s=I_b` para flexión.

La contribución uniforme es:

$$
\mathbf{K}_0=A\,s\,\mathbf{B}_0^T\mathbf{D}\mathbf{B}_0.
$$

### 7.1. Modo bilineal

Se parte del vector:

$$
\mathbf{h}_0=[1,-1,1,-1]^T.
$$

Si `G0` contiene los gradientes físicos de las funciones de forma en el centro
y `X` las coordenadas locales de los nodos, se elimina la parte contenida en
los campos lineales:

$$
\mathbf{h}=\mathbf{h}_0-\mathbf{G}_0^T\mathbf{X}^T\mathbf{h}_0.
$$

El operador modal `L` coloca `h` en las dos direcciones del plano. La energía
de los modos debe corregirse cuando las diagonales dejan de bisecarse. La
normalización usada es:

$$
\chi=\left(1-\frac{\|\mathbf{h}-\mathbf{h}_0\|^2}{12}\right)^{-1}.
$$

El denominador debe ser positivo; en caso contrario la geometría no admite la
normalización y se rechaza.

### 7.2. Métrica de líneas medias

Se definen los vectores entre puntos medios de lados opuestos:

$$
\mathbf{a}=\frac12(\mathbf{x}_2+\mathbf{x}_3-
                         \mathbf{x}_1-\mathbf{x}_4),
$$

$$
\mathbf{b}=\frac12(\mathbf{x}_3+\mathbf{x}_4-
                         \mathbf{x}_1-\mathbf{x}_2).
$$

La métrica geométrica es:

$$
\mathbf{M}=A\left(
\frac{\mathbf{a}\otimes\mathbf{a}}
     {(\mathbf{a}\cdot\mathbf{a})^2}
+
\frac{\mathbf{b}\otimes\mathbf{b}}
     {(\mathbf{b}\cdot\mathbf{b})^2}
\right).
$$

La rigidez modal queda:

$$
\mathbf{K}_{hg}=s\,\mathbf{L}
\left(\frac{E\chi}{12}\mathbf{M}\right)\mathbf{L}^T.
$$

El bloque final es:

$$
\mathbf{K}_{mb}=\mathbf{K}_0+\mathbf{K}_{hg}.
$$

Esta misma construcción se emplea para membrana y flexión, cambiando el
material, el factor de sección y el mapeo de grados de libertad. Los modos
bilineales evitan mecanismos *hourglass* sin penalizar los campos uniformes.

## 8. Cortante transversal de MacNeal

La clase Reissner–Mindlin de referencia usa las deformaciones convencionales:

$$
\gamma_{xz}=w_{,x}+\theta_y,
\qquad
\gamma_{yz}=w_{,y}-\theta_x.
$$

La clase Nastran activa no interpola directamente estas deformaciones sobre un
cuadrilátero irregular. Construye un operador dual desde las cuatro aristas,
siguiendo la estructura `UEV/ANGLEI/EDGSHR` del código público NASTRAN-95.

### 8.1. Vectores de arista

Sea `e_i` el vector unitario de la arista que parte del nodo `i`, `l_i` su
longitud y:

$$
a_i=-\mathbf{e}_i\cdot\mathbf{e}_{i-1}.
$$

El vector dual interpolado para la arista `i` es:

$$
\mathbf{s}_i=
\frac{\mathbf{e}_i+a_i\mathbf{e}_{i-1}}{1-a_i^2}N_i
+
\frac{\mathbf{e}_i+a_{i+1}\mathbf{e}_{i+1}}{1-a_{i+1}^2}N_{i+1}.
$$

Las diferencias de `w` y el promedio de las rotaciones en las aristas
incidentes generan un operador `S` de dimensión `4x12`.

### 8.2. Puntos de muestreo

Se toman dos muestras de cada dirección en las líneas naturales centrales:

```text
gamma_x: (0, -1/sqrt(3)), (0, +1/sqrt(3))
gamma_y: (-1/sqrt(3), 0), (+1/sqrt(3), 0)
```

### 8.3. Formulación en flexibilidad

La rigidez se expresa mediante:

$$
\mathbf{K}_s=\mathbf{S}^T\mathbf{Z}^{-1}\mathbf{S},
\qquad
\mathbf{Z}=\mathbf{Z}_{r}+\mathbf{Z}_{MID3}.
$$

El código no forma `Z^-1` explícitamente. Resuelve:

$$
\mathbf{Z}\mathbf{Y}=\mathbf{S},
\qquad
\mathbf{K}_s=\mathbf{S}^T\mathbf{Y}.
$$

Esto mejora la estabilidad numérica y conserva la simetría.

### 8.4. Flexibilidad residual de flexión

Si `Delta_x` y `Delta_y` son las proyecciones de las líneas medias y
`D_b,11`, `D_b,22` las rigideces de flexión:

$$
f_x=\frac{\Delta x^2}{12A D_{b,11}},
\qquad
f_y=\frac{\Delta y^2}{12A D_{b,22}}.
$$

La corrección de aspecto usa `epsilon=0.025`:

$$
\alpha=\frac{\varepsilon}
{\varepsilon+(1-\varepsilon)(\Delta x/\Delta y)^2},
$$

$$
\beta=\frac{\varepsilon}
{\varepsilon+(1-\varepsilon)(\Delta y/\Delta x)^2}.
$$

### 8.5. Corrección torsional para elementos oblicuos

Sea `c` el coseno entre las dos líneas medias:

$$
c=\frac{\mathbf{a}\cdot\mathbf{b}}
        {\|\mathbf{a}\|\,\|\mathbf{b}\|}.
$$

El factor utilizado es:

$$
m(c^2)=1+
\frac{0.4100259495463487\,c^2}
     {1-0.9156556465380681\,c^2}.
$$

Propiedades:

- `m(0)=1` para líneas medias ortogonales;
- el denominador es positivo para `0 <= c² <= 1`;
- la función es positiva, monótona y finita en el dominio;
- reproduce los controles DMIG moderados 81/89 y el control de alta
  oblicuidad 90;
- el caso cargado 60 comprueba una oblicuidad intermedia.

Esta ley es una identificación numérica propia del proyecto. No se presenta
como una constante publicada por MSC.

Los bloques residuales son:

$$
\mathbf{Z}_{r,x}=f_x
\begin{bmatrix}
1+\alpha m&1-\alpha m\\
1-\alpha m&1+\alpha m
\end{bmatrix},
$$

$$
\mathbf{Z}_{r,y}=f_y
\begin{bmatrix}
1+\beta m&1-\beta m\\
1-\beta m&1+\beta m
\end{bmatrix}.
$$

### 8.6. Flexibilidad física MID3

Cuando existe `MID3`, para cada muestra `j`:

$$
f_{s,j}=\frac{1}{2J_j t_s G_3}.
$$

El modo de pareja usa:

$$
\mathbf{P}=\frac12
\begin{bmatrix}
1+m&1-m\\
1-m&1+m
\end{bmatrix}.
$$

La congruencia con las raíces de `f_s,j` mantiene la simetría y positividad
cuando el jacobiano cambia entre las dos muestras. MID3 se añade **en serie**
como flexibilidad; por eso activarlo hace más flexible el elemento, en lugar de
añadir una rigidez Reissner–Mindlin en paralelo.

## 9. Estabilización de la rotación drilling: K6ROT

La teoría plana no aporta rigidez física independiente a `theta_z`. La
implementación extiende a cuatro nodos la penalización ya validada en CTRIA3.

En cada vértice se calcula la macrorrotación de membrana:

$$
\omega_i=\frac12(v_{,x}-u_{,y})_i.
$$

El residuo nodal es:

$$
r_i=\theta_{zi}-\omega_i.
$$

La contribución completa es:

$$
\mathbf{K}_{6ROT}=\sum_{i=1}^{4}
\left(K6ROT\,10^{-6}G_1t\right)
\left(4J_i\right)\mathbf{r}_i\mathbf{r}_i^T.
$$

Las derivadas y `J_i` se evalúan en el vértice natural correspondiente. No se
trata de añadir cuatro números a la diagonal: el residuo acopla `theta_z` con
los desplazamientos `u,v` que forman la macrorrotación.

Esta forma:

- conserva los movimientos rígidos;
- coincide con el principio usado en CTRIA3;
- estabiliza `theta_z` sin crear una rotación física independiente;
- reproduce el acoplamiento drilling–membrana de las DMIG validadas.

Si `K6ROT=0`, no se genera esta contribución. El solver no ejecuta un AUTOSPC
general ni elimina automáticamente grados de libertad inactivos. Cualquier
rotación sin rigidez debe restringirse explícitamente en el modelo; una fila
libre nula o un mecanismo estructural producen un error de matriz singular.

## 10. Ensamblaje local y transformación global

Los bloques se dispersan en la matriz local `24x24`:

$$
\mathbf{K}_{loc}=
\mathbf{K}_m+\mathbf{K}_b+\mathbf{K}_s+\mathbf{K}_{6ROT}.
$$

Los índices locales son:

- membrana: `[T1,T2]` de cada nodo;
- placa: `[T3,R1,R2]` de cada nodo;
- drilling: `[R3]`, acoplado con `[T1,T2]` por `K6ROT`.

La matriz de rotación local se forma con las filas `[e_x, e_y, e_z]` y se
repite en traslaciones y rotaciones de los cuatro nodos. La transformación a
coordenadas básicas es:

$$
\mathbf{K}_{basic}=\mathbf{T}^T\mathbf{K}_{loc}\mathbf{T}.
$$

El resultado se simetriza para eliminar ruido de redondeo:

$$
\mathbf{K}\leftarrow\frac12(\mathbf{K}+\mathbf{K}^T).
$$

Durante el ensamblaje global se aplica además la transformación asociada al
campo `CD` de cada `GRID`. Por ello los grados de libertad nodales pueden estar
definidos en sistemas distintos sin cambiar la formulación del elemento.

## 11. Recuperación de fuerzas y tensiones

La recuperación se realiza en el centro `(xi,eta)=(0,0)`.

### 11.1. Membrana

$$
\boldsymbol{\varepsilon}_m=\mathbf{B}_m(0,0)\mathbf{u}_m,
$$

$$
\boldsymbol{\sigma}_m=\mathbf{D}_m\boldsymbol{\varepsilon}_m,
\qquad
\mathbf{N}=t\boldsymbol{\sigma}_m.
$$

Se devuelven `Nx`, `Ny` y `Nxy`.

### 11.2. Momentos

$$
\boldsymbol{\kappa}=\mathbf{B}_b(0,0)\mathbf{q}_p,
$$

$$
\mathbf{M}=-I_b\mathbf{D}_b\boldsymbol{\kappa}.
$$

El signo negativo adapta la convención interna de curvaturas al criterio de
momentos usado por Nastran.

### 11.3. Tensiones de fibra

Suponiendo distribución lineal a través del espesor:

$$
\boldsymbol{\sigma}_{bot}=
\boldsymbol{\sigma}_m+\frac{6\mathbf{M}}{t^2},
$$

$$
\boldsymbol{\sigma}_{top}=
\boldsymbol{\sigma}_m-\frac{6\mathbf{M}}{t^2}.
$$

Para cada fibra se calculan tensiones principales:

$$
\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}
\pm\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2},
$$

y von Mises en tensión plana:

$$
\sigma_{VM}=\sqrt{
\sigma_x^2-\sigma_x\sigma_y+\sigma_y^2+3\tau_{xy}^2}.
$$

### 11.4. Resultantes de cortante Qx/Qy

Se calculan primero las fuerzas nodales producidas por `K_s` y el vector de
desplazamientos de placa. Después se reconstruyen por mínimos cuadrados las
cuatro fuerzas generalizadas del operador `S`. Se elimina la medida de
integración `2J_j` y se promedian las dos muestras de cada dirección.

Esta recuperación es consistente con la rigidez implementada, pero no replica
exactamente la convención propietaria que MSC Nastran usa para imprimir el
valor centroidal. Puede haber diferencias relativas grandes en componentes
pequeñas de `Qx/Qy` sin una diferencia equivalente en la matriz, los
desplazamientos, las reacciones o los momentos.

## 12. Estrategia de validación

La implementación se valida en varios niveles.

### 12.1. Pruebas unitarias

- dimensión `24x24`;
- simetría de Maxwell–Betti;
- seis modos de cuerpo rígido;
- ejes locales ortonormales;
- deformación constante de membrana;
- desacoplamiento membrana–placa;
- uso independiente de materiales de membrana y flexión;
- forma exacta de la penalización `K6ROT`;
- puntos de control y monotonía de la métrica oblicua;
- rechazo de colinealidad, jacobiano no positivo y alabeo.

### 12.2. Comparación de matrices DMIG

El error relativo de Frobenius se define como:

$$
e_F=\frac{\|\mathbf{K}_{solver}-\mathbf{K}_{Nastran}\|_F}
          {\|\mathbf{K}_{Nastran}\|_F},
$$

con:

$$
\|\mathbf{A}\|_F=\sqrt{\sum_i\sum_j A_{ij}^2}.
$$

Esta medida resume el error global, pero podría ocultar un término pequeño si
otro bloque domina la norma. Por ello cada caso aplica también:

- comparación término a término con `rtol` y `atol`;
- comparación del bloque de membrana `[T1,T2]`;
- comparación del bloque placa/cortante `[T3,R1,R2]`;
- comparación del bloque drilling `[R3]`;
- comparación del acoplamiento drilling–membrana cuando procede.

Los bloques se rotan previamente a los ejes del elemento. De este modo una
membrana mucho más rígida no puede ocultar un error de cortante o `K6ROT`.

El test permite exportar la matriz calculada en CSV mediante:

```python
"debug_export": True,
"output_k_file": "draft_nombre_matriz_python.csv",
```

### 12.3. Casos cubiertos

- 51–55: membrana, flexión, combinación, orientación 3D y `K6ROT`.
- 65–75: cortante residual y físico, `G3`, `TS/T`, ausencia de MID3,
  corrección de aspecto y sesgo.
- 81, 83 y 89: cuadriláteros afines oblicuos, MID3 activo/vacío y matriz
  completa.
- 90: control DMIG de alta oblicuidad.
- 56–58, 60–61, 66, 68 y 76–91: modelos cargados, multipropiedad, paneles
  inclinados, estructuras mixtas y desplazamientos impuestos.

Los valores de desplazamientos, cargas, reacciones y resultados de Nastran se
escriben directamente en el test E2E. El test no lee el archivo F06 durante la
ejecución.

## 13. Resultado de los casos de cierre

### 13.1. Caso 90

El caso 90 tiene mayor oblicuidad que los controles anteriores. Demostró que
la ley escalar previa era demasiado flexible. La investigación aisló el error
en la continuación geométrica del modo torsional; no fue necesario modificar:

- el operador por aristas;
- la flexibilidad residual base;
- la contribución MID3;
- la membrana o flexión bilineal;
- `K6ROT`.

Con la ley racional final:

- la matriz completa cumple error de Frobenius menor que `1e-6`;
- el bloque placa/cortante queda por debajo de `4e-6`;
- los términos del bloque cumplen `rtol=1e-5`, `atol=2e-3`.

### 13.2. Caso 91

El caso 91 contiene dos CQUAD4 coplanares de alta oblicuidad, propiedades
distintas, fuerzas y momentos. Los resultados son:

- error vectorial de desplazamientos: `2.07e-5`;
- error vectorial de reacciones: `6.30e-5`;
- fuerzas de membrana, momentos y tensiones de fibra: menos del `0.1 %` por
  componente;
- mayor diferencia de `Qx/Qy`: `47.1 %`, equivalente a `0.0216` en valor
  absoluto.

La última diferencia se clasifica como limitación de recuperación centroidal,
no como error de rigidez.

### 13.3. Casos 57 y 58

No deben emplearse como controles principales para volver a calibrar la
rigidez. Contienen un `MID3` con módulo del orden de `1e-20`:

- caso 57: `cond(K)=2.90e25`; errores vectoriales `U=1.72e-3`, `R=4.07e-3` y
  resultados `9.28e-4`;
- caso 58: `cond(K)=8.52e26`; traslaciones próximas al `0.3 %`, pero pequeñas
  rotaciones, momentos secundarios y `Qy` son sensibles; error vectorial de
  resultados `2.94e-2`.

Una diferencia relativa grande en una magnitud casi nula no implica que la
respuesta estructural dominante sea incorrecta.

### 13.4. Casos fuera de alcance

Los seis `xfail` estrictos son:

| Casos | Motivo | Resultado esperado |
|---|---|---|
| 54, 59, 62, 63, 64 | Uno o varios CQUAD4 alabeados | El constructor rechaza la geometría |
| 78 | Suavizado de normales `SNORM` predeterminado | No coincide sin implementar SNORM |

Los casos facetados equivalentes 82/88 con `SNORM,0` pasan. Estos `xfail` no
representan defectos pendientes del CQUAD4 plano; protegen la frontera del
alcance.

La colección completa observada termina con:

```text
335 passed, 6 xfailed
```

## 14. Limitaciones explícitas

La implementación no incluye:

- CQUAD4 alabeados o proyección sobre plano medio;
- suavizado de normales `SNORM` entre elementos;
- materiales ortótropos, `MAT8` o laminados;
- acoplamiento membrana–flexión de composites;
- `MID4`;
- `ZOFFS` o excentricidad de la superficie media;
- espesores nodales variables y `TFLAG`;
- cargas distribuidas o presión sobre la cara;
- recuperación exacta de la convención centroidal propietaria de `Qx/Qy`;
- análisis no lineal, grandes desplazamientos o material no lineal.

La formulación tampoco intenta reproducir todas las opciones del CQUAD4
comercial. Sólo se consideran soportados los campos cubiertos por el lector,
el modelo de datos y los tests.

## 15. Criterio de cierre

La rigidez CQUAD4 plana se considera terminada porque:

1. se verifican simetría, modos rígidos, ejes y campos uniformes;
2. las DMIG activas se comparan globalmente, por bloques y término a término;
3. se han aislado MID1, MID2, MID3, `12I/T^3`, `TS/T` y `K6ROT`;
4. existen controles rectangulares, trapezoidales, sesgados y de alta
   oblicuidad;
5. los modelos multipropiedad y mixtos reproducen desplazamientos, reacciones,
   fuerzas, momentos y tensiones;
6. las geometrías no soportadas se rechazan explícitamente.

Las mejoras futuras no bloqueantes serían:

- reproducir con mayor precisión la impresión de `Qx/Qy`;
- realizar estudios de convergencia con refinamiento de malla;
- separar el postproceso de la solución global;
- optimizar el ensamblaje disperso;
- implementar `SNORM` únicamente si lo requiere el alcance futuro.

No se recomienda ampliar la tolerancia para aceptar alabeo. Una futura
formulación no coplanar debería tratarse como un elemento nuevo, con cinemática
y validación DMIG propias.

## 16. Referencias

- R. H. MacNeal, *A Simple Quadrilateral Shell Element*, Computers &
  Structures, volumen 8, páginas 175–183, 1978.
- R. H. MacNeal, [*Derivation of Element Stiffness Matrices by Assumed Strain
  Distributions*](https://doi.org/10.1016/0029-5493(82)90262-X), Nuclear
  Engineering and Design 70(1), 3–12, 1982.
- NASA, [NASTRAN-95: rutina `quad4d.f`](https://github.com/nasa/NASTRAN-95/blob/master/mis/quad4d.f).
- Hexagon Manufacturing Intelligence, *MSC Nastran Reference Guide*.
