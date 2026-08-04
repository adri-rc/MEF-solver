# Desarrollo y validación de CQUAD4

## Objetivo

El CQUAD4 pretende representar una lámina cuadrilateral isoparamétrica de cuatro nodos con membrana, flexión Reissner-Mindlin y cortante transversal. La clase activa es [`CQUAD4ElementSolver`](../src/solver/elements/cquad4_element.py) y produce una matriz `24 x 24`.

El elemento está en fase de desarrollo: su ensamblaje existe y varias propiedades algebraicas pasan las pruebas, pero la equivalencia con Nastran, la estabilidad y el postproceso no están cerrados.

## Algoritmo implementado

1. Calcula el centro geométrico y un triedro local basado en diagonales y puntos medios.
2. Proyecta los nodos sobre el plano local `(x, y)`.
3. Evalúa las cuatro funciones de forma bilineales y sus derivadas naturales.
4. Construye el Jacobiano `2 x 2`; rechaza determinantes no positivos.
5. Construye:
   - `B_m` de tamaño `3 x 8` para membrana;
   - `B_b` de tamaño `3 x 12` para curvaturas;
   - `B_s` de tamaño `2 x 12` para cortante transversal.
6. Integra membrana y flexión con `2 x 2` puntos.
7. Integra cortante transversal con `1 x 1` puntos.
8. Mapea `[w, theta_x, theta_y]` a `[Tz, Rx, Ry]` y dispersa los bloques en `24 x 24`.
9. Transforma la matriz local a global con `Tᵀ K T`.

Existe una opción `use_membrane_SRI` que integra las deformaciones normales con `2 x 2` y el cortante en el plano con `1 x 1`. El solver global no expone esta opción y utiliza siempre `False`.

## Qué está verificado

- tamaño de matriz;
- simetría dentro de las tolerancias actuales;
- triedro local en geometría canónica;
- ausencia de acoplamiento artificial entre los bloques de membrana y placa en coordenadas locales;
- uso de material MID1/MID2 separado;
- rechazo de geometría degenerada y numeración con Jacobiano no positivo;
- ejecución de la rigidez en elementos trasladados, distorsionados y ligeramente alabeados.

## Bloqueos actuales

### 1. Recuperación de tensiones desactivada

`calculate_stresses()` está dentro de un bloque de cadena multilínea y no forma parte de la clase ejecutable. `LinearStaticSolver.reactions_and_stresses()` permite solo CROD, CBAR y CTRIA3, por lo que cualquier modelo con CQUAD4 falla después de resolver los desplazamientos.

El código comentado puede servir como borrador, pero antes de activarlo hay que hacer coincidir sus convenciones de giro y signos con las utilizadas realmente en `_build_local_stiffness()`.

### 2. Rotación drilling sin rigidez

Los índices locales `[5, 11, 17, 23]` no reciben rigidez. Una penalización diagonal antigua permanece comentada y no usa `PARAM,K6ROT`. Esto introduce mecanismos cuando `Rz` no está restringido y no reproduce la estrategia usada en CTRIA3.

La solución no debe ser simplemente imponer una diagonal arbitraria. Debe elegirse y documentarse una formulación consistente —penalización basada en macrorrotación, *drilling* compatible o restricción controlada— y estudiar su influencia en energía y condicionamiento.

### 3. Modos de energía cero

En el cuadrado canónico se midieron 12 autovalores menores que `max(abs(K))*1e-9`. La prueba actual solo comprueba que al menos los seis menores sean pequeños. Además de los GDL *drilling*, la integración reducida de cortante puede introducir modos espurios (*hourglass*).

Se necesita inspeccionar los autovectores y clasificar:

- 6 movimientos rígidos;
- modos asociados a `Rz` sin rigidez;
- modos espurios de cortante/flexión;
- posibles duplicidades dependientes de la tolerancia.

### 4. Diferencia respecto a Nastran

La comparación exploratoria dio un error relativo de Frobenius del 12.56 % en membrana para el caso 51 y del 78.02 % en flexión para el caso 52. Deben aislarse las causas por bloque, no ajustarse coeficientes globales a ciegas.

Hipótesis a contrastar:

- definición exacta de ejes y plano medio en cuadriláteros distorsionados;
- signos y orden de `Rx/Ry` frente a `theta_x/theta_y`;
- formulación de curvaturas `B_b`;
- semántica PSHELL de MID2/MID3 en blanco;
- integración y formulación específica que emplea el CQUAD4 de Nastran;
- ausencia de términos de estabilización o mejora contra *locking/hourglass*;
- efecto de alabeo y proyección a un único plano.

Nastran no tiene por qué usar el mismo elemento bilineal Mindlin básico. La meta del TFM debe expresarse como reproducir casos y tendencias dentro de tolerancias definidas, o documentar conscientemente las diferencias de formulación.

### 5. Semántica PSHELL

El lector convierte MID2/MID3 vacíos a `-1`, pero `ShellProperty.__post_init__()` los reemplaza por otros materiales. En Nastran, un MID vacío puede desactivar el comportamiento correspondiente; no equivale necesariamente a reutilizar MID1. Los casos actuales recurren a materiales con `E` casi cero para aislar bloques.

Debe representarse explícitamente “sin material/sin rigidez” mediante `None`, sin sobrecargar un ID negativo ni crear materiales artificiales.

## Plan de validación propuesto

### Fase 1: pruebas cinemáticas y constitutivas

1. Verificar partición de la unidad y suma nula de derivadas de `N`.
2. Comparar Jacobiano y derivadas físicas con cálculo manual en cuadrado y trapecio.
3. Aplicar campos rígidos y comprobar deformaciones exactamente nulas en cada punto de Gauss.
4. Aplicar campos lineales de membrana y comprobar deformación constante.
5. Aplicar campos de curvatura y cortante controlados, verificando signos por separado.

### Fase 2: matrices por bloque

1. Exponer para test `K_mem`, `K_bend`, `K_shear` y `K_drill`, además de la matriz compuesta.
2. Comparar membrana del cuadrado con la solución Q4 clásica.
3. Comparar flexión y cortante con una implementación de referencia publicada.
4. Verificar rango, autovectores y energía de cada bloque.
5. Ensayar espesores desde placa delgada hasta moderadamente gruesa para detectar *locking*.

### Fase 3: referencia Nastran

1. Automatizar los casos 51 y 52 sin rutas absolutas.
2. Añadir casos puros: membrana `x`, membrana `y`, cortante `xy`, flexión `x`, flexión `y`, torsión y cortante transversal.
3. Comparar geometría cuadrada antes de introducir distorsión.
4. Incorporar rotaciones rígidas 3D y gran *offset* de coordenadas.
5. Comparar matriz, desplazamientos, reacciones y resultantes por separado.

### Fase 4: solución de malla

1. *Patch test* con varios CQUAD4.
2. Voladizo de placa con refinamientos `1x1`, `2x2`, `4x4`, `8x8`.
3. Placa simplemente apoyada con solución analítica.
4. Mallas distorsionadas y alabeadas con métricas geométricas.
5. Equilibrio de resultantes y convergencia de energía.

## Criterio de terminado

CQUAD4 puede declararse integrado cuando:

- no genera mecanismos no físicos en modelos correctamente restringidos;
- supera *patch tests* de membrana y flexión;
- tiene estrategia *drilling/hourglass* documentada y estable;
- reproduce los casos Nastran aceptados con tolerancias justificadas;
- recupera resultantes/tensiones con convenciones verificadas;
- completa `LinearStaticSolver.solve()` y el informe sin excepciones;
- sus casos forman parte de una suite pytest limpia y portable.
