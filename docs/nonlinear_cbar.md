# Análisis geométricamente no lineal de CBAR

## Estado y alcance

`NonlinearStaticSolver` admite modelos `MAT1/PBAR/PBARN1/CBAR` con grandes
desplazamientos y rotaciones. La implementación vive en
`src/solver/elements/cbar_non_linear.py` y no modifica la ruta lineal de
`BarElementSolver`. La campaña externa CBAR comprende los casos 128--166 de MSC
Nastran SOL 400.

La referencia funcional principal es la documentación oficial de
[MSC Nastran SOL 400 2022.4](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.4_SOL_400_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.4_SOL_400_Getting_Started_Guide.pdf?save_local=true).
Para CBAR/PBAR, SOL 400 admite grandes desplazamientos y rotaciones, pero no
grandes deformaciones de sección. Las fuerzas y tensiones se expresan respecto
a ejes que corrotan con el elemento.

El alcance implementado es:

- `MAT1` elástico isótropo y sección `PBAR` constante;
- `PBARN1,SECT=S,BEH2=BAR`, con `INT2=LC` e `INT2=LS` validados;
- axial, torsión y flexión en ambos planos, incluido `I12`;
- interpolación lineal/cúbica Euler--Bernoulli para `LC` e interpolación
  lineal-cortante de un punto para `LS`;
- factores de área transversal `K1/K2` en `LS`, incluido el valor unitario
  que SOL 400 aplica cuando el campo es cero;
- offsets `WA/WB`, releases `PA/PB`, orientación por vector o `G0`;
- `FORCE`, `MOMENT`, `SPC`, `SPC1`, `SPCADD` y valores impuestos;
- cargas conservativas (`FOLLOWK=NO`) y Newton completo (`PFNT`).

Quedan fuera plasticidad, grandes deformaciones de sección, alabeo, cargas
distribuidas o seguidoras, PBEAM, CBARAO y control por longitud de arco.

## Cinemática corrotacional

Para cada elemento se guardan longitud y terna de referencia. Las rotaciones
nodales totales se interpretan como vectores de rotación y se convierten a
matrices mediante el mapa exponencial. Los offsets rotan con sus nodos:

```text
xa = Xa + ua + Ra*WA
xb = Xb + ub + Rb*WB
l  = ||xb-xa||
ex = (xb-xa)/l
```

La terna actual se construye en dos pasos. Primero se aplica la rotación más
corta que lleva el eje de referencia a la cuerda actual. Después se añade el
promedio del giro físico de ambas secciones alrededor de la cuerda. Esta
descomposición *swing--twist* evita contaminar el giro de sección con las
rotaciones de flexión.

Las orientaciones nodales relativas a esa terna proporcionan los giros locales
de extremo. El vector deformacional básico es:

```text
d = [L0*log(l/L0), theta_xb-theta_xa,
     theta_ya, theta_za, theta_yb, theta_zb]
```

El primer término es la deformación logarítmica integrada en la longitud de
referencia. Los casos axiales 130/131 la identifican directamente: para una
tensión `sigma`, Nastran devuelve `l/L0=exp(sigma/E)`. `LGDISP=-1` conserva,
en cambio, la formulación lineal exacta de referencia.

## Interpolación PBARN1

La documentación oficial de SOL 400 distingue dos elementos CBAR de dos
nodos. `INT2=LC` usa interpolación lineal a lo largo del eje y cúbica para el
desplazamiento transversal, con tres puntos de integración. Es la rama
Euler--Bernoulli: los factores `PBAR K1/K2` no cambian su rigidez. Los casos
154 y 160, idénticos salvo por esos factores, producen exactamente el mismo
resultado MSC.

`INT2=LS` usa desplazamiento transversal y giro de sección lineales y un punto
de integración. En el plano local XY, para los GDL
`[vA, theta_zA, vB, theta_zB]`, se emplean:

```text
B_curvature = [0, -1/L, 0, 1/L]
B_shear     = [-1/L, -1/2, 1/L, -1/2]

K_LS = E*I*L*(B_curvature outer B_curvature)
     + K*A*G*L*(B_shear outer B_shear)
```

El segundo plano aplica el cambio de signo de `theta_y` propio del CBAR. El
acoplamiento `I12` actúa entre las dos curvaturas constantes. Releases y
offsets se aplican después con las mismas convenciones que en `LC`.

Para un voladizo con fuerza transversal en el extremo, las dos interpolaciones
dan, en su límite lineal:

```text
delta_LC = F*L^3/(3*E*I)
delta_LS = F*L^3/(4*E*I) + F*L/(K*A*G)
theta    = F*L^2/(2*E*I)
```

La diferencia `1/4` frente a `1/3` explica por qué sustituir `LS` por la matriz
Timoshenko cerrada del CBAR lineal no reproducía el caso 142. Los pares
155/156 y 157/158 aíslan ambos planos. El caso 159 demuestra que SOL 400 usa
`K=1` cuando `K1/K2` valen cero en `LS`; en `LC` esos campos se ignoran.

## Ley de sección, fuerza interna y tangente

La rama `LC` reutiliza desde `BarElementSolver` la matriz Euler con `K1/K2`
anulados expresamente. La rama `LS` construye la matriz anterior dentro de
`CBARNonlinear`. Las seis medidas básicas se embeben en un vector local de 12
GDL y se calculan las acciones:

```text
q_local = K_local d_local
```

Fuerzas y momentos se rotan después a coordenadas básicas. Los momentos de
offset añaden `r x f`. Este vector espacial es el que se ensambla y equilibra
contra las cargas nodales.

La tangente es la derivada central de la fuerza interna completa respecto a
los 12 desplazamientos y vectores de rotación totales. Incluye por tanto la
variación de cuerda, terna, offsets y acciones geométricas. No se fuerza su
simetría: con momentos espaciales y vectores de rotación totales esas
coordenadas incrementales no forman un par conjugado que garantice una matriz
simétrica. `spsolve` resuelve directamente el sistema no simétrico.

Esta elección prioriza una fuerza física auditable y una derivada consistente.
El coste es alto: una evaluación de tangente requiere 24 evaluaciones de
fuerza por CBAR. La sustitución futura por una linealización analítica debe
reproducir antes los tests de diferencia finita y de sonda Nastran.

## Recuperación

`CBARNonlinearState` contiene:

- fuerza interna y tangente 12x12;
- longitud y terna actuales;
- deformaciones básicas/locales y acciones locales A/B;
- deformación/tensión/fuerza axial y giro/torsor;
- tensiones normales en los puntos C--F de ambos extremos;
- extremos de tensión y energía elástica de estado.

La tensión axial usa `E*log(l/L0)`. La flexión reutiliza la recuperación PBAR
lineal en la terna corrotada. Las tensiones centroidales axiales se comparan en
el E2E. La convención exacta de MSC para tensiones de fibra bajo rotación 3D
combinada no está cerrada todavía; no debe inferirse un esfuerzo centroidal
promediando C--F si los puntos de recuperación no son simétricos.

### Presentación en el informe `.rpt`

El estado final convergido de cada CBAR se presenta en cuatro bloques:

1. `SECTION RESPONSE`: fuerza y tensión axial, torsor, tensión torsional
   máxima en los puntos configurados, deformación axial logarítmica y giro por
   unidad de longitud;
2. `NONLINEAR STATE`: longitud actual entre extremos con offset, energía
   elástica elemental y tensiones normales extrema máxima/mínima entre ambos
   extremos;
3. `END FORCES AND MOMENTS`: acciones locales `Fx/Fy/Fz/Mx/My/Mz` en A y B;
4. `STRESSES AT RECOVERY POINTS`: contribuciones axial, de flexión, total y
   torsional en C--F para los dos extremos.

Las claves del diccionario `P0`--`P3` se convierten únicamente durante la
presentación a C--F; el resultado numérico no se modifica. La tabla
`NONLINEAR STATE` se activa por la presencia conjunta de `current_length` y
`strain_energy`, de modo que un CBAR lineal conserva las otras tres tablas sin
mostrar un estado no lineal ficticio. La correspondencia detallada con el F06
está en [`result_report.md`](result_report.md).

Los decks 161--166 aíslan esa convención de forma escalonada y solicitan los
resultados C--F en estados finales e intermedios. Sus F06/PCH terminan sin
`FATAL`, alcanzan factor 1.0 y ya forman parte del E2E CBAR. La tensión axial
centroidal ajustada sobre C--F coincide aproximadamente entre 0.01 % y 1.7 %;
se comparan además fibras alejadas del eje neutro. Algunas fibras casi neutras
y el reparto de momentos del segundo tramo de 165/tercer tramo de 166 siguen
siendo diagnósticos: no se han ocultado con tolerancias amplias ni `xfail`.

Los modelos mixtos 167--176 comprueban que la recuperación y el ensamblaje de
CBAR conviven con la respuesta axial corrotacional de CROD. Sus diez referencias
se mantienen en un E2E separado y alcanzan carga completa. Véase
[`NONLINEAR_FINAL_CLOSURE_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_FINAL_CLOSURE_CASES.md).

## Validación

Las pruebas se reparten en tres niveles:

- `test_cbar_non_linear.py`: once pruebas de límite lineal, matrices `LC/LS`,
  valor unitario de cortante, objetividad ante rotación rígida, deformación
  logarítmica, torsión, releases, offsets y derivada de tangente;
- `test_nonlinear_CBAR_intermediate_validation.py`: trayectorias de los casos
  132/137, flexibilidad tangente de 138/139 y estados 0.2/0.5/0.8 de los
  diagnósticos `LC/LS` 155--159;
- `test_nonlinear_CBAR_validation_E2E_.py`: diccionario explícito de los 39
  casos con desplazamientos, reacciones, esfuerzos axiales, torsor y controles
  de incremento seleccionados, incluidas fibras C--F de 161--166;
- `test_nonlinear_mixed_CROD_CBAR_validation_E2E.py`: diez referencias
  integradas 167--176 con desplazamientos, reacciones y tensión axial de un
  CBAR y un CROD por modelo.

La flexibilidad `dT2/dFy` del caso 137 difiere de la sonda MSC entre 0.1 % y
0.8 % en factores 0.2--1.0. Los 39 E2E CBAR y los diez E2E mixtos alcanzan
carga completa. Las tolerancias reflejan la complejidad: los casos básicos
quedan por debajo de
0.1 %, los ensamblajes grandes se comparan entre 0.3 % y 4 %, y el marco con
release 147 usa 7 % por una reacción transversal secundaria.

Los casos diagnósticos 154--160 aíslan `LC/LS`, cada plano de cortante y la
presencia o ausencia de `K1/K2`. Todos convergen sin `FATAL` y forman parte del
E2E final. Frente a MSC, los máximos errores vectoriales entre esos siete casos
son `8.04e-5` en desplazamientos, `6.03e-5` en reacciones y `1.52e-4` en
energía. El caso 142 pasa con `rtol=5e-4` y deja de ser `xfail`.

Los DAT y la lectura de sus salidas se describen en
[`NONLINEAR_CBAR_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CBAR_CASES.md)
y [`validacion.md`](validacion.md).
