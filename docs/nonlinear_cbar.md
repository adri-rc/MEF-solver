# Análisis geométricamente no lineal de CBAR

## Estado y alcance

`NonlinearStaticSolver` admite modelos `MAT1/PBAR/PBARN1/CBAR` con grandes
desplazamientos y rotaciones. La implementación vive en
`src/solver/elements/cbar_non_linear.py` y no modifica la ruta lineal de
`BarElementSolver`. La campaña externa comprende los casos 128--153 de MSC
Nastran SOL 400.

La referencia funcional principal es la documentación oficial de
[MSC Nastran SOL 400 2022.4](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.4_SOL_400_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.4_SOL_400_Getting_Started_Guide.pdf?save_local=true).
Para CBAR/PBAR, SOL 400 admite grandes desplazamientos y rotaciones, pero no
grandes deformaciones de sección. Las fuerzas y tensiones se expresan respecto
a ejes que corrotan con el elemento.

El alcance implementado es:

- `MAT1` elástico isótropo y sección `PBAR` constante;
- `PBARN1,SECT=S,BEH2=BAR`, con `INT2=LC` como variante validada;
- axial, torsión y flexión en ambos planos, incluido `I12`;
- flexibilidad de cortante `K1/K2` mediante la ley Timoshenko ya usada por el
  CBAR lineal;
- offsets `WA/WB`, releases `PA/PB`, orientación por vector o `G0`;
- `FORCE`, `MOMENT`, `SPC`, `SPC1`, `SPCADD` y valores impuestos;
- cargas conservativas (`FOLLOWK=NO`) y Newton completo (`PFNT`).

Quedan fuera plasticidad, grandes deformaciones de sección, alabeo, cargas
distribuidas o seguidoras, PBEAM, CBARAO y control por longitud de arco. La
interpolación avanzada `PBARN1,INT2=LS` se aproxima con `PBAR K1/K2`, pero no
coincide todavía con MSC en el caso corto 142; el solver emite una advertencia
y el E2E lo conserva como `xfail` estricto.

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

## Ley de sección, fuerza interna y tangente

La matriz local `MAT1/PBAR` se reutiliza desde `BarElementSolver`; así no se
duplican las convenciones de `I1/I2/I12`, `K1/K2` ni la condensación de
releases. Las seis medidas básicas se embeben en un vector local de 12 GDL y
se calculan las acciones:

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

## Validación

Las pruebas se reparten en tres niveles:

- `test_cbar_non_linear.py`: límite lineal, objetividad ante rotación rígida,
  deformación logarítmica, torsión, releases, offsets y derivada de tangente;
- `test_nonlinear_CBAR_intermediate_validation.py`: trayectorias de los casos
  132/137 y flexibilidad tangente de las sondas 138/139;
- `test_nonlinear_CBAR_validation_E2E_.py`: diccionario explícito de los 26
  casos con desplazamientos, reacciones, esfuerzos axiales, torsor y controles
  de incremento seleccionados.

La flexibilidad `dT2/dFy` del caso 137 difiere de la sonda MSC entre 0.1 % y
0.8 % en factores 0.2--1.0. Veinticinco E2E alcanzan carga completa. Las
tolerancias reflejan la complejidad: los casos básicos quedan por debajo de
0.1 %, los ensamblajes grandes se comparan entre 0.3 % y 4 %, y el marco con
release 147 usa 7 % por una reacción transversal secundaria. El caso 142
documenta la única discrepancia física abierta de `INT2=LS`.

Los casos diagnósticos 154--160 aíslan `LC/LS`, cada plano de cortante y la
presencia o ausencia de `K1/K2`. Están preparados pero requieren ejecución MSC
antes de proponer una nueva ley; no forman todavía parte del E2E final.

Los DAT y la lectura de sus salidas se describen en
[`NONLINEAR_CBAR_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CBAR_CASES.md)
y [`validacion.md`](validacion.md).
