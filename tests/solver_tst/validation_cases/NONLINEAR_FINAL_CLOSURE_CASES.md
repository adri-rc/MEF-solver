# Cierre SOL 400: fibras CBAR y modelos mixtos CBAR/CROD

## Objetivo y estado

Los casos 161--176 forman la campaña final de validación geométrica no lineal.
Los DAT pasan el preprocesado, la comprobación de rango y el solver local. Sus
F06/PCH de MSC Nastran terminan sin `FATAL`, alcanzan factor de carga 1.0 y sus
referencias finales ya están transcritas a los E2E CBAR y mixto.

Las referencias 161--166 están en
`test_nonlinear_CBAR_validation_E2E_.py`. Las 167--176 están separadas en
`test_nonlinear_mixed_CROD_CBAR_validation_E2E.py`; cada una compara al menos
un resultado de elemento CBAR y otro CROD además de la respuesta nodal.

Todos los modelos usan unidades coherentes N-mm-MPa, `SOL 400`,
`ANALYSIS=NLSTAT`, `PARAM,LGDISP,1`, cargas nodales conservativas y pasos fijos.
Cada tarjeta Bulk Data va precedida por un comentario con sus campos. Los decks
solicitan `OUTCTRL=(SOLUTION,INTERM)` y las tablas `NLSTRESS`, `BOUTPUT`,
`ELFORCE`, `STRAIN`, desplazamientos, reacciones y energía.

## Diagnóstico de tensiones de fibra CBAR: 161--166

Esta primera familia separa las operaciones que pueden cambiar la convención
de signo o el eje en que MSC evalúa los puntos C--F.

| Caso | Modelo | Qué permite identificar |
|---:|---|---|
| 161 | un CBAR sobre X, puntos C--F simétricos | referencia de signos con axial, flexión biaxial y torsión combinados |
| 162 | un CBAR inclinado en 3D | transporte de las fibras con la terna corrotada |
| 163 | sección asimétrica e `I12 != 0` | acoplamiento biaxial sin inferir el centroide promediando C--F |
| 164 | orientación `G0` y cargas en `CORD2R` | composición de orientación, transformación de carga y rotación finita |
| 165 | codo espacial de dos CBAR, dos MAT1/PBAR | convención de extremos A/B y continuidad entre ternas distintas |
| 166 | tres CBAR, `LROFF` y giro impuesto | fibras con offsets, varios elementos y SPC no nulo de rotación |

Para cerrar la recuperación no basta con el máximo/mínimo. Del último estado y
de, al menos, dos estados intermedios hay que conservar para cada CBAR:

- las tensiones normales C, D, E y F en los extremos A y B;
- fuerza axial, momentos locales 1/2, cortantes y torsor;
- deformación axial y curvaturas o deformaciones de sección disponibles;
- desplazamientos y rotaciones nodales globales;
- identificación inequívoca del factor de carga de cada bloque.

Los casos 161 y 162 fijan primero los signos. El 163 decide cómo actúa `I12` y
los puntos asimétricos. Sólo después deben interpretarse 164--166; así un error
de transformación no puede compensar por casualidad un error de signo.

## Modelos mixtos CBAR/CROD: 167--176

| Caso | Nodos | CBAR | CROD | Configuración principal |
|---:|---:|---:|---:|---|
| 167 | 3 | 1 | 1 | voladizo con tirante, FORCE y MOMENT |
| 168 | 4 | 3 | 1 | pórtico arriostrado, dos materiales |
| 169 | 6 | 2 | 3 | mástil espacial con tres vientos y propiedades CROD distintas |
| 170 | 6 | 5 | 4 | pórtico de dos vanos con cruces, varias cargas nodales |
| 171 | 8 | 8 | 6 | plataforma espacial con columnas, anillo y diagonales |
| 172 | 9 | 10 | 8 | pórtico de dos plantas y asiento vertical impuesto |
| 173 | 9 | 12 | 10 | cubierta piramidal, varias secciones e integración `LC/LS` |
| 174 | 12 | 16 | 8 | marco espacial de dos niveles, cargas en `CORD2R` y desplazamiento impuesto |
| 175 | 10 | 13 | 15 | nave espacial con offsets, tres materiales y desplazamiento impuesto |
| 176 | 12 | 20 | 10 | torre arriostrada de tres niveles, varias FORCE/MOMENT y asiento impuesto |

La secuencia incrementa deliberadamente el número de elementos hasta 30 y
combina `MAT1`, `PBAR/PBARN1`, `PROD`, `CBAR`, `CROD`, `FORCE`, `MOMENT`,
`LOAD`, `SPC1`, `SPC`, `SPCADD`, sistemas de carga, offsets y propiedades
distintas. No se usa `K6ROT`: ese parámetro estabiliza los giros de láminas y no
forma parte de la formulación de CROD/CBAR.

## Orden de ejecución y auditoría

Ejecutar primero 161--163. Si aparecen `NLSTRESS/BOUTPUT` para CBAR y todos
alcanzan factor 1.0, continuar con 164--166. Después lanzar 167--171 y, por
último, 172--176. No cambiar el nombre del job: carpeta, número de caso y SID
principal coinciden.

Devolver al menos `.f06` y `.pch` junto al DAT. Conservar también `.itr` y
`.sts` durante la auditoría inicial, porque pueden aclarar un fallo incremental;
se decidirá después si merecen versionarse.

En cada F06 comprobar:

1. ausencia de `USER FATAL MESSAGE` y `SYSTEM FATAL MESSAGE`;
2. terminación normal y factor de carga final 1.0;
3. número real de incrementos convergidos y posibles bisecciones;
4. `DISPLACEMENT VECTOR` y `FORCES OF SINGLE-POINT CONSTRAINT`;
5. tablas no lineales de CBAR y CROD, incluidos C--F en 161--166;
6. equilibrio entre cargas externas, reacciones y acciones de elemento;
7. resultados intermedios de `OUTCTRL=(SOLUTION,INTERM)`.

Un F06 que termina normalmente pero no contiene resultados CBAR no valida el
modelo: debe comprobarse que MSC reconoce `PBARN1`. Si un caso complejo falla,
no relanzar toda la campaña; devolver el F06 fallido para aislar si se trata de
convergencia, formato de una tarjeta o una combinación no admitida.
