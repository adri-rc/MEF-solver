# Limitaciones y hoja de ruta

## Diagnóstico

El proyecto tiene una arquitectura válida para un TFM y una base de validación inusualmente rica. El riesgo principal no es la falta de funcionalidades, sino confundir compatibilidad de lectura con equivalencia física. La prioridad debe ser cerrar, automatizar y documentar el alcance lineal actual antes de añadir sólidos 3D.

## Limitaciones funcionales

### Entrada y casos de carga

- no se seleccionan SUBCASE, LOAD o SPC por SID;
- no hay `LOAD`, `SPCADD`, cargas distribuidas, gravedad ni térmicas;
- FORCE/MOMENT solo admiten `CID=0`;
- no se aplican SPC permanentes de `GRID.PS`;
- MPC, RBE y elementos rígidos se rechazan;
- la semántica de campos PSHELL vacíos no coincide con Nastran.

### Elementos

- CROD: formulación estable, pero requiere restringir todos los GDL inactivos.
- CBAR: sin `G0`, offsets, liberaciones, `I12` ni cortante Timoshenko.
- CTRIA3: formulación de cortante calibrada de forma limitada; dependencia importante de K6ROT/MID3.
- CQUAD4: postproceso desactivado, *drilling* ausente, modos de energía cero y discrepancia de rigidez.
- no existen elementos sólidos.

### Núcleo numérico

- ensamblaje mediante bucles escala mal;
- no se comprueba condicionamiento ni se diagnostican GDL causantes de singularidad;
- no hay AUTOSPC controlado;
- no existe escalado de matrices ni selección de solver;
- una llamada repetida a `assemble()` acumula otra vez K y F;
- no se cachean los solvers de elemento entre ensamblaje y recuperación.

### Resultados y calidad de software

- informe orientado a CROD y con contador CROD incorrecto;
- errores de análisis se capturan sin código de salida fallido;
- directorio de postproceso vacío;
- suite global rota por imports y ruta de datos;
- faltan metadatos de dependencias y configuración pytest;
- existen código comentado extenso, borradores ignorados y texto con problemas de codificación;
- no hay CI ni una especificación versionada del formato de resultados.

## Hoja de ruta priorizada

### P0 — Línea base reproducible

Objetivo: que cualquier revisor pueda clonar, instalar y ejecutar exactamente la misma suite.

1. Crear `pyproject.toml` con versiones compatibles y comandos de test.
2. Configurar pytest para recoger solo `tests/`.
3. Corregir imports obsoletos y ubicación del caso 18.
4. Eliminar rutas absolutas de scripts de validación.
5. Añadir CI y publicar el resumen de pruebas.
6. Corregir README raíz, que actualmente describe por error un solver “implícito no lineal”.

### P1 — Cerrar CQUAD4

Objetivo: elemento estable y validado de extremo a extremo.

1. Separar y probar bloques `K_mem/K_bend/K_shear/K_drill`.
2. Resolver convenciones de giro y discrepancias de casos 51/52.
3. Implementar tratamiento *drilling* y control de modos espurios.
4. Restaurar y validar recuperación de resultantes/tensiones.
5. Integrarlo en `reactions_and_stresses()` y el informe.
6. Añadir convergencia de malla y *patch tests* multielemento.

El detalle está en [cquad4.md](cquad4.md).

### P2 — Semántica Nastran y casos múltiples

Objetivo: definir con precisión qué significa “compatible con Nastran”.

1. Representar MID1/MID2/MID3 opcionales sin materiales ficticios.
2. Seleccionar un SUBCASE y solo sus LOAD/SPC asociados.
3. Implementar combinaciones `LOAD` y `SPCADD`.
4. Transformar cargas desde CID y aplicar `GRID.PS`.
5. Validar referencias de propiedad/material antes del ensamblaje.
6. Documentar tarjetas y campos soportados a nivel de campo, no solo de tarjeta.

### P3 — Robustez y rendimiento

Objetivo: pasar de ejemplos docentes a mallas medianas con diagnósticos útiles.

1. Ensamblar en formato COO mediante arrays de filas, columnas y valores.
2. Hacer `assemble()` idempotente o impedir llamadas duplicadas.
3. Cachear instancias/estado de elemento.
4. Diagnosticar filas nulas, componentes desconectadas y modos de mecanismo.
5. Incorporar chequeos de calidad geométrica y condicionamiento.
6. Añadir balance global de fuerzas y momentos al resultado.

### P4 — Completar CBAR y postproceso

Objetivo: una capacidad 1D coherente y presentable.

1. Implementar `G0`, `OFFT`, `WA/WB`, `PA/PB` e `I12`, o rechazarlos explícitamente si aparecen.
2. Decidir Euler-Bernoulli frente a Timoshenko y tratar `K1/K2` en consecuencia.
3. Separar recuperación de resultados en `src/postprocessing/`.
4. Crear escritores específicos por elemento y un formato estructurado (JSON/HDF5 además de texto).
5. Corregir resumen, contadores y propagación de errores de `main.py`.

### P5 — Elementos sólidos 3D

Objetivo: incorporar sólidos sin arrastrar ambigüedades del núcleo.

Orden recomendado:

1. `CTETRA4`: formulación de deformación constante, útil como primera implementación y para verificar el ensamblaje 3D.
2. `CHEXA8`: isoparamétrico trilineal con integración `2x2x2`, pruebas de *locking* y distorsión.
3. después, elementos de orden superior si el alcance del TFM lo permite.

Antes de esta fase deben existir:

- ley constitutiva 3D bien probada;
- utilidades generales de Jacobiano e integración multidimensional;
- *patch tests* 3D;
- casos analíticos de tracción, cortante y flexión;
- comparación independiente y estudios de convergencia.

## Entregables recomendados para el TFM

- especificación de alcance y subconjunto Nastran;
- derivación de cada formulación y convenciones de signos;
- arquitectura del software y decisiones de diseño;
- matriz de verificación requisito-prueba-caso;
- resultados de *patch tests*, comparación Nastran y convergencia;
- estudio de sensibilidad a distorsión, espesor y K6ROT;
- limitaciones conocidas, sin ocultar diferencias respecto a Nastran;
- caso demostrador reproducible desde entrada hasta informe.

La aportación académica será más sólida si cada afirmación de compatibilidad se vincula a un caso y una tolerancia, y si las diferencias de formulación se explican en lugar de ajustarse únicamente por coincidencia numérica.
