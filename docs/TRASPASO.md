# Documento de traspaso

Fecha de referencia: **5 de agosto de 2026**
Repositorio: **MEF-solver**  
Rama: **`feat/ctria3`**  
Commit verificado: **`a374c52` (`feat(review_fixes): main workflow refactors and fixes`)**

Este documento permite continuar el proyecto desde otra instalación sin depender del historial de la conversación que originó los cambios. Antes de modificar el código, se debe recuperar el commit indicado, comprobar que el árbol de trabajo está limpio y ejecutar la validación descrita al final.

## 1. Objetivo del proyecto

Desarrollar, como Trabajo de Fin de Máster, un solver de elementos finitos estructural inspirado en Nastran, con estas prioridades:

1. Análisis estático lineal de pequeñas deformaciones.
2. Lectura directa de modelos Nastran BDF/DAT para ejecutar el mismo archivo tanto en Nastran como en este proyecto.
3. Elementos 1D (`CROD`, `CBAR`), 2D (`CTRIA3`, `CQUAD4`) y, en una fase posterior, sólidos 3D.
4. Validación numérica trazable frente a las entradas y salidas originales de Nastran guardadas en `tests/solver_tst/validation_cases`.
5. Mantener separadas la capacidad de leer una tarjeta, la formulación implementada y el grado de validación conseguido.

El trabajo actual está centrado en consolidar el flujo completo del solver y validar la formulación de `CQUAD4`, sin presentar todavía el programa como sustituto general de Nastran.

## 2. Decisiones que deben conservarse

Estas decisiones fueron expresamente fijadas durante la revisión y no deben revertirse sin una nueva decisión del autor:

- **PSHELL con campos MID vacíos:** se conserva la sustitución automática implementada en el modelo de datos. No se debe cambiar para imitar la semántica exacta de Nastran en esta fase.
- **Tolerancias de pruebas:** no se han endurecido ni deben modificarse todavía. Algunos casos contienen magnitudes muy próximas a cero y sus diferencias se aceptan provisionalmente; se analizarán con métricas específicas más adelante.
- **Entradas Nastran no implementadas:** no se deben rechazar de forma general. El objetivo es poder ejecutar el mismo `.dat` sin editarlo en Nastran y en MEF-solver mientras se amplía la cobertura. Las características aún no soportadas deben documentarse y tratarse de forma compatible cuando sea posible.
- **`assemble()` no idempotente:** no se ha cambiado. Cada instancia de `LinearStaticSolver` debe ensamblarse una sola vez antes de resolver.

### Qué significaría hacer `assemble()` idempotente

Una operación idempotente produce el mismo estado aunque se invoque varias veces. En el solver actual, una segunda llamada a `assemble()` vuelve a sumar las matrices y los vectores elementales, por lo que duplica `K` y `F`. Hacerla idempotente requeriría una de estas políticas:

- reinicializar siempre la matriz global y los vectores antes de ensamblar;
- detectar que ya se ensambló y no repetir la operación; o
- invalidar y reconstruir explícitamente el estado cuando cambie el modelo.

No se aplicó ninguna de ellas porque afecta al ciclo de vida del solver y no era necesaria para la validación actual. El contrato vigente está advertido en el código: crear la instancia, llamar una vez a `assemble()` y después a `solve()`.

## 3. Estado actual verificado

El commit de referencia está en `origin/feat/ctria3` y el árbol estaba limpio al crear este documento. Las comprobaciones realizadas fueron:

```text
python -m pytest -q
298 passed in 1.88s

python -m compileall -q src tests main.py
sin errores
```

También se comprobó previamente la ejecución completa de:

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_1\Ej1_v1.dat
```

La ejecución generó el informe `output/Ej1_v1.rpt` y recuperó correctamente desplazamientos, reacciones y resultados de `CROD`.

### Capacidades principales

| Área | Estado actual |
|---|---|
| Flujo BDF/DAT → modelo → ensamblaje → solución → resultados | Operativo para el subconjunto implementado |
| CROD | Rigidez, ensamblaje y recuperación de tensión validados |
| CBAR | Orientación, offsets, releases, acoplamiento `I12`, corrección de Timoshenko y resultados de extremos implementados |
| CTRIA3 | Membrana, flexión, cortante y recuperación de tensiones integrados; quedan discrepancias numéricas concretas por estudiar |
| CQUAD4 | Formulaciones Nastran y Reissner-Mindlin separadas; postproceso integrado; casos DMIG 51-54 automatizados |
| Case Control | Selección de un subcase, `LOAD`, `SPC`, `LOAD` combinado y `SPCADD` implementados |
| Sistemas de coordenadas | Transformaciones de `FORCE`/`MOMENT`, `GRID.CD`, `GRID.PS`, `THETA/MCID` y orientación de `CBAR` implementadas |
| Reacciones | Incluido el caso de modelos totalmente restringidos |
| Sólidos 3D | No implementados; fase futura |

### Situación específica de CQUAD4

La implementación activa usa deformación uniforme y modos bilineales con métrica de líneas medias, integración reducida estabilizada para cortante, `K6ROT` y una corrección de equilibrio para elementos alabeados. La formulación Reissner-Mindlin clásica vive en una clase independiente.

La comparación automatizada obtiene **0,152382 %** en los casos 51 y 52, **0,372560 %** en el caso 53 y **2,150270 %** en la matriz completa del caso 54 alabeado. Los casos planos exigen menos del 1 % y el alabeado menos del 3 %. No se afirma igualdad término a término con el código propietario de Nastran.

La matriz completa del caso 51 presenta una discrepancia aproximada del **50,26 %** porque el archivo representa una configuración de membrana pura para Nastran, mientras que la decisión vigente de sustituir automáticamente los MID vacíos de `PSHELL` activa contribuciones adicionales en MEF-solver. Por esa razón, mientras se mantenga esa decisión, la comparación correcta es el bloque de membrana activo y no la matriz completa.

### Limitaciones conocidas que no bloquean la suite

- Los casos CTRIA3 36, 38, 44, 45, 48 y 49 conservan diferencias en magnitudes de flexión o reacciones que requieren un estudio numérico dedicado.
- La reconstrucción CQUAD4 está contrastada con cuatro DMIG, pero el caso alabeado aún conserva un error del 2,15 % y necesita más geometrías independientes.
- El parser resuelve un único subcase seleccionado por ejecución; todavía no devuelve ni resuelve todos los subcases de un modelo en una sola llamada.
- Quedan campos de propiedades y elementos shell pendientes de tratamiento físico completo, entre ellos `12I/T³`, `TS/T`, `Z1/Z2`, `MID4`, `NSM`, `ZOFFS`, `TFLAG` y espesores nodales `Ti`.
- El informe de salida puede mejorar la nomenclatura y separación de resultados por tipo de elemento.
- `cquad4_element.py` conserva un bloque antiguo entre comillas triples (`_build_local_stiffness_old`) que no se ejecuta y puede retirarse en una limpieza posterior.

## 4. Cambios realizados por archivo

El commit `a374c52` modifica 27 archivos, con 1797 inserciones y 189 eliminaciones.

### Ejecución, dependencias y documentación

| Archivo | Cambio relevante |
|---|---|
| `README.md` | Corrige el alcance a análisis estático lineal y actualiza el uso del programa. |
| `main.py` | Añade `--subcase`, corrige conteos y presentación de elementos, y propaga excepciones para devolver un código de salida no nulo. |
| `requirements.txt` | Declara NumPy, SciPy, pandas, pyNastran y pytest. |
| `pytest.ini` | Limita el descubrimiento a `tests` y evita recoger pruebas antiguas dentro de `src`. |
| `docs/README.md` | Índice y resumen funcional de la documentación. |
| `docs/arquitectura.md` | Arquitectura, modelo de datos y flujo del análisis. |
| `docs/formulacion.md` | Ecuaciones, convenciones y formulaciones implementadas. |
| `docs/uso.md` | Instalación, ejecución, entradas y salidas. |
| `docs/validacion.md` | Cobertura y estrategia de comparación con Nastran. |
| `docs/cquad4.md` | Estado técnico y validación específica del cuadrilátero. |
| `docs/roadmap.md` | Limitaciones y evolución propuesta. |

### Preproceso y modelo de datos

| Archivo | Cambio relevante |
|---|---|
| `src/preprocessing/data_model.py` | Añade transformaciones nodales, metadatos de subcase y validaciones físicas/de referencias. Mantiene la sustitución automática de MID en `PSHELL`. |
| `src/preprocessing/nastran_reader_v2.py` | Selección de subcase, cargas combinadas, `SPCADD`, coordenadas, `GRID.PS`, orientación shell y ejes/offsets de `CBAR`; sustituye la terminación del proceso por excepciones. |

### Solver y elementos

| Archivo | Cambio relevante |
|---|---|
| `src/solver/linear_static_solver.py` | Ensambla en los sistemas `GRID.CD`, transforma cargas y desplazamientos, integra CQUAD4 en el postproceso, calcula reacciones totalmente restringidas y mejora diagnósticos de singularidad/conflictos. |
| `src/solver/transformations.py` | Valida vectores de orientación nulos o no finitos. |
| `src/solver/elements/bar_element.py` | Implementa ejes resueltos, offsets rígidos, releases, `I12`, corrección de Timoshenko y recuperación ampliada de tensiones. |
| `src/solver/elements/ctria3_element.py` | Mejora validaciones geométricas, simetría, uso de `G`, rotación de resultados y compatibilidad de nombres `Z1/Z2`. |
| `src/solver/elements/cquad4_element.py` | Base compartida de geometría, constitutivas, `K6ROT` y recuperación de resultados. |
| `src/solver/elements/cquad4_nastran_element.py` | Rigidez Nastran reconstruida, modos bilineales y corrección de alabeo. |
| `src/solver/elements/cquad4_reissner_mindlin_element.py` | Formulación clásica Reissner-Mindlin con integración completa. |

### Pruebas

| Archivo | Cambio relevante |
|---|---|
| `tests/preprocessing_tst/test_case_control_and_coordinates.py` | Cubre subcases, combinaciones `LOAD`, `SPCADD`, `CID/CD` y `GRID.PS`. |
| `tests/solver_tst/test_solver_regressions.py` | Cubre reacciones con restricción total, conflictos SPC y postproceso global CQUAD4. |
| `tests/solver_tst/elements/test_cquad4_nastran_matrix.py` | Compara automáticamente las DMIG Nastran de los casos 51 a 54. |
| `tests/solver_tst/elements/CQUAD4Element_test.py` | Comprueba seis modos rígidos y estabilización *drilling*. |
| `tests/solver_tst/elements/CTRIA3ElementSolver_test.py` | Repara importaciones y compatibilidad de la prueba. |
| `tests/solver_tst/elements/run_cquad4_validation.py` | Mantiene compatibilidad y delega en la comparación portable de los casos 51-54. |
| `tests/solver_tst/elements/test_E2E_CTRIA_K_solution.py` | Corrige rutas/literales copiados y evita exportaciones de depuración. |
| `tests/solver_tst/test_validation_end_to_end.py` | Corrige un literal de referencia copiado. |
| `tests/solver_tst/test_validation_end_to_end_ctria3.py` | Corrige literales de referencia copiados sin cambiar tolerancias. |

Los literales corregidos se contrastaron con los resultados originales almacenados junto a cada caso: caso 36 (`Rz`), caso 43 (`T3`) y la ruta del caso 18.

## 5. Archivos que requieren atención al mover el proyecto

El archivo `tests/solver_tst/draft_ctria_validationcases.py` coincide con la regla `draft_*` de `.gitignore`. Se corrigió localmente su sintaxis envolviendo su contenido como texto para que no rompiera `compileall`, pero **no forma parte del commit ni se descargará en otro ordenador**.

No es necesario para ejecutar pytest. Si se quiere conservar ese borrador, hay que copiarlo expresamente o renombrarlo y versionarlo en un commit separado. No debe asumirse que su presencia forma parte del estado reproducible descrito aquí.

Los informes generados bajo `output`, como `output/Ej1_v1.rpt`, tampoco son fuente del programa y pueden regenerarse.

## 6. Próximos pasos recomendados

Orden propuesto para continuar:

1. **Auditoría numérica de CTRIA3:** estudiar los casos 36, 38, 44, 45, 48 y 49 sin cambiar aún sus tolerancias. Añadir métricas diagnósticas separadas para error absoluto, relativo y magnitudes próximas a cero.
2. **Consolidar CQUAD4:** añadir más casos alabeados y de cortante, pruebas de parche multielemento y convergencia, y estudiar la diferencia residual del caso 54.
3. **Ampliar validación de CBAR:** crear casos combinados con `I12`, `K1/K2`, offsets y releases, comparados directamente con salidas originales de Nastran.
4. **Completar parámetros shell:** incorporar progresivamente el efecto físico de los campos pendientes de `PSHELL` y `CQUAD4`, manteniendo la entrada `.dat` compatible y sin convertir lo no implementado en un rechazo general.
5. **Arquitectura multicaso:** separar modelo base, configuración de cada subcase y resultados para resolver varios subcases en una misma ejecución.
6. **Resultados e informes:** usar nombres específicos por elemento y añadir formatos estructurados que faciliten comparaciones automáticas con F06/OP2.
7. **Rendimiento:** perfilar antes de optimizar; después considerar ensamblaje disperso vectorizado y caché de transformaciones/formulaciones elementales.
8. **Sólidos 3D:** comenzar sólo cuando las formulaciones 1D/2D y su estrategia de validación estén consolidadas.

La idempotencia de `assemble()` queda como decisión futura de arquitectura, no como corrección pendiente inmediata.

## 7. Recuperación en el otro ordenador

Con Git disponible, desde el clon del mismo repositorio:

```powershell
git fetch origin
git switch feat/ctria3
git pull --ff-only origin feat/ctria3
git rev-parse --short HEAD
git status --short
```

El primer comando de comprobación debe mostrar `a374c52` o un descendiente que incluya ese commit. `git status --short` debe quedar vacío antes de empezar a modificar archivos.

Crear el entorno y validar:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m compileall -q src tests main.py
```

Si el commit que añade este propio documento todavía no está en remoto al cambiar de ordenador, deberá confirmarse y publicarse primero; `a374c52` contiene las correcciones del solver, pero no este documento posterior.

## 8. Mensaje sugerido para reanudar la conversación

En una nueva tarea de Codex abierta sobre el mismo proyecto:

> Lee completamente `docs/TRASPASO.md`, comprueba la rama, el commit y `git status`, y ejecuta la suite antes de modificar nada. Continúa con el primer próximo paso pendiente respetando las decisiones de la sección 2. Usa como referencia principal los inputs y resultados originales de Nastran en `tests/solver_tst/validation_cases` y no cambies las tolerancias existentes salvo que yo lo autorice.

Después conviene indicar explícitamente qué próximo paso se quiere abordar, ya que este documento transfiere el contexto técnico pero no inicia trabajo adicional por sí solo.

## 9. Actualización CQUAD4 de 13 de agosto de 2026

Esta sección prevalece sobre las cifras y nombres CQUAD4 anteriores de este
documento. El estado reproducible actual está descrito en
`docs/cquad4_tesis.md`, `docs/cquad4.md` y `docs/validacion.md`.

- La suite completa queda en `333 passed, 4 xfailed`. Los `xfail` son 62-64
  (alabeo) y 78 (`SNORM` predeterminado no implementado), no fallos de CQUAD4
  coplanar con normal propia.
- Se incorporaron como literales de test los F06/OP2 84-88. Cubren paneles
  planos multipropiedad, mezclas CQUAD4/CTRIA3/CBAR/CROD, `K6ROT=50/100`, FORCE,
  MOMENT y desplazamientos impuestos. El test nunca lee F06/OP2 en ejecución.
- Las DMIG 83 y 89, junto con 81, demostraron que la flexibilidad aumentada de
  un paralelogramo oblicuo acopla las dos muestras de cada dirección. La
  métrica positiva/acotada está en
  `CQUAD4NastranElementSolver._build_augmented_shear_stiffness`; sus constantes
  son una identificación DMIG, no una constante publicada por MSC.
- Los errores completos de 81/83/89 son `1.03e-7`, `3.14e-6` y `8.89e-7` y se
  comprueban también término a término. No volver al polinomio de cuarto orden:
  cambiaba de signo para la oblicuidad del caso 60 y rompía una regresión E2E.
- `MID3` vacío significa `SHRFLX=False` y ninguna flexibilidad física `Z_s`.
  Aun así, NASTRAN-95 reutiliza internamente `MID2` al formar la corrección
  residual/torsional; por ello no debe documentarse simplemente como “MID3 no
  hereda MID2” sin esta precisión.
- `SNORM` crea normales comunes en GRID compartidos dentro de una tolerancia
  angular (20 grados por defecto en MSC). El solver usa normales independientes
  y para compararlo de igual a igual se emplea `PARAM,SNORM,0.0`.
- La discrepancia coplanar abierta está en la recuperación centroidal `QX/QY`,
  no en la rigidez. En el caso 84, `U`, reacciones, `N/M` y tensiones coinciden
  estrechamente, pero un `QX` pequeño difiere `0.00977` (`72 %` relativo).
  Mantener tolerancias exclusivas para `QX/QY`; no alterar la rigidez cerrada
  para ajustar esta salida.
- Los casos 90 (DMIG) y 91 (F06) están preparados para una comprobación
  predictiva de alta oblicuidad. El siguiente paso es ejecutarlos en Nastran e
  incorporar resultados sin recalibrar primero la formulación.
