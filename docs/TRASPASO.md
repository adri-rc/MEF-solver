# Traspaso técnico del proyecto

Actualizado el **21 de agosto de 2026**. Este documento es el
punto de entrada para continuar MEF-solver en otra cuenta, sesión u ordenador.
El punto de referencia histórico es la rama `feat/ctria3`, commit `0cc154c`.
La fase actual se desarrolla en `feat/non-linear-CBAR`. El commit `df88ee4`
contiene la campaña 128--153 y la primera implementación CBAR no lineal. El
cierre `LC/LS` basado en 154--160 permanece como cambio local. Los 33 casos
alcanzan carga completa en MSC y coinciden con el solver, incluido 142.

## 1. Incorporación rápida para una IA

Antes de modificar nada:

```powershell
git status --short --branch
git log -12 --date=short --oneline --decorate
git lfs pull
.\.venv\Scripts\python.exe -m pytest -q
```

Después, leer en este orden:

1. [`AGENTS.md`](../AGENTS.md): invariantes y forma de trabajo obligatoria.
2. Este documento completo: estado, decisiones e historial reciente.
3. [`arquitectura.md`](arquitectura.md) y [`uso.md`](uso.md): flujo y ejecución.
4. [`roadmap.md`](roadmap.md): prioridades abiertas.
5. Para trabajar en CQUAD4, [`cquad4.md`](cquad4.md), después
   [`cquad4_tesis.md`](cquad4_tesis.md) y [`validacion.md`](validacion.md).
6. Los tests próximos al código que se vaya a cambiar.

Si una afirmación documental contradice el código o los tests, comprobarla y
corregir el `.md` en el mismo cambio. No asumir que una idea de documentos
antiguos o del material de tesis pertenece a la ruta activa.

## 2. Estado verificable

Comando de referencia:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Resultado observado el 21 de agosto de 2026:

```text
587 passed, 6 xfailed, 2 warnings
```

Las dos advertencias se deben a que pytest no puede escribir `.pytest_cache`.
Los seis `xfail` son estrictos y deliberados:

- 54, 59, 62, 63 y 64: uno o varios CQUAD4 alabeados;
- 78: ensamblaje facetado con suavizado `SNORM` predeterminado.

No hay fallos activos de CQUAD4 plano. Los casos facetados 82/88 usan normales
independientes (`SNORM,0`) y pasan.

### Estado de Git que debe conocerse

- Remoto: `https://github.com/adri-rc/MEF-solver.git`.
- La rama de trabajo está muy por delante de `main`; continuar desde
  `feat/ctria3` o desde una rama nueva creada a partir de ella, no desde
  `main`, salvo que antes se integre el historial.
- Los `.dat`, `.bdf`, `.f06`, PDF e imágenes están administrados por Git LFS.
  Un clon necesita Git LFS y `git lfs pull` para materializar 231 referencias
  observadas en este estado.
- Al preparar este traspaso, `docs/TFM/` era una carpeta local no versionada
  con fuentes LaTeX y dos artículos PDF. No forma parte del clon remoto. Debe
  decidirse expresamente si se incorpora, se transfiere por separado o se
  mantiene privada; una IA no debe editarla ni borrarla por iniciativa propia.
- Las modificaciones documentales de este traspaso tampoco llegarán a otra
  máquina hasta que el usuario las confirme, haga commit y push, o copie el
  directorio de trabajo completo.

## 3. Alcance funcional

- Estática lineal, pequeñas deformaciones y `MAT1` isótropo.
- Elementos: CROD, CBAR, CTRIA3 y CQUAD4 coplanar.
- Propiedades: PROD, PBAR y subconjunto PSHELL.
- Cargas: FORCE, MOMENT y combinaciones LOAD.
- Restricciones: SPC, SPC1, SPCADD, campo GRID.PS y valores impuestos.
- Selección explícita de SUBCASE cuando hay más de uno.
- Posición GRID.CP, grados de libertad GRID.CD y cargas en sistemas locales se
  transforman a la base adecuada.
- `PARAM,K6ROT` se aplica a CTRIA3 y CQUAD4.
- No hay `AUTOSPC` general: los mecanismos reales deben fallar.
- El informe `.rpt` distingue análisis lineal/no lineal y contiene resumen del
  modelo, configuración, desplazamientos, reacciones y tablas físicas
  específicas para CROD, CBAR, CTRIA3 y CQUAD4. En SOL 400 incluye además los
  historiales de incrementos convergidos y correcciones de Newton al final. La
  correspondencia con F06 está en `docs/result_report.md`.
- El preproceso distingue `SOL 101` y estática no lineal `SOL 400`; para esta
  última lee `ANALYSIS=NLSTAT`, `NLSTEP`/`NLPARM`, `LGDISP` y `FOLLOWK`.
- La campaña CROD comprende 36 referencias, casos 92--127, ejecutadas con MSC
  Nastran 2026.1/NLPERF. Treinta y cinco alcanzan carga completa y el 120
  documenta un punto límite bajo control de fuerza. Los modelos 123--127
  validan estructuras de 10--30 barras, propiedades mixtas y SPC no nulos.
- `NonlinearStaticSolver` resuelve SOL 400 para MAT1 con CROD o CBAR, cargas
  nodales conservativas, `FIXED`/`ADAPT`, Newton completo y recortes.
- `CBARNonlinear` implementa una formulación corrotacional con deformación
  axial logarítmica, rotaciones finitas, `I12`, `K1/K2`, offsets y releases.
- La campaña CBAR comprende 33 decks, casos 128--160, todos con soluciones MSC
  y referencias E2E explícitas. Los 33 pasan.
- Los casos 154--160 identifican que `LC` es Euler--Bernoulli e ignora
  `K1/K2`, mientras `LS` usa interpolación lineal-cortante de un punto y un
  factor unitario si `K1/K2=0`. El antiguo `xfail` 142 queda cerrado.

No existe todavía solución dinámica, no linealidad material, contacto, cargas
seguidoras, cargas distribuidas, MPC/RBE, sólidos, CQUAD4
alabeado ni suavizado `SNORM`.

## 4. Arquitectura y flujo activo

- `main.py`: CLI y orquestación del trabajo.
- `src/preprocessing/data_model.py`: entidades independientes del formato.
- `src/preprocessing/nastran_reader_v2.py`: pyNastran → `FEMModel`, selección
  de subcaso, configuración lineal/no lineal y reducción de cargas/restricciones.
- `src/preprocessing/data_model.py`: añade `AnalysisConfiguration` y
  `NonlinearStepControl`, independientes de pyNastran.
- `src/solver/linear_static_solver.py`: ensamblaje LIL, partición, `spsolve`,
  reacciones y recuperación numérica.
- `src/solver/nonlinear_analysis.py`: incrementos, Newton-Raphson, residuo,
  tangente global, recortes e historial CROD/CBAR.
- `src/solver/elements/crod_non_linear.py`: fuerza interna, energía y tangente
  material/geométrica consistente del CROD.
- `src/solver/elements/cbar_non_linear.py`: cinemática corrotacional, matrices
  `PBARN1 LC/LS`, acciones, tangente numérica y recuperación del CBAR.
- `src/solver/elements/`: formulaciones de CROD, CBAR, CTRIA3 y CQUAD4.
- `src/postprocessing/report_writer.py`: presentación del informe; ya no vive
  dentro de `main.py`.
- `tests/solver_tst/validation_cases/`: DAT, F06, PCH y CSV de referencia.

El flujo es:

```text
DAT/BDF -> pyNastran -> FEMModel -> K/F -> contorno -> spsolve
        -> recuperación por elemento -> report_writer -> RPT
```

Ese diagrama describe la ruta lineal. Para `SOL 400`, `main.py` selecciona el
módulo hermano `nonlinear_analysis.py`, propietario de incrementos, residuo,
tangente, recortes y recuperación incremental.
pyNastran 1.4.1 rechaza sintácticamente la tarjeta `NLSTEP`, pero conserva sus
líneas; el lector interpreta de forma estricta `GENERAL`, `FIXED`/`ADAPT` y
`MECH`.

La matriz elemental se expresa primero en coordenadas básicas. En el
ensamblaje se transforma después al sistema `CD` de cada GRID. La recuperación
convierte los desplazamientos de nuevo a coordenadas básicas antes de aplicar
el sistema local del elemento. Véase [`arquitectura.md`](arquitectura.md) antes
de cambiar transformaciones o el orden de los GDL.

`run_analysis_job()` vuelve a lanzar las excepciones después de imprimir el
contexto, de modo que la CLI termina con error. En `--all`, el primer fallo
detiene el lote; la política de continuar y acumular fallos sigue abierta.

## 5. Decisiones cerradas de CQUAD4

La clase activa es `CQUAD4NastranElementSolver`. Su dominio son cuatro nodos
coplanares, ordenados consecutivamente y con jacobiano positivo. La tolerancia
de coplanaridad es `max(1e-10, 1e-8*L)`. La ruta antigua de
proyección/corrección de alabeo se eliminó.

La formulación contiene:

1. sistema local por diagonales y bisectriz de diagonales unitarias;
2. funciones de forma bilineales;
3. deformación uniforme y modos bilineales normalizados para membrana/flexión;
4. métrica geométrica de líneas medias;
5. cortante por `edge strings` de MacNeal;
6. flexibilidad residual más MID3 en serie;
7. corrección de aspecto con `epsilon=0.025`;
8. métrica torsional oblicua racional;
9. `K6ROT` como penalización `theta_z - omega`, igual que CTRIA3 y extendida a
   cuatro vértices;
10. recuperación central de N, M, Q y tensiones de fibra.

Semántica PSHELL que no debe revertirse:

- MID2 vacío desactiva flexión; no hereda MID1.
- MID3 vacío desactiva la flexibilidad física de cortante; no hereda MID2.
- MID3 vacío conserva la flexibilidad residual de flexión del QUAD4.
- `12I/T^3` y `TS/T` se leen y aplican.

### Ajuste de los casos 90/91

La ley oblicua anterior era demasiado blanda para el caso 90. Se mantuvo el
operador por aristas y se modificó un único escalar torsional:

```text
m(x) = 1 + 0.4100259495463487*x / (1 - 0.9156556465380681*x)
x = cos(midline_angle)^2
```

Controles principales:

- caso 90: DMIG completa menor que `1e-6` de Frobenius y bloque
  placa/cortante menor que `4e-6`;
- caso 91: errores vectoriales `U=2.07e-5`, `R=6.30e-5`, y N, M y fibras por
  debajo de 0.1 % por componente;
- Qx/Qy llega a 47.1 % en una componente pequeña, 0.0216 absoluto.

La discrepancia Q es de recuperación centroidal, no de rigidez. No ajustar
`_build_augmented_shear_stiffness` para corregirla sin nueva evidencia DMIG.

### Casos 57–59

- 57: `cond(K)=2.90e25`; errores vectoriales `U=1.72e-3`, `R=4.07e-3` y
  resultados `9.28e-4`.
- 58: `cond(K)=8.52e26`; traslaciones próximas, pero pequeñas rotaciones,
  momentos secundarios y Q son sensibles; resultados `2.94e-2`.
- 59: alabeado, fuera de alcance y `xfail` estricto.

No usar 57/58 para recalibrar la rigidez: el MID3 casi nulo hace que sean
controles extremadamente mal condicionados.

## 6. Validaciones que orientan el trabajo

- `tests/solver_tst/elements/CQUAD4Element_test.py`: geometría, simetría,
  modos rígidos, hourglass, K6ROT, ley oblicua y rechazo no coplanar.
- `tests/solver_tst/elements/test_cquad4_pshell_material_semantics.py`: MID2,
  MID3, `12I/T^3` y `TS/T`.
- `tests/solver_tst/elements/test_E2E_CQUAD_Nastran_K_solution.py`: DMIG,
  Frobenius, término a término y bloques.
- `tests/solver_tst/test_validation_end_to_end_cquad4.py`: referencias E2E
  escritas directamente; no lee F06 durante el test.
- El caso 80 de ese último módulo mezcla CBAR, CTRIA3 y CQUAD4 y comprueba
  desplazamientos, cargas, reacciones y resultados de las tres familias.
- `tests/postprocessing_tst/test_report_writer.py`: contrato del informe RPT,
  incluida configuración, convergencia iterativa y estado final no lineal.
- `tests/preprocessing_tst/test_nonlinear_analysis_control.py`: normalización
  SOL 400, `NLSTEP`/`NLPARM`, enrutado al solver no lineal y lectura de los 36
  decks CROD; también comprueba estructura, campos avanzados y cabeceros de
  los 33 decks CBAR.
- `tests/solver_tst/elements/test_crod_non_linear.py`: cinemática, fuerzas,
  torsión, objetividad y derivada numérica de la tangente.
- `tests/solver_tst/test_nonlinear_intermediate_validation.py`: estados de 112
  y flexibilidad tangente de 114/115.
- `tests/solver_tst/test_nonlinear_validation_end_to_end.py`: diccionarios
  explícitos de los 35 casos convergidos; 123--127 documentan además metadatos,
  tolerancias y descripciones legibles por caso.
- `tests/solver_tst/validation_cases/NONLINEAR_CROD_CASES.md`: orden de
  ejecución y resultados de MSC que deben conservarse.
- `docs/nonlinear_cbar.md`: formulación corrotacional implementada, validación
  y límites conocidos.
- `tests/solver_tst/validation_cases/NONLINEAR_CBAR_CASES.md`: matriz y orden
  de ejecución de los casos 128--153 y auditoría de la primera tanda.

Los parámetros `debug_export` y `output_k_file` controlan la escritura de la
matriz calculada en CSV. Los ficheros `draft_*` están ignorados por Git.

## 7. Histórico reciente y por qué importa

| Commit | Fecha | Decisión o resultado |
|---|---|---|
| `0cc154c` | 2026-08-16 | Caso mixto 80 revisado; referencia de ensamblaje CBAR/CTRIA3/CQUAD4. |
| `fd96a8d` | 2026-08-14 | Escritura RPT separada en `postprocessing/report_writer.py` y cubierta por tests. |
| `624549f` | 2026-08-14 | Nombres físicos de resultados de lámina corregidos y homogeneizados. |
| `0bc6cd9` | 2026-08-14 | Cierre de la implementación y batería CQUAD4 coplanar. |
| `13a9ca8`–`dd2c362` | 2026-08-13/14 | Ajustes de flexibilidad transversal, ley oblicua y tolerancias con evidencia E2E/DMIG. |
| `eaf23d8` | 2026-08-05 | Separación entre formulación CQUAD4 tipo Nastran y Reissner-Mindlin. |
| `a374c52`–`5d7e65e` | 2026-08-04 | Revisión general del flujo y anterior punto de traspaso. |

Para el detalle exacto, usar `git show <commit>`; la tabla expresa intención y
no sustituye al diff. El nombre histórico de la rama (`feat/ctria3`) ya no
describe todo su contenido: incluye también el cierre de CQUAD4 y el nuevo
postproceso.

## 8. Forma de continuar

Prioridad recomendada:

1. investigar la recuperación de tensión de fibra CBAR bajo rotación 3D;
2. integración continua y entorno bloqueado/reproducible;
3. estudio de convergencia, ensamblaje COO y benchmark de escalabilidad;
4. política de errores del modo por lotes.

No es necesario crear más casos CQUAD4 de rigidez plana para cerrar la
implementación actual. Cualquier ampliación debe comenzar con un deck mínimo,
una referencia Nastran u otra fuente auditable y una prueba que falle antes de
modificar el código.

Al terminar cada bloque de trabajo, dejar constancia de:

- problema y alcance decidido;
- archivos y decisiones técnicas;
- pruebas focalizadas y suite completa;
- limitaciones o incertidumbres restantes;
- actualización de este documento si cambia el estado global.

## 9. Mensaje recomendado para una cuenta nueva

```text
Quiero continuar el desarrollo de este repositorio MEF-solver. Antes de hacer
cambios, lee por completo AGENTS.md y docs/TRASPASO.md y sigue el orden de
lectura que allí se indica. Revisa también git status, la rama activa, los 12
commits más recientes y el estado de Git LFS. Ejecuta la suite completa si el
entorno está disponible. Después, sin modificar todavía ningún archivo,
explícame: (1) alcance y arquitectura actuales, (2) decisiones técnicas que no
deben revertirse, (3) validación y límites conocidos, (4) histórico reciente,
(5) cambios locales o material no versionado y (6) próximos trabajos de la
hoja de ruta. Distingue hechos verificados de inferencias y avísame si falta
algún archivo LFS o si no estás en la rama que contiene 0cc154c. Espera mi
siguiente objetivo antes de implementar nada.
```
