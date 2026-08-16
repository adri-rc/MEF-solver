# Guía de trabajo para asistentes de IA

Estas instrucciones se aplican a todo el repositorio. El objetivo es conservar
la trazabilidad numérica de un solver académico comparado con MSC Nastran.

## Incorporación al proyecto

Antes de proponer o modificar código:

1. ejecutar `git status --short --branch` y no sobrescribir cambios locales;
2. leer `README.md`, `docs/TRASPASO.md` y `docs/roadmap.md`;
3. consultar el documento técnico específico de la zona que se vaya a tocar;
4. revisar los tests y casos de validación relacionados;
5. ejecutar como línea base
   `.\.venv\Scripts\python.exe -m pytest -q` en Windows, si existe el entorno.

El estado verificable, el historial reciente y el orden de lectura están en
`docs/TRASPASO.md`. El código y los tests tienen prioridad si un documento se
ha quedado desactualizado; en ese caso, actualizar también la documentación.

## Invariantes que deben preservarse

- El alcance actual es estática lineal de pequeñas deformaciones con `MAT1`
  isótropo y seis GDL por nodo.
- La ruta activa admite CROD, CBAR, CTRIA3 y CQUAD4 coplanar. No presentar el
  proyecto como un sustituto general de Nastran.
- No existe `AUTOSPC` general: un mecanismo estructural real debe producir un
  fallo, no estabilizarse silenciosamente.
- Las coordenadas `GRID.CP`, los GDL en `GRID.CD` y los vectores de carga en su
  `CID` deben transformarse en las etapas documentadas en
  `docs/arquitectura.md`.
- Para PSHELL, `MID2` vacío desactiva flexión y `MID3` vacío desactiva la
  flexibilidad física de cortante; no heredan automáticamente otros MID.
- CQUAD4 rechaza geometrías alabeadas. Los `xfail` estrictos de alabeo y
  `SNORM` son límites explícitos, no fallos pendientes del elemento coplanar.
- No recalibrar la rigidez CQUAD4 para ajustar únicamente `Qx/Qy` centroidales
  sin una nueva referencia DMIG: la discrepancia conocida es de recuperación.
- La presentación de resultados pertenece a
  `src/postprocessing/report_writer.py`; el solver conserva la recuperación
  numérica y `main.py` sólo orquesta el flujo.

## Forma de trabajo

- Hacer cambios pequeños y apoyarlos primero con el test focalizado. Antes de
  entregar, ejecutar la colección completa cuando el entorno lo permita.
- Toda nueva capacidad física necesita un caso mínimo, una referencia externa
  auditable y una prueba E2E que falle antes de implementar el cambio.
- No relajar tolerancias ni convertir un fallo en `xfail` sin documentar la
  causa numérica y la evidencia de referencia.
- Mantener nombres de código e identificadores en inglés. La documentación de
  usuario y el resumen de entrega se mantienen en español.
- Los modelos `.dat/.bdf`, resultados `.f06` y otros binarios usan Git LFS.
  Comprobar que no sean punteros de texto antes de diagnosticar una lectura.
- `output/`, cachés y ficheros `draft_*` son artefactos locales; no son fuente
  de verdad ni deben incorporarse al repositorio.
- `docs/TFM/` puede contener material de tesis local no versionado. No editarlo,
  borrarlo ni añadirlo a Git sin una petición explícita del usuario.
- No crear commits, hacer push ni reescribir el historial salvo petición
  explícita. Al terminar, informar de tests, archivos cambiados y riesgos.

## Documentación que acompaña a los cambios

- Cambios de uso, dependencias o CLI: `README.md` y `docs/uso.md`.
- Cambios de módulos o flujo: `docs/arquitectura.md`.
- Cambios de formulación: `docs/formulacion.md` y el documento del elemento.
- Nuevas referencias o tolerancias: `docs/validacion.md`.
- Cambio del estado general, prioridades o decisiones: `docs/TRASPASO.md` y
  `docs/roadmap.md`.
