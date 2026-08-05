# Limitaciones y hoja de ruta

## Estado de partida

El flujo estático lineal y la suite de 298 pruebas están operativos. CROD,
CBAR, CTRIA3 y CQUAD4 se leen, ensamblan y recuperan resultados. CQUAD4 dispone
de formulaciones Nastran y Reissner-Mindlin separadas; los casos DMIG 51-54
cumplen sus tolerancias documentadas.

El proyecto sigue siendo un solver académico de alcance controlado, no un
sustituto general de Nastran.

## Limitaciones principales

- solo estática lineal, pequeñas deformaciones y material `MAT1` isótropo;
- sin cargas distribuidas, gravedad, temperatura, MPC/RBE ni contacto;
- una ejecución resuelve un único subcaso seleccionado;
- campos shell como `MID4`, `Z1/Z2`, `ZOFFS`, `TFLAG` y espesores nodales no
  tienen todavía tratamiento físico completo;
- CQUAD4 alabeado coincide con la DMIG del caso 54 al `2.15 %`, no término a
  término, y necesita más referencias independientes;
- no hay sólidos 3D, dinámica, plasticidad ni no linealidad geométrica;
- el informe de texto aún puede especializarse mejor por tipo de elemento;
- `assemble()` no es idempotente y no existe CI publicada.

## Prioridades

### P0 - Consolidar la validación 2D

1. Añadir casos CQUAD4 alabeados, de cortante y con diferentes métricas de
   calidad geométrica.
2. Ejecutar *patch tests* multielemento y convergencia de placas.
3. Cerrar las discrepancias CTRIA3 conocidas con métricas absolutas y
   relativas separadas.
4. Registrar versión, deck y procedencia de cada DMIG/F06.

### P1 - Completar semántica shell

1. Implementar `MID4`, `12I/T^3`, `TS/T`, `Z1/Z2`, `ZOFFS`, `TFLAG` y `Ti`.
2. Validar materiales ortótropos y laminados antes de generalizar la
   estabilización isotrópica actual.
3. Separar recuperación y escritores en `src/postprocessing/`.

### P2 - Robustez y automatización

1. Añadir CI con pytest, `compileall` y resumen numérico de validación.
2. Hacer `assemble()` idempotente o impedir explícitamente la segunda llamada.
3. Mejorar diagnósticos de condicionamiento, mecanismos y componentes
   desconectadas.
4. Perfilar el ensamblaje antes de optimizarlo.

### P3 - Nuevas capacidades

1. Cargas distribuidas, gravedad, temperatura y restricciones multipunto.
2. Arquitectura para resolver varios subcasos en una ejecución.
3. `CTETRA4` como primer sólido, seguido de `CHEXA8` con *patch tests* 3D.
4. Solo después, dinámica y no linealidad si entran en el alcance del TFM.

## Entregables recomendados para la tesis

- alcance y subconjunto Nastran soportado a nivel de campo;
- derivación y comparación de las dos formulaciones CQUAD4;
- matriz de trazabilidad requisito-prueba-caso;
- resultados de DMIG, *patch tests*, convergencia y sensibilidad geométrica;
- explicación explícita de las diferencias con Nastran;
- caso demostrador reproducible desde el deck hasta el informe.
