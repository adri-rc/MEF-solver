# Estado y validación

## Principio de validación

Para un solver de elementos finitos no basta con obtener desplazamientos plausibles. La validación debe separar:

1. geometría y transformación;
2. funciones de forma, Jacobiano y matrices B;
3. ley constitutiva;
4. matriz local;
5. transformación y ensamblaje global;
6. imposición de contorno y solución;
7. reacciones y recuperación de tensiones;
8. convergencia de malla.

El repositorio ya sigue parcialmente esta estructura mediante pruebas unitarias, *patch tests* y casos comparados con archivos F06/PCH/CSV de MSC Nastran.

## Resultado reproducido en esta auditoría

Fecha: **4 de agosto de 2026**. Entorno: `.venv` incluido en el repositorio.

La ejecución global:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

se detiene durante la colección con tres errores:

1. dos pruebas de `src/preprocessing/deprecated/` importan `src.preprocessing.nastran_reader`, aunque el fichero está dentro de `deprecated/`;
2. `CTRIA3ElementSolver_test.py` importa `ctria3_element_FINALL`, módulo que ya no existe en fuentes, aunque queda un `.pyc` antiguo.

Al limitar la ejecución a `tests/` e ignorar ese test CTRIA obsoleto:

```powershell
.\.venv\Scripts\python.exe -m pytest tests -q --ignore=tests/solver_tst/elements/CTRIA3ElementSolver_test.py
```

el resultado observado fue:

```text
193 passed, 8 failed
```

Los ocho fallos son explicables y reproducibles:

- 7 fallos CQUAD4: `calculate_stresses()` no existe en la clase activa porque su implementación está comentada (seis casos parametrizados y un *patch test*).
- 1 fallo de datos: el caso 18 se espera en `validation_cases/case_18`, pero está almacenado en `validation_cases/case_22/case_18`.

Por tanto, las 193 pruebas correctas son una señal sólida para el núcleo activo, pero no equivalen a una suite global limpia.

## Cobertura por componente

### Preproceso y modelo de datos

Las pruebas cubren:

- creación y consulta de entidades;
- cierre de constantes elásticas y errores de consistencia;
- lectura de deck completo y *punch*;
- `INCLUDE` simples y múltiples;
- mapeo CROD, CBAR, CTRIA3 y CQUAD4;
- cargas, restricciones y rechazo de entidades no soportadas;
- puntos de recuperación PBAR.

Faltan pruebas sobre selección real de SUBCASE/SID, combinaciones de carga, `GRID.PS`, sistemas de coordenadas de cargas y semántica completa de PSHELL.

### Utilidades y solver global

Se prueban:

- mapeo de IDs arbitrarios a GDL;
- reglas de Gauss de orden 1 a 3;
- ortogonalidad, conservación de energía y regla de la mano derecha en transformaciones;
- partición libre/restringida;
- desplazamientos impuestos;
- detección de matrices singulares;
- reacciones y tensiones de CROD/CBAR;
- ensamblajes rotados y multimaterial.

### CROD y CBAR

Hay ensayos analíticos de rigidez, transformaciones 3D y recuperación. Además, `test_validation_end_to_end.py` declara 18 casos Nastran para barras y estructuras mixtas. Uno de esos casos no se ejecuta por la ubicación incorrecta del caso 18; los restantes pasaron en la ejecución auditada.

### CTRIA3

`test_validation_end_to_end_ctria3.py` declara 22 casos que incluyen:

- membrana, flexión y cortante;
- planos XY, XZ, YZ y planos 3D arbitrarios;
- mallas de varios elementos;
- materiales múltiples y elementos distorsionados;
- distintos valores de K6ROT;
- una estructura mixta CTRIA3/CBAR.

Hay además 10 comparaciones explícitas de `K_global` con matrices CSV exportadas de Nastran. Estos casos pasaron en la ejecución auditada. No obstante, el test unitario histórico más extenso de CTRIA3 no se recoge por su import obsoleto; debe actualizarse para que valide la clase que realmente usa el solver.

## Estado de CQUAD4

Las pruebas activas que no dependen de recuperación de tensiones sí pasan para seis geometrías: cuadrado, trasladado, distorsionado, con materiales desacoplados, con gran *offset* de coordenadas y ligeramente alabeado. Comprueban forma `24 x 24`, simetría, ejes locales, desacoplamiento membrana/flexión y errores geométricos.

Sin embargo, la prueba denominada `test_rigid_body_modes` solo verifica que los seis autovalores menores sean pequeños; no verifica que existan **exactamente** seis. Para el cuadrado canónico auditado se observaron 12 autovalores bajo el umbral usado por el propio test. La prueba debe comprobar nulidad, modos físicos y modos de energía cero por separado.

### Comparación exploratoria con matrices Nastran

Se regeneró la matriz del código actual con la geometría y propiedades de los scripts de los casos 51 y 52, y se comparó directamente con los CSV Nastran de igual orden de GDL.

| Caso | Bloque aislado | Error relativo de Frobenius | Mediana del error relativo en entradas Nastran no nulas |
|---|---|---:|---:|
| 51 | Membrana | 12.56 % | 12.80 % |
| 52 | Flexión | 78.02 % | 61.42 % |

Estas cifras son diagnósticas, no criterios finales de aceptación. Confirman que la matriz CQUAD4 actual todavía no reproduce la referencia, especialmente en flexión. El script `run_cquad4_validation.py` exporta CSV manualmente y contiene rutas absolutas; la comparación debe convertirse en pruebas automáticas, portables y sin escritura sobre los datos de referencia.

## Criterios recomendados

Para cada elemento nuevo o corregido:

- simetría: `||K-Kᵀ|| / ||K||` cercana a precisión de máquina;
- invariancia ante traslación y rotación del elemento;
- exactamente los modos de cuerpo rígido esperados después de tratar *drilling*;
- energía no negativa para desplazamientos arbitrarios;
- *patch test* de deformación constante;
- reproducción analítica en geometrías regulares;
- sensibilidad controlada a distorsión y espesor;
- comparación independiente con Nastran u otra referencia;
- convergencia de desplazamientos y energía al refinar malla;
- equilibrio global: suma de cargas y reacciones dentro de tolerancia.

Para comparar coeficientes de matrices conviene usar simultáneamente tolerancia absoluta y relativa:

```python
np.testing.assert_allclose(K_python, K_ref, rtol=rtol, atol=atol)
```

Los valores cercanos a cero deben evaluarse con `atol`; aplicar solo error relativo produce cocientes sin significado. También conviene comparar bloques físicos separados antes de comparar la matriz completa.

## Acciones inmediatas sobre la suite

1. Añadir configuración pytest para recoger únicamente `tests/`.
2. Actualizar o retirar imports de módulos obsoletos.
3. Mover el caso 18 a la ruta declarada o corregir su `folder`.
4. Convertir casos 51/52 en tests CQUAD4 automáticos.
5. Evitar que los tests escriban CSV dentro de `validation_cases`.
6. Registrar tolerancias y origen/versionado de cada referencia Nastran.
