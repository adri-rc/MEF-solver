# Uso del programa

## Entorno observado

El repositorio contiene un entorno virtual local con las siguientes versiones verificadas:

| Componente | Versión |
|---|---:|
| Python | 3.10.4 |
| NumPy | 1.26.4 |
| SciPy | 1.15.3 |
| pyNastran | 1.4.1 |
| pytest | 9.0.1 |

No existe todavía un `requirements.txt` ni un `pyproject.toml`. Por tanto, el entorno `.venv` permite reproducir el estado local, pero las dependencias no están declaradas de una forma portable. La creación de esa especificación es una tarea prioritaria de la [hoja de ruta](roadmap.md).

En PowerShell, desde la raíz del repositorio:

```powershell
.\.venv\Scripts\Activate.ps1
python main.py --file tests\solver_tst\validation_cases\case_1\Ej1_v1.dat
```

También puede usarse el intérprete del entorno sin activarlo:

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_1\Ej1_v1.dat
```

## Modos de ejecución

```text
python main.py --file RUTA_MODELO.dat
python main.py --all
python main.py
```

- `--file` analiza el archivo indicado. El lector acepta un deck completo o un fichero *punch* de Bulk Data.
- `--all` procesa todos los `*.dat` situados directamente en `input/`.
- Sin argumentos, procesa el primer `*.dat` encontrado en `input/`.

Aunque `parse_nastran_input()` puede leer `.bdf`, la búsqueda automática solo considera la extensión `.dat`. Los `INCLUDE` son resueltos por pyNastran respecto al deck de entrada.

El flujo de `main.py` es:

1. leer y mapear el deck;
2. construir el solver;
3. ensamblar `K` y `F`;
4. resolver desplazamientos;
5. escribir `output/<nombre>.rpt`.

Actualmente `run_analysis_job()` captura las excepciones y muestra un mensaje, pero no devuelve un código de error distinto de cero. En automatización conviene ejecutar pytest o invocar las clases directamente hasta que se corrija este comportamiento.

## Subconjunto de entrada admitido

| Grupo | Tarjetas/parámetros admitidos | Observaciones |
|---|---|---|
| Nodos | `GRID` | pyNastran transforma `CP` a coordenadas cartesianas globales. |
| Materiales | `MAT1` | Se necesitan al menos dos valores de `E`, `G` y `NU`; el tercero se calcula. |
| Propiedades | `PROD`, `PBAR`, `PSHELL` | Solo se mapea un subconjunto de los campos Nastran. |
| Elementos | `CROD`, `CBAR`, `CTRIA3`, `CQUAD4` | CQUAD4 usa por defecto la reconstrucción Nastran y completa el postproceso. |
| Cargas | `FORCE`, `MOMENT` | Solo se aplican si `CID=0`. |
| Restricciones | `SPC`, `SPC1` | Se admiten desplazamientos impuestos no nulos. |
| Parámetros | `PARAM,K6ROT` | Se utiliza en CTRIA3 y CQUAD4. |

El lector falla de forma explícita si encuentra elementos, cargas o restricciones no soportados. En particular, rechaza MPC, RBE y elementos rígidos.

### Limitaciones de selección de casos

Los identificadores `SID` se almacenan, pero el solver no interpreta `SUBCASE`, `LOAD`, `SPC`, combinaciones `LOAD` ni conjuntos `SPCADD`. Ensambla **todas** las cargas y restricciones mapeadas del deck en un único caso. Para obtener un resultado inequívoco, el fichero de entrada debe contener un solo conjunto efectivo de cargas y condiciones de contorno.

Tampoco se aplican las restricciones permanentes del campo `PS` de `GRID`.

### Convención de unidades

El programa no convierte unidades. El usuario debe emplear un sistema coherente en todo el modelo. Por ejemplo:

- SI: m, N, Pa;
- sistema habitual en modelos Nastran: mm, N, MPa (`N/mm²`).

Las dimensiones de áreas, inercias y momentos deben ser coherentes con el sistema elegido.

## Ejemplo mínimo CROD

El caso completo [Ej1_v1.dat](../tests/solver_tst/validation_cases/case_1/Ej1_v1.dat) incluye el modelo estructural mediante [model_ej1.bdf](../tests/solver_tst/validation_cases/case_1/model_ej1.bdf). Contiene:

- un `CROD`;
- una propiedad `PROD` y un material `MAT1`;
- restricciones suficientes para eliminar los grados de libertad sin rigidez;
- una fuerza axial global.

Es importante restringir los grados de libertad inactivos. El solver reserva seis GDL para cada nodo incluso cuando el elemento solo aporta rigidez a algunos de ellos. Un CROD aislado, por ejemplo, solo aporta rigidez axial y torsional respecto a su eje local.

## Resultados

`LinearStaticSolver.solve()` devuelve un diccionario:

```python
{
    "displacements": numpy.ndarray,
    "reactions": {node_id: numpy.ndarray},
    "stresses": {element_id: dict}
}
```

- `displacements` usa seis posiciones consecutivas por nodo.
- `reactions` contiene vectores de seis componentes para nodos restringidos.
- `stresses` contiene magnitudes diferentes según el tipo de elemento.

El informe de texto de `main.py` es todavía provisional:

- su tabla de tensiones tiene formato CROD, aunque pueda haber CBAR o CTRIA3;
- el contador `CROD` usa por error el número de fuerzas nodales;
- los resultados detallados de CBAR y CTRIA3 no se imprimen;
- el informe puede resumir mejor las magnitudes específicas de CQUAD4.

Para validación y desarrollo resulta más fiable inspeccionar el diccionario devuelto por `solve()`.

## Ejecución de pruebas

La colección completa se ejecuta con:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

El resultado observado y sus fallos conocidos están documentados en [validacion.md](validacion.md).
