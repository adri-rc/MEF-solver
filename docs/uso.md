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

El repositorio declara dependencias mínimas en `requirements.txt`, pero todavía
no dispone de un fichero de bloqueo con versiones exactas ni de
`pyproject.toml`. En un clon nuevo de Windows puede prepararse el entorno así:

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Los decks y referencias Nastran se almacenan con Git LFS. Después de clonar,
conviene ejecutar `git lfs pull` y comprobar que un `.dat` contiene el deck
real en lugar del texto `version https://git-lfs.github.com/spec/v1`.

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

Para modelos lineales, el flujo de `main.py` es:

1. leer y mapear el deck;
2. construir el solver;
3. ensamblar `K` y `F`;
4. resolver desplazamientos;
5. escribir `output/<nombre>.rpt`.

Un deck `SOL 400` con `ANALYSIS=NLSTAT` se dirige a
`NonlinearStaticSolver`. Esta ruta está separada de `LinearStaticSolver` y se
limita a MAT1 con PROD/CROD o PBAR/CBAR, cargas nodales conservativas y
Newton-Raphson completo.
El diccionario de resultados añade historial de incrementos/iteraciones,
factor de carga, energía elástica y número de recortes.

`run_analysis_job()` muestra el contexto del fallo y vuelve a lanzar la
excepción. Por tanto, la ejecución desde `main.py` termina con código distinto
de cero ante un error; en modo `--all`, el primer trabajo fallido detiene el
lote.

## Subconjunto de entrada admitido

| Grupo | Tarjetas/parámetros admitidos | Observaciones |
|---|---|---|
| Nodos | `GRID` | pyNastran transforma `CP` a coordenadas cartesianas globales. |
| Materiales | `MAT1` | Se necesitan al menos dos valores de `E`, `G` y `NU`; el tercero se calcula. |
| Propiedades | `PROD`, `PBAR`, `PBARN1`, `PSHELL` | PBARN1 conserva `SECT`, `BEH2` e `INT2` para SOL 400. |
| Elementos | `CROD`, `CBAR`, `CTRIA3`, `CQUAD4` | CQUAD4 usa por defecto la reconstrucción Nastran y completa el postproceso. |
| Cargas | `FORCE`, `MOMENT`, `LOAD` | Los vectores se transforman desde su `CID`; se reducen combinaciones. |
| Restricciones | `SPC`, `SPC1`, `SPCADD`, `GRID.PS` | Se admiten desplazamientos impuestos no nulos. |
| Parámetros | `PARAM,K6ROT`, `PARAM,LGDISP`, `PARAM,FOLLOWK` | `K6ROT` se utiliza en láminas; `LGDISP=-1/1/2` controla CROD/CBAR y `FOLLOWK` debe ser `NO`. |
| Control de análisis | `SOL 101`; `SOL 400`, `ANALYSIS=NLSTAT`, `NLSTEP` y `NLPARM` | SOL 400 se resuelve para el subconjunto MAT1 con CROD/CBAR. |

El lector falla de forma explícita si encuentra elementos, cargas o restricciones no soportados. En particular, rechaza MPC, RBE y elementos rígidos.

### Selección de casos

Si el deck contiene un solo `SUBCASE`, se selecciona automáticamente. Si
contiene varios, debe indicarse `--subcase ID`; omitirlo genera un error para
evitar mezclar casos. El lector reduce los `LOAD` y `SPCADD` asociados. En un
fichero Bulk Data sin Case Control se mapean las tarjetas base disponibles.

Las restricciones permanentes del campo `PS` de `GRID` también se aplican.

### Entrada no lineal CROD/CBAR

El lector normaliza `NLSTEP` con continuaciones `GENERAL`, `FIXED` o `ADAPT`
y `MECH`. Admite también `NLPARM` heredado. Guarda el resultado en
`model.analysis`, junto con `LGDISP` y `FOLLOWK`. El solver ejecuta `FIXED` o
`ADAPT`, realiza recortes por bisección y actualiza la tangente en cada
iteración. No admite `ARCLN`, controles térmicos/acoplados, cargas seguidoras,
actualizaciones parciales de rigidez ni más de una estrategia incremental por
entrada.

Los desplazamientos impuestos no nulos se introducen con el campo `D` de
`SPC`; el solver no lineal los escala con el factor de carga hasta alcanzar el
valor completo. La tarjeta de carga `SPCD` no está implementada. `K6ROT` sólo
se aplica a las láminas y no estabiliza los GDL rotacionales de CROD.

Los ejemplos ejecutables por MSC Nastran están en `case_92` a `case_127` y
todos conservan las salidas auditadas de MSC Nastran 2026.1. Los 35 casos que
alcanzan carga completa se comparan con valores explícitos; el 120 documenta
la terminación esperada cerca del punto límite. Los modelos 123--127 amplían
la campaña a estructuras de 10--30 CROD. La guía está en
[`NONLINEAR_CROD_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CROD_CASES.md).
El significado de las tarjetas está desarrollado en
[`nonlinear_crod.md`](nonlinear_crod.md).

Para CBAR, la combinación validada es `PBAR` más
`PBARN1,SECT=S,BEH2=BAR`, con `INT2=LC` o `INT2=LS`. Se admiten `I12`,
`K1/K2`, offsets, releases y SPC no nulos. `LC` selecciona flexión
Euler--Bernoulli lineal/cúbica y no usa `K1/K2`; `LS` selecciona el elemento
lineal-cortante de un punto y sí usa esos factores. Los casos 128--160 y sus
referencias se documentan en [`nonlinear_cbar.md`](nonlinear_cbar.md).

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

`NonlinearStaticSolver.solve()` conserva esas tres claves para el informe y
añade:

```python
{
    "increments": list[NonlinearIncrementResult],
    "iterations": list[NonlinearIterationResult],
    "load_factor": float,
    "strain_energy": float,
    "cutbacks": int,
}
```

- `displacements` usa seis posiciones consecutivas por nodo.
- `reactions` contiene vectores de seis componentes para nodos restringidos.
- `stresses` contiene magnitudes diferentes según el tipo de elemento.

El informe `output/<nombre>.rpt` se genera desde
`src/postprocessing/report_writer.py`. El cabecero identifica explícitamente
`LINEAR STATIC` o `NONLINEAR STATIC`, y tras el resumen del modelo muestra la
secuencia `SOL`, el tipo de análisis y los parámetros globales aplicables.

En SOL 400 se imprimen además:

- los parámetros principales de `NLSTEP`/`NLPARM`, incluidos método y tamaños
  de incremento, límites de Newton/recortes, criterios y tolerancias;
- un resumen con factor final, incrementos, iteraciones, recortes y energía;
- después del estado final, el historial de incrementos convergidos y el de
  cada corrección de Newton con factor de carga, residuos y norma de la
  corrección.

El bloque `FINAL CONVERGED STATE RESULTS` contiene el estado final. Sus
resultados se organizan por familia:

- CROD: fuerzas y torsor; tensiones, deformación y giro unitario;
- CBAR: fuerzas y momentos en los extremos A/B; tensiones en los puntos PBAR;
- CTRIA3 y CQUAD4: una tabla de resultantes `Nx`, `Ny`, `Nxy`, `Mx`, `My`,
  `Mxy`, `Qx`, `Qy` y otra de tensiones en las fibras `TOP` y `BOTTOM`;
- desplazamientos y reacciones conservan tablas nodales independientes.

Todas las componentes se identifican por su nombre físico; ya no se utilizan
cabeceras genéricas `PRIMARY 1`, `PRIMARY 2`, etc. Las magnitudes se expresan
en el sistema coherente de unidades elegido en el deck de entrada. Para placas,
`TOP` corresponde a `z=+t/2` y `BOTTOM` a `z=-t/2` en los ejes locales.
El significado de las tablas, las razones por las que cambia el número de
iteraciones y su correspondencia con el F06 están en
[`result_report.md`](result_report.md).

## Ejecución de pruebas

La colección completa se ejecuta con:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

El resultado observado y sus fallos conocidos están documentados en [validacion.md](validacion.md).
