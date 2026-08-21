# Caso lineal final: cupón de frame aeronáutico

Este caso es un modelo de integración `SOL 101` pensado para la defensa del
TFM. Representa un cupón plano de frame/web aeronáutico de `1200 x 600 mm`,
con dos huecos circulares de radio `150 mm`, una fila inferior de seis
remaches idealizados, caps y nervios y dos horquillas superiores de
introducción de carga.

No debe presentarse como una correlación experimental ni como una validación
general de Nastran. El mismo DAT se ejecutó sin modificar con MSC Nastran
Student Edition 2026.1 el 21 de agosto de 2026: terminó normalmente, sin
mensajes `FATAL`, y sus resultados cierran la validación cruzada numérica de
este modelo. El [F06 auditado](tfm_final_airframe_frame.f06) se conserva junto
al deck y las referencias seleccionadas están transcritas en el test E2E. La
ejecución corresponde al build `2026.1-2a19f742`; el SHA-256 del F06 es
`80aea144d3a084c886353e56d3e88e4290d625b6a894ec9203ffaad6ba721a5e`.

## Dimensionado y contenido

Se usa el sistema coherente `mm - N - MPa` y pequeñas deformaciones:

| Parte | Idealización | Propiedad |
|---|---|---|
| Web | 256 CQUAD4 | PSHELL 10, aluminio, `t=2.0 mm` |
| Contorno de los dos huecos | 128 CTRIA3 | PSHELL 11, `t=2.0 mm`, MID3 casi nulo |
| Cap perimetral y montante central | 56 CBAR | PBAR 20, T equivalente `40x4 + 25x3 mm` |
| Nervios superior/inferior de cada hueco | 20 CBAR | PBAR 21, sección equivalente secundaria |
| Horquillas de tiro | 4 CROD | PROD 30, `A=50 mm2`, `J=20 mm4` |

El modelo contiene 377 GRID y 464 elementos. La web es coplanar. Los GDL
`T3/R1/R2` de sus nodos están guiados por el útil para que el ensayo sea de
membrana en el plano. Esto evita atribuir precisión física a la flexibilidad
de cortante transversal de CTRIA3, que es una limitación declarada del
solver. `R3` permanece activo para el drilling de las láminas y la flexión en
plano de los CBAR.

La sección PBAR 20 representa el cap de una T formada por una ala de
`40 x 4 mm` y un alma de `25 x 3 mm`: `A=235 mm2`,
`I1=21389 mm4`, `I2=14854 mm4`, `J=1078 mm4`. Es una idealización centrada;
no se aplican offsets porque la excentricidad real cap-web requeriría una
decisión de modelización adicional.

## Cargas y condiciones de contorno

- SPC 100: guía plana del útil; además fija los GDL inactivos de los nodos de
  horquilla.
- SPC 110: los GRID `181, 184, 187, 366, 369, 372` de la fila inferior
  quedan empotrados; junto con SPC 100 equivalen a `123456`.
- SPC 120: la horquilla derecha recibe `T2=+0.8 mm`.
- SPCADD 130: combinación aplicada en el subcaso.
- LOAD 201: `Fy=+12000 N` en la horquilla izquierda y `Fx=+2500 N` en el
  encuentro superior del montante central.
- LOAD 202: `Mz=+250000 Nmm` en ese mismo encuentro.
- LOAD 200: suma de 201 y 202.

La carga por fuerza, el momento y el desplazamiento impuesto se combinan en
un único subcaso. Es útil para una diapositiva final de integración; para una
correlación más diagnóstica conviene duplicar después el deck en tres casos
independientes, sin cambiar malla ni propiedades.

## Generación y ejecución con MEF-solver

Desde la raíz del repositorio:

```powershell
.\.venv\Scripts\python.exe tests\solver_tst\validation_cases\case_linear_test_final\section_frame_model.py
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_linear_test_final\tfm_final_airframe_frame.dat
```

El informe se escribe en `output/tfm_final_airframe_frame.rpt`. El generador
es determinista y debe conservarse como fuente editable de la malla.

### Correlación con MSC Nastran 2026.1

Los puntos más claros para la defensa son:

| Magnitud | MEF-solver | MSC Nastran (F06) |
|---|---:|---:|
| Horquilla izquierda, GRID 376, `T2` | `0.744011933 mm` | `0.7440119 mm` |
| Horquilla derecha, GRID 377, `T1` | `0.055691117 mm` | `0.05569107 mm` |
| Reacción del actuador en GRID 377, `F2` | `12993.822458 N` | `12993.82 N` |
| Máxima tensión axial CROD, EID 8003/8004 | `152.976475 MPa` | `152.9765 MPa` |
| Máxima tensión normal CBAR, EID 7056 | `131.623096 MPa` | `131.6232 MPa` |
| Máximo von Mises de shell, EID 5030 | `43.074543 MPa` | `43.07457 MPa` |

La comparación del OP2 completo, antes de seleccionar puntos, produce los
siguientes errores relativos de norma L2. Así se evita apoyar la conclusión
sólo en máximos aislados:

| Campo completo | Error relativo L2 |
|---|---:|
| Desplazamientos | `1.052e-6` |
| Reacciones SPC | `1.450e-5` |
| Fuerzas / tensiones CROD | `4.945e-8 / 5.701e-8` |
| Fuerzas / tensiones CBAR | `5.138e-6 / 5.361e-6` |
| Resultantes / tensiones CTRIA3 | `6.798e-6 / 5.701e-6` |
| Resultantes / tensiones CQUAD4 | `1.073e-5 / 8.485e-6` |

El equilibrio global de fuerzas cierra a menos de `1e-6 N` y el de momentos a
menos de `1e-3 Nmm` en la prueba automatizada. Los valores MSC de varios
GRID/EID de cada familia se escriben explícitamente en
`test_linear_final_validation_case.py`; el test no depende del F06 ni del OP2
durante su ejecución. El F06 redondea sus tablas impresas, mientras que los
errores vectoriales anteriores se calcularon con la precisión del OP2.

## Apertura y ejecución en Patran

La secuencia corresponde al formulario documentado por Hexagon en el capítulo
**Read Input File** de la guía de interfaz Patran-MSC Nastran:

1. Crear una base nueva (`File > New`) y seleccionar
   `Preferences > Analysis`; elegir `MSC Nastran` como código de análisis.
   Patran no asigna unidades, por lo que debe mantenerse `mm-N-MPa`.
2. Abrir `Analysis`, seleccionar `Action = Read Input File`, pulsar
   `Select Input File` y escoger `tfm_final_airframe_frame.dat`.
3. En `Entity Selection`, mantener activados nodos, elementos, propiedades,
   materiales, cargas, restricciones y casos de carga. Aplicar la importación.
4. Revisar `Summary Data`, `Reject Cards` y `View Log File` antes de aceptar.
   No deben faltar GRID, CQUAD4, CTRIA3, CBAR, CROD, PSHELL, PBAR, PROD,
   MAT1, FORCE, MOMENT, SPC/SPC1/SPCADD, LOAD ni el subcaso.
5. Para ejecutar desde la base importada: `Analysis`,
   `Action = Analyze`, solución lineal estática, modelo completo y ejecución
   completa. Conservar el BDF que regenere Patran para comprobar que no ha
   modificado la intención del deck original.
6. La vía más auditable es ejecutar directamente el DAT original con MSC
   Nastran y leer después el OP2: `Analysis > Access Results`, objeto
   `Read Output2`, seleccionando el OP2 generado. `PARAM,POST,-1` y las
   peticiones `PLOT` ya están incluidas.

Guía oficial de referencia: [Patran Interface to MSC Nastran Preference Guide,
Volume 1, Structural Analysis](https://documentation-be.hexagon.com/bundle/Patran_2024.1_Interface_to_MSC_NASTRAN_Preference_Guide_Volume_1_Structural_Analysis/raw/resource/enus/Patran_2024.1_Interface_to_MSC_NASTRAN_Preference_Guide_Volume_1_Structural_Analysis.pdf?save_local=true),
capítulos 3 a 5.

Para la presentación conviene crear grupos visuales por PID/familia, mostrar
el espesor de PSHELL, renderizar las secciones de viga y superponer símbolos
de FORCE, MOMENT y SPC. Una segunda vista con deformada amplificada y contorno
de von Mises en TOP/BOTTOM permite explicar qué magnitudes son directamente
comparables entre MSC y el informe del proyecto.
