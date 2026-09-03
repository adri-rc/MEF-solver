# MEF-solver

Solver de elementos finitos lineal estático para análisis estructural,
desarrollado como Trabajo de Fin de Máster.

Para retomar el proyecto con otra cuenta o asistente de IA, comenzar por
[docs/TRASPASO.md](docs/TRASPASO.md). Las reglas de trabajo del repositorio
están en [AGENTS.md](AGENTS.md).

El solver activo implementa estática lineal para CROD, CBAR, CTRIA3 y CQUAD4
coplanar, además de estática no lineal geométrica `SOL 400` para modelos
`MAT1` con CROD o CBAR. La ruta no lineal usa incrementos fijos o adaptativos
y Newton-Raphson completo. CROD se valida con los casos 92--127 y CBAR con la
campaña 128--166 de MSC Nastran; los casos mixtos CROD/CBAR 167--176 tienen un
E2E independiente. El alcance y sus límites están en
[docs/nonlinear_crod.md](docs/nonlinear_crod.md) y
[docs/nonlinear_cbar.md](docs/nonlinear_cbar.md); el detalle de CQUAD4 está en
[docs/cquad4.md](docs/cquad4.md). La lectura del informe `.rpt` y su
correspondencia con F06 se explican en
[docs/result_report.md](docs/result_report.md). La revisión obtiene
`696 passed, 6 xfailed`; todos los `xfail` documentan geometrías
CQUAD4 alabeadas o suavizado `SNORM` fuera de alcance.

La asignación de materiales PSHELL del CTRIA3 es explícita: `MID1` activa
membrana, `MID2` activa flexión y `MID3` sólo aporta el término físico de
cortante de la formulación implementada. Los campos vacíos no heredan otro
material y `MID3` sin `MID2` se rechaza como MSC Nastran. La campaña mínima
auditable está en
[`case_ctria3_pshell_material_semantics`](tests/solver_tst/validation_cases/case_ctria3_pshell_material_semantics/README.md).

El caso lineal de integración para la defensa está en
[`case_linear_test_final`](tests/solver_tst/validation_cases/case_linear_test_final/README.md):
combina 384 shells, 76 CBAR, 4 CROD, fuerzas, momento, empotramientos y un
desplazamiento impuesto. El mismo DAT terminó sin `FATAL` en MSC Nastran
2026.1; la referencia F06 está transcrita en el E2E y los errores relativos
globales de desplazamientos, reacciones y resultados de elemento son menores
de `1.5e-5`.

La evolución tridimensional está en
[`case_linear_test_final_T`](tests/solver_tst/validation_cases/case_linear_test_final_T/README.md):
forma una sección de frame en T con 417 GRID, 288 CQUAD4, 128 CTRIA3, 88
CBAR y 2 CROD, cuatro uniones puntuales tipo remache, dos fuerzas, un momento
y movimiento diferencial impuesto. El E2E transcribe el F06 de MSC Nastran
2026.1 para ambos CROD y 15 elementos de cada familia restante. La respuesta
global, la membrana de láminas y las tensiones de fibra correlacionan; los
momentos y cortantes centroidales CTRIA3 quedan documentados como límite de
recuperación, no como resultado validado.

El diagnóstico complementario
[`case_linear_test_final_T_no_CTRIA_MID3`](tests/solver_tst/validation_cases/case_linear_test_final_T_no_CTRIA_MID3/README.md)
cambia únicamente `PSHELL 11.MID3` a un `MAT1 99` con `E=1e-20`. Con esa
flexibilidad transversal degenerada, los momentos CTRIA3 pasan a un error L2
completo de `7.106e-5` y ambos códigos recuperan cortantes de cero ingenieril.
Las demás familias y los GDL no asociados a `T3` conservan la correlación.
Esto localiza la discrepancia en la formulación de cortante/recuperación, pero
no valida físicamente la variante: su condición es `5.660208e28` y `T3` se
convierte en un casi-mecanismo que no correlaciona. La recuperación física de
`M/Q` CTRIA3 como Nastran permanece fuera del alcance implementado.

El caso no lineal final está en
[`case_nonlinear_test_final`](tests/solver_tst/validation_cases/case_nonlinear_test_final/README.md):
es un puente arco atirantado espacial de 32 GRID, 48 CBAR y 30 CROD, con 15
fuerzas, 3 momentos, empotramiento y asiento diferencial impuesto. Converge
en 50 incrementos sin recortes tanto en MSC Nastran 2026.1 como localmente y
ya forma parte del E2E mixto. Los errores L2 son `8.24e-5` en desplazamientos
y `1.62e-4` en reacciones.

La ampliación de escala está en
[`case_nonlinear_test_final_updated`](tests/solver_tst/validation_cases/case_nonlinear_test_final_updated/README.md):
un puente arco espacial de tres vanos y terreno escalonado con 202 GRID, 144
CBAR y 428 CROD. Combina 44 fuerzas, 3 momentos y asientos de cimentación de
`-2/-6/-10 mm`, con `AUTOSPC=NO`. MSC Nastran 2026.1 y el solver local
completan 20 incrementos sin recortes. Los errores L2 son `3.28e-5` en
traslaciones, `2.90e-5` en reacciones de fuerza y `2.13e-5` en tensiones CROD;
la recuperación completa de fibras CBAR queda cuantificada por separado
(`1.67 %` L2 y `11.57 %` máximo significativo en el larguero de tablero).

El CBAR no lineal valida las dos variantes `PBARN1,SECT=S,BEH2=BAR`:
`INT2=LC` usa interpolación lineal/cúbica Euler--Bernoulli y `INT2=LS` usa
interpolación lineal con cortante en un punto y factores `PBAR K1/K2`.
El `.rpt` presenta su respuesta axial/torsional, longitud y energía del estado
convergido, acciones locales en ambos extremos y tensiones C--F.

```powershell
python -m pip install -r requirements.txt
python main.py --file input/modelo.dat
python main.py --file input/modelo.dat --subcase 2
python -m pytest -q
```
