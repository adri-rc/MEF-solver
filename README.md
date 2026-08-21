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
campaña 128--153 de MSC Nastran. El alcance y sus límites están en
[docs/nonlinear_crod.md](docs/nonlinear_crod.md) y
[docs/nonlinear_cbar.md](docs/nonlinear_cbar.md); el detalle de CQUAD4 está en
[docs/cquad4.md](docs/cquad4.md). La lectura del informe `.rpt` y su
correspondencia con F06 se explican en
[docs/result_report.md](docs/result_report.md). La revisión obtiene
`570 passed, 7 xfailed`; seis `xfail` documentan geometrías alabeadas/SNORM y
el séptimo la interpolación CBAR `INT2=LS` aún no reproducida.

La variante CBAR validada usa `PBARN1,SECT=S,BEH2=BAR,INT2=LC`. `INT2=LS` se
aproxima mediante `PBAR K1/K2`, pero mantiene una discrepancia conocida en la
viga corta 142 y se declara mediante advertencia y `xfail` estricto.

```powershell
python -m pip install -r requirements.txt
python main.py --file input/modelo.dat
python main.py --file input/modelo.dat --subcase 2
python -m pytest -q
```
