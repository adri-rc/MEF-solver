# MEF-solver

Solver de elementos finitos lineal estático para análisis estructural,
desarrollado como Trabajo de Fin de Máster.

Para retomar el proyecto con otra cuenta o asistente de IA, comenzar por
[docs/TRASPASO.md](docs/TRASPASO.md). Las reglas de trabajo del repositorio
están en [AGENTS.md](AGENTS.md).

El solver activo implementa estática lineal para CROD, CBAR, CTRIA3 y CQUAD4
coplanar, además de estática no lineal geométrica `SOL 400` para modelos
`MAT1`/`PROD`/`CROD`. La ruta no lineal usa incrementos fijos o adaptativos,
Newton-Raphson completo y una tangente material más geométrica; está validada
con los 36 casos 92--127 de MSC Nastran 2026.1. Treinta y cinco alcanzan carga
completa y el 120 conserva el fallo de punto límite. El alcance y sus límites están
en [docs/nonlinear_crod.md](docs/nonlinear_crod.md), y el detalle de CQUAD4 en
[docs/cquad4.md](docs/cquad4.md). La lectura del informe `.rpt` y su
correspondencia con F06 se explican en
[docs/result_report.md](docs/result_report.md). La revisión obtiene
`430 passed, 6 xfailed`; los `xfail` documentan sólo
geometrías alabeadas y el suavizado `SNORM`, ambos fuera de alcance.

```powershell
python -m pip install -r requirements.txt
python main.py --file input/modelo.dat
python main.py --file input/modelo.dat --subcase 2
python -m pytest -q
```
