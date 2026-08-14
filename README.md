# MEF-solver

Solver de elementos finitos lineal estático para análisis estructural,
desarrollado como Trabajo de Fin de Máster.

Actualmente implementa CROD, CBAR, CTRIA3 y una formulación CQUAD4 validada
frente a matrices DMIG y resultados F06 de Nastran para geometrías coplanares.
Lee modelos Nastran mediante `pyNastran` y resuelve `K u = F`. El detalle y
las limitaciones del CQUAD4 están en [docs/cquad4.md](docs/cquad4.md). La
revisión final obtiene `337 passed, 6 xfailed`; los `xfail` documentan sólo
geometrías alabeadas y el suavizado `SNORM`, ambos fuera de alcance.

```powershell
python main.py --file input/modelo.dat
python main.py --file input/modelo.dat --subcase 2
python -m pytest -q
```
