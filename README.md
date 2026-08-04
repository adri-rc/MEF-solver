# MEF-solver

Solver de elementos finitos lineal estático para análisis estructural,
desarrollado como Trabajo de Fin de Máster.

Actualmente implementa CROD, CBAR, CTRIA3 y una formulación CQUAD4 en
desarrollo. Lee modelos Nastran mediante `pyNastran` y resuelve `K u = F`.

```powershell
python main.py --file input/modelo.dat
python main.py --file input/modelo.dat --subcase 2
python -m pytest -q
```
