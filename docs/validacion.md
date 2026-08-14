# Validación

## Estado final observado

Comando:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Resultado del 14 de agosto de 2026:

```text
337 passed, 6 xfailed
```

Las advertencias visibles proceden de que pytest no puede escribir su caché en
el entorno de ejecución; no afectan a los resultados.

## Capas de comprobación

1. Modelo de datos, lector, referencias y sistemas de coordenadas.
2. Funciones de forma, jacobianos, transformaciones y cuadratura.
3. Propiedades matemáticas de cada elemento: simetría, modos rígidos,
   orientación, desacoplamientos y geometrías inválidas.
4. Matrices DMIG frente a Nastran.
5. Validación E2E de desplazamientos, cargas, reacciones y resultados de
   elemento con referencias transcritas en los tests.

## Matrices CQUAD4

El test `test_E2E_CQUAD_Nastran_K_solution.py` aplica dos criterios
complementarios:

- error relativo de Frobenius para medir la discrepancia global;
- comparación término a término mediante `allclose`.

Además repite la comparación en ejes del elemento para membrana,
placa/cortante, drilling y sus acoplamientos. Esta descomposición evita que la
mayor magnitud de la membrana oculte un término incorrecto de cortante.

Los casos activos cubren:

- membrana, flexión y combinación de ambas;
- `K6ROT=0/100`;
- `MID3` vacío y activo;
- cambios de `G3` y `TS/T`;
- rectángulos, trapecios y cuadriláteros afines sesgados;
- orientación arbitraria de un plano en 3D;
- DMIG completas de oblicuidad moderada y alta.

El caso 90 reveló que la ley oblicua anterior era demasiado blanda en alta
oblicuidad. Se mantuvieron el operador por aristas, la flexibilidad residual y
la contribución MID3 y se ajustó únicamente la métrica torsional racional. El
resultado cumple:

- matriz completa: error de Frobenius menor que `1e-6`;
- placa/cortante: menor que `4e-6`;
- términos del bloque: `rtol=1e-5`, `atol=2e-3`.

## E2E CQUAD4

`test_validation_end_to_end_cquad4.py` no lee F06 durante su ejecución. Los
valores de Nastran están escritos directamente en los casos, lo que permite
auditar cada referencia.

El caso 91, pareja cargada del DMIG 90, confirma el ajuste con dos paneles y
propiedades distintas:

- error vectorial de desplazamientos: `2.07e-5`;
- error vectorial de reacciones: `6.30e-5`;
- `N`, `M` y tensiones de fibra: menos de 0.1 % por componente;
- mayor discrepancia `Qx/Qy`: 47.1 %, con 0.0216 absoluto.

La diferencia de `Q` pertenece a la recuperación centroidal, no a la matriz de
rigidez. No se observa una discrepancia asociada en desplazamientos, reacciones
o momentos.

### Casos 57–59

- **57:** `cond(K)=2.90e25`; errores vectoriales `U=1.72e-3`, `R=4.07e-3` y
  resultados `9.28e-4`. La mayor diferencia relativa de desplazamiento está
  en una rotación casi nula.
- **58:** `cond(K)=8.52e26`; traslaciones próximas al 0.3 %, pero rotaciones,
  momentos secundarios y `Qy` sensibles. Error vectorial de resultados
  `2.94e-2`.
- **59:** geometría alabeada, fuera de alcance. El test es `xfail` estricto y
  espera que el solver la rechace.

## Significado de los seis xfail

| Casos | Motivo | Interpretación |
|---|---|---|
| 54, 59, 62, 63, 64 | CQUAD4 alabeados | Fuera del alcance plano; rechazo esperado |
| 78 | suavizado `SNORM` por defecto | No implementado; usar la variante `SNORM,0` |

Son límites declarados, no errores pendientes del CQUAD4 plano. Los casos 82 y
88 confirman ensamblajes facetados con normales independientes.

## Criterio de cierre

La implementación CQUAD4 se considera cerrada para elementos coplanares porque:

- sus invariantes y modos rígidos están cubiertos;
- las DMIG activas coinciden globalmente, por bloques y término a término;
- los modelos con varios elementos, materiales y orientaciones reproducen la
  respuesta global de Nastran;
- el dominio no soportado se rechaza en código y se documenta en tests.

La única mejora específica no bloqueante es reproducir con mayor precisión la
convención de impresión centroidal de `Qx/Qy`.
