# Factores PSHELL del CTRIA3: matrices de rigidez

Estos cuatro modelos de un único CTRIA3 aíslan los dos factores escalares de
`PSHELL`. Los casos K01 y K02 cambian `TS/T` manteniendo `12I/T^3=1`; los
casos K03 y K04 cambian `12I/T^3` manteniendo `TS/T=0.833333`.

Todos los decks declaran `AUTOSPC=NO`, `K6ROT=0` y solicitan la matriz `KAAX`
mediante `EXTOUT=DMIGPCH`. MSC Nastran conserva las referencias en los
archivos PCH y los CSV derivados de ellos.

MSC Nastran Student Edition 2026.1 ejecutó los cuatro decks el 3 de septiembre
de 2026 sin mensajes `FATAL` y finalizó con `END OF JOB`. Los errores relativos
de Frobenius solver/MSC son `9.96e-5` y `4.98e-4` para K01/K02, y menores que
`1.9e-10` para K03/K04. Los términos representativos de K01/K02 permanecen
por debajo de `1.23e-3` de error relativo.

El material `MID3` de K01/K02 se hace deliberadamente flexible para que su
aportación pueda aislarse sin convertir este control en una recalibración de
la formulación de placa. K03/K04 emplean un `MID3` casi nulo para aislar el
escalado de flexión. La campaña valida el uso de los factores leídos del DAT,
no la equivalencia general de la placa CTRIA3 con la formulación nativa MSC.
