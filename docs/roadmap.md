# Hoja de ruta

## Estado de base

El flujo DAT/BDF → modelo → ensamblaje disperso → solución → resultados está
operativo para CROD, CBAR, CTRIA3 y CQUAD4 plano. La colección completa termina
con `688 passed, 6 xfailed`; los seis límites corresponden a CQUAD4.

La rigidez CQUAD4 plana se considera cerrada tras los casos 90/91. No se
planifica extenderla a elementos alabeados en la siguiente fase.

Desde ese cierre se separó la escritura del informe en
`src/postprocessing/report_writer.py`, se añadieron sus pruebas y se verificó
el caso mixto 80 con CBAR, CTRIA3 y CQUAD4. `requirements.txt` declara
dependencias mínimas, aunque todavía falta un bloqueo reproducible.

El modelo `case_linear_test_final` añade una demostración lineal a escala de
integración para la defensa: 384 shells, 76 CBAR, 4 CROD, FORCE/MOMENT, seis
empotramientos nodales y un SPC impuesto. Su control automatizado cubre
integridad, solución y equilibrio. La validación cruzada quedó cerrada con una
ejecución sin `FATAL` de MSC Nastran 2026.1 y referencias transcritas de todas
las familias; el mayor error relativo L2 de los campos completos es `1.450e-5`
en las reacciones SPC.

`case_linear_test_final_T` extiende esa demostración a una sección espacial de
frame en T: 417 GRID y 506 elementos, incluidos un ala CQUAD4 perpendicular,
tres nervios simétricos formados por 12 CBAR y cuatro apoyos puntuales que
idealizan remaches. MSC Nastran 2026.1 termina el deck con `AUTOSPC=NO`; el E2E
compara ambos CROD y 15 elementos de cada otra familia. La respuesta global,
la membrana y las tensiones de fibra se aceptan. Los momentos/cortantes
centroidales CTRIA3 no correlacionan (`0.682/17.30` L2), sin alterar la
rigidez para forzar el ajuste.

La investigación causal se considera cerrada con
`case_linear_test_final_T_no_CTRIA_MID3`, que cambia únicamente
`PSHELL 11.MID3` a `MAT1 99`, `E=1e-20`. MSC alcanza `END OF JOB` con
`AUTOSPC=NO`; el error L2 completo de momentos CTRIA3 cae a `7.106e-5` y los
cortantes son cero ingenieril (normas menores que `4.1e-11` local y
`8.7e-24` MSC). La sonda no valida una estructura física: tiene
`cond(K)=5.660208e28` y un casi-mecanismo `T3` que no correlaciona
(`7.3586e10 mm` local frente a `-0.3421 mm` de MSC en GRID 203). Reproducir
físicamente la recuperación `M/Q` de Nastran permanece fuera de alcance.

El modelo `case_nonlinear_test_final` aporta la demostración equivalente para
`SOL 400`: un puente espacial de 48 CBAR y 30 CROD sometido a 15 fuerzas, 3
momentos y asientos diferenciales. La regresión local converge en 50
incrementos sin recortes y muestra una diferencia del `1.97 %` frente a
`LGDISP=-1`. MSC Nastran 2026.1 reproduce el estado completo con errores L2
`8.238e-5` en desplazamientos y `1.619e-4` en reacciones; la referencia ya
forma parte del E2E mixto con una envolvente del `2 %` gobernada por una fibra
CBAR al `1.36 %`.

`case_nonlinear_test_final_updated` amplía ese cierre a una geometría de tres
vanos, 202 GRID y 572 elementos con bases de pilas/estribos explícitas. Ambos
códigos completan 20 incrementos sin recortes. La respuesta global y axial
mantiene errores L2 del orden de `10^-5`; la recuperación completa de fibras
CBAR da `1.673 %` L2 y un máximo significativo del `11.57 %` en los largueros
continuos. El E2E separa las tolerancias de desplazamientos, reacciones,
estados centroidales y fibras para no ocultar esa diferencia localizada.

El preproceso distingue estática lineal y `SOL 400`, normaliza `NLSTEP` o
`NLPARM`, `LGDISP`, `FOLLOWK` y los campos relevantes de `PBARN1`.
`NonlinearStaticSolver` resuelve MAT1 con CROD o CBAR mediante Newton completo,
stepping fijo/adaptativo y recortes. CROD dispone de 36 casos (92--127), CBAR
de 39 (128--166) y la integración mixta de doce referencias (167--176 más los
dos puentes finales), todas con resultados MSC auditables.

## Prioridad 0 completada: no linealidad geométrica CROD

1. Completada la campaña MSC y conservadas sus salidas auditables.
2. Transcritos estados intermedios, tangente por sonda y resultados finales.
3. Implementados `CRODNonlinear` y `NonlinearStaticSolver` sin modificar la
   formulación lineal.
4. Integrados CLI, recuperación e informe con configuración, historial de
   incrementos, correcciones de Newton y estado final diferenciado.
5. Documentado el punto límite 120 y el alcance sin arco de longitud.

## Prioridad 0 completada: no linealidad geométrica CBAR

1. Implementado `CBARNonlinear` corrotacional sin modificar la ruta lineal.
2. Integrados `I12`, `K1/K2`, offsets, releases, orientación `G0`, cargas y
   SPC impuestos.
3. Identificada la deformación axial logarítmica con 130/131 y validada la
   tangente mediante las sondas 137--139.
4. Transcritos 39 E2E explícitos; todos pasan, incluido 142.
5. Identificadas y validadas las interpolaciones `INT2=LC` Euler--Bernoulli e
   `INT2=LS` lineal-cortante mediante los diagnósticos 154--160. Queda abierta
   la convención de tensión de fibra bajo rotación 3D combinada.
6. Integrada la recuperación CBAR en el `.rpt`: respuesta axial/torsional,
   longitud y energía del estado no lineal, acciones A/B y tensiones C--F, sin
   alterar las tablas lineales existentes.

### Cierre externo ejecutado

Los seis diagnósticos 161--166, los diez modelos mixtos CROD/CBAR 167--176
(de 2 a 30 elementos) y los puentes finales de 78 y 572 elementos terminan sin
`FATAL` y alcanzan factor 1.0 en MSC. Las referencias finales están
transcritas: 161--166 amplían el E2E CBAR; 167--176 y ambos puentes usan un
módulo E2E mixto independiente. La respuesta global y la tensión axial
centroidal quedan cerradas; sigue abierta la reproducción estricta de algunas
fibras casi neutras y del reparto de momentos en los tramos más rotados de
165/166. El informe presenta esas fibras con los nombres C--F de Nastran y
documenta qué magnitudes tienen correspondencia directa o reconstruida en F06.

## Prioridad 1: consolidación

1. Añadir integración continua para ejecutar pytest en cada cambio.
2. Añadir `pyproject.toml` o un fichero de bloqueo con versiones reproducibles.
3. Añadir pruebas de convergencia sobre mallas refinadas, además de DMIG/E2E.
4. Definir una política de ejecución por lotes: detenerse o continuar y
   acumular errores por modelo.
5. Revisar y aislar los módulos y scripts antiguos o exploratorios antes de
   ampliar el alcance físico.

El desacuerdo centroidal `M/Q` CTRIA3 del frame T ya está aislado como límite
de formulación mediante la sonda sin MID3 efectivo. No queda como tarea de
consolidación ni justifica recalibrar la rigidez validada.

## Prioridad 2: rendimiento

1. Ensamblar mediante tripletas COO vectorizadas.
2. Cachear transformaciones y objetos de elemento entre ensamblaje y
   recuperación.
3. Medir lectura, rigidez elemental, ensamblaje, solución y postproceso sobre
   familias de mallas crecientes.
4. Estudiar reordenamientos y solvers iterativos cuando el tamaño lo justifique.

## Prioridad 3: nuevas capacidades

Cada ampliación debe incorporar un caso mínimo, una referencia externa y una
prueba E2E antes de anunciarse como soportada. Candidatos:

- presión de lámina, gravedad y otras cargas distribuidas;
- MPC y elementos rígidos;
- propiedades avanzadas de vigas;
- materiales ortótropos y laminados;
- elementos de mayor orden y sólidos;
- exportación a formatos de visualización.

## Mejoras opcionales de CQUAD4

- Investigar la recuperación centroidal propietaria de `Qx/Qy` sin tocar la
  rigidez validada.
- Añadir resultados en más estaciones o vértices.
- Implementar `SNORM` sólo si aparecen modelos facetados que lo requieran.
- Mantener el rechazo de nodos no coplanares. Una futura formulación alabeada
  sería un elemento nuevo, no una relajación de la tolerancia actual.
