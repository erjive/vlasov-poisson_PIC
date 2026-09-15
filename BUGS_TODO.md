# Lista de bugs y mejoras — seguimiento

Este repo (`vlasov-poisson_PIC`) es, en el código fuente, idéntico byte a
byte a `old_VlasovPoisson_PIC_sp` — la versión histórica, previa al
refactor que introdujo `l_part` (arreglo, dispersión en momento angular
$L$) en `VlasovPoisson_PIC_sp`. Este repo sigue usando `Lfix` (escalar,
$L$ fijo para todas las partículas).

Esta lista rastrea el port a esta rama (`fix/bugs-mejoras`) de los bugs y
mejoras ya encontrados y corregidos en `VlasovPoisson_PIC_sp` (rama
`fix/bugs-mejoras`, ver su propio `BUGS_TODO.md`) que son
**independientes de la arquitectura** `Lfix` vs `l_part` — es decir, el
mismo bug con el mismo síntoma, verificado leyendo el código de este
repo antes de portar cada fix (no asumido por analogía). Un commit por
ítem resuelto.

## Resueltos

- [x] **Makefile: `FLAGS` vacío para gfortran.** Idéntico al bug
  original: las 4 líneas candidatas de `FLAGS` del bloque `gfortran`
  estaban comentadas, rompiendo el build con gfortran (faltaba
  `-Jobjs`, se perdía `-fopenmp`). — commit `fix(build): activate
  gfortran FLAGS in Makefile`

- [x] **`density.f90`: ghost zones mal reflejadas** en `density` y
  `avg_density` (`rho(i-1)=rho(i)`/`avg_rho(i-1)=avg_rho(i)` en vez de
  `rho(1-i)=rho(i)`), mismo bug exacto que en el otro repo. — commit
  `fix(density): correct ghost-zone mirroring in density/avg_density`

## Pendientes de portar (confirmados presentes en este repo)

- [x] **`density.f90`/`poisson_rk.f90`: `collapse(2)` sin protección
  (condición de carrera OpenMP).** Confirmado en el loop combinado
  `rho`/`curr`/`avg_rho` de `density()` y en el loop de interpolación
  `pot_part`/`force_part` de `poisson_rk()`. Arreglado paralelizando
  solo en el índice externo (mismo patrón que ya usa correctamente
  `avg_density()`, que además ya tenía un parche con `!$OMP ATOMIC`
  para el mismo problema — correcto pero subóptimo, sigue como está
  por ahora; se puede alinear al mismo patrón junto con el cell-list
  de abajo). Verificado: `density`/`avg_density`/`force`/`potential`/
  `energy` dan salida idéntica byte a byte entre 1 y 8 hilos
  (autogravitante, `spatial_output=5`). — commit `fix(openmp): remove
  collapse(2) data race in density() and poisson_rk()`
- [x] **`utils.f90` `deallocate_mem`: código muerto y roto.** Mismos
  bugs que en el otro repo: `deallocate(p_part_hp)` — typo por
  `p_part_h` (`p_part_hp` nunca se allocatea); `deallocate(res)` bajo
  `conv_test=="on"` pero `res` nunca se allocatea en
  `alloc_mem_set0`; `force`/`pot`/`dev_pot` desallocateados sin
  comprobar `autointeraction` y luego otra vez si
  `autointeraction=.true.` (doble free). No se llamaba desde ningún
  lado. Arreglado para reflejar exactamente `alloc_mem_set0`, y se
  agregó la llamada real al final de `main.f90`. Probado con
  `gaussian1` (2500 partículas) en ambas configuraciones
  (`autointeraction` `.true.`/`.false.`): terminan limpio con "Memory
  deallocated". — commit `fix(memory): repair deallocate_mem and
  actually call it`
- [x] **`utils.f90` `set_timestep`: `dt = dtr` anulaba
  `min(dtr,dtp)`.** Misma línea suelta después del `if/else` que
  pisaba el resultado incondicionalmente. Portado junto con el mismo
  paquete de mejoras del otro repo: criterio de aceleración
  `dtp=courant*sqrt(2*drc/Fmax)` (menos restrictivo que
  `courant*dpc/Fmax`), guarda contra `Fmax=0`, y la condición de
  entrada `BGtype/="null" .or. autointeraction` (antes solo
  `BGtype/="null"`, por lo que con fondo nulo y autogravedad `Fmax`
  nunca se recalculaba). Probado con `gaussian1`, ambas
  configuraciones de `autointeraction`. — commit `fix(timestep):
  restore force-based dt bound with a less restrictive, physically
  motivated criterion`
- [x] **`utils.f90` `save_data`: `save1Ddata` recibe arreglos con
  ghost cells directo** (`r`, `rho`, `avg_rho`, `curr`, `force`,
  `pot`, todos `(1-ghost:Nr)`) en vez de recortarlos a `(1:Nr)` antes
  de pasarlos al dummy explícito `(1:Nr)` — mismo corrimiento de
  índice por asociación de secuencia. Verificado el `.rl` de salida
  antes/después: antes arrancaba en `r` negativo (punto fantasma);
  después arranca en `r=dr/2` (primer punto físico), como debe ser.
  — commit `fix(io): correct index shift when saving ghost-augmented
  grid arrays`
- [x] **`utils.f90` `reduce_arrays`: operadores de comparación
  inconsistentes** (`r_part(i)<=rmax` en el conteo vs
  `r_aux(i)<rmax` en la copia). — commit `fix(reduce_arrays,grid):
  consistent <=rmax comparison; allocate r(0:Nr) for rmin>0`
- [x] **`functions.f90` `Sn`/`Wn`: no rechazan `n<1`** (`else if
  (n>4)`/`else if (n>3)`, sin `else` genérico) — mismo bug. — commit
  `fix(functions): Sn/Wn now reject any invalid order, not just
  n>4/n>3`
- [x] **`utils.f90` `alloc_mem_set0`: `r` siempre se allocatea como
  `(1-ghost:Nr)`,** pero `construct_grid` llena `r(0)` cuando
  `rmin>0` (donde `ghost=0`, o sea el rango allocateado es `(1:Nr)`)
  — mismo out-of-bounds write a `r(0)`. Probado con `rmin=1.0`:
  corre limpio, termina con "Memory deallocated". — commit
  `fix(reduce_arrays,grid): consistent <=rmax comparison; allocate
  r(0:Nr) for rmin>0` (bundleado con el ítem de `reduce_arrays` de
  arriba, ambos tocan `utils.f90` y se probaron juntos)

## Mejoras de rendimiento a portar

- [x] **`density()`/`avg_density()`/`poisson_rk()`: búsqueda de
  vecinos por fuerza bruta `O(Nr×Npart)`.** Portado el cell-list de
  `utils.f90` (`build_cell_list`, `O(Nr+Npart)`), reemplazando
  también el parche `!$OMP ATOMIC` de `avg_density()` por el mismo
  patrón "paralelizar solo en el índice externo" que ya usa
  `density()`. Validado bit a bit contra la versión sin optimizar
  (autogravitante, `spatial_output=5`, 50 pasos): `density`,
  `avg_density`, `force`, `potential` y `energy` idénticos. También
  determinismo 1 vs 8 hilos verificado de nuevo sobre esta versión.
  Benchmark autogravitante (90000 partículas, 200 pasos,
  `gaussian1`): 3m46s → 2m57s (**~1.28× más rápido**, 1 hilo) — mucho
  menos que el 2.9× visto en el otro repo, esperable: esta corrida
  gasta buena parte del tiempo en la cuadratura de `phik` en
  `analysish.f90` (Simpson de 512 puntos sin optimizar, ver más
  abajo — fuera de alcance de este port), que diluye la ganancia del
  cell-list. — commit `perf(density,poisson_rk): replace
  O(Nr*Npart) brute-force deposit/interpolation with a cell list`

## Mejoras de rendimiento adicionales (esta sesión, más allá del port)

- [x] **`analysish.f90`: el loop de `phik` no estaba paralelizado
  (`!!$OMP` comentado), y recalculaba la parte del integrando que no
  depende del modo 5 veces por partícula (una por cada modo 0-4) ×2
  (una por cada función de prueba $\Phi_1,\Phi_2$).** El `!$OMP` que
  sí estaba activo (`Qr`/`Jr`) paraleliza sobre `Npart`; el que
  faltaba paralelizaba (comentado) sobre el loop *externo* de modos
  (`do i=0,mode`, solo 5 iteraciones) — mal ajuste para 8 hilos.
  Reescrito completo: loop externo ahora es sobre partículas
  (paraleliza con `Npart`, no con 5), calculando una sola vez por
  partícula la parte del integrando independiente del modo
  ($g(Q)=e^{-\sin^2(Q/2)/\sigma_Q^2}e^{-(J-J_0)^2/\sigma_J^2}J^2$,
  igual que la optimización ya hecha en `VlasovPoisson_PIC_sp`) y
  reusándola para los 5 modos, para ambas $\Phi_1,\Phi_2$ a la vez.
  Misma precisión de cuadratura que antes (512 subintervalos, sin
  reducir). De paso, mismo hallazgo que en el otro repo: el resultado
  de `phik` siempre fue real (`real(auxsum)*2.0`, la parte imaginaria
  se descartaba en silencio) — reemplazado `exp(-ilQ)` por
  `cos(lQ)` directo, sin cambiar el resultado (la reducción por
  simetría a $[0,\pi]$ ya usada es válida — verificado analíticamente
  que $\mathrm{Re}[\phi]$ es par y $\mathrm{Im}[\phi]$ impar respecto
  a $Q=\pi$).
  Validado: `hk1.tl`/`hk2.tl` para `state="aa"` coinciden con la
  versión anterior a la precisión de máquina esperada al reordenar
  una suma en coma flotante (7-8 cifras significativas, no bit a bit
  — la suma pasa de secuencial a reducción paralela). Benchmark
  (~10072 partículas tras el corte, 2000 pasos, `spatial_output=100`,
  8 hilos): **73.2s → 4.3s, ~17× más rápido** — la versión anterior
  corría 100% serial (73s de usuario ≈ 73s real), confirmando que el
  `!!$OMP` deshabilitado nunca se ejecutaba. — commit
  `perf(analysish): parallelize over particles and reuse the
  mode-independent quadrature weight`
- [x] **Intentado y revertido: paralelizar `grav_force.f90` (rama
  `Isochrone` + término centrífugo) con `!$OMP PARALLEL DO`.**
  Convertida la sintaxis de arreglos completos a loops explícitos con
  `!$OMP`, mismo patrón que ya usa la rama `sphere` en el mismo
  archivo. **Empeoró el rendimiento en vez de mejorarlo**: medido en
  aislamiento (`autointeraction=.false.`, `Isochrone`, 50000 pasos,
  ~10072 partículas, `spatial_output` grande para que
  `density`/`analysish`/`save_data` casi no se llamen), 8.5s → 26.2s
  (**~3× más lento**), con el tiempo de usuario saltando de ~9s a
  2m39s. Causa: `grav_force()` se llama 2 veces por paso de leapfrog
  (kick-drift-kick) — 100000 llamadas en 50000 pasos, cada una con
  al menos 2 regiones paralelas nuevas (rama `Isochrone` + término
  centrífugo) — y el trabajo real por partícula es mínimo (un puñado
  de `sqrt`/divisiones). El overhead de crear y sincronizar equipos
  de hilos de OpenMP decenas de miles de veces domina por completo
  sobre el trabajo que se paraleliza. Revertido por completo
  (`git checkout -- src/grav_force.f90`); no vale la pena perseguir
  esta variante de la idea — a diferencia de `analysish()` (que se
  llama solo cada `spatial_output` pasos, con mucho trabajo por
  llamada), `grav_force()` se ejecuta demasiado seguido con
  demasiado poco trabajo por llamada para que valga la pena
  paralelizarlo así.
- [x] **Output en ASCII de texto plano, pesado y lento de escribir**
  (`vlasov_fdist.2D` — el hallazgo que motivó esta ronda de mejoras:
  una corrida N≈10⁴ con `spatial_output=100` llegó a 110 MB con solo
  5.5% de progreso, 100% CPU en formateo `ES16.8`, ~5h proyectadas).
  Portado el módulo `hdf5_io.f90` de `VlasovPoisson_PIC_sp`
  (adaptado: sin `l_part`, el dataset de partículas es `f` directo en
  vez de `l_part*f`), agregado el parámetro `output_format`
  (`ascii`/`hdf5`, nuevo campo al final de `input_parameters` —
  extiende la lista ya incompleta, ver el bug de campos faltantes de
  arriba) y el linking de HDF5 al `Makefile` (mismos flags/paths que
  el otro repo, ya verificados en esta máquina:
  `h5fc`/`libhdf5-dev` disponibles). Un archivo `.h5` por corrida, un
  grupo por snapshot, comprimido con gzip.
  Validado: estructura del `.h5` correcta (`h5dump -H`, grupo
  `/grid/r`, grupos `/step_<l>` con los atributos/datasets
  esperados, `force`/`potential` solo si `autointeraction`), y
  valores numéricos verificados idénticos a la salida ASCII
  (`avg_rho` en varios puntos de la malla, comparado con `h5dump -d`
  contra `vlasov_avg_density.rl` de la misma corrida). Benchmark
  (~10072 partículas, 2000 pasos, `spatial_output=100`): **7.0s →
  4.1s (~42% más rápido), 11 MB → 4.3 MB (~61% más chico)** —
  consistente con el 34%/62% visto en `VlasovPoisson_PIC_sp`. `hk1.tl`,
  `hk2.tl` y `vlasov_rhomix.tl` se quedan en ASCII (fuera de alcance,
  chicos), igual que en el otro repo. — commit `feat(io): add optional
  HDF5 output, selected via output_format parameter`

## Autogravedad: dónde se va el tiempo realmente (investigado, no una mejora en sí)

Instrumentado temporalmente `poisson_rk.f90` con `cpu_time()` alrededor
de `avg_density()`, el shooting RK2, y la interpolación de
`pot_part`/`force_part`, en una corrida autogravitante (~10072
partículas, 3000 pasos, 8 hilos). Resultado (CPU-segundos acumulados,
suma entre los 8 hilos):

| Parte | CPU-s | % |
|---|---|---|
| `avg_density()` (depósito) | 35.4s | 56% |
| interpolación a partículas | 27.2s | 43% |
| **shooting RK2 de Poisson** | **0.41s** | **<1%** |

**El propio solver de Poisson no es el cuello de botella** — en
simetría esférica se reduce a una EDO 1D, $O(N_r)$, ya
algorítmicamente óptimo (un método tipo árbol/FMM, pensado para 3D sin
simetría, sería un paso *atrás* acá, no una mejora). El costo real
está en el depósito/interpolación partícula↔malla — que ya tenían el
cell-list y la paralelización de esta sesión. Instrumentación
revertida después de medir (no es un cambio permanente).

- [x] **`poisson_rk.f90`/`density.f90`: `Wn`/`Sn` evaluados dos veces
  por par (partícula, punto de malla) con los mismos argumentos** —
  una vez para `pot_part`/`rho`, otra para `force_part`/`curr`.
  Corregido calculando el valor una sola vez y reusándolo. Validado
  bit a bit idéntico contra la versión anterior. **Medido sin
  ganancia real** (10.60s → 10.69s en el mismo benchmark autogravitante,
  dentro del ruido) — casi seguro `gfortran -O3` ya eliminaba la
  llamada redundante por sí solo (common subexpression elimination,
  dado que `Wn`/`Sn` son funciones simples y sin efectos secundarios).
  Se deja el cambio de todos modos: es correcto, más claro, y no
  depende de que el compilador siga optimizándolo así con otros
  flags. — commit `perf(poisson_rk,density): deduplicate repeated
  Wn/Sn evaluations`

- [ ] **Idea más grande, no implementada: ordenar las partículas por
  posición/celda periódicamente**, para mejorar la localidad de
  caché tanto en el depósito (`avg_density`) como en la interpolación
  (`poisson_rk`) — hoy el depósito ya recorre partículas agrupadas
  por celda (vía `build_cell_list`), pero la interpolación sigue
  iterando en el orden original del arreglo `r_part`, sin esa
  localidad. Reordenar físicamente los arreglos de partículas (no
  solo un índice auxiliar) podría mejorar el uso de caché en ambos
  loops a la vez, pero es un cambio de mayor riesgo/alcance (toca
  cualquier lugar que indexe partículas por su posición original) —
  no evaluado en esta sesión.

## Rendimiento con `OMP_NUM_THREADS`: 8 hilos no es óptimo en esta laptop

Máquina de referencia: Intel i7-4710HQ, **4 núcleos físicos, 8 hilos
lógicos por hyperthreading** (confirmado con
`lscpu -p=CORE,CPU`: cores 0-3, cada uno con dos CPUs lógicas). Medido
con `OMP_PLACES=cores`/`OMP_PROC_BIND=close` (para que a 2/4 hilos cada
uno caiga en un núcleo físico distinto, no en el par de hyperthreads
de un mismo núcleo) en dos benchmarks:

**A — `analysish()`, sin autogravedad** (`state="aa"`, ~10072
partículas, 2000 pasos, `spatial_output=100`):

| Hilos | Tiempo | Speedup | Eficiencia |
|---|---|---|---|
| 1 | 12.94s | 1.00× | 100% |
| 2 | 7.15s | 1.81× | 90.5% |
| 4 | 5.32s | 2.43× | 60.8% |
| 8 | 4.14s | 3.13× | 39.1% |

**B — autogravedad** (`avg_density`+interpolación en `poisson_rk`,
llamado una vez por paso; mismo N, 3000 pasos,
`spatial_output` grande para aislar este camino):

| Hilos | Tiempo | Speedup | Eficiencia |
|---|---|---|---|
| 1 | 7.14s | 1.00× | 100% |
| 2 | 6.26s | 1.14× | 57.0% |
| 4 | 4.37s | 1.63× | 40.8% |
| 8 | **8.93s** | **0.80×** | **10.0%** |

**En el caso autogravitante, 8 hilos es directamente más lento que 1
solo hilo** (no solo "subóptimo") — `poisson_rk()` se llama una vez
por paso, cada llamada abre 2 regiones paralelas
(`avg_density`+interpolación), y con hyperthreading la contención
entre hilos que comparten núcleo físico termina costando más que lo
que aportan los 4 "hilos" extra falsos. El caso sin autogravedad
(`analysish`, llamado solo cada `spatial_output` pasos, con mucho más
trabajo por llamada) sí sigue mejorando hasta 8, pero con rendimientos
muy decrecientes.

**Acción tomada**: agregado el target `make run` al `Makefile`, que
corre `exe/VP_PIC` con `OMP_NUM_THREADS=4` y
`OMP_PLACES=cores`/`OMP_PROC_BIND=close` por default (override con
`make run OMP_THREADS=N`), en vez de dejar que OpenMP use las 8 CPUs
lógicas por defecto. 4 es un default razonable para máquinas de este
tipo (4 núcleos físicos); en una máquina con más núcleos reales el
override es directo. — commit `build: add "make run" defaulting to
OMP_NUM_THREADS=4, document the 8-thread regression`

## Decaimiento de $h_k$: primer intento (jitter Weyl/Kronecker) — no funcionó

Herramienta nueva: `paper_runs/notebooks/hk_exact.ipynb`, calcula
$h_k(t)$ exacto (semi-analítico, sin partículas) para el test de la
Tabla 1 y lo compara contra corridas PIC reales. Con eso medimos por
primera vez, con números y no solo intuición, qué tan lejos está la
simulación del decaimiento verdadero:

- [x] **Convergencia con $N_c$ (sin jitter)**: corridas reales
  $N_c\sim10^3$ y $N_c\sim10^4$ (mismos parámetros de la Tabla 1)
  contra la curva exacta. $h_0$ ya converge sin depender de $N_c$.
  Para $h_k$ ($k>0$) el patrón es "sigue a la curva exacta un tramo,
  se estabiliza en un piso, y **vuelve a subir** ('repunte')" — y ese
  repunte se corre más tarde con más partículas ($N_c\sim10^3$:
  repunta ~t=2500-3000; $N_c\sim10^4$: ~t=5000-7000), la firma
  clásica de recurrencia por rejilla (Birdsall & Langdon), ahora
  vista directamente en vez de solo sospechada. Ver
  `hk_convergence.png`.

- [ ] **Intentado y no mejora: jitter sub-celda (Weyl/Kronecker,
  razón áurea/$\sqrt2-1$) en el estado `"aa"`.** Implementado en
  `vlasov-poisson_PIC/src/initial_data.f90` (mismo offset que el
  experimento uncommitted equivalente en `VlasovPoisson_PIC_sp`) y
  medido con el mismo par de corridas $N_c\sim10^3/10^4$. **Resultado
  mixto y, para $N_c\sim10^4$, claramente peor**: durante un tramo
  largo ($t\sim1500$-$5000$) $h_1$-$h_4$ con jitter quedan casi un
  orden de magnitud *por arriba* de la versión sin jitter, con un
  patrón de "muescas" periódicas — el offset determinístico de Weyl
  parece introducir su propia estructura correlacionada con la
  rejilla $(i,j)$ en vez de romperla limpiamente. Para $N_c\sim10^3$
  el resultado es más parejo (mejor en $h_1$/$h_2$, peor en $h_3$),
  tampoco una mejora clara. Revertido
  (`git checkout -- src/initial_data.f90`, nunca comiteado). Ver
  `hk_jitter_N~1e3.png`, `hk_jitter_N~1e4.png` y la Sec. 7 del
  notebook para el detalle completo.

  **Candidatos para probar después** (no evaluados todavía): una
  secuencia de baja discrepancia genuinamente 2D (Halton/Sobol en
  lugar de dos secuencias de Weyl 1D independientes en $r$ y $p$, que
  quedan acopladas a través del mismo índice `indx`), o jitter con
  amplitud distinta a "una celda completa", o perturbar directamente
  en $(Q_3,J_3)$ en vez de en $(r,p_r)$ (la regularidad problemática
  es la de $J_3$, no la de $r$ — perturbar $r,p_r$ solo la rompe de
  forma indirecta, a través del mapeo no lineal).

## Decaimiento de $h_k$: segundo intento — Halton 2D (mejora modesta, mixta) y rejilla directa en $(Q_3,J_3)$ (peor, revertido)

Los dos candidatos de la sección anterior, probados por separado, cada
uno con el mismo par $N_c\sim10^3/10^4$ (recalibrando `Nrc,Npc` para
que el número de partículas *después* del corte coincida con el
baseline: 988-1009 y 10072-10470 según la variante, contra 1001/10072
del baseline — si no, la comparación queda contaminada por tener más
partículas). Herramienta: `paper_runs/notebooks/hk_exact.ipynb`
Sec. 8.

### Halton 2D (`aa_halton`, mantenido en el código)

Mismo lugar que el jitter anterior — perturbar $(r,p_r)$ — pero con
una secuencia de Halton genuinamente 2D (bases 2 y 3, inversa radical)
en vez de dos Weyl 1D acopladas por el mismo índice `indx`.

Media/σ en $t\in[8000,10000]$, baseline vs. halton:

| | $N_c\sim10^3$ h1 | h2 | h3 | h4 | $N_c\sim10^4$ h1 | h2 | h3 | h4 |
|---|---|---|---|---|---|---|---|---|
| baseline | 1.46e-9 | 6.99e-10 | 7.28e-10 | 8.19e-10 | 5.07e-12 | 4.74e-12 | 3.50e-10 | 4.26e-10 |
| halton | 3.37e-10 | 3.45e-10 | 5.56e-10 | 6.09e-10 | 3.47e-11 | 7.28e-11 | 1.33e-10 | 1.21e-10 |

A $N_c\sim10^3$: **mejora consistente en los 4 modos** (piso 2-4×
más bajo que el baseline). A $N_c\sim10^4$: **mixto** — peor en h1/h2
(el baseline ya había alcanzado un piso muy bajo, ~5e-12, que el
jitter de Halton eleva a ~3-7e-11), pero mejor en h3/h4 (2.6-3.5× más
bajo). $h_0$ no cambia (1.510e-8 en ambos, como se espera — el modo
$k=0$ ya convergía sin depender del método de muestreo).

**Conclusión**: mejora real pero no uniforme, y notablemente mejor
que el intento anterior (Weyl/Kronecker 1D), que a $N_c\sim10^4$ era
peor en *todos* los modos por casi un orden de magnitud. No resuelve
la recurrencia (sigue habiendo un piso de ruido, no decae a cero),
pero la reduce en la mayoría de los casos. Mantenido en el código
(estado `aa_halton`, no reemplaza `aa`) por ser una mejora neta, útil
como opción, aunque no concluyente. No commiteado como reemplazo del
estado por defecto.

### Rejilla directa en $(Q_3,J_3)$ con $Q_3$ aleatorio (`aa_qj`, revertido)

Implementado y depurado (tres bugs encontrados y corregidos en el
camino, ver commit para el detalle completo):

1. **Bug de normalización**: copiaba la línea de normalización de
   `aa` usando `drc*dpc` (paso de rejilla en $(r,p)$) como peso de
   cuadratura, cuando las partículas viven en una rejilla uniforme en
   $(J,Q)$ con paso propio `(dJc,dQc)` — corregido a `dJc*dQc`.
2. **Bug de rango de $J$ (el que de verdad rompía todo)**: el rango
   de $J_3$ para la rejilla se determinaba evaluando $J_r$ en las 4
   esquinas de la caja $(r_{minc}..r_{maxc},p_{minc}..p_{maxc})$ — pero
   esas esquinas son órbitas *no ligadas* para esta combinación
   $L_0=2$/$r_{minc}=1$ (energía positiva: $r$ chico con $L$ grande
   está por encima de la barrera centrífuga), así que $J_r$ salía NaN
   en las 4, y como las comparaciones con NaN son siempre falsas,
   `Jminc`/`Jmaxc` quedaban congelados en sus centinelas
   $\pm10^{30}$ — cascada a `Jgrid~1e30`, $r\sim10^{60}$,
   $p\sim10^{-30}$, y $h_k(t=0)\sim10^{-48}$ en vez de $\sim10^{-8}$
   (exactamente el síntoma reportado: "difiere por muchos órdenes de
   magnitud", idéntico entre la versión con rejilla determinista en
   $Q$ y la primera versión con $Q$ aleatorio, porque en ambas el
   rango de $J$ era la misma basura).
3. **Defecto de diseño, encontrado al corregir el (2)**: barrer la
   caja *completa* para el rango de $J$ tampoco sirve — la caja
   contiene puntos arbitrariamente cerca del borde $E=0$ (escape),
   donde $J_r=1/\sqrt{-2E}-\ldots$ diverge (se midió
   $J_{maxc}\sim10^2$ a partir de una sola celda casi marginal, contra
   un ancho real del soporte de $\sigma_r\sim0.1$ — >99.9% de una
   rejilla uniforme en ese rango cae donde $F$ es
   astronómicamente pequeña). Solución: usar el soporte conocido de
   la propia $F$ ($J_r\in[10^{-4}\sigma_r,\,6\sigma_r]$) en vez de
   barrer la caja.

Con los tres bugs corregidos, $h_k(t=0)$ da el orden de magnitud
correcto ($\sim4\times10^{-8}$ vs. $1.5\times10^{-8}$ de referencia).
Pero el resultado final (media/σ en $t\in[8000,10000]$) es **peor que
el baseline en casi todos los modos, en ambos $N_c$**:

| | $N_c\sim10^3$ h1 | h2 | h3 | h4 | $N_c\sim10^4$ h1 | h2 | h3 | h4 |
|---|---|---|---|---|---|---|---|---|
| baseline | 1.46e-9 | 6.99e-10 | 7.28e-10 | 8.19e-10 | 5.07e-12 | 4.74e-12 | 3.50e-10 | 4.26e-10 |
| $(Q_3,J_3)$ | 3.15e-9 | 2.02e-9 | 2.04e-9 | 1.65e-9 | 5.33e-10 | 3.59e-10 | 6.37e-10 | 4.28e-10 |

Además $h_0$ sale sistemáticamente ~2.87× por arriba del valor de
referencia (4.32-4.34e-8 vs. 1.510e-8) en ambos $N_c$ — un sesgo de
cuadratura consistente, no ruido, de origen no identificado (posible
candidato: el recorte `Jgrid<=0 → 1e-6` cerca del borde $J=0$, donde
$F\propto J^2$ se anula pero la rejilla uniforme en $J$ no capta bien
la caída parabólica). No investigado más a fondo dado que el
resultado ya no es prometedor.

**Conclusión**: la hipótesis original ("la regularidad problemática
es la de $J_3$, perturbar $r,p_r$ solo la rompe indirectamente") no
se confirma en la práctica — mover la regularidad a $J_3$
directamente, incluso con $Q_3$ aleatorio (sin el problema de
cancelación tipo DFT del intento con rejilla uniforme en $Q$), da un
piso de recurrencia *más alto*, no más bajo, que el baseline.
Revertido (`git diff`/eliminado el bloque `aa_qj` de
`src/initial_data.f90`, nunca quedó en `main` ni en ninguna corrida
"oficial" — los tres bugs y la implementación completa quedan en el
historial de commits de esta rama para referencia).

**Estado de los candidatos de la sección anterior**: Halton 2D →
mejora modesta y no uniforme, mantenida (`aa_halton`). Rejilla directa
$(Q_3,J_3)$ → peor, revertida. Ninguna de las dos resuelve la
recurrencia de fondo; sigue pendiente si se quiere profundizar más
(p. ej. combinar Halton con un $N_c$ bastante mayor, o investigar por
qué el piso de $N_c\sim10^4$ del baseline ya es tan bajo en h1/h2
específicamente — posible artefacto de esa corrida particular más que
una propiedad general).

## Decaimiento de $h_k$: tercer intento — Monte Carlo puro (`aa_random`), sin rejilla en absoluto

Idea de fondo distinta a las dos secciones anteriores: en vez de
perturbar una rejilla (Weyl, Halton) o reubicarla en otras variables
($Q_3,J_3$), eliminar la rejilla por completo. Motivación (Birdsall &
Langdon): toda variante basada en rejilla, por más baja discrepancia
que sea la perturbación, conserva *alguna* correlación entre vecinos
en frecuencia $\omega(J)$ — eso es lo que produce el refasamiento
coherente. El muestreo Monte Carlo puro la elimina: el precio es un
piso de ruido estadístico más alto ($\sim1/\sqrt{N}$), pero sin
repunte coherente ni comportamiento errático por modo.

El estado `aa_random` (aceptación-rechazo directo sobre $F(r,p_r)$) ya
existía en el código, sin usar ni verificar. Tenía **dos bugs
acoplados**, encontrados al verificarlo contra el valor exacto de
$h_0$ (que por teoría de mezcla debe ser constante en el tiempo,
$\approx1.4998\times10^{-8}$ — ver Sec. 1-5 de `hk_exact.ipynb`):

1. **Doble conteo**: pesaba cada partícula aceptada por
   `f(i)=w=F(Qr,Jr)` otra vez — pero el muestreo por rechazo ya
   distribuye la densidad de partículas proporcional a $F$, así que
   pesarlas también por $F$ sesga la forma reconstruida hacia las
   regiones ya densas (efectivamente $\sim F^2$, renormalizado) en vez
   de representar $F$. Fix: peso uniforme `f(i)=1.0d0` (partículas de
   igual masa, la representación MC correcta).
2. **Casi elimino `drc*dpc` sin necesidad**: mi primer intento de
   arreglo también quitó el factor `drc*dpc` de la normalización final,
   razonando que partículas ya-aleatorias no necesitan un peso de
   tamaño de celda. Esto dio $h_0(t=0)\sim1.56\times10^{-10}$ — cien
   veces menor que el valor exacto. La razón: `density.f90`,
   `energy.f90` y `analysish.f90` tratan `f()` como "valor de $F$ en
   ese punto" y multiplican por `drc*dpc` de forma **incondicional**,
   sin importar cómo se colocó la partícula — así que `drc*dpc` sí
   debe seguir apareciendo en la normalización de `aa_random`
   (compensándolo), no eliminarse. Con esto restaurado,
   $h_k(t=0)$ calza con la curva exacta al ~1-2% en los 5 modos — la
   mejor concordancia inicial de las cuatro variantes probadas en
   esta sesión (jitter, halton, $(Q_3,J_3)$, random).

**Resultado final** (media/σ en $t\in[8000,10000]$), mismo par
$N_c\sim10^3/10^4$ (`Nrc=Npc=32` → 1024 partículas, `Nrc=Npc=100` →
10000 partículas — `aa_random` no tiene el paso de corte `r0`, así
que `Npart=Nrc*Npc` directamente, sin recalibrar):

| | $N_c\sim10^3$ h1 | h2 | h3 | h4 | $N_c\sim10^4$ h1 | h2 | h3 | h4 |
|---|---|---|---|---|---|---|---|---|
| baseline | 1.46e-9 | 6.99e-10 | 7.28e-10 | 8.19e-10 | 5.07e-12 | 4.74e-12 | 3.50e-10 | 4.26e-10 |
| MC puro | 2.96e-10 | 3.86e-10 | 3.98e-10 | 2.53e-10 | 2.17e-10 | 1.37e-10 | 1.40e-10 | 1.23e-10 |

A $N_c\sim10^3$: **mejora en los 4 modos** (1.8-5× más bajo),
comparable o mejor que Halton en el mismo régimen. A $N_c\sim10^4$:
mismo patrón mixto que Halton — peor en h1/h2 (el baseline ya tenía
ahí un piso inusualmente bajo, ~5e-12, sin recurrir todavía), mejor en
h3/h4 (2.5-3.5× más bajo, el baseline ya había empezado a recurrir).

**Lo distinto de MC puro, y la razón para preferirlo pese a no ganar
en todos los modos**: el piso que da es **plano y predecible** —
todos los modos h1-h4 caen en la misma banda estrecha,
$\sim1.2$-$2.2\times10^{-10}$, en ambos $N_c$. El piso del baseline
(y de Halton) es **errático por modo**, abarcando casi dos órdenes de
magnitud (5e-12 a 4.5e-10) porque la recurrencia por rejilla golpea a
cada modo en un momento distinto — cuando un modo "no ha recurrido
todavía" en una corrida particular, eso es suerte de esa corrida, no
una propiedad confiable del método. El piso de MC, en cambio, es una
consecuencia directa y predecible de $N_c$ (ruido estadístico
$\sim1/\sqrt{N_c}$), sin la lotería de qué modo recurre cuándo.

**Conclusión**: bug corregido y commiteado (afecta a cualquier uso
futuro de `aa_random`, no solo a este experimento). Como alternativa
a la recurrencia, MC puro no elimina el piso de ruido pero sí elimina
el repunte coherente y su dependencia errática del modo — es la
opción más alineada con la teoría clásica de PIC (Birdsall & Langdon)
para este problema específico. No reemplaza `aa` por defecto (el
piso de ruido sigue siendo real, solo que predecible), queda como
estado adicional disponible.

## `output_format="raw"`: binario crudo, alternativa a HDF5

- [x] **Agregado un tercer `output_format="raw"`** (`src/raw_io.f90`),
  motivado por medir directamente el overhead propio de HDF5: escribir
  el mismo volumen de datos vía Fortran `stream`/`unformatted` plano
  en vez de a través de la API de grupos/datasets/atributos de HDF5
  fue **~20× más rápido** en una prueba aislada (2000 registros
  seguidos, sin cómputo de por medio: 0.18s crudo vs 3.7-4.0s HDF5) —
  el overhead de metadatos por objeto de HDF5 domina, no la
  compresión (gzip apenas ayuda acá: el 98% del volumen es
  `r_part`/`p_part`/`f`, ruido de partículas de alta entropía que no
  comprime; solo el 2% que vive en la malla sí comprime bien).

  **Corrección importante, medida después de implementarlo**: esa
  cifra de ~20× es para el caso sintético de guardar *cada paso*. En
  uso real, con `spatial_output=100` (1 de cada 100 pasos), el
  cómputo (`grav_force`+leapfrog+`analysish` periódico) domina sobre
  el I/O casi siempre, y la diferencia total de la corrida completa
  es mucho más chica: **~6% más rápido** que HDF5+gzip con
  `spatial_output=100` (3.74s vs 3.98s, N~10072, 2000 pasos), **~13%**
  con `spatial_output=10` (32.8s vs 37.6s). El archivo crudo además
  sale *más grande* que HDF5+gzip (sin comprimir: 5.1-5.2 MB vs 4.3
  MB en el mismo benchmark) — gzip sí gana en tamaño, aunque no en
  velocidad.

  **Trade-off real**: `raw` no es autodescriptivo — el layout exacto
  de bytes está documentado a mano en el comentario de cabecera de
  `raw_io.f90` (y debe mantenerse sincronizado ahí si cambia), y
  necesita el lector a medida `paper_runs/scripts/rawgraph_io.py` en
  vez de `h5py`/`h5dump` gratis. Vale la pena para corridas con
  `spatial_output` muy chico (guardado muy frecuente) donde el I/O sí
  llega a dominar; para el uso típico (`spatial_output=100`) la
  diferencia con HDF5 es marginal y probablemente no justifica perder
  las herramientas de HDF5.

  Validado: `rawgraph_io.RawRun` leyendo un `.raw` coincide bit a bit
  (`np.allclose`) contra `h5py` leyendo el `.h5` de la *misma* corrida
  (grid, `rho`, `avg_rho`, `r_part`, `p_part`, `f`, tiempo, energía,
  en el primer/décimo/último snapshot). `hygraph.py` extendido para
  detectar el formato por extensión (`.h5` vs `.raw`) y usar el mismo
  visor para ambos. — commit `feat(io): add output_format="raw",
  ~20x faster than HDF5 in isolation but only ~6-13% in realistic
  runs`

## No aplica / ya está bien en este repo

- **Paralelización de `initial_data`**: en `VlasovPoisson_PIC_sp` el
  bug era un índice `indx` incrementado a mano en un loop triple
  `(k,i,j)` con `l_part` como tercera dimensión — dependencia
  secuencial que bloqueaba el `!$OMP`. Acá la malla es 2D (`Lfix`
  escalar, sin dimensión `L`), así que el índice cerrado
  `(i-1)*Npc+j` siempre estuvo disponible y los loops de
  `gaussian1`/`aa`/`Plummer`/`compact`/`compact2` ya están
  correctamente paralelizados con `!$OMP PARALLEL DO`. Nada que
  portar acá.
- **`analysish.f90` bug de índice `mode` vs `i` en `phik`**: no
  aplica — esta versión llama `phik(Jr(j), i, ...)` (con `i`, el modo
  actual de la iteración), no con la constante `mode`. Ese bug se
  introdujo más tarde, durante el refactor a `l_part`, no está
  presente en esta base histórica.

## Bugs nuevos encontrados en este repo (no relacionados con el port, sin corregir)

Encontrados al intentar correr el código para validar los fixes de
arriba. Ninguno existe en `VlasovPoisson_PIC_sp` (arquitectura y/o
formato de parámetros distintos), así que no hay nada que "portar" —
son bugs propios de este repo.

- [x] **`analysish.f90` `hk1`/`hk2`: falta el factor $8\pi^2 L_0$
  de la Ec. 44 del paper.** Confirmado con datos reales (no ya
  hipótesis): en el respaldo externo
  (`.../paper/l=2.0/m=0.0001_c=0.01/Untitled.ipynb`) el propio autor
  grafica `8.0*np.pi**2*2*hk[0]` directo desde `hk1.tl` — es decir,
  aplicaba el factor $8\pi^2 L_0$ ($L_0=2$) a mano en post-proceso
  porque el Fortran no lo hace. Verificado además que la reducción de
  la Ec. 44 vía $\mathcal F=\mathcal F_0(r,p_r)\delta(L-L_0)$ da
  exactamente $8\pi^2L_0\int\mathcal F_0\hat\Phi_k^*e^{-ikQ^3}\,dr\,dp_r$,
  y que el resto del bucle (`Qr`/`Jr` a $L_0$ fijo, `phik`, el peso
  `drc·dpc`) ya calculaba correctamente esa integral doble — solo
  faltaba el prefactor. Agregado `hk1 = 8π²·Lfix·hk1` (e igual para
  `hk2`) al final de cada bucle en `analysish.f90`. Probado: con
  `state="aa"`, $L_0=2$, el $h_0(t=0)$ resultante
  ($1.5076\times10^{-8}$) es exactamente el valor crudo anterior
  ($9.547\times10^{-11}$) × $8\pi^2\times2$. `VlasovPoisson_PIC_sp`
  ya tenía este fix (generalizado con `l_part`); acá se portó a la
  versión `Lfix` escalar. Ver
  `VlasovPoisson_PIC_sp/Vlasov_Poisson_evolutions/h0_normalization_check.md`
  §§9-10 para la derivación completa y la confirmación numérica
  (acuerdo al 0.01% con el valor "Analytical" de la Tabla 1 una vez
  aplicado el factor, más un factor de masa efectiva ×100 — ver
  ítem separado más abajo, todavía sin resolver). — commit
  `fix(analysish): add the missing 8*pi^2*L0 factor to hk1/hk2`
- [ ] **$h_k$ ($k>0$) no decae como en la Tabla 2 del paper — ni en
  mis corridas de reconstrucción ni en los datos originales
  archivados.** Usando los datos reales de
  `.../paper/l=2.0/m=0.0001_c=0.01/vlasov_fdist.2D` (1601 snapshots,
  20112 partículas, la corrida que muy probablemente produjo la
  Fig. 4 del paper), los modos $h_k$ con $k>0$ en $t\in[8000,10000]$
  quedan en $\log_{10}h_k\approx-9.2$ a $-9.8$, contra $-11.3$ a
  $-11.6$ publicado en la Tabla 2 para $N_c\sim10^4$ — 2 órdenes de
  magnitud menos mezclados de lo esperado, y esto es en los datos
  *originales*, no en una reconstrucción con parámetros adivinados.
  Candidato más probable: el artefacto de "recurrencia" por colocar
  cada partícula computacional exactamente en el centro de su celda
  $(r,p)$ — un reticulado casi regular en $J_3$ hace que la mezcla en
  $Q_3$ se revierta coherentemente antes de alcanzar el piso de ruido
  esperado. `VlasovPoisson_PIC_sp` tiene un experimento en curso
  (jitter sub-celda tipo Weyl/Kronecker en `initial_data.f90`, sin
  commitear/validar todavía) apuntado a este mismo problema. **Foco
  actual de trabajo**: investigar/resolver esto para el caso
  $L$ fijo específicamente (este repo), antes de decidir la cuestión
  de normalización/masa de arriba.
- [ ] **Valor de masa efectiva de la Tabla 1/Fig. 4: 0.01, no 0.0001
  como dice el texto.** Confirmado con los datos reales archivados
  (ver nota §10 arriba): con $a_0=0.0001$ (el valor del texto) y el
  factor $8\pi^2 L_0$ aplicado, $h_0$ da $1.5068\times10^{-8}$; con
  $a_0=0.01$ da $1.5068\times10^{-6}$, que coincide al 0.01% con el
  valor "Analytical" publicado ($1.507\times10^{-6}$). Decisión
  editorial pendiente (no técnica): corregir el texto del artículo o
  re-correr las Tablas 1-2 con $a_0=0.0001$ real.

- [ ] **`input_parameters`: desincronizado con `read_initial_param`
  (`utils.f90`), le faltan 6 campos.** `read_initial_param` lee, en
  orden, ... `state`, `j1`, `j2`, `sj1`, `sj2`, `sq1`, `sq2`,
  `bsplineorder`, ... (líneas 39-46) — `j1`/`j2` son los centros en
  $J$ de las dos funciones de prueba $\Phi_1$/$\Phi_2$ que usa
  `analysish.f90` (`phik(Jr(j), i, j1, sq1, sj1)` y análogo con
  `j2,sq2,sj2`), y `sj1`/`sj2`/`sq1`/`sq2` sus anchos en $J$/$Q$. El
  `input_parameters` versionado en el repo **no tiene esas 6
  líneas**: pasa directo de `state` (línea 23, `aa`) a
  `bsplineorder` (línea 24, `1`). Al leerlo con `read(*,*)`, Fortran
  simplemente seguía consumiendo lo que encontraba en las líneas
  siguientes como si fueran `j1..sq2` — `j1` terminaba leyendo el
  valor de `bsplineorder` (`1`, que sí parsea como real), pero `j2`
  intentaba leer la palabra `leapfrog` (línea 25, pensada para
  `integrator`) como número real, y el programa abortaba con
  `Fortran runtime error: Bad real number in item 1 of list input`
  en la línea 41 de `utils.f90`. Es decir: **el `input_parameters`
  de este repo, tal como está commiteado, no corre — ninguna corrida
  documentada pudo haberse hecho con él sin antes agregarle esas 6
  líneas.** Reproducido en esta sesión al intentar un smoke test con
  él directamente. Corregido *solo localmente* (no commiteado, fuera
  del alcance del port) agregando 6 líneas con valores por defecto
  razonables (`0.0`/`0.15`/`0.1`/`0.05`/`0.1`/`0.2`, tomados de los
  dos casos de prueba $\Phi_1$/$\Phi_2$ descritos en
  `VlasovPoisson_PIC_sp/Vlasov_Poisson_evolutions/main.md`) solo
  para poder ejecutar los smoke tests de esta sesión.
- [ ] **`initial_data.f90`, estado `"aa"`: el filtro de corte
  reutiliza `r0` (pensado como "centro de la gaussiana en r" para
  `gaussian1`/`gaussian2`/`compact`/`compact2`) como si fuera la
  fracción de corte del máximo de `f`.** El bloque (líneas ~239-250)
  hace, en esencia, `if (f(...)<=r0*f_max) then r_part(...)=100000
  end if` (la rama previa, `(r0-0.00)*f_max<=f<=r0*f_max`, es código
  muerto: intervalo de medida cero, nunca dispara en la práctica) —
  esto es *matemáticamente correcto* como criterio de corte (`f <=
  fracción·f_max` → descartar) **siempre que `r0` sea una fracción
  pequeña en $[0,1)$**, y de hecho coincide exactamente con el
  parámetro `cutoff` que existe como campo propio, correctamente
  nombrado, en `VlasovPoisson_PIC_sp` (`if (f(indx)<=
  cutoff*f_max)`). El problema es que en este repo no hay un campo
  `cutoff` separado — se reutiliza `r0`, sin renombrar ni
  documentar el cambio de significado para el estado `"aa"` — y el
  `input_parameters` commiteado tiene `r0=3.0` (razonable como
  centro de gaussiana para otros estados, pero no como fracción de
  corte). Con `r0=3.0`, como `f<=f_max` siempre, la condición
  `f<=3·f_max` es **siempre verdadera** — el filtro marca el 100% de
  las partículas para descarte, y `reduce_arrays` las elimina a
  todas: `state="aa"` con el `input_parameters` del repo produce una
  distribución inicial vacía en silencio (sin ningún error o aviso).
  Reproducido en esta sesión. No corregido — evitado para las
  corridas de comparación con el artículo usando explícitamente
  `r0` = el valor real de `cutoff` del artículo (0.01), no el `r0`
  (centro de gaussiana) del archivo de parámetros del artículo (que
  vale 0.0 y no aplica al estado `"aa"` de todos modos).

## Fuera de alcance de este port (decisión pendiente, no arquitectura-independiente o ambigua)

- **`eps` fijo/no fijo**: en este repo `eps = Lfix/(10*pmax)` está
  **activo** (no es el bug — es exactamente el código que confirmó la
  regresión en el otro repo). No tocar.
- **Normalización de `hk1`/`hk2` en `analysish.f90`** (sin factor
  $8\pi^2 L$): ver `VlasovPoisson_PIC_sp/Vlasov_Poisson_evolutions/h0_normalization_check.md`
  — es la misma ambigüedad de interpretación del paper que ese
  documento dejó abierta, no un bug claro. No tocar sin decidir esto
  primero.
- **Optimización de la cuadratura de `phik`** (acá es un Simpson de
  512 puntos recalculado por partícula y por modo, no la versión con
  cuadratura reusada): es una mejora de rendimiento real y análoga a
  la del otro repo, pero no fue parte de lo pedido ("~10 bugs y
  optimizaciones" ya identificadas); queda para una pasada aparte si
  se quiere.
- **Bugs de `grav_force.f90`** (`sphere` sin `abs()`, fondos
  `iso`/`isotrun`/`nfw`/`burkert` sin actualizar
  `pot_part`/`force_part`): siguen **pendientes también en el otro
  repo** (nunca se corrigieron ahí), así que no hay nada que "portar"
  — si se quieren corregir, es trabajo nuevo en ambos repos a la vez.
