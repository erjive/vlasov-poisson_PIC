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

- [x] **`build_cell_list` (`utils.f90`), la única pieza serial que
  quedaba en el paso autogravitante, paralelizada** (era ~15% del
  tiempo de `avg_density()`, medido antes de tocar nada). Es un
  counting sort (cuenta partículas por celda, prefix-sum, reordena) —
  no es un `!$OMP PARALLEL DO` trivial porque el paso de reordenamiento
  final tiene una dependencia de escritura por celda
  (`cursor(c) = cursor(c) + 1`). Implementado como counting sort
  paralelo en dos pasadas (técnica estándar, la misma familia que los
  radix sorts paralelos):
  1. Cada partícula obtiene su índice de celda; cada hilo acumula su
     **propio** histograma por celda (`local_count`) — sin atómicos,
     cada hilo solo toca su propia fila.
  2. Reducción de los histogramas por hilo a `cell_count`/`cell_start`
     (igual que antes) más el prefix-sum *entre hilos* para el offset
     de escritura propio de cada hilo en cada celda
     (`local_offset(c,th)`), partiendo el rango de cada celda en
     sub-rangos disjuntos por hilo.
  3. Segunda pasada (mismo *schedule* que la primera, para que cada
     partícula la procese el mismo hilo que la contó): cada hilo
     coloca sus propias partículas en su propio sub-rango — de nuevo
     sin atómicos, sin choques de escritura entre hilos por
     construcción.

  **El orden del arreglo `local_count`/`local_offset` importó
  muchísimo**: la primera versión, indexada `(hilo,celda)`, dio solo
  ~1.6× de aceleración con 4 hilos — Fortran es *column-major*, así
  que con ese orden las entradas de hilos distintos para la *misma*
  celda quedan a solo `nth` elementos de distancia en memoria, lo
  bastante cerca para caer en la misma línea de caché (false sharing:
  cada hilo solo escribe su propia fila, pero la línea de caché salta
  entre núcleos igual, sin que haya una carrera de datos real).
  Invertido a `(celda,hilo)` — cada hilo pasa a tener un bloque
  contiguo de `Nr` elementos, completamente separado del de los
  demás — la aceleración subió a **3.08×** con 4 hilos.

  Medido con el mismo benchmark de la sección anterior (~10072
  partículas, autogravedad, 4 hilos, 2000 pasos):

  | | antes (serial) | después (paralelo, layout `(celda,hilo)`) |
  |---|---|---|
  | `build_cell_list` (2000 llamadas, wall) | 0.184s | **0.0598s** (3.08×) |
  | `avg_density()` total (2000 llamadas, wall) | 1.263s | 1.126s (10.9%) |
  | corrida completa (wall) | 2.364s | 2.281s (~3.5%) |

  Consistente con el ~15% de `avg_density()` medido al principio: ya
  no queda margen grande ahí, la parte de depósito/interpolación
  (ya paralela) domina el costo restante.

  **Verificado, no solo medido**: comparado contra la versión serial
  original en una corrida autogravitante completa (mismo estado `aa`,
  11 snapshots HDF5) — `avg_rho`, `force`, `potential`, `rho`,
  `r_part`, `p_part`, `f` salen **bit a bit idénticos** (diferencia
  relativa exactamente 0.0 en los 11 snapshots), energía total a nivel
  de redondeo de punto flotante (~1e-15, ruido de orden de suma
  distinto entre hilos, no un error). Repetido dos veces (una por cada
  versión del layout) con el mismo resultado. — commit `perf(utils):
  parallelize build_cell_list (2-pass counting sort), 3.08x`

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

## $h_k$: guardar la fase (complejo, no solo magnitud) + intento de predecir el piso de ruido teórico

Dos tareas relacionadas con "cómo saber si $h_k$ realmente tiende a
cero pese al piso de ruido": (a) `analysish.f90` solo guardaba
$|h_k|$, lo cual hace inútil promediar sobre corridas independientes
para cancelar ruido; (b) intentar predecir el piso de ruido esperado
desde primeros principios, sin correr nada, para tener una referencia
independiente de la corrida real.

- [x] **`analysish.f90` ahora también guarda $h_k$ complejo**
  (`hk1_complex.tl`/`hk2_complex.tl`, columnas `t, Re(h_0),Im(h_0),
  ..., Re(h_4),Im(h_4)`), sin tocar `hk1.tl`/`hk2.tl` (se siguen
  guardando igual, compatibilidad con todo el análisis anterior). La
  razón: el ruido de discretización tiene fase aleatoria entre
  corridas con distinta semilla, pero una señal física real tiene fase
  consistente — promediar $|h_k|$ (lo único que había) sobre varias
  corridas **no cancela el ruido** (es como promediar una magnitud
  Rayleigh: nunca baja de su propia escala, sin importar cuántas
  corridas se promedien). Promediar el $h_k$ **complejo** primero, y
  recién ahí tomar magnitud, sí cancela el ruido de fase aleatoria y
  deja sobrevivir una señal física real. Verificado: a $t=0$,
  $\mathrm{Re}(h_k)\approx|h_k|$ guardado en `hk1.tl` (con
  $\mathrm{Im}(h_k)\sim10^{-17}$, ruido de redondeo), consistente. —
  commit `feat(analysish): save complex h_k (Re/Im), not just
  magnitude`

- [ ] **Intento de predecir el piso de ruido de `aa_random` desde
  primeros principios — validado el método, pero con una discrepancia
  real de ~150-230x contra lo medido, sin explicar del todo.**
  Para partículas de igual peso $m_p=a_0/N$ distribuidas
  $\propto F(Q,J)$, la varianza de un estimador MC estándar da
  $$\mathrm{Var}[h_k] \approx (8\pi^2L_0)^2\frac{a_0^2}{N}\,
  \mathbb E_{q_J}[|\hat\Phi_k(J)|^2],\quad
  q_J(J)=\hat F_0(J)\big/\!\int\hat F_0\,dJ'$$
  con media de Rayleigh $\mathbb E[|h_k|]\approx\sqrt{\pi/4}\sqrt{\mathrm{Var}[h_k]}$.

  Un primer intento de esta fórmula tenía un error de normalización
  (a $q_J(J)$ le faltaba exactamente el factor $8\pi^2L_0$ que sí
  aparece en la propia definición de $h_k$ — la masa total es
  $a_0=8\pi^2L_0\int\!\!\int F\,dQ\,dJ$, no $\int\!\!\int F$ a secas).
  Detectado y corregido comparando contra una simulación Monte Carlo
  idealizada en Python (sorteo por rechazo de $(Q_0,J)\sim F$, streaming
  libre exacto $Q(t)=Q_0+\omega(J)t$, evaluado en la misma ventana
  $t\in[8000,10000]$ medida en la corrida real) — dos construcciones
  MC independientes (fase uniforme al azar; streaming exacto en la
  ventana real) concuerdan entre sí y con la fórmula ya corregida
  (dentro de ~10%), confirmando que el método está bien planteado
  **para lo que modela**.

  Pero comparado contra el piso medido en la corrida PIC real
  (`aa_random`, sección anterior), la predicción sale
  **~150-230x más alta** que lo medido, en los 4 modos y ambos $N_c$.
  Pista encontrada al investigar: la serie temporal real en
  $t\in[8000,10000]$ es **suave y lentamente decreciente** (cambia
  ~0.7% entre muestras consecutivas, $\Delta t=2.5$), no el salto
  brusco muestra-a-muestra que se esperaría de un ensamble ya
  completamente aleatorizado en fase (que es lo que asume la fórmula
  de ruido). Lectura más probable: a $t\sim10^4$ el ensamble de
  `aa_random` **todavía no llegó** al régimen "completamente mezclado
  en fase" — sigue en una relajación más lenta y específica de cómo
  se dispersan las frecuencias $\omega(J)$ entre partículas cercanas
  en $J$, no capturada por una teoría de ruido genérica tipo
  shot-noise. No investigado más a fondo (p. ej. medir explícitamente
  la escala de tiempo de esa relajación, o repetir la comparación a
  $t$ mucho mayor) — detalle completo, con las cuentas y las tres
  comparaciones lado a lado, en la Sec. 10 de `hk_exact.ipynb`. —
  commit `docs: derive and cross-validate the aa_random noise-floor
  prediction, document the unexplained gap against measured data`

- [x] **Prueba directa: ¿sube el piso al doblar el tiempo de la
  corrida? No.** Se corrió `aa_random` a $N_c\sim10^4$ hasta
  $t=20000$ (el doble de antes), usando el nuevo `field_output` para
  guardar el volumen pesado (`r_part`/`p_part`/`f`) 100x más espaciado
  (18MB en vez de ~1.3GB) sin perder resolución en `hk1.tl` (sigue
  cada 100 pasos). Resultado: la media por ventana de 2000 unidades de
  tiempo de $h_1$ ronda 1.2e-10 a 1.8e-10 en las 6 ventanas cubiertas
  — **sin tendencia sistemática**, ni rastro de acercarse al piso
  teórico (~2.2e-8, casi dos órdenes de magnitud arriba). La curva
  sigue igual de lisa que antes (~0.5% de cambio entre muestras
  consecutivas en toda la ventana $t\in[8000,20000]$, prácticamente
  igual al 0.67% medido con la mitad del tiempo) — ninguna señal de
  volverse "ruido blanco erizado" al acercarse al régimen que asume la
  fórmula de la Sec. 10.

  **Lectura**: duplicar el tiempo no acercó nada la curva medida a la
  predicción estadística. Dado que tampoco cambió cualitativamente
  (sigue lisa, sigue en la misma banda de amplitud), la lectura más
  probable ya no es "todavía no llegó, hay que esperar más" sino que
  el supuesto de fondo — que con tiempo suficiente las fases de
  partículas con $J$ en un rango angosto $\Delta J\sim\sigma_J$
  terminan pareciendo completamente al azar — **simplemente no aplica
  bien acá**: con $N_c\sim10^4$ partículas concentradas en un rango
  angosto de $J$, el conjunto discreto de frecuencias $\{\omega(J_p)\}$
  nunca se ve "suficientemente denso y genérico" para que la suma se
  comporte como ruido blanco genuino — en cambio oscila de forma lisa
  y acotada, más parecido a un batido (*beat pattern*) de pocas
  frecuencias dominantes que a ruido estadístico de muchos grados de
  libertad independientes.

  Punto a favor, de todos modos: el piso no crece con el tiempo — la
  simulación es estable, no acumula ruido, se queda oscilando en la
  misma banda ($\sim10^{-10}$) tanto a $t=10000$ como a $t=20000$.
  Sigue pendiente, si se quiere cerrar la pregunta del todo, correr a
  un $t$ sustancialmente mayor o medir directamente el espectro de
  $\{\omega(J_p)\}$ de una corrida real para estimar la escala de
  batido esperada — detalle completo, con la tabla por ventanas y el
  gráfico, en la Sec. 11 de `hk_exact.ipynb`. — commit `docs: run
  aa_random to t=20000, show the measured floor does not rise toward
  the theoretical prediction`

- [x] **Medido el espectro $\{\omega(J_p)\}$ directamente — explica
  por completo por qué la curva sigue lisa.** En vez de seguir
  estirando el tiempo a ciegas, se usaron las $N=10^4$ partículas
  *reales* de `articleN1e4_rand2x` (posiciones en $t=0$, del snapshot
  HDF5 que guardó `field_output`) para calcular $\omega(J_p)$ de cada
  una y medir dos escalas de tiempo distintas:
  - $T_1\sim2\pi/\sigma_\omega$ (dispersión global de
    $\{\omega(J_p)\}$): la escala en que decae la señal física
    coherente. Medido: $T_1\approx1.8\times10^3$ — coincide con que la
    señal ya decae bien antes de $t=10^4$, consistente con todo lo
    visto hasta ahora.
  - $T_2\sim2\pi/\langle\Delta\omega_{vecino}\rangle$ (espaciamiento
    típico entre frecuencias de partículas *vecinas*, ordenadas): la
    escala en la que el carácter discreto de tener "solo" $N$
    frecuencias distintas (no un continuo) debería empezar a notarse
    como grano/ruido. Medido: $T_2\approx3$-$7\times10^6$ — **más de
    150 veces mayor que los $t=2\times10^4$ ya simulados**
    ($t_{simulado}/T_2\approx0.007$).

  Esto responde la pregunta abierta de la sección anterior: la curva
  sigue lisa simplemente porque **ni de cerca** se llegó a la escala
  de tiempo donde el grano de tener $10^4$ frecuencias distintas
  debería manifestarse — no es que el supuesto de fondo esté mal, es
  que $T_2$ es astronómicamente más grande que cualquier tiempo de
  simulación práctico para este $N$. Además $T_2$ escala
  aproximadamente con $N$ (una estimación simple
  $T_2\sim2\pi N/\sigma_\omega\approx1.8\times10^7$ da el mismo orden
  que el valor medido directamente) — confirma por qué subir $N_c$
  nunca mostró señales de "grano" en las secciones anteriores: va en
  la dirección que hace $T_2$ **más grande todavía**, nunca menor.

  Llegar a $t\sim T_2$ por fuerza bruta (~100-300× más tiempo de
  simulación) no es práctico. Si se quisiera verificar esta predicción
  directamente, tendría más sentido *bajar* $N$ deliberadamente (para
  bajar $T_2$ a un rango simulable) y ver si ahí sí aparece el
  comportamiento "ruidoso" esperado — no probado en esta sesión.
  Detalle completo, con histogramas del espectro y del espaciamiento,
  en la Sec. 12 de `hk_exact.ipynb`. — commit `docs: measure the
  omega(J) spectrum directly, explaining the smooth h_k(t) via a beat
  timescale T2 >> simulated t`

## RESUELTO: el piso de $h_k$ era el suavizado `eps` del termino centrifugo (inconsistencia dinamica/analisis)

**Causa raiz encontrada.** `grav_force.f90:131-132` integra el termino
centrifugo **suavizado**:
$0.5L_0^2/(r^2+\epsilon^2)$, fuerza $L_0^2 r/(r^2+\epsilon^2)^2$;
pero `analysish.f90:38` reconstruye $E$, $J_3$ y $Q_3$ con las formulas
**sin suavizar** ($0.5L_0^2/r^2$). Las particulas se mueven en un
Hamiltoniano y se analizan con otro.

Y $\epsilon$ **no es chico**: `utils.f90:141` lo define como
`eps = Lfix/(10*pmax)` donde `pmax` es un **default hardcodeado**
(`parameters.f90:24`, `pmax=2.0`) que **nunca se lee del archivo de
entrada** -- lo que se lee es `pmaxc`, otra variable. Con $L_0=2$ queda
$\epsilon=0.1$, del mismo orden que los anchos de la propia DF
($\sigma_J=\sigma_Q=0.1$). `pmax` existe legitimamente para la condicion
CFL (`dtr = courant*dr/pmax`); el bug es haber derivado de ahi una
longitud de suavizado con significado fisico.

**Evidencia** (diagnostico offline sobre snapshots ya guardados de
`articleN1e4_quad`, sin correr simulaciones nuevas):

| cantidad | deriva relativa $t=0 \to 10^4$ |
|---|---|
| $E$ **con** suavizado (la que integra el codigo) | rms **1.78e-07** |
| $E$ **sin** suavizado (la que usa `analysish`) | rms **4.67e-04** |
| $\delta J$ reconstruido | **6.590e-04** |

y el diagnostico de cuatro rutas sobre el mismo snapshot:

| ruta | $|h_1|$ en $t=10^4$ |
|---|---|
| A reconstruido $(Q_{rec},J_{rec})$ | 5.0078e-12 (reproduce `hk1.tl`) |
| B analitico $(Q_0+\omega t, J_{grid})$ | **5.6315e-16** (= valor exacto) |
| A' $Q_{rec}+J_{grid}$ | 9.6177e-12 |
| A'' $Q_{anal}+J_{rec}$ | 5.0915e-12 |

Con los **mismos pesos** $f_p$, la ruta analitica da el valor exacto:
la colocacion y los pesos siempre estuvieron bien. En $t=0$ la inversion
Newton-Raphson de `aa_quad` es exacta a precision de maquina
($\delta J$ rms $=2.6\times10^{-16}$, $\delta Q$ rms $=8.8\times10^{-15}$),
asi que el error se acumula en la evolucion -- pero de forma
independiente de $\Delta t$ (test `courant` 0.5 vs 0.125: piso identico
a 4-5 cifras), consistente con una inconsistencia de **modelo**, no de
truncamiento. La contribucion del leapfrog se estima en $\sim6\times10^{-5}$
rad, 3600x menor que los 0.227 rad de $\delta Q$ medidos.

**Esto explica de golpe por que TODOS los esquemas de muestreo chocaban
con el mismo piso** ($\sim5\times10^{-12}$, identico entre `aa` y `aa_quad`
al 0.3%, plano en el tiempo): nunca fue un problema de muestreo.

**Fix propuesto (sin aplicar todavia)**: poner $\epsilon=0$ cuando
$L_0\neq0$ -- el bloque ya esta guardado por `if(Lfix /= 0.0d0)`, y con
$L_0\neq0$ la barrera centrifuga impide $r\to0$, asi que el suavizado no
protege de nada ahi; y ademas exponer `eps` como parametro de entrada en
vez de derivarlo del default de CFL. Impacto esperado: el piso deberia
caer de $5\times10^{-12}$ hacia $\sim5.6\times10^{-16}$ (~4 ordenes) y
**mejora todos los estados por igual**, no solo `aa_quad`.

**FIX APLICADO Y VERIFICADO.** `eps` ya no se deriva de `pmax`: se lee del
archivo de entrada (nuevo campo, 0.0 en todos los `paper_runs/input_*`),
con un warning en tiempo de ejecucion si alguien pone `eps/=0` con
`Lfix/=0`. El suavizado sigue disponible para el caso $L_0\to0$, donde no
hay barrera centrifuga y si hace falta.

Resultado con `aa_quad` ($N_Q{=}40\times N_J{=}800$), $\epsilon=0$ vs el
$\epsilon=0.1$ anterior, contra el exacto:

| $t$ | exacto | $\epsilon=0$ | err.rel | $\epsilon=0.1$ | err.rel |
|---|---|---|---|---|---|
| 1000 | 4.301e-10 | 4.301e-10 | **1.3e-06** | 4.402e-10 | 2.4e-02 |
| 2000 | 2.138e-12 | 2.138e-12 | **2.5e-04** | 4.110e-12 | 9.2e-01 |
| 3000 | 2.507e-13 | 2.512e-13 | **1.9e-03** | 4.896e-12 | 1.9e+01 |
| 5000 | 1.844e-14 | 1.873e-14 | **1.6e-02** | 5.000e-12 | 2.7e+02 |

El PIC ahora **sigue la curva exacta hasta $\sim10^{-14}$**; a $t=5000$ el
error relativo pasa de 270x a 1.6%. Grafica: `hk_epsfix.png`.

**Advertencia**: cambia el modelo de fuerza, asi que invalida
cuantitativamente todas las corridas previas de esta rama (todas usaron
$\epsilon=0.1$) -- incluidas las tablas de $h_k$ de las secciones
anteriores. Las conclusiones *cualitativas* sobre muestreo (recurrencia
de rejilla, batido de MC) siguen valiendo, pero los pisos numericos
reportados ahi estaban dominados por este bug, no por el muestreo.

## Piso de $h_k$ con `aa_quad`: NO es el muestreo, NO es el integrador (diagnostico que llevo a la causa raiz)

`aa_quad` (rejilla de cuadratura limpia en $(Q_3,J_3)$, $N_Q{=}40\times N_J{=}800$)
reproduce $h_0..h_4$ exactos **a 8 cifras** en $t=0$ y baja el piso tardio a
$\sim5\times10^{-12}$, **30x** por debajo de Monte Carlo puro ($1.5\times10^{-10}$).
Pero no alcanza el $\sim5.6\times10^{-16}$ que la *misma* rejilla logra offline
con streaming exacto $Q(t)=Q_0+\omega(J)t$. Tres mediciones acotan el culpable:

- **No es el muestreo ni el layout.** El baseline `aa` (rejilla en $(r,p_r)$,
  distribucion de particulas completamente distinta) cae en el **mismo** piso:
  4.992e-12 vs 5.008e-12, 0.3% de diferencia, y plano en el tiempo.
- **No es la cuadratura de la colocacion.** Esa misma rejilla, evaluada offline
  con streaming analitico, da 5.6e-16.
- **No es el integrador temporal.** Corriendo con `courant` 0.5 vs 0.125
  ($\Delta t/4$), el piso es identico a 4-5 cifras:

  | $t$ | courant=0.5 | courant=0.125 | razon |
  |---|---|---|---|
  | 1500 | 9.6884e-12 | 9.6884e-12 | 1.00 |
  | 2000 | 4.1095e-12 | 4.1104e-12 | 1.00 |
  | 3000 | 4.8827e-12 | 4.8836e-12 | 1.00 |

  Error de fase del leapfrog seria $O(\Delta t^2)$: habria bajado ~16x. No se movio.

Por descarte, el piso se origina en la evaluacion por snapshot dentro de
`analysish.f90` (la reconstruccion de $(Q,J)$ desde $(r,p)$, o la cuadratura
de Simpson de $\hat\Phi_k$). **El mecanismo concreto no esta identificado.**
Dos hipotesis revisadas y descartadas por analisis estatico:

- *Cancelacion catastrofica en* `argaux = (s1+s2-2s)/(s2-s1)`: para una
  particula tipica ($E\approx-0.078$, $L_0=2$) sale $s_1\approx5.5$,
  $s_2\approx9.3$, $s_2-s_1\approx3.8$ -- bien condicionado. Solo degenera
  para particulas a $\lesssim10^{-8}$ de sus puntos de retorno, fraccion
  despreciable.
- *Error de la cuadratura de Simpson de* $\hat\Phi_k$ (nquad=512): su error
  relativo es una **constante** independiente de $J$ (la dependencia en $J$
  es un prefactor exacto), asi que solo reescala $h_k$ -- decae con la senal,
  no produce piso.

**Siguiente paso propuesto (sin correr simulaciones):** tomar el snapshot
$(r_p,p_p)$ ya guardado en HDF5 y recomputar $h_1$ offline por dos rutas --
la cadena de formulas de `analysish` en doble precision vs. una ruta de mayor
precision (mpmath) o el angulo analitico $Q_0+\omega t$ conocido de la rejilla
inicial. La diferencia aisla el mecanismo en minutos.

## Con `eps=0`, el limite pasa a ser el leapfrog: $O(\Delta t^2)$ confirmado

Antes del fix de `eps`, el test de Courant era decisivamente negativo
(piso identico a $\Delta t/4$), justamente porque $\epsilon$ dominaba y es
independiente de $\Delta t$. Repetido con $\epsilon=0$ (`aa_quad`,
$N_J{=}800$), el resultado se invierte: error **absoluto** contra el exacto,

| $t$ | c=0.5 | c=0.125 | razon |
|---|---|---|---|
| 1500 | 4.21e-16 | 2.51e-17 | 16.8x |
| 2000 | 5.36e-16 | 3.33e-17 | 16.1x |
| 2500 | 9.57e-16 | 5.99e-17 | 16.0x |
| 3000 | 4.82e-16 | 3.00e-17 | 16.1x |

$\Delta t/4 \Rightarrow$ error$/16$: exactamente $O(\Delta t^2)$, en cuatro
tiempos independientes. Esto **descarta** la hipotesis alternativa de que
ya estuvieramos midiendo el piso de la propia curva de referencia -- la
referencia es buena al menos hasta 3e-17. El limite ahora es el error de
fase del leapfrog.

**Consecuencia practica para las producciones**: a $t=10^4$ el valor
exacto de $h_1$ es $5.6\times10^{-16}$, asi que

| Courant | err. absoluto | err. relativo a $t=10^4$ | costo |
|---|---|---|---|
| 0.5 | ~5e-16 | ~100% (no resuelve) | 1x |
| 0.25 | ~1.3e-16 | ~23% | 2x |
| 0.125 | ~3e-17 | ~5% | 4x |

Con el courant=0.5 usado en toda la sesion **no se puede llegar a
$t=10^4$**: el error igualaria la senal. Las re-corridas usan 0.25 como
compromiso.

**Nota sobre $N_J$ y Nyquist**: $n_{osc}(t)=k\,\Delta\omega\,t/2\pi
\approx k\cdot5.5\times10^{-3}t$, y hace falta $N_J\gtrsim2n_{osc}$. A
$t=10^4$: $k{=}1$ necesita $N_J\gtrsim110$ (holgado con 800), pero
$k{=}4$ necesita $N_J\gtrsim440$ -- o sea $N_J{=}800$ esta **al limite**
para los modos altos; para $k{=}4$ conviene $N_J\approx1600$. Y al reves:
con $N=10^3$ ($N_J{=}40$) el esquema de cuadratura solo sirve hasta
$t\sim1500$ para $k{=}1$ -- no por ruido, sino por aliasing.

## Corrida de produccion a $t=10^4$ con `eps=0`: que limita cada modo

`aa_quad` ($N_Q{=}40\times N_J{=}800$, $\epsilon=0$, courant=0.25) hasta
$t=10^4$, contra el exacto:

| $t$ | $h_1$ exacto | $h_1$ PIC | err.rel | antes ($\epsilon{=}0.1$) |
|---|---|---|---|---|
| 2000 | 2.138e-12 | 2.138e-12 | **6.3e-05** | 9.2e-01 |
| 4000 | 5.725e-14 | 5.701e-14 | **4.2e-03** | 8.7e+01 |
| 6000 | 7.340e-15 | 7.529e-15 | **2.6e-02** | 6.8e+02 |
| 10000 | 5.631e-16 | 4.613e-16 | **1.8e-01** | 8.9e+03 |

(el 18% a $t=10^4$ calza con los ~23% predichos desde el test de Courant).

**Por modo a $t=10^4$** -- y aca hay que corregir una atribucion previa:

| modo | exacto | PIC | err.rel | $n_{osc}$ vs $N_J/2$ |
|---|---|---|---|---|
| $k$=1 | 5.631e-16 | 4.613e-16 | 1.8e-01 | 55 vs 400 |
| $k$=2 | 1.648e-17 | 1.866e-17 | 1.3e-01 | 110 vs 400 |
| $k$=3 | 1.971e-18 | 1.969e-16 | 99x | 165 vs 400 |
| $k$=4 | 4.303e-19 | 2.782e-16 | 650x | 220 vs 400 |

En la seccion anterior se atribuyo el problema de los modos altos a
Nyquist ("$N_J{=}800$ esta al limite para $k{=}4$"). **Eso era
incorrecto**: $k{=}3,4$ cumplen Nyquist con holgura (165 y 220 contra
$N_J/2=400$) y aun asi fallan. Discriminado con los datos de Courant ya
existentes a $t=3000$, mismo $N_J$: el error absoluto de los modos altos
escala como $\Delta t^2$ (razones 16.1x, 16.1x, 13.2x, 11.2x para
$k=1..4$), o sea **es el leapfrog, no aliasing**.

Cuadro correcto: el error de fase del leapfrog impone un **piso
absoluto** que crece con $t$ y varia algo por modo; a $t=10^4$ con
courant=0.25 vale $\sim2$-$3\times10^{-16}$. Como $h_3^{exacto}=2\times10^{-18}$
y $h_4^{exacto}=4\times10^{-19}$ estan por debajo de ese piso, esos modos
son irresolubles en esa corrida -- no por falta de particulas sino por
$\Delta t$. Para resolverlos haria falta courant mas chico (el piso baja
como $\Delta t^2$) o un integrador de orden mayor.

## Resumen: los tres limites de $h_k$, separados y medidos

Con $\epsilon=0$ y courant=0.25, a $t=10^4$, cada esquema queda limitado
por un mecanismo **distinto** -- y cada uno se midio por separado:

| Esquema | piso a $t=10^4$ | lo limita |
|---|---|---|
| `aa_random` (MC puro) | $\sim10^{-10}$ | ruido de muestreo $1/\sqrt{N}$ |
| baseline `aa` (rejilla en $r,p$) | $\sim2.2\times10^{-13}$ | el *layout*: espaciado irregular en $J$ |
| `aa_quad` (cuadratura en $Q,J$) | $\sim4.6\times10^{-16}$ | fase del leapfrog, $O(\Delta t^2)$ |

Baseline `aa` con $\epsilon=0$ contra el exacto:

| $t$ | $h_1$ exacto | `aa` $\epsilon{=}0$ | err.rel | `aa_quad` | err.rel |
|---|---|---|---|---|---|
| 2000 | 2.138e-12 | 2.039e-12 | 4.6e-02 | 2.138e-12 | 6.3e-05 |
| 6000 | 7.340e-15 | 2.188e-13 | 2.9e+01 | 7.529e-15 | 2.6e-02 |
| 10000 | 5.631e-16 | 2.475e-13 | 4.4e+02 | 4.613e-16 | 1.8e-01 |

El fix de $\epsilon$ mejoro el baseline 20x (piso 5e-12 -> 2.2e-13), pero
sigue 480x por encima de la cuadratura. O sea: **los dos hallazgos
cuentan**, el bug de $\epsilon$ y el layout. Mientras $\epsilon$ estuvo
presente, su piso de 5e-12 tapaba por completo la diferencia entre
esquemas -- por eso durante toda la sesion todos parecian chocar contra
el mismo muro.

## Integradores nuevos: `yoshida4` (simplectico de orden 4) y `analytic` (exacto)

Con el leapfrog convertido en el limite dominante de `aa_quad`, se
agregaron dos integradores (`main.f90`):

- **`yoshida4`**: composicion de Yoshida (1990) de tres pasos leapfrog
  con coeficientes $w_1,w_0,w_1$, $w_1=1/(2-2^{1/3})$,
  $w_0=-2^{1/3}/(2-2^{1/3})$. El sub-paso central va **hacia atras** en
  el tiempo ($w_0<0$), que es lo que cancela el termino $O(\Delta t^2)$.
  Error de fase $O(\Delta t^2)\to O(\Delta t^4)$ al costo de 3
  evaluaciones de fuerza por paso (3 resolvedores de Poisson por paso si
  `autointeraction=.true.`, donde ese camino ya domina). Sigue siendo
  simplectico: no reintroduce deriva secular.

- **`analytic`**: avance exacto en forma cerrada. Sin autointeraccion y a
  $L$ fijo el movimiento radial es integrable, $J_3$ se conserva exacto y
  $Q_3(t)=Q_3(0)+\omega(J_3)t$, asi que se avanza al tiempo **absoluto**
  $t$ (no incrementalmente) y **no hay error de integracion de ningun
  tipo**. Pensado como camino de validacion: aisla todo lo que viene
  despues (cuadratura, `analysish`, normalizaciones) de cualquier error
  del integrador. Su costo no depende de $\Delta t$.
  Requiere setup integrable (`autointeraction=.false.`, `forcetype="bg"`,
  `BGtype="Isochrone"`, `Lfix/=0`, `eps=0`) y es incompatible con
  `reduceparticles=.true.` -- ambas cosas se verifican y abortan con
  mensaje claro. La inversion $(Q_3,J_3)\to(r,p_r)$ se factorizo en
  `utils.f90:invert_QJ_to_rp`, compartida con el estado `aa_quad`.

**`rk4` sigue deliberadamente sin implementar**: RK4 no es simplectico,
reintroduciria la deriva secular en energia y $J_3$ que justamente
acabamos de eliminar. El mensaje de aborto ahora lo explica y remite a
`yoshida4`.

Verificado (smoke test, 400 pasos): ambos reproducen los $h_k(t{=}0)$
exactos, y a $t=10$ `yoshida4` y `analytic` coinciden **a las 9 cifras
impresas** -- como debe ser cuando el error de Yoshida es despreciable.
Benchmarks de convergencia y costo, pendientes.

## Benchmark de `yoshida4`: domina al leapfrog en los dos ejes

Barrido con `aa_quad` ($N_Q{=}40\times N_J{=}800$, $\epsilon=0$, $t=3000$),
mismo $\Delta t$ para ambos integradores. Error absoluto en $h_1$ contra
el exacto ($2.506924\times10^{-13}$):

| Integrador | courant | $\Delta t$ | error abs | razon al duplicar $\Delta t$ |
|---|---|---|---|---|
| leapfrog | 1.0 | 0.050 | 1.946e-15 | -- |
| leapfrog | 2.0 | 0.100 | 8.033e-15 | 4.13x |
| leapfrog | 4.0 | 0.200 | 3.568e-14 | 4.44x |
| yoshida4 | 1.0 | 0.050 | **1.150e-19** | -- |
| yoshida4 | 2.0 | 0.100 | 1.436e-18 | 12.5x |
| yoshida4 | 4.0 | 0.200 | 2.255e-17 | 15.7x |

**Orden efectivo medido**: leapfrog $p=2.05$ y $2.15$; yoshida4 $p=3.64$
y $3.97$. Segundo y cuarto orden confirmados.

Al **mismo** $\Delta t$, yoshida4 es ~17000x mas preciso. Y la comparacion
practica, con tiempos re-medidos **limpios** (sin nada mas corriendo):

| | $\Delta t$ | tiempo | error abs |
|---|---|---|---|
| leapfrog | 0.05 | 333.2 s | 1.946e-15 |
| yoshida4 | 0.20 | **106.2 s** | **2.255e-17** |

**3.1x mas rapido y 86x mas preciso** -- dominacion estricta, sin
compromiso. La razon: el multiplicador de 3x evaluaciones de fuerza de
Yoshida cae sobre una fraccion chica del costo total, porque **`analysish`
domina** (se llama cada `spatial_output` pasos y cada llamada cuesta
$O(N_{part}\times n_{quad}\times n_{modos})$). Al mismo $\Delta t$ yoshida4
solo costo 1.1-1.24x mas que leapfrog, no 3x.

**Recomendacion**: usar `yoshida4` por defecto para este tipo de corrida.
Nota: con `autointeraction=.true.` el balance cambia -- ahi las 3
evaluaciones son 3 resolvedores de Poisson por paso, y ese camino si
domina; habria que re-medir en ese regimen antes de generalizar.

Graficas: `bench_yoshida.png` (convergencia y precision por costo),
`hk_yoshida_modos.png` ($h_k(t)$ de los 4 modos contra el exacto) y
`hk_yoshida_errores.png` (error relativo por modo, yoshida4 vs leapfrog
al mismo $\Delta t$), todas en `paper_runs/notebooks/`.

**Los 4 modos a $t=3000$, mismo $\Delta t=0.05$:**

| modo | exacto | yoshida4 | err.rel | leapfrog | err.rel |
|---|---|---|---|---|---|
| $h_1$ | 2.5069e-13 | 2.5069e-13 | **4.6e-07** | 2.5264e-13 | 7.8e-03 |
| $h_2$ | 6.9088e-15 | 6.9087e-15 | **1.3e-05** | 7.3761e-15 | 6.8e-02 |
| $h_3$ | 8.1284e-16 | 8.1278e-16 | **6.7e-05** | 2.2333e-15 | 1.7e+00 |
| $h_4$ | 1.6679e-16 | 1.6678e-16 | **6.7e-05** | 4.0506e-15 | 2.3e+01 |

El contraste crece con el modo, como corresponde: el error de fase entra
como $e^{-ikQ}$, asi que pesa $\propto k$. En $h_4$ el leapfrog se equivoca
por un factor 23 mientras yoshida4 acierta a 7 cifras. Es decir: **los
modos altos, que con leapfrog eran irresolubles, con yoshida4 si se
resuelven** -- y sin pagar mas tiempo de computo.

## `analysish` factorizado: la simulacion completa 10.4x mas rapida

El integrando de `phik` **se factoriza exactamente**:

$$g(Q,J)=\underbrace{e^{-\sin^2(Q/2)/\sigma_q^2}}_{A(Q)}\cdot\underbrace{e^{-(J-J_0)^2/\sigma_j^2}J^2}_{B(J)}$$

asi que la cuadratura entera en $Q$ se separa de la particula:

$$\phi_k(J,i)=B(J)\cdot\underbrace{\Big[\texttt{quadnorm}\sum_k w_k A(Q_k)\cos(iQ_k)\Big]}_{C(i),\ \textbf{constante de toda la corrida}}$$

$C(i)$ **no depende de la particula**, pero se recalculaba dentro del
bucle de particulas: cada una pagaba una cuadratura completa de
$(n_{quad}+1)\times(\text{modos}+1)$ -- $2\times513$ `exp`, $2\times513$
`sin` y $2\times5\times513$ `cos`/mult **cada una**. Sacado fuera, cada
particula cuesta ahora 2 `exp` y $2\times5$ multiplicaciones.

**Medido** (mismo benchmark, `aa_quad` $N_c=3.2\times10^4$, $N_t=60000$,
4 hilos, cronometrado limpio):

| | tiempo |
|---|---|
| antes | 333.2 s |
| despues | **32.0 s** |

**10.4x mas rapida la corrida completa** -- o sea que `analysish`
representaba ~90% del costo total. Es la misma clase de redundancia que
la duplicacion de `Wn`/`Sn` documentada arriba, pero aquella midio ~1x y
esta 10x sobre el total.

**Fisica identica**: $h_0$ y $h_1$ salen bit a bit iguales; $h_2$-$h_4$
difieren a lo sumo $10^{-22}$ en valor **absoluto**, que es el ultimo
digito del formato `ES16.8` con que se escribe el archivo. Las
diferencias *relativas* de hasta 7e-9 caen todas en minimos locales de
las curvas (p. ej. $h_3$ en $t=2585$ vale 1.42e-15 con vecinos de
3.08e-15), donde el ultimo digito se amplifica.

**Consecuencia para el benchmark de Yoshida**: aquella medicion ("solo
1.1-1.24x mas caro que leapfrog") se hizo cuando `analysish` dominaba el
costo y el empuje de particulas era despreciable. Ahora que `analysish`
es 10x mas barato, el empuje pesa mucho mas, asi que el 3x de
evaluaciones de fuerza de Yoshida **va a notarse mas**. Estimacion
gruesa: su costo relativo subiria de ~1.24x a ~1.9x. Sigue valiendo
ampliamente la pena (17000x de precision al mismo $\Delta t$), pero la
frase "practicamente gratis" ya no aplica y habria que re-medirlo.

**Y sobre calcular $h_k$ en post-proceso**: con esto la pregunta queda
casi sin objeto. El post-proceso nunca podia ser mas rapido en total
(la aritmetica es la misma, mas escribir y releer ~4 GB por corrida);
ahora ademas el costo en linea es marginal. La recomendacion es $h_k$ en
linea + snapshots ralos vía `field_output` para diagnosticos posteriores
-- que es justo lo que permitio encontrar el bug de `eps`.

## Orden de los integradores, medido: `yoshida4` confirmado, `yoshida6` parcialmente

Agregado `yoshida6` (Yoshida 1990 "Solucion A", 7 etapas, coeficientes
$w_3,w_2,w_1,w_0,w_1,w_2,w_3$ con $w_0=1-2(w_1+w_2+w_3)$) via una
composicion **manejada por tabla**, de modo que agregar ordenes mas
altos sea solo agregar coeficientes.

**Metodologia -- dos intentos fallidos antes de uno valido**, vale la
pena registrarlos porque cada uno falla de una forma distinta:

1. Comparando contra `hk_exact` con $\Delta t$ chico: los errores salieron
   ~2.7e-20 y **planos** (orden medido 0.00). No se estaba midiendo el
   integrador sino un piso: a ese nivel dominan la cuadratura en $J$ del
   PIC y la de la propia referencia semi-analitica.
2. Subiendo $\Delta t$ a 0.8-3.2: los errores (1e-14 a 7e-14) llegaron al
   4-30% del propio $|h_1|$ -- **fuera del regimen asintotico**, ordenes
   medidos sin sentido (5.86, luego 1.19).
3. Valido: referencia = corrida con `integrator="analytic"` (comparte
   particulas, rejilla, cuadratura y `analysish`, asi que la diferencia
   es **solo** el integrador), ventana $\Delta t=0.05$-$0.4$, y salida
   subida de `ES16.8` a `ES24.16` (resolucion 1e-29 en vez de 1e-21).

**Resultado** (error absoluto en $h_1$ a $t=3000$ contra `analytic`):

| Integrador | $\Delta t$ | error abs | razon | orden $p$ |
|---|---|---|---|---|
| leapfrog | 0.05 / 0.10 / 0.20 | 1.95e-15 / 8.03e-15 / 3.57e-14 | 4.13x, 4.44x | **2.05, 2.15** |
| yoshida4 | 0.05 / 0.10 / 0.20 / 0.40 | 8.70e-20 / 1.41e-18 / 2.25e-17 / 3.60e-16 | 16.2x, 16.0x, 16.0x | **4.02, 4.00, 4.00** |
| yoshida6 | 0.40 / 0.20 | 7.31e-20 / 6.07e-22 | 120x | **6.91** |
| yoshida6 | 0.20 / 0.10 | 6.07e-22 / 1.27e-21 | 0.48x | -1.06 (piso) |

`yoshida4` queda **confirmado sin ambiguedad**: tres razones consecutivas
de 16.00 sobre tres decadas de $\Delta t$, lo que valida de paso que sus
coeficientes estan bien.

`yoshida6` es evidencia **mas debil**: solo hay **un** intervalo utilizable
(0.4->0.2), $p=6.91$. Debajo de $\Delta t=0.2$ su error ya esta en ~1e-21 y
choca contra un piso **no identificado** -- no es el formato de salida
(resuelve 1e-29), ni el redondeo de la suma (~2e-26 estimado), ni la
tolerancia del Newton-Raphson de `advance_analytic` (~8e-25 estimado).
Queda abierto. Se puede afirmar que es de orden alto ($\geq$6, netamente
mejor que 4), pero no con la solidez del 4.00,4.00,4.00 de `yoshida4`.

**Conclusion practica: no hace falta sexto orden para este problema.**
`yoshida4` a $\Delta t=0.05$ da error 8.7e-20, muy por debajo de
$h_4^{exacto}=4.3\times10^{-19}$ a $t=10^4$: los modos altos ya se
resuelven con cuarto orden. El sexto cuesta 7 evaluaciones de fuerza por
paso en vez de 3 y resuelve un problema que no tenemos. Coincide con la
estimacion a priori: para reducir el error un factor $R$, orden $p$
cuesta $\sim n_{etapas}R^{1/p}$, asi que el 4to gana hasta $R\sim10^4$ y
el 6to recien paga pasado $R\sim10^6$.

Queda implementado igual, por si hiciera falta, y la tabla hace trivial
agregar el 8vo orden (15 etapas) mas adelante.

## Barridos seriales eliminados: la corrida pasa de 333 s a 14 s (23x en total)

Los barridos seriales sobre `Npart` que estaban diagnosticados pero sin
corregir (ver la seccion de autogravedad, items A y B) finalmente se
arreglaron. Se volvieron **el** cuello de botella despues de factorizar
`analysish`: al desaparecer el 90% del costo (que era paralelo), la
fraccion paralela medida se derrumbo de $p=0.78$ a $p=0.11$ -- techo de
Amdahl 1.1x, o sea el codigo practicamente habia dejado de escalar.

**Que se hizo**, en pasos medidos por separado:

1. `main.f90`: leapfrog y Yoshida fusionados de **6 barridos seriales a 2
   pasadas paralelas**, eliminando de paso las copias `r_part_p`/
   `p_part_p` (el *drift* se hace in-place, los valores viejos nunca
   hacian falta). El bucle de simetria en el origen: la condicion
   `rmin == 0` es constante de corrida, asi que se evalua una vez en vez
   de una por particula por paso, y si `rmin > 0` la pasada se saltea.
2. `grav_force.f90`: el bloque isocrono y el centrifugo fusionados en una
   pasada cada uno, con `sqrt(1+r^2)` evaluada **una** vez en lugar de
   2-3, y el denominador centrifugo formado una vez en vez de dos.
3. Recien **despues** de fusionar, se agregaron las directivas OMP.

**Ese orden fue deliberado y resulto ser la clave.** Un intento anterior
de paralelizar `grav_force` habia dado una **regresion de 3x**
(documentada arriba): se habian puesto regiones paralelas alrededor de
operaciones de arreglo diminutas, y el costo de lanzar hilos dominaba.
Con los bucles ya fusionados, cada region hace trabajo sustancial y el
overhead se amortiza -- esta vez no hubo regresion sino 2.11x.

**Medido** (mismo benchmark de siempre: `aa_quad`, $N_c=3.2\times10^4$,
$N_t=60000$, cronometrado limpio):

| etapa | 1 hilo | 4 hilos |
|---|---|---|
| inicio del dia | -- | 333.2 s |
| tras factorizar `analysish` | -- | 32.0 s |
| tras fusionar el push | 34.2 s | 31.3 s |
| tras fusionar `grav_force` (sin OMP) | 30.2 s | 25.1 s |
| **tras agregar OMP** | 30.3 s | **14.3 s** |

**23x mas rapida la corrida completa** respecto del inicio del dia, y la
fraccion paralela se recupero de $p=0.11$ a $p\approx0.70$ (techo 3.3x).

**Fisica intacta**: el error contra la referencia exacta es 1.9460e-15
antes y despues (identico a 5 cifras); la diferencia entre versiones es
3.8e-27, y el maximo sobre todos los modos y tiempos es 6.6e-24 -- doce
ordenes por debajo del error propio del integrador.

**Nota sobre hilos**: con 8 hilos el codigo seguia siendo mas lento que
con 1 antes de este cambio (47.9 s vs 34.2 s) por hyperthreading. El
default de `make run` (4 hilos) sigue siendo el correcto. Y para la
serie de convergencia conviene igual **paralelismo a nivel de trabajos**
(4 corridas concurrentes de 1 hilo) antes que una sola corrida a 4
hilos: 4x de throughput contra 2.1x de latencia.

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
