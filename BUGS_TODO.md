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

- [ ] **`density()`/`avg_density()`/`poisson_rk()`: búsqueda de
  vecinos por fuerza bruta `O(Nr×Npart)`.** Portar el cell-list de
  `utils.f90` (`build_cell_list`, `O(Nr+Npart)`).

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
