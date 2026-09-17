# Preguntas abiertas

Lo que queda por entender o por hacer después del estudio de phase mixing con
autogravedad. Cada entrada dice qué se sabe (con los números medidos), qué falta y
cuál sería el siguiente paso concreto. Las conclusiones establecidas están en
`docs/introduccion/vlasov_intro.tex` y en `paper_runs/notebooks/selfgrav_mixing.ipynb`;
los bugs, en `BUGS_TODO.md`.

Convención: corridas en `exe/sg/<nombre>`, configuración base
`paper_runs/base_selfgrav.par` (gaussiana, `a0=1e-3`, cuadratura `Nrc=400 Npc=25`,
`yoshida4`, `dt=0.1`).

---

## 1. La componente de h_1 que crece con a0=1e-2

**Qué se sabe.** Además de la parte estática S (efecto de coordenadas, ver
*Resueltas*), h_1 tiene una componente que gira a omega = 0.055-0.062:

| a0 | t in [2000,15000] | t in [15000,20000] |
|---|---|---|
| 1e-4 | 1.4% de S | 1.3% de S |
| 1e-3 | 8.7% de S | 10.2% de S |
| 1e-2 | 18% de S | **86% de S** (rodea el origen: saltos de fase) |

Con a0=1e-3 quedó establecido que, después de t=2000, es **ruido de discreción**
(ver *Resueltas*). Con a0=1e-2 no se ha comprobado, y ahí la componente **crece** en
el tiempo, lo que un piso de ruido fijo no explica por sí solo: podría ser ruido
amplificado por la autogravedad más fuerte, o una inestabilidad genuina.

**Siguiente paso.** Repetir la prueba de resolución con a0=1e-2: corrida hasta
t=20000 con Nrc=800 (~22 min) y comparar con `long_a0_1e-2` usando
`paper_runs/scripts/giro_resolucion.py`. Si baja y se descorrelaciona, es numérica;
si no cambia, es física y merece un análisis de estabilidad.

---

## 2. Un montaje limpio para amortiguamiento de Landau

**Problema.** El estado inicial es equilibrio del isócrono *solo*; al activar la
autogravedad toda la componente es a la vez "perturbación" y fuente del potencial.
Eso mezcla el reajuste del equilibrio (el cambio de h_0 del 1.6% en t < 600, las
amplitudes de fila que cambian entre 0.83 y 1.005) con la respuesta a la perturbación.

**Siguiente paso.** Construir un equilibrio autoconsistente F_eq(J) del potencial
total (iterando Poisson con el mapa numérico, `paper_runs/scripts/aa_numerico.py`) y agregar encima una
perturbación pequeña y separada. Solo entonces tiene sentido medir omega_r y gamma.

---

## 3. Régimen a0 >= 1e-2

Con a0=1e-2 el cambio de h_0 es 14% y la componente que gira crece hasta dominar
(pregunta 1).
El mapa analítico ya no sirve; hay que usar el numérico. Depende de 1 y 2. Para perturbaciones grandes,
relajación violenta (Lynden-Bell 1967).

---

## 4. Modos discretos

Resolver el problema lineal de autovalores (método matricial de Kalnajs) para este
equilibrio y comparar omega_r, gamma con las corridas. Un modo dentro de la banda
orbital es resonante y se amortigua por Landau; fuera, no.

---

## 5. Levantar L fijo

Con dispersión en L hay dos frecuencias y resonancias entre ellas. Existe una rama
con `l_part` en `VlasovPoisson_PIC_sp`. Es otro proyecto.

---

## Resueltas

### La parte estática de la meseta es un efecto de coordenadas (confirmado)

Antes era una interpretación consistente con todo (insensible a la discretización,
lineal en a0, fase congelada, J isócrona oscilando a la frecuencia orbital), pero no
probada. **Prueba directa** (`paper_runs/scripts/aa_meseta.py`, sobre
`sg/long20k_snap`): h_1 recalculado con el mapa ángulo-acción numérico del potencial
real (isócrono + potencial propio promediado en t >= 2000; el mapa,
`aa_numerico.py`, reproduce el analítico del isócrono a 1e-15 en J y 1e-12 en Q).

| mapa | parte estática, t in [2000,15000] | parte que gira (rms) | oscilación de J (t in [1500,2000]) |
|---|---|---|---|
| isócrono | 1.767e-3 | 1.540e-4 | 0.435% |
| numérico, potencial promedio | **7.9e-7** | 1.541e-4 | **0.001%** |
| numérico, potencial instantáneo | 1.4e-6 | 1.539e-4 | |

La parte estática baja 2200 veces y la acción de cada partícula queda constante. La
diferencia entre los dos h_1 es un desplazamiento fijo de 1.77e-3 con fase 0,
presente desde t=0. El cambio de h_0 durante la mezcla inicial baja de 1.60% a 1.35%
(potencial promedio) o 0.76% (instantáneo): la mayor parte es un cambio real de las
acciones mientras el potencial propio se reajusta (cambia 21% entre t=0 y t=2000).

### Con a0=1e-3, la componente que gira es física hasta t=2000 y ruido de discreción después

La componente no cambia con el mapa (correlación 0.99988), así que no es un efecto de
coordenadas. Pruebas:

- **Sin autogravedad** (`sg/long20k_nosg`), la misma cuadratura cancela h_1 hasta
  ~2e-10 en t > 12000: la rejilla sola no tiene piso.
- **Barridos existentes, t in [1600,2000]**: 5.675e-5 con N=1e4 y 5.680e-5 con N=1e5;
  igual con Npc, dr y orden del spline. Hasta t=2000 es el continuo resuelto.
- **Streaming libre** (`rotadores.py`): en variables verdaderas las partículas tienen J
  constante a 6e-4 pero la fase se desvía ~0.01 rad de una recta, y h_1 en t=2000 es
  una cancelación de 6e4 (sum|G| = 0.79, |sum G| = 1.4e-5). La reconstrucción no es
  concluyente (correlación 0.60-0.90): basta una modulación de fase coherente de
  ~2e-4 rad para producir la componente.
- **Prueba de resolución hasta t=20000** (`giro_resolucion.py`):

| ventana | Nrc=800 (doble en J) | Npc=50 (doble en Q) |
|---|---|---|
| [1600,2000] | 1.001x, correlación 1.000 | 1.000x, 1.000 |
| [2000,6000] | 0.31x, 0.76 | 0.65x, 0.61 |
| [6000,12000] | 0.53x, 0.84 | 0.35x, 0.61 |
| [12000,20000] | 0.58x, 0.26 | 0.55x, 0.16 |

Al duplicar las partículas en cualquiera de las dos direcciones la componente tardía
baja a 0.3-0.65x y deja de parecerse a la de la corrida base: depende de la
discretización, luego es numérica. La parte estática vale 1.767e-3 en las tres.

## Pendientes técnicos menores

- `main.f90`: `call system('mkdir -p')` y `cp` impiden correr en Windows nativo.
- `yoshida6`: un piso de ~1e-21 impide medir su orden limpiamente.
- La rama `clean/comentarios` no está integrada a `main`.
