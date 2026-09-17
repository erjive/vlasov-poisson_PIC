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

## 1. Mecanismo del ruido de discreción que crece en el tiempo

**Qué se sabe.** La componente de h_1 que gira a frecuencia orbital es, después de
t=2000, ruido de discreción tanto con a0=1e-3 como con a0=1e-2 (ver *Resueltas*).
Pero crece con el tiempo:

| a0 | Nrc | rms por ventana de 2000, de t=2000 a 20000 | ley |
|---|---|---|---|
| 1e-3 | 400 | 3.4e-5 ... 1.9e-4 (se satura hacia t~10000) | ~t^0.9 |
| 1e-3 | 800 | 1.4e-5 ... 1.4e-4 | ~t^0.9 |
| 1e-2 | 400 | 2.7e-4 ... 1.2e-2 | ~t^1.9 (R2 0.96; exponencial 0.93) |
| 1e-2 | 800 | 1.0e-4 ... 3.4e-3 | ~t^1.9 (R2 0.99; exponencial 0.92) |

Con a0=1e-2, duplicar Nrc multiplica la amplitud por ~0.35 en todas las ventanas
(0.23-0.50), cerca de 2^{-3/2}, y deja igual la ley de crecimiento.

**Interpretación plausible, no demostrada.** El potencial propio calculado con N
partículas finitas fluctúa; esas fluctuaciones, proporcionales a a0, desvían las
fases orbitales de cada partícula, y el desfase acumulado crece con t hasta
saturarse. Encaja con: amplitud dependiente de la resolución, crecimiento algebraico
y no exponencial, y crecimiento más rápido con más masa. No encaja con una
inestabilidad física, que crecería exponencialmente con una tasa propia.

**Siguiente paso, si importa.** Medir el escalamiento completo: Nrc=1600 con a0=1e-2
(~45 min), predicción ~0.35 otra vez; variar dr a tiempos largos; medir la difusión
de las acciones verdaderas por partícula (`rotadores.py`) y compararla con la
amplitud. Para el objetivo físico (Landau) basta saber que es numérico y escoger N y
t_max para que quede por debajo de la señal.

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
(ruido de discreción, ver *Resueltas*).
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

### Con a0=1e-2, la componente que crece también es ruido de discreción

`sg/long20k_a0_1e-2_nrc800` frente a `sg/long_a0_1e-2` (`giro_resolucion.py`):

| ventana | cociente Nrc=800/Nrc=400 | correlación |
|---|---|---|
| [1600,2000] | 1.01 | 0.999 |
| [2000,6000] | 0.24 | 0.78 |
| [6000,12000] | 0.41 | 0.48 |
| [12000,20000] | 0.29 | 0.46 |
| [15000,20000] | 0.25 | 0.56 |

Hasta t=2000 es física resuelta; después, la amplitud depende de la resolución. La
parte estática vale 1.233e-2 en ambas. El crecimiento en el tiempo tiene la misma ley
en las dos resoluciones (pregunta 1).

## Pendientes técnicos menores

- `main.f90`: `call system('mkdir -p')` y `cp` impiden correr en Windows nativo.
- `yoshida6`: un piso de ~1e-21 impide medir su orden limpiamente.
- La rama `clean/comentarios` no está integrada a `main`.
