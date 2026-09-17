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

## 1. La componente de h_1 que gira a frecuencia orbital

**Qué se sabe.** Además de la parte estática S (que resultó ser un efecto de
coordenadas, ver *Resueltas*), h_1 tiene una componente que gira a
omega = 0.055-0.062 (dentro de la banda orbital [0.047, 0.071]):

| a0 | t in [2000,15000] | t in [15000,20000] |
|---|---|---|
| 1e-4 | 1.4% de S | 1.3% de S |
| 1e-3 | 8.7% de S | 10.2% de S |
| 1e-2 | 18% de S | **86% de S** (rodea el origen: saltos de fase) |

Por filas (`long_fino`, t <= 2000): hasta t=2000 es el continuo de filas que aún no
termina de cancelarse (el ajuste S_i + A_i exp(-i w_i t) por fila reproduce 14.5 de
los 15.4% medidos en [1300,2000]). Los residuos de cada fila están en su propia
omega_i y en 2 omega_i y son incoherentes entre filas (coherencia 0.017): **no hay
un modo con frecuencia común en t <= 2000**.

**Lo que no cuadra.** Ese continuo, extrapolado con los A_i y w_i ajustados, bajaría a
0.5% de S en [2000,15000] y 0.1% en [15000,20000]. Se mide ~9-10%. Algo después de
t=2000 detiene la cancelación.

**Corrida hasta t=20000 con instantáneas** (`sg/long20k_snap`, reproduce
`long_a0_1e-3` a 2.5e-21), ajuste por filas en ventanas de 2000
(`filas_J.por_ventanas`):

| ventana | abs(S)/h0 | arg S | rms(h-S)/abs(S) | rms C/abs(S) | rms residuos/abs(S) | amplitud de fila relativa a t=0 |
|---|---|---|---|---|---|---|
| [2600,4600] | 1.7663e-3 | +0.000 | 0.032 | 0.036 | 0.009 | 0.91 |
| [8600,10600] | 1.7659e-3 | -0.000 | 0.083 | 0.074 | 0.017 | 0.68 |
| [14600,16600] | 1.7657e-3 | -0.001 | 0.086 | 0.092 | 0.025 | 0.37 |
| [16600,18600] | 1.7655e-3 | -0.001 | 0.104 | 0.112 | 0.024 | 0.26 |

Tres cosas quedan establecidas:

1. **La parte estática no se mueve** en 16 000 unidades de tiempo (cuarta cifra), y
   es la suma de contribuciones de fila casi todas en fase.
2. **La componente que gira es el continuo de filas en cada ventana**: el rms de
   C = sum A_i exp(-i w_i t) coincide con el de h - S, y la suma de residuos es 3-10
   veces menor. No hay una oscilación coherente a frecuencia común.
3. **Las filas no son rígidas**: su amplitud cae de 0.91 a 0.26. Las 25 partículas de
   una fila comparten J *isócrona* inicial pero no J verdadera (dispersión de J_iso
   dentro de la fila ~2.9e-3), así que giran a frecuencias distintas y la fila se
   desfasa internamente en ~1/(abs(omega') dJ) ~ 1e4. Las frecuencias de fila también
   cambian (rms 2e-4 entre la primera y la última ventana, con una parte no suave en J
   de 2.2e-5). Por eso falló la extrapolación desde t <= 2000: suponía filas rígidas.

**Lo que sigue abierto.** Por qué ese continuo no se cancela por debajo de ~10% de S.
Hipótesis a probar:
- **Discreción de la rejilla**: con 400 filas cuyas fases evolucionan de forma no
  suave en J, la suma no se cancela por debajo de un piso incoherente ~ 1/sqrt(Nrc).
  Predicción: con Nrc=1600 el componente baja a la mitad. Costo: t=20000 con
  N=4e4 son ~45 min.
- ~~**Mismo efecto de coordenadas**~~ **Descartada.** Con el mapa numérico del
  potencial real (ver *Resueltas*), la parte que gira no cambia: rms 1.540e-4 frente
  a 1.541e-4, correlación 0.99988 entre los dos cálculos de h_1. No es un efecto del
  mapa, es una propiedad de la distribución.
- **Dinámica colectiva genuina**: es lo que queda si se descarta la discreción.
- Para a0=1e-2, donde la componente crece a 86% de S, repetir `por_ventanas` con
  instantáneas (no hecho).

**Preguntas.** ¿Qué la sostiene con a0=1e-3? ¿Por qué crece con a0=1e-2? ¿Es dinámica
colectiva genuina (candidata a lo que se quería medir) o un piso de discreción?
Con el mapa numérico, repetir el análisis por filas agrupando por J verdadera: si
las filas vuelven a ser rígidas y la componente sigue ahí, la discreción pierde peso.

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

Con a0=1e-2 el cambio de h_0 es 14% y la componente que gira crece hasta dominar.
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

## Pendientes técnicos menores

- `main.f90`: `call system('mkdir -p')` y `cp` impiden correr en Windows nativo.
- `yoshida6`: un piso de ~1e-21 impide medir su orden limpiamente.
- La rama `clean/comentarios` no está integrada a `main`.
