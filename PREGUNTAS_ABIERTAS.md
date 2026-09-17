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

## 1. Confirmar que la parte estática de la meseta es un efecto de coordenadas

**Qué se sabe.** Con autogravedad, `|h_1|/h_0` se detiene en `1.77e-3` (a0=1e-3).
- No cambia con dr (0.2-0.025), orden de B-spline (1-3), nodos en Q (25-200) ni
  N (1e3-1e5): tres cifras iguales.
- Escala con a0: pendiente log-log +0.96 en t in [15000,20000].
- Su fase es 0.000 y constante de t=2000 a t=20000: no es un modo.
- La J isócrona de partículas individuales oscila ~0.45% pico a pico a su
  frecuencia orbital (`long_fino`).
- Por filas de acción (`paper_runs/scripts/filas_J.py`): cada fila lleva una parte
  estática de 0.2-2% de su amplitud, todas en fase (sum|S_i|/|sum S_i| = 1.19);
  sin autogravedad son 1e-9.

**Interpretación.** El diagnóstico usa el mapa ángulo-acción del isócrono cuando el
potencial es isócrono + autogravedad. Es consistente con todo lo anterior, **pero no
está confirmada de forma directa.**

**Siguiente paso.** Construir el mapa ángulo-acción numérico del potencial real
(cuadraturas con la sustitución `r = rm + ra sin(theta)`, ya validada contra el
isócrono a 1e-15 en `docs/introduccion/figuras/generar_figuras.py`), con el potencial
promediado de las instantáneas tardías, y recalcular h_1. **Predicción: la parte
estática desaparece.** Si no desaparece, la interpretación es falsa.

---

## 2. La componente de h_1 que gira a frecuencia orbital

**Qué se sabe.** Además de la parte estática S, h_1 tiene una componente que gira a
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
- **Mismo efecto de coordenadas**: la agrupación en filas de J isócrona es la
  equivocada. Con el mapa numérico de la pregunta 1, agrupar por J verdadera
  debería devolver filas rígidas y la cancelación suave.
- Para a0=1e-2, donde la componente crece a 86% de S, repetir `por_ventanas` con
  instantáneas (no hecho).

**Preguntas.** ¿Qué la sostiene con a0=1e-3? ¿Por qué crece con a0=1e-2? ¿Es dinámica
colectiva genuina (candidata a lo que se quería medir) o un efecto de coordenadas de
orden superior?

---

## 3. Un montaje limpio para amortiguamiento de Landau

**Problema.** El estado inicial es equilibrio del isócrono *solo*; al activar la
autogravedad toda la componente es a la vez "perturbación" y fuente del potencial.
Eso mezcla el reajuste del equilibrio (el cambio de h_0 del 1.6% en t < 600, las
amplitudes de fila que cambian entre 0.83 y 1.005) con la respuesta a la perturbación.

**Siguiente paso.** Construir un equilibrio autoconsistente F_eq(J) del potencial
total (iterando Poisson con el mapa numérico de la pregunta 1) y agregar encima una
perturbación pequeña y separada. Solo entonces tiene sentido medir omega_r y gamma.

---

## 4. Régimen a0 >= 1e-2

Con a0=1e-2 el cambio de h_0 es 14% y la componente que gira crece hasta dominar.
El mapa analítico ya no sirve. Depende de 1 y 3. Para perturbaciones grandes,
relajación violenta (Lynden-Bell 1967).

---

## 5. Modos discretos

Resolver el problema lineal de autovalores (método matricial de Kalnajs) para este
equilibrio y comparar omega_r, gamma con las corridas. Un modo dentro de la banda
orbital es resonante y se amortigua por Landau; fuera, no.

---

## 6. Levantar L fijo

Con dispersión en L hay dos frecuencias y resonancias entre ellas. Existe una rama
con `l_part` en `VlasovPoisson_PIC_sp`. Es otro proyecto.

---

## Pendientes técnicos menores

- `main.f90`: `call system('mkdir -p')` y `cp` impiden correr en Windows nativo.
- `yoshida6`: un piso de ~1e-21 impide medir su orden limpiamente.
- La rama `clean/comentarios` no está integrada a `main`.
