# El bulto como escenario de estudio propio

Fecha: 2026-09-27. Resultados en `docs/demo_eta/demo_eta.pdf` (sección "El
bulto: un dato que no es un equilibrio"); guion en `reproducir/scripts/bulto.py`;
salidas en `exe/demo_eta/bulto/`; videos en `reproducir/videos/demo_eta/bulto_*.mp4`.

## El dato

Es el dato de las corridas originales de autogravedad
(`reproducir/corridas/09_autogravedad`), generado por el propio código
(`state = aa_quad`, `dftype = gauss`) en una malla regular de las variables
ángulo-acción del isócrono desnudo:

    F0(Q, J) = exp(-sin²(Q/2)/s_Q²) · J² exp(-J²/s_J²),   s_Q = s_J = 0.1

Toda la masa está concentrada en una fase orbital, con un ancho de ~0.2 rad. No
es un equilibrio más una perturbación: en el lenguaje del estudio η equivale a
ε ≫ 1 en todos los armónicos hasta k ~ 1/s_Q ≈ 10. Las corridas originales usaban
a₀ = 1e-4 a 1e-2.

## Por qué es otro escenario

Comparte la pregunta de fondo con el estudio η (¿puede la autogravedad impedir
que la mezcla de fases borre una estructura?), pero el régimen y el mecanismo
son distintos.

| | Estudio η (D1–D10) | Bulto (B0–B5) |
|---|---|---|
| Dato inicial | equilibrio autoconsistente + perturbación chica en ángulo | toda la masa en una fase orbital; no es un equilibrio |
| Parámetro que ordena | η (acoplamiento / ancho de banda) y ν = ω_b/γ | μ = ω_b/ΔΩ (pozo propio / cizalla) |
| Por qué h_k no decae | un modo colectivo sale de la banda: efecto **lineal** | el bulto atrapa a sus propias partículas y forma un núcleo autoligado: efecto **no lineal** desde el principio |
| Dónde ocurre | η > η_c ≈ 1.1 | a₀ ≳ 0.03 (μ ≈ 0.3–0.5), con η de solo 0.05–0.16 |
| Referencia teórica | la teoría lineal lo predice (`lineal.py`) | no hay teoría lineal que lo describa |

Se conectan en las islas de atrapamiento. En D6 y D7 la oscilación atrapa una
fracción de las partículas; en el bulto las atrapa casi todas. El bulto es el
extremo no lineal del mismo fenómeno.

## Lo que se midió (2026-09-27)

Siete corridas: N = 400 × 25 como el original, Δt = 0.05, t_fin = 10 000 ≈ 31 τ₁
(τ₁ = 327, banda del isócrono desnudo [0.0513, 0.0705]).

| corrida | a₀ | N_pc | η | μ | ⟨\|h₁\|⟩ 10–20 τ₁, potencial del instante | ídem, isócrono desnudo | masa en el núcleo |
|---|---|---|---|---|---|---|---|
| B0 | sin autogravedad | 25 | 0 | 0 | 7.0e-5 | 7.0e-5 | — |
| B1 | 1e-3 | 25 | 0.005 | 0.10 | 1.1e-4 | 1.0e-3 | — |
| B2 | 1e-2 | 25 | 0.05 | 0.31 | 2.6e-3 | 1.1e-2 | — |
| B3 | 0.03 | 25 | 0.16 | 0.52 | 0.22 | 0.22 | 40 % |
| B4 | 0.1 | 25 | 0.60 | 0.89 | 0.80 | 0.81 | 79 % |
| B5 | 0.3 | 25 | 2.5 | 1.16 | 0.69 | 0.83 | 99 % |
| B4q | 0.1 | 100 | 0.60 | 0.89 | 0.85 | 0.86 | 79 % |

- η es el de la distribución promediada en Q (`eta.py`). μ se estima con la
  oscilación del potencial propio en la primera vuelta radial. "Núcleo" es la
  fracción de la masa con J < 0.06 al final, en el potencial del instante; solo
  tiene sentido cuando hay núcleo (B3–B5).
- **Hay una transición entre a₀ = 1e-2 y 0.03**, con μ entre 0.3 y 0.5. Por
  encima, el bulto colapsa en las primeras vueltas en un núcleo denso y h_k deja de
  tender a cero. En este modelo radial el núcleo es una cáscara esférica delgada,
  atrapada en el pozo que crea su propio potencial. Se mueve como un objeto: en B4
  su radio medio oscila ±0.8 alrededor de 5.7, con frecuencia 0.075. Alrededor
  queda un halo de partículas arrancadas que sí se mezclan.
- **La intuición original era correcta:** al subir la masa, el bulto atrae a las
  partículas y los h_k dejan de tender a cero. Pero con las masas originales
  (a₀ ≤ 1e-2) el efecto es chico.

## Cuidados de medición

- **El mapa importa.** El dato no es estacionario y no hay un equilibrio fijo que
  sirva de marco. Con masa chica, el h₁ medido con el mapa del isócrono desnudo (el
  que usa el código) se estanca diez veces por encima del medido con el potencial
  real. Es en buena parte el efecto de coordenadas ya documentado en
  `PREGUNTAS_ABIERTAS.md`. Con masa grande los dos mapas coinciden.
- **Ruido tardío.** En B2 los armónicos altos crecen al final (|h₄| ≈ 3e-2 en
  20–31 τ₁). Es compatible con el crecimiento del ruido de discreción que
  `PREGUNTAS_ABIERTAS.md` mide para a₀ = 1e-2 (∝ t^1.9). No se interpreta como
  física.
- **Resolución.** Cuadruplicar N_pc (B4q) cambia el |h₁| tardío un 6 % y no cambia
  la masa en el núcleo. No se han variado N_rc ni Δt.

## Qué falta para que sea un resultado

1. **Afinar la transición:** más masas entre 0.01 y 0.1, para fijar la μ crítica
   con incertidumbre.
2. **Ver de qué depende:** variar el ancho del bulto en ángulo (s_Q; un bulto más
   ancho tiene menos potencial propio por unidad de masa) y en acción (s_J), y
   comprobar si la transición cae siempre en el mismo μ.
3. **Convergencia** en N_rc y Δt, y en N_pc más allá de 100.
4. **Describir el estado final:** masa en el núcleo, amplitud y frecuencia de su
   oscilación, y perfil del halo, en función de la masa. La fracción atrapada es
   el observable más directo; conviene definirla con un criterio dinámico (por
   ejemplo, partículas cuya J oscila con el núcleo) y no con un corte en J.
5. **Medir siempre con el potencial real,** o con un marco que siga al núcleo, y no
   con el isócrono desnudo.
6. **Buscar la teoría:** la escala natural es μ = ω_b/ΔΩ; falta una estimación
   analítica de la μ crítica, por ejemplo comparando el tiempo de cizalla del bulto
   con su tiempo de caída libre en su propio pozo.

## Prioridad sugerida

Terminar primero los puntos de la auditoría sobre el estudio η
(`auditoria_deepseekpro/respuesta_revision_cientifica.md`), que es el más
avanzado. El bulto puede seguir como segunda línea, o como una sección aparte
del artículo.
