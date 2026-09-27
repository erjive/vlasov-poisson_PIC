# Respuesta a la revisión científica (fase 2)

Fecha: 2026-09-27. Responde a `revision_cientifica.md` punto por punto. Para
cada uno: **estado** (de acuerdo / de acuerdo en parte / en desacuerdo),
**evidencia** (con datos o código ya existentes, verificados hoy) y **qué falta**.
Las cifras citadas se recalcularon para esta respuesta; los comandos están en
`reproducir/scripts/demo_eta.py` y `bulto.py`.

Resumen: la revisión es útil y en lo esencial correcta en su diagnóstico de
qué falta (convergencia, barrido en ε, energía, separar lo demostrado de lo
interpretado). Tres afirmaciones concretas no se sostienen al contrastarlas con
el código y los datos: B.1 en su premisa, B.9 y la fórmula de F.4.

---

## B. Debilidades científicas

### B.1 La pérdida del 19 % de D3 "no tiene cota γ lineal" — en desacuerdo con la premisa, de acuerdo con la acción

La premisa es que el resultado depende de que el matrix pencil resuelva
γ_lin < 1e-5. No depende: la solución lineal de L5 da una cota directa, sin
ajuste alguno. Su envolvente de ‖δΦ‖ pasa de 0.687 (5 τ₁) a 0.690 (20 τ₁); la
tasa efectiva entre ventanas es −6e-6 a +2e-6. La PIC (D3) da 3.9e-5 con el
mismo cálculo, y el matrix pencil con incertidumbre (commit 2e3e767) da
γ_PIC = (5.2 ± 1.0)e-5 contra γ_lin = (−0.4 ± 0.3)e-5. Una γ lineal de 4–5e-5
implicaría que la solución lineal pierde un 20 % entre 5 y 20 τ₁, y no pierde
nada. **La pérdida no es lineal.**

Lo que la cota no decide es si es no lineal o numérica. Para eso sigue haciendo
falta el barrido en ε (D.1) y la convergencia (D.2). El test sintético del
matrix pencil (K.1) vale igual como documentación del método, y es barato.

**Hecho (K.1):** curva de resolución sintética (`demo_eta.py resolucion`):
señales con γ conocida de 0 a 1e-4, a ω = 0.0725, muestreadas cada 10 y
ajustadas en 10–20 τ₁ con la misma rutina. Limpias, γ se recupera a 1e-16; con un
batido del 5 % (ω = 0.0735, γ = 1e-4) y ruido complejo de 1e-3, a ~1e-6 en todo el
rango (p. ej. 3e-5 → 2.89e-5 ± 0.1e-5). La diferencia de D3 está cincuenta veces
por encima de esa resolución.

**Pendiente:** D.1 (barrido ε = 0.003, 0.01, 0.03, 0.1 en L5).

### B.2 Cero tests de convergencia en lo no lineal — de acuerdo

Correcto y ya señalado en `observaciones.md` §3.3. Única excepción parcial: el
bloque del bulto (B4 contra B4q, N_pc 25 → 100) da la misma masa en el núcleo
(79 %) y un |h₁| tardío 6 % mayor. No cubre Δt ni N_rc, ni a D3/D5/D6.

**Pendiente:** D.2 (N_rc = 800, N_pc = 50 y Δt/2 para D3, D5, D6).

### B.3 η_c empírico y de una familia — de acuerdo en parte

De acuerdo en que η_c es una medición y no una predicción. Pero el trabajo no
afirma universalidad: `bateria_eta.tex` (§ "El barrido en η", resultado 2) dice
explícitamente "η ordena, pero no colapsa del todo" y muestra, en teoría lineal,
que a η = 1 el resultado depende del ancho (A7, A4, A9: amortiguado, borde,
discreto), de la forma (E1, E2, E3) y del borde (G075…G3), y no de L₀ con el
ancho escalado (D2 ≈ A4). Lo que falta es un **barrido en η** para cada variante
que dé η_c(g, J_t, W₀), no solo el punto η = 1.

**Pendiente:** D.4 en el solver lineal (barato); revisar que ningún texto sugiera
universalidad.

### B.4 D8 al límite del rango dinámico — de acuerdo

Coincide con lo que dice el propio documento (§ "Qué limita la medida"): con
ε = 1 la parte estática de segundo orden fija un piso de ~7e-4, y la amplitud
óptima estimada es ε ≈ 0.3–0.5.

**Pendiente:** D8 con ε ≈ 0.3–0.5 y la parte estática restada.

### B.5 Límites del modelo de L fijo — de acuerdo

**Pendiente:** una sección dedicada en el documento (no una nota), con la lista
de J.2.

### B.6 Separación lisa/fina ad hoc — de acuerdo

**Pendiente:** repetir la separación con escalas 0.25, 0.5 y 1.0 en r, y
definirla en el texto de forma reproducible. Es barato: no requiere correr nada.

### B.7 Ajustes sin incertidumbre — resuelto en parte

Desde el commit 2e3e767 cada polo lleva su incertidumbre: la mayor desviación
entre matrix pencil M = 2, 3 sobre tres recortes de la ventana. Cambió una
lectura: en D1 el polo solo está determinado a ±3 % en ω y ±15 % en γ, aunque
PIC y teoría lineal coincidan a 1e-5. Falta aplicarlo al ω tardío de D5
(30–44 τ₁), que la revisión señala con razón.

### B.8 Energía por corrida — de acuerdo

Los HDF5 guardan `kinetic_energy`, `potential_energy` y `total_energy` en cada
instantánea: es extraerlos.

**Pendiente:** F.1, tabla y figura de ΔE/E(t) para las 24 corridas.

### B.9 El factor en Q del h_k del Fortran — en desacuerdo

La revisión dice que e^{−sin²(Q/2)/s_Q²} convoluciona los armónicos. No en el
código. `src/analysish.f90` (líneas ~125–175) calcula primero el coeficiente de
Fourier de ese factor, C(k) = ⟨A(Q) cos kQ⟩, con una cuadratura de Simpson que
no toca partículas, y después suma f_j·B(J_j)·C(k)·e^{−ikQ_j}. Es decir, h_k es la
proyección de f sobre el armónico k de la función de prueba Φ(Q,J) = A(Q)B(J),
y cada armónico queda multiplicado por la constante C(k). Habría convolución si
se multiplicara cada partícula por A(Q_j)·e^{−ikQ_j}, y el código no hace eso. La
frase de la batería es correcta; se puede precisar diciendo que h_k es la
proyección sobre Φ_k(J) = B(J)C(k).

---

## C. El manuscrito

- **C.1, C.2, C.6, C.7, C.8:** de acuerdo. La separación
  demostrado / interpretado / hipótesis de C.8 es la reorganización más útil y se
  puede hacer ya. Pendiente.
- **C.3:** Δt, Δr, N e integrador sí están, pero en una sola sección
  (§ configuración numérica), no junto a cada figura. De acuerdo en repetirlos en
  los pies de figura.
- **C.4:** de acuerdo.
- **C.5:** de acuerdo; ver B.1, B.2 y B.4.
- **C.10:** de acuerdo; hay que verificar la cita de Hadžić et al. (ARMA 2025,
  arXiv:2301.07662) y el enunciado exacto que se usa.
- **C.11:** en parte. Los datos iniciales de la demo son deterministas (malla
  regular, `state = checkpoint` o `aa_quad`); la semilla no interviene. Falta
  registrar el commit del código por corrida: `params_usados.par` no lo guarda.

---

## D–G. Experimentos y análisis propuestos

De acuerdo con el orden de prioridad de K. Observaciones puntuales:

- **F.4, ancho de la isla:** la fórmula de la revisión, 2√(2|δΦ|/|Ω'|), tiene un
  √2 de más. Para H = ½|Ω'|p² + A cos q la separatriz tiene semiancho
  2√(A/|Ω'|) (verificado numéricamente), que es la que usa `bateria_eta.tex`.
  La propuesta de ajustar la isla de D6 con el péndulo sigue siendo buena.
- **F.5, fracción de partículas atrapadas:** de acuerdo; es el observable más
  directo, y el bloque del bulto ya lo sugiere (masa en el núcleo: 40 %, 79 %,
  99 % para a₀ = 0.03, 0.1, 0.3).
- **G.2, virial con L fijo:** de acuerdo en derivarla antes de usarla.
- **E, barrido en L₀:** de acuerdo con E.4: primero en el solver lineal. Ya hay
  un punto (D1, D2 con L₀ = 1) que coincide con L₀ = 2 con el ancho escalado.

## J. Limitaciones

- **J.1, "a₀ ≥ 1e-2: la partida silenciosa ya no es exacta":** aplica al dato
  construido en el mapa del isócrono desnudo (las corridas originales y el
  bulto), no a la demo η, cuyos datos se construyen en el mapa del equilibrio
  autoconsistente (`equilibrio.py`, inversión a 1e-9 en J). En el bulto el efecto
  se ve: con a₀ = 1e-2 el |h₁| medido con el mapa desnudo se estanca diez veces
  por encima del medido con el potencial real.

---

## Plan de trabajo propuesto (orden)

| # | Punto | Tipo | Costo |
|---|---|---|---|
| 1 | K.1 curva de resolución sintética del matrix pencil | análisis | **hecho** |
| 2 | F.1 / B.8 energía de las 24 corridas | análisis | minutos |
| 3 | B.6 separación lisa/fina a tres escalas | análisis | minutos |
| 4 | B.9 precisar la frase en la batería; C.8 reorganizar conclusiones | texto | — |
| 5 | D.1 barrido en ε de D3 (4 corridas) | PIC | ~40 min |
| 6 | D.2 convergencia N y Δt de D3, D5, D6 | PIC | ~6 h |
| 7 | F.4 / F.5 isla con péndulo y fracción atrapada | análisis | horas |
| 8 | B.4 D8 con ε óptimo | PIC | ~30 min |
| 9 | D.4 robustez de η_c (solver lineal) | lineal | ~1 h |
| 10 | D.3 barrido en η a ε fija | PIC | ~2 h |
| 11 | B.5 sección de límites de L fijo | texto | — |
| 12 | E.4 η_c(L₀) en el solver lineal | lineal | ~1 h |
