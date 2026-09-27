# Respuesta a la revisión científica (fase 2)

Fecha: 2026-09-27. Responde a `revision_cientifica.md` punto por punto. Para
cada uno: **estado** (de acuerdo / de acuerdo en parte / en desacuerdo),
**evidencia** (con datos o código ya existentes, verificados hoy) y **qué falta**.
Las cifras citadas se recalcularon para esta respuesta; los comandos están en
`reproducir/scripts/demo_eta.py` y `bulto.py`.

Resumen: la revisión es útil y en lo esencial correcta en su diagnóstico de
qué falta (convergencia, barrido en ε, energía, separar lo demostrado de lo
interpretado). Varias afirmaciones concretas no se sostienen al contrastarlas
con el código, los datos o la literatura: B.1 en su premisa, B.9, la fórmula de
F.4, la escala N^−1/2 de B.2, E.2 y, sobre todo, B.5. El teorema de Hadžić et al.
está formulado en nuestra misma reducción de |L| fijo, lo que corrige también
nuestros documentos y mejora el encuadre del trabajo (segunda lectura, §1).

---

## Segunda lectura (2026-09-27): lo que la primera respuesta omitió

Una relectura completa de la revisión, sección por sección (A–M), encontró
puntos que la primera respuesta no atendió, uno en que la revisión y nuestros
propios documentos estaban equivocados, y varios que se pudieron resolver con los
datos existentes. Todo lo de esta sección se verificó hoy.

### 1. El teorema de Hadžić et al. está formulado en nuestro mismo modelo (corrige B.5, C.10 y nuestros documentos)

La revisión (B.5, J.3) y nuestros documentos decían que la comparación con
Hadžić, Rein, Schrecker y Straub es "motivación, no equivalencia", porque su
resultado sería para sistemas isótropos. **Es falso.** El artículo (verificado:
*Damping versus oscillations for a gravitational Vlasov-Poisson system*, Arch.
Ration. Mech. Anal. **249**, 45 (2025), arXiv:2301.07662) trabaja en **la misma
reducción**: el sistema de Vlasov–Poisson radial con todas las partículas con el
mismo |L|, en el potencial de una masa puntual central fija, con equilibrios
pequeños f = ε(E₀−E)₊^k. Su teorema 1.2:

- (a) si 1/2 < k ≤ 1, existe ε₀(k) tal que para 0 < ε < ε₀ no hay amortiguamiento:
  hay un autovalor por debajo del fondo del espectro esencial (un modo discreto
  bajo Ω_min);
- (b) si k > 1, existe ε₀(k) tal que para 0 < ε < ε₀ hay amortiguamiento de Landau
  (débil, sin tasa) y el espectro puntual es vacío.

Solo importa la regularidad en el borde de vacío, y los autores señalan que King
(k = 1) no se amortigua. La prueba usa la ausencia de autovalores embebidos y un
principio de Birman–Schwinger para el hueco principal: la misma divergencia de
∫F′/(Ω−ω) en Ω_min para g ≤ 1 que la batería usaba como interpretación.

**Consecuencia para la contribución (A.3, A.6).** El encuadre mejora: el estudio η
es la exploración numérica, en el mismo modelo reducido y con otro fondo
(isócrono en vez de masa puntual), de lo que el teorema deja abierto. Primero, el
valor de ε₀(k) para k > 1: η_c ≈ 1.1–1.2 con g = 2 es una estimación de ε₀(2) para
el fondo isócrono. Segundo, qué pasa por encima de ε₀. Tercero, la dinámica no
lineal. Los tres regímenes lineales que encontramos son los del teorema (g = 2
amortigua con masa chica; King oscila justo por debajo de Ω_min, en PIC
ω = 0.06088 ± 0.00005 contra Ω_min = 0.06101; g = 2 deja de amortiguar sobre η_c).

**Referencias que faltaban** (de la bibliografía del artículo): Ramming y Rein,
Phys. D **365**, 72 (2018), estudio numérico de soluciones oscilantes del problema
radial, el antecedente numérico más cercano; Rioseco y Sarbach, Class. Quantum
Grav. **37**, 195027 (2020), mezcla de fases en un potencial central externo;
Hadžić, Rein y Straub, Arch. Ration. Mech. Anal. **243**, 611 (2022), galaxias que
oscilan linealmente; Weinberg, ApJ **421**, 481 (1994), modos débilmente
amortiguados, pertinente para A4. **Hecho:** corregidos `bateria_eta.tex` (reducción, § borde de vacío, §
estado de la teoría) y `demo_eta.tex` (conceptos, límites, conclusiones).

### 2. Afirmaciones de la revisión que no se sostienen (además de B.1, B.9 y F.4)

- **B.2, "la parte fina debe escalar como N^−1/2".** Eso vale para muestreo
  aleatorio. Con partida silenciosa (malla regular en (Q,J)) el ruido no es de
  Monte Carlo: `PREGUNTAS_ABIERTAS.md` mide que pasar de N_rc = 400 a 800 baja el
  ruido de discreción de a₀ = 1e-2 de 1.2e-2 a 3.4e-3, un factor 3.5, no √2. La prueba
  correcta es que la parte fina cambie con N y la lisa no; D5N lo mide.
- **B.3, "robustez frente al tipo de perturbación s(J)".** η_c es la posición de
  un polo del operador linealizado y no depende del dato inicial; s(J) cambia
  cuánto se excita el modo y la forma de las colas, no η_c.
- **E.2, "la anisotropía radial favorece inestabilidades".** La inestabilidad de
  órbitas radiales es no esférica (ℓ ≥ 1); este modelo solo tiene perturbaciones
  radiales (ℓ = 0), y con F monótona en E el criterio de Antonov excluye modos
  radiales inestables. No puede aparecer aquí.
- **E.3, "dependencia en κ/Ω".** Con L fijo la frecuencia acimutal no entra en la
  dinámica; lo que cambia con L₀ es la forma de Ω(J) (el ancho y la pendiente en
  el borde). La pregunta bien formulada es si η_c depende de esa forma; D2 (L₀ = 1)
  ≈ A4 (L₀ = 2) con el ancho escalado es un primer punto.
- **G.1, "volumen de Liouville vía 8π²L₀∫∫F".** En un PIC de pesos fijos esa
  integral es Σw y se conserva por construcción: no prueba nada. El control útil
  del integrador es la comparación con la mezcla libre exacta (G.3, abajo).
- **C.6 y D.5 remiten a una sección "E.5"** que no existe (E termina en E.4). Por
  el contexto, se refieren a F.4–F.6.
- **D.2 pide "D5/4"**: ambiguo (¿Δt/4?). Se hace Δt/2; si D5 cambia, se agrega Δt/4.

### 3. Resuelto hoy con los datos existentes

- **B.7, ω tardío de D5 con incertidumbre.** En h₁, 30–43.5 τ₁:
  ω_PIC = 0.07056 ± 1e-4 contra ω_lin = 0.07022 ± 4e-6; la diferencia,
  (3.4 ± 1)e-4, es real a ~3σ (la revisión citaba 5.7e-4).
- **B.8 / F.1, energía:** ≤ 4e-6 en las 17 corridas de la demo; ≤ 1.7e-4 en el
  bulto (`demo_eta.py energia`, figura en el documento).
- **B.6, lisa/fina a tres escalas:** el exceso de D5 sobre la teoría lineal es
  4.44, 4.40 y 4.31 con suavizados 0.25, 0.5 y 1.0, y 4.45 sin suavizar
  (proyección sobre el perfil del modo).
- **G.2 / F.8, virial de L fijo** (⟨p²⟩ + ⟨L²/r²⟩ − ⟨r ∂_rΦ⟩ = 0 en un estado
  estacionario): residuo estacionario de 4e-5 en Z_A4 y Z_L5 y 4e-4 en Z_M5, del
  orden del error del gradiente de Φ en la malla; en las perturbadas fluctúa con
  la respuesta (2e-4 en D5, 4e-3 en D6, 9e-3 en D10).
- **G.3, solución libre exacta con residuos:** B0 (sin autogravedad) contra
  Σ w e^{−ik(Q₀+Ω(J₀)t)}: error absoluto ≤ 1.6e-7 en h₁ durante 31 τ₁, mientras h₁
  baja de 0.8 a 1e-5; en k = 2, 3, ≤ 9e-7.
- **H.1, paso de tiempo en el pericentro:** Ω_p·Δt ≤ 7.8e-3 en las 22 corridas;
  error de energía por pericentro ≤ 1e-10.
- **J.1 / H.1(ii), resolución de la malla en J cerca del borde:** el modo de L5
  está a 2.5e-4 de Ω_min, que son **solo 5.7 espaciados** de la malla con N_rc = 400
  (11.5 con 800; 23 en `lineal.py`). Es una explicación numérica posible de la
  pérdida de D3 que la primera respuesta no mencionaba; D3N la decide. L6 y M5
  están a 75 y 234 espaciados.
- **C.11, versión del código por corrida:** `demo_eta.py` y `bulto.py` registran
  ahora sha256 del ejecutable y del `.par`, commit y fecha de cada corrida en
  `reproducir/corridas/12_demo_eta/METADATOS.txt`; las anteriores, de forma
  retroactiva (todas con el ejecutable de f649049, sin cambios en `src/` desde el
  23 de septiembre).

### 4. Puntos que siguen pendientes y no estaban en el plan

- **C.2:** derivar η (por qué el corrimiento de la frecuencia media mide el
  acoplamiento) desde la relación de dispersión en acoplamiento débil. Con el
  teorema de Hadžić, la pregunta se precisa: estimar ε₀(k) analíticamente.
- **C.7 / G.5:** comparar la saturación de D5 y D6 con las predicciones de
  O'Neil (amplitud, periodo de rebote). Para una onda en plasma hay resultados
  del umbral de atrapamiento; para un continuo de frecuencias gravitatorio no
  conozco fórmula cerrada. Hay que revisarlo en la literatura antes de afirmarlo.
- **G.6:** agregar g = 1/2 (el borde del teorema) al barrido del borde.
- **G.7:** continuación analítica (Fouvry–Prunet) para separar polo de
  transitorio.
- **H.2:** que el ancho de la isla escale con √ε: hace falta D6 a dos o tres
  amplitudes con isla.
- **H.4:** apéndice de control con el ruido de discreción de `PREGUNTAS_ABIERTAS.md`.
- **I.1:** panel de residuo PIC − lineal en las figuras de D1 y D2.
- **F.6 y F.7:** frecuencia de rebote medida con la FFT de J(t) de las partículas
  del borde (con D6L), y quiebre de pendiente en las colas.
- **F.2, F.3:** perfiles ρ(r,t) y distribución de p_r por capa; baja prioridad
  (δΦ ya da ρ por Poisson).
- **L:** la estructura de artículo propuesta es razonable; la demo no es el
  artículo.
- **K.11 y K.12** no estaban en mi plan: G.7 y levantar L fijo (proyecto aparte).

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

### B.5 Límites del modelo de L fijo — de acuerdo en los límites, en desacuerdo con la comparación

De acuerdo con los límites frente a un sistema esférico real (lista de J.2).
**Hecho:** sección dedicada en `demo_eta.tex`. En desacuerdo con que la
comparación con Hadžić et al. sea "motivación, no equivalencia": su teorema está
en esta misma reducción (segunda lectura, §1).

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

## Plan de trabajo (estado al 2026-09-27)

| # | Punto | Tipo | Estado |
|---|---|---|---|
| 1 | K.1 curva de resolución sintética del matrix pencil | análisis | hecho |
| 2 | F.1 / B.8 energía de las 24 corridas | análisis | hecho |
| 3 | B.6 separación lisa/fina a tres escalas | análisis | hecho |
| 4 | B.9 frase en la batería; C.8 conclusiones; C.1 predicciones; C.3 pies | texto | hecho |
| 5 | B.5 límites de L fijo y corrección de la comparación con Hadžić | texto | hecho |
| 6 | B.7, G.2, G.3, H.1, J.1, C.11 | análisis | hecho |
| 7 | D.1 barrido en ε de D3 (ε = 0.003, 0.01, 0.03) | PIC | corridas hechas; falta analizar |
| 8 | D.2 convergencia N (800×50) y Δt/2 de D3, D5, D6 | PIC | en curso |
| 9 | F.4 / F.5 / F.6 isla, fracción atrapada, rebote (con D6L) | PIC + análisis | D6L en cola |
| 10 | B.4 D8 con ε óptimo | PIC | pendiente |
| 11 | D.4 / G.6 robustez de η_c: g (incl. 1/2), J_t, W₀ | lineal | pendiente |
| 12 | D.3 barrido en η a ε fija; H.2 isla contra √ε | PIC | pendiente |
| 13 | C.7 / G.5 comparación con O'Neil | teoría | pendiente |
| 14 | C.2 derivación de η; estimar ε₀(k) | teoría | pendiente |
| 15 | E.4 η_c(L₀) en el solver lineal | lineal | pendiente |
| 16 | G.7 continuación analítica | teoría + lineal | pendiente |
| 17 | fondo de masa puntual (BGtype = sphere) para comparar con el teorema tal cual | PIC + lineal | pendiente |
