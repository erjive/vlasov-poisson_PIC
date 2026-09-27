# Revisión científica del estudio eta: segunda fase

Fecha: 2026-09-27
Alcance: investigación científica y manuscrito (no solo código)
Contexto usado: código `src/`, auditoría previa
(`auditoria_deepseekpro/observaciones.md`), documentos
`docs/demo_eta/demo_eta.tex` y `docs/experimento_eta/bateria_eta.tex`,
notas `docs/introduccion/vlasov_intro.tex`, archivo de preguntas abiertas
(`PREGUNTAS_ABIERTAS.md`), resultados y figuras citadas en los documentos.

**Nota de método.** El acceso a web no funcionó en esta sesión, de modo que
las referencias externas (Hadžić–Rein–Schrecker–Straub, O'Neil, Antonov,
Kalnajs, Fouvry–Prunet) se toman de lo citado dentro del propio repositorio.
No se verificaron en la literatura, y **no se afirma novedad frente a
resultados externos** más allá de lo que los documentos internos demuestran.

---

## A. Contribución científica real

### A.1 Qué problema aborda

El problema físico es antiguo y preciso: **en un sistema autogravitante sin
colisiones, ¿la autogravedad de una componente de masa pequeña puede sostener
una oscilación colectiva que la mezcla de fases, por sí sola, destruiría?**
La mezcla de fases (cizalla de ángulos por dispersión de frecuencias) hace
decaer cualquier perturbación; la autogravedad puede acoplarla en un modo
colectivo y, si el modo sale de la banda de frecuencias orbitales, la
perturbación sobrevive sin amortiguar.

Es el problema de **amortiguamiento de Landau frente a modos discretos** en
sistemas esféricos, en el régimen de masa pequeña, estudiado con una
reducción extrema que permite un fondo exacto y una teoría lineal de
referencia numéricamente limpia.

### A.2 Pregunta de investigación precisa

La pregunta operativa que el trabajo construye es:

> Supongamos que la autogravedad corre las frecuencias orbitales en una
> fracción `δΩ` y que la banda de frecuencias tiene ancho `ΔΩ`. ¿Existe un
> valor crítico `η_c` del cociente `η = δΩ / (ΔΩ/Ω̄)` tal que, para `η < η_c`
> toda perturbación se mezcla (amortiguamiento), y para `η > η_c` aparece un
> modo discreto sin amortiguar? Y cerca de esa transición, ¿qué hace la
> dinámica **no lineal** a amplitud finita?

La segunda pregunta (no lineal) es la parte genuinamente nueva: la teoría
lineal predice `γ→0` al acercarse a `η_c` desde abajo, de modo que el
parámetro de O'Neil `ν = ω_b/γ` diverge para **cualquier** `ε`. Por tanto,
justo donde el modo está a punto de sobrevivir, la teoría lineal deja de
valer por atrapamiento, no por amplitud.

### A.3 Qué es genuinamente nuevo

A juicio de esta revisión, hay **tres** contribuciones potencialmente nuevas,
con distinto grado de madurez:

1. **El criterio `η` como parámetro ordenador de la transición** (nuevo como
   organización conceptual; ver B.3 sobre su estado actual). Une dos
   cantidades que suelen tratarse por separado: el corrimiento de frecuencia
   por autogravedad y el ancho de banda que fija la mezcla. Es una idea
   sensata y bien motivada, pero hoy es una **medición empírica**, no una
   predicción teórica.

2. **La región no lineal del borde de la transición** (D5, D6). La
   observación de que a `ε=0.075` la PIC deja de decaer después de ~20 τ₁,
   con saturación y rebote, mientras la teoría lineal sigue cayendo, es el
   resultado más interesante del trabajo. Es exactamente el régimen
   `ν≳1` con `γ→0`, que la literatura de plasmas (O'Neil) describe pero que
   en autogravedad esférica con mezcla de un continuo de frecuencias no
   está, a juicio de esta revisión, probado con esta limpieza.

3. **La fragilidad de un modo discreto pegado al borde** (D3 contra D4).
   Que un modo discreto a `x=-0.01` pierda 19% de amplitud a `ε=0.1`
   mientras otro a `x=-0.19` pierde solo un 5% extra, con una explicación
   en términos de partículas empujadas hacia la resonancia, es una
   observación física no trivial **si** se confirma que la pérdida no es
   lineal (ver B.1).

### A.4 Qué es metodológico, no científico

Buena parte del trabajo documentado es método, no resultado:

- La construcción del equilibrio autoconsistente (`equilibrio.py`).
- La partida silenciosa en malla `(Q,J)` con perturbación en pesos.
- El diagnóstico `h_k` factorizado con fase compleja (exportado a
  `analysish.f90`).
- La validación PIC contra `lineal.py` al 1% en régimen lineal.
- La validación del mapa ángulo-acción a 1e-15 en J y 1e-12 en Q.
- Las pruebas de que la componente "estática" de `h₁` es un efecto de
  coordenadas (documentado en `PREGUNTAS_ABIERTAS.md`).

Esto es infraestructura científica de calidad, pero no es el aporte. Un
manuscrito que venda el método como contribución diluye el trabajo.

### A.5 Qué podría interesar a la comunidad

- A la comunidad de **dinámica estelar/galáctica**: la conexión entre la
  dicotomía del borde tipo Hadžić (exponente del vacío) y la transición
  controlada por `η`; y el régimen no lineal cerca de `γ→0`.
- A la de **plasmas**: una realización gravitatoria del régimen de O'Neil
  cerca de una transición marginal, con la diferencia clave de que aquí hay
  un continuo de frecuencias orbitales, no una sola resonancia.
- A la de **métodos PIC**: un ejemplo limpio de cómo usar un solver lineal
  como referencia y reservar el PIC para lo no lineal.

### A.6 Si la contribución es débil, decirlo

**Lo es en su estado actual para los resultados no lineales.** La demo
demuestra que el aparato funciona; no demuestra todavía los resultados no
lineales con el estándar de publicación. Las tres contribuciones de A.3
están hoy en estado de "observación sugestiva" (D5/D6/D8), "medición
empírica" (η_c) y "afirmación que depende de una cota no demostrada" (D3).
La sección B detalla por qué.

---

## B. Principales debilidades científicas

### B.1 [CRÍTICA] La pérdida del 19% de D3 no tiene cota γ lineal

El resultado "la teoría lineal no predice la pérdida del 19%" depende por
completo de que `γ_lineal < 10⁻⁵`. Si la γ lineal verdadera es 4–5×10⁻⁵, la
pérdida es lineal y el resultado no es nuevo. **No se ha demostrado la
resolución del matrix pencil** para γ en ese rango, ni se ha hecho el
barrido `ε→0.003` que cerraría el argumento. Es la debilidad más grave del
trabajo.

**Cómo resolverla.** (i) Publicar la curva de resolución: polo sintético
con γ conocida de 1e-6 a 1e-4, extraída con el mismo matrix pencil y las
mismas ventanas; (ii) barrer ε en D3 y mostrar que la pérdida de amplitud
normalizada por ε tiende a cero al bajar ε; (iii) si no tiende a cero, el
mecanismo "el modo empuja partículas a la resonancia" queda descartado y hay
que reformular la interpretación.

### B.2 [CRÍTICA] Cero tests de convergencia en los resultados no lineales

D3 (19%), D5 (saturación a 44 τ₁), D6 (isla) y D8 (dicotomía) descansan en
**una sola resolución** (N=10⁴, Δt no citado por corrida, malla r∈[0,20],
Δr=0.1). La "parte fina" de D5 (hasta 2×10⁻²) podría ser ruido de
discretización amplificado por 1/ε. Sin:

- partir Δt,
- cuadruplicar N (Nrc=800, Npc=50),
- mostrar que la parte fina escala como N⁻¹ᐟ² y lo liso no cambia,

ninguno de los resultados no lineales es publicable. Esto **ya está
identificado** en la auditoría previa (`observaciones.md`, §3.3) y no se ha
hecho.

### B.3 [MAJOR] `η_c` es empírico y de una sola familia

`η_c≈1.1–1.2` sale del solver lineal para una familia (Wilson g=2,
J_t=0.138, W₀=3, L₀=2). No hay evidencia de robustez frente a:

- forma del borde (g=1, g=1/2);
- ancho de la banda (variar J_t);
- forma del potencial (variar L₀, o campo externo);
- tipo de perturbación (s(J) diferente).

Sin el bloque V de la batería (confirmar robustez de η_c), la frase "η_c es
un parámetro universal" no está justificada. Hoy es **una medición**, no una
predicción.

### B.4 [MAJOR] D8 está al límite del rango dinámico

La diferencia de D8 (g=1) es 2.6×10⁻³, ~13× el piso de 2×10⁻⁴, y se pierde
a ~15 τ₁. Con ε≈0.3–0.5 (no ε=1) y sustracción de la componente estática
podría ampliarse la ventana. Hoy es sugestivo, no concluyente.

### B.5 [MAJOR] El modelo de L fijo limita la analogía con sistemas esféricos

El documento es honesto al llamar al modelo "extremadamente anisótropo",
pero conviene más explicitud: de un sistema esférico isótropo, **solo el
borde Ω_min** es representativo. No hay dispersión en L, no hay resonancias
entre frecuencias radiales y acimutales, y la medida de fase `8π²L` fija
una estructura orbital rígida. La comparación con Hadžić et al. es
motivación, no equivalencia.

### B.6 [MODERADA] La separación "lisa/fina" es ad hoc

El suavizado sobre 0.5 en r (5 celdas) separa lo "liso" de lo "fino" sin
justificación. Falta robustez frente a otras escalas (0.25, 1.0) y una
definición que el lector pueda reproducir.

### B.7 [MODERADA] Ajustes de polo sobre piso, sin incertidumbre

Extraer ω de D5 en 30–44 τ₁ (Δω=5.7×10⁻⁴ sobre señal 10⁻²) está al límite.
No se reportan incertidumbres de ajuste ni la sensibilidad a la ventana.

### B.8 [MODERADA] Conservación de energía por corrida no reportada

Con islas y atrapamiento a tiempos largos, la deriva de energía es el
control mínimo del integrador. Debe ser una figura o tabla estándar.

### B.9 [MENOR] El factor en Q del h_k del Fortran está mal descrito

La batería dice que `e^{-sin²(Q/2)/s_Q²}` "solo multiplica cada armónico por
una constante". Falso: es una función periódica de Q que **convoluciona** los
armónicos; solo es constante por armónico en `s_Q→∞`. No afecta al estudio
(donde se usa `s_Q` grande), pero es un error conceptual del manuscrito.

---

## C. Debilidades del manuscrito como tal

Clasificación referee (CRÍTICA / MAJOR / MODERADA / MENOR), científica, no
de estilo.

### C.1 Motivación

**MODERADA.** La motivación es correcta pero incompleta: no deja claro
desde el principio cuál es la **predicción falsable** del trabajo. Una
motivación fuerte diría "predecimos que la PIC debe (i) seguir a la teoría
lineal donde ν≪1; (ii) saturar donde ν≳1; (iii) perder el modo discreto
pegado al borde por empuje a resonancia". La demo prueba (i) y sugiere
(ii)/(iii); la motivación debería anunciarlo.

### C.2 Formulación teórica

**MODERADA.** La reducción a L fijo, las variables ángulo-acción y la
ecuación linealizada están bien expuestas. Lo que falta es una derivación
**explícita** del `η` (por qué `δΩ` medido así y no de otra forma) y una
discusión de en qué sentido `η_c` debería ser universal. Esa es la pieza
teórica que separaría el trabajo de una medición.

### C.3 Método numérico

**MENOR.** Suficientemente descrito para reproducir; la parte PIC está
validada. Lo que falta es el estándar de estatus: N, Δt, Δr, semilla,
integrador por corrida en las figuras de resultados, no solo en una tabla.

### C.4 Validación

**FUERTE** (esto es lo mejor del trabajo). La validación contra `lineal.py`
(1% en ω y γ; 1e-4 en ω de modos discretos) es exactamente lo que un árbitro
quiere. No es débil.

### C.5 Resultados

**MAJOR.** Los resultados lineales (D1, D2, D4 parcial, D8 parcial) son
sólidos. Los no lineales (D3, D5, D6, D8 en lo no lineal) son presentados
con más confianza de la que la evidencia soporta (ver B.1, B.2, B.4).

### C.6 Interpretación física

**MAJOR.** La interpretación de D3 ("empuje a resonancia") es plausible pero
**no es la única**: también podría ser pérdida lineal mal acotada (B.1) o
discretización (B.2). La interpretación de D5/D6 ("atrapamiento tipo
O'Neil") es correcta en concepto, pero no se da un criterio cuantitativo
independiente del atrapamiento (ver E.5).

### C.7 Discusión

**MODERADA.** Faltan dos piezas: (i) comparación explícita con predicciones
analíticas del régimen de O'Neil (amplitud de saturación, frecuencia de
rebote) más allá del orden de magnitud; (ii) límites del modelo de L fijo
como sección dedicada, no como nota.

### C.8 Conclusiones

**MAJOR.** Las conclusiones actuales mezclan lo demostrado con lo
interpretado. Hay que separar:

- **Demostrado**: la PIC reproduce la teoría lineal donde ν≪1; la transición
  existe en el solver lineal; la PIC se aparta de la teoría lineal cerca de
  η_c.
- **Interpretado**: atrapamiento tipo O'Neil; fragilidad del modo discreto
  por empuje a resonancia.
- **Hipótesis**: universalidad de η_c; conexión con la dicotomía de borde
  en sistemas esféricos completos.

### C.9 Figuras

Ver sección F.

### C.10 Referencias

**MODERADA.** Las referencias citadas (O'Neil, Antonov, Kalnajs, Hadžić et
al., Fouvry–Prunet) son pertinentes. No se puede verificar en esta sesión si
la atribución exacta es correcta; el manuscrito **debe** verificar la cita
de Hadžić et al. (ARMA 2025) y el resultado exacto que se usa.

### C.11 Reproducibilidad

**FUERTE.** Los `.par` de configuración, los scripts (`demo_eta.py`,
`lineal.py`, `equilibrio.py`) y el apéndice de reproducción están en orden.
Falta: semillas, hashes de corridas, y versión exacta del código (commit)
por corrida.

---

## D. Experimentos adicionales más valiosos

Cada uno responde una pregunta, con método, observable, resultado esperado,
qué se aprende y dificultad.

### D.1 Barrido en ε del modo discreto pegado al borde (obligatorio)

- **Pregunta**: ¿la pérdida del 19% de D3 es no lineal (empuje a resonancia)
  o lineal (γ no resuelta)?
- **Método**: repetir D3 para ε = 0.003, 0.01, 0.03, 0.1; misma malla; restar
  Z_L5.
- **Observable**: pérdida de amplitud normalizada por ε en 15 τ₁.
- **Esperado**: si el mecanismo es no lineal, la pérdida/ε tiende a cero al
  bajar ε. Si es lineal, tiende a un valor fijo.
- **Se aprende**: si el fenómeno "modo discreto frágil" existe o es un
  artefacto de cota.
- **Dificultad**: baja (4 corridas cortas, infraestructura existente).

### D.2 Convergencia en Δt y N para D3/D5/D6 (obligatorio)

- **Pregunta**: ¿la parte fina de D5/D6 es física o ruido de discretización?
- **Método**: Nrc=800 y Npc=50 (N=4×10⁴); Δt y Δt/2; para D5 también D5/4.
- **Observable**: h₁(t) en ventanas; separar componente lisa y fina.
- **Esperado**: fina ∝ N⁻¹ᐟ²; lisa invariante.
- **Se aprende**: si los resultados no lineales sobreviven.
- **Dificultad**: media (corridas 4× más costosas).

### D.3 Barrido η a ε fija pequeña (importante)

- **Pregunta**: ¿dónde, en función de η, la PIC deja de seguir a la teoría
  lineal para un ε fijo y pequeño, y coincide ese lugar con ν≈1?
- **Método**: η = 0.1, 0.5, 0.8, 1.0, 1.1 con ε=0.075.
- **Observable**: tiempo de despegue de la PIC respecto a `lineal.py` frente
  a η; comparación con ν(η).
- **Se aprende**: si ν, y no ε, es el parámetro que decide la no linealidad,
  en presencia de mezcla de un continuo.
- **Dificultad**: media (5 corridas).

### D.4 Robustez de η_c frente a la forma de la distribución (importante)

- **Pregunta**: ¿η_c≈1.1–1.2 es universal o depende del borde?
- **Método**: g=1 (King), g=1/2, g=2 con W₀=1.5 y 6; medir η_c en el solver
  lineal. Bloque V de la batería.
- **Observable**: η_c frente a g y W₀.
- **Se aprende**: si la transición está controlada solo por η o también por
  la regularidad del borde.
- **Dificultad**: media (más lineal que PIC).

### D.5 Estudio sistemático de L₀ (importante, con precaución)

Ver sección F (análisis L₀) y la crítica en E.5 de abajo. **Respuesta corta:
sí, un barrido en L₀ es científicamente significativo, pero hay que
formularlo bien.**

- **Pregunta**: ¿cómo dependen η_c, γ y el régimen no lineal del parámetro
  de forma L₀?
- **Método**: L₀ = 0.5, 1, 2, 4 a J_t ajustado para mantener ancho de banda
  y forma comparables.
- **Observable**: curva η_c(L₀) del solver lineal; una corrida PIC en el
  régimen no lineal por valor de L₀.
- **Se aprende**: si la barrera centrífuga cambia la transición, o si la
  reducción en (r,p_r) es universal.
- **Dificultad**: alta (requiere reequilibrar y revalidar).

---

## E. Análisis de la suposición de L fijo: ¿vale la pena variar L₀?

### E.1 Qué cambia con L₀

Con L₀ fijo, el potencial efectivo es `V_eff = Φ + L₀²/(2r²)`. Variar L₀
cambia:

- la posición del mínimo (órbita circular r_c),
- la frecuencia epicíclica κ y su relación con Ω,
- el ancho de banda ΔΩ para un dado rango de acciones,
- la forma de la barrera centrífuga y, por tanto, la regularidad efectiva
  del borde de vacío.

### E.2 Física accesible

Un barrido en L₀ exploraría, en el mismo formalismo 2D:

- **Transición L₀→0**: el límite L₀=0 (órbitas radiales puras) es singular en
  la convención actual (el código aborta en `Lfix=0`); acercarse a L₀→0
  cruza hacia un régimen donde la barrera centrífuga desaparece y el
  equilibrio cambia de forma.
- **Dependencia de η_c con la forma del potencial**: a L₀ distinto, ΔΩ y δΩ
  no escalan igual; eso prueba si η es un buen parámetro o si esconde una
  dependencia en κ/Ω.
- **Estabilidad y colapso**: a L₀ pequeño, el sistema se vuelve más
  radialmente dominado; es sabido que la anisotropía radial favorece
  inestabilidades. Con L₀ variable, uno puede ver si la transición
  amortiguado→discreto se adelanta o se atrasa con la anisotropía.
- **Relajación y mezcla**: la dispersión de L en un sistema real mezcla las
  dos frecuencias; el límite L₀→0 es el régimen donde la medida de fase
  `8π²L` pesa distinto las órbitas, y el observable h_k cambia de peso.

### E.3 Veredicto crítico: ¿produce resultados significativos?

**Sí, pero no es la prioridad inmediata.** Un barrido en L₀ es
científicamente significativo **si** se formula como "dependencia de la
transición y del régimen no lineal con la forma del potencial efectivo", no
como "más corridas variando un botón". Sin una predicción previa de cómo
debería moverse η_c con L₀, el barrido es descriptivo y su interés
disminuye.

**La pregunta correcta es**: ¿η_c depende solo de η definido como en la
sección A.2, o aparece una dependencia adicional en κ/Ω (que varía con L₀)?
Si η_c es constante frente a L₀, η queda validado como parámetro universal;
si no, el trabajo gana una ley de escala η_c(κ/Ω) que **sí** sería un
resultado nuevo y teóricamente interesante.

### E.4 Recomendación

Hacer el barrido **en el solver lineal primero** (barato): medir η_c(L₀)
para 4–5 valores, a ancho de banda comparable. Solo si aparece una
dependencia no trivial, invertir en corridas PIC. Esto evita el gasto
inútil.

---

## F. Análisis adicional posible con los datos existentes

Sin correr nada nuevo, los datos actuales pueden dar resultados extra:

### F.1 Conservación de energía por corrida (obligatorio)

Extraer de los snapshots existentes la serie temporal de energía (T, W,
Φ_ext, Φ_self separados) y reportar la deriva relativa. Control mínimo del
integrador; probablemente ya calculada por `energy.f90`, falta reportarla.

### F.2 Perfiles de densidad y su evolución

ρ(r,t) en ventanas: ¿la componente "respira" con amplitud proporcional a ε
solamente, o cambia la forma? Un perfil ρ(r,t) a t fijo separaría
deformación global de ruido local.

### F.3 Distribuciones de p_r por capa radial

p_r en anillos: ¿la perturbación cambia la anisotropía local? Con L₀ fijo,
la anisotropía en r es fija por construcción, pero la cola de p_r cerca del
borde muestra si hay partículas aceleradas por atrapamiento.

### F.4 Estructura del espacio fase (Q,J)

Más allá de figuras: ajustar la isla de D6 con el modelo de péndulo.
Ancho de la isla en J predicho `ΔJ_isla ≈ 2√(2|δΦ|/|Ω'|)`; comparar con la
medida. Si cierra, es una **demostración cuantitativa** del atrapamiento, no
una impresión visual.

### F.5 Clasificación de partículas atrapadas vs. libres por órbita

Integrar J(t) por partícula y clasificar por desviación cuadrática media
⟨(ΔJ)²⟩: atrapadas (J oscila con ω_b) vs. libres (J constante). Esto da un
**criterio independiente** de la figura y un número: fracción de partículas
atrapadas.

### F.6 Análisis de frecuencia de las partículas del borde

FFT de J(t) para partículas en J≈J_t: picos en ω_b. Comparar ω_b medida con
la teórica √(|Ω'||δΦ|). Más evidencia cuantitativa del régimen de O'Neil.

### F.7 Escalamiento de colas

Extender la regla de colas t⁻⁽α⁺¹⁾ a las ventanas largas ya disponibles:
¿la cola de D5 obedece t⁻³ (g=2, α=2) y la de D8 t⁻² (g=1, α=1) antes de
saturar? Si se distingue un quiebre de pendiente, eso marca la transición
mezcla→modo.

### F.8 Virialización

Calcular el cociente virial 2T/|W| en el tiempo: ¿se mantiene ≈1 (equilibrio)
o la perturbación lo aparta? Con L₀ fijo hay que derivar la forma virial
correcta (ver G.2); si no se deriva, reportar el cociente tal cual como
diagnóstico, no como ley.

---

## G. Predicciones analíticas comprobables

### G.1 Leyes de conservación (ya disponibles)

- Masa total: trivial pero reportable como control.
- Energía total: idem (F.1).
- Volumen de fase de Liouville: con integrador simpléctico debería
  conservarse a la precisión del método; medir la desviación de la medida
  `8π²L₀ ∫∫ F dr dp` da un control del PIC.

**Cómo convertirlo en test**: reportar la serie temporal de cada una en las
14 corridas y una cota.

### G.2 Relación virial

Para L fijo, la relación virial se obtiene integrando `d(rp_r)/dt`. Hay que
derivarla con el término centrífugo y con el potencial propio; no citar la
forma isótropa. **Test**: en el equilibrio sin perturbar (corridas Z),
2T + W_centrífugo + W_ext + W_self ≈ 0 a precisión del integrador; en las
perturbadas, la desviación mide la respuesta dinámica.

### G.3 Solución libre exacta (ya usada)

La envolvente gaussiana `e^{-σ_Ω²t²/2}` (auditoría previa, §2.2) se reproduce
bien; hay que mostrarla como figura cuantitativa con residuos, no solo
citar el acuerdo.

### G.4 Frecuencias del isócrono (ya usada)

Ω(J)=(J+c)⁻³ con c=2.414; κ²=4πρ(r_c)+GM(<r_c)/r_c³ verificada. Test
nuevo posible: derivar Ω'(J) analíticamente y usarla en el ancho de isla
(F.4) y la frecuencia de rebote ω_b.

### G.5 Amplitud de saturación de O'Neil

El modelo de O'Neil predice la amplitud de saturación en función de ν. Hay
que **escribir la fórmula** y compararla con la meseta de D5/D6. Si no hay
fórmula cerrada en la literatura para el caso con continuo de frecuencias,
decirlo y comparar con el valor asintótico de `lineal.py` como referencia.

### G.6 Criterio de estabilidad del borde

La analogía con Hadžić: `F ∝ (E_t−E)^g`; predecir para g=1/2, 1, 2 si hay
modo discreto a η fijo. Ya está en la batería; convertir en tabla de
predicción vs. medición para tres g.

### G.7 Continuación analítica

La batería menciona Fouvry–Prunet (2022) para decidir si ω−iγ es polo de la
relación de dispersión. Es un test analítico-numérico valioso; hoy solo es
una mención. Si se implementa, separa "polo" de "respuesta transitoria".

---

## H. Distinguir física de artefactos numéricos

Para cada resultado clave:

### H.1 D3: pérdida del 19%

**Posibles artefactos**: (i) γ lineal no resuelta (B.1); (ii) N finito y
discretización en J que difuminan el borde; (iii) Δt demasiado grande cerca
del pericentro. **Test**: barrido ε (D.1) + convergencia (D.2) + medir el
pericentro mínimo y el criterio de Courant. Nota: la auditoría previa
documenta que el criterio de Courant en pericentro está implementado (E11).

### H.2 D5/D6: saturación e isla

**Posibles artefactos**: la parte fina ∝ N⁻¹ᐟ² (B.2); la isla puede ser una
estructura de la rejilla en (Q,J) más que del continuo. **Test**: convergencia
en N; además, ajustar la isla con el modelo de péndulo (F.4): si el ancho
medido escala con √ε y √|δΦ| como predice, es física.

### H.3 D8: dicotomía de borde

**Posible artefacto**: rango dinámico insuficiente (B.4); la respuesta
persistente puede ser una cola algebraica mal distinguida de un modo.
**Test**: ε óptimo y sustracción de la parte estática; y análisis de
frecuencia para ver si la parte tardía tiene frecuencia Ω_min limpia.

### H.4 Ruido de discreción creciente (documentado en PREGUNTAS)

Ya está identificado como numérico con pruebas de resolución. Para el
manuscrito: incluirlo como **apéndice de control**, no como resultado.

---

## I. Figuras: evaluación y propuestas

### I.1 Figuras actuales

| Figura | Pregunta científica | ¿Apoya la afirmación? | Falta | Veredicto |
|---|---|---|---|---|
| Plano (η,ν) con corridas | ¿Dónde cae cada corrida? | Sí | — | Conservar |
| h₁(t) D1/D2 | ¿PIC sigue a lineal? | Sí | incertidumbre del polo | Conservar; añadir residuo |
| h₁(t) D3/D4 | ¿El modo pierde amplitud? | Parcial: solo una resolución | curva ε, γ lineal acotada | Conservar; reforzar |
| h₁(t) D5/D6 | ¿Satura/rebota? | Parcial: parte fina no verificada | convergencia N | Conservar; añadir N |
| h₁(t) D8 | ¿Dicotomía de borde? | Parcial: ventana estrecha | ε óptimo | Rehacer con ε=0.3–0.5 |
| Fase (r,p_r), (Q,J) por corrida | ¿Estructura? | Cualitativo | ajuste de isla | Conservar como apoyo visual |
| Videos | ¿Evolución? | Cualitativo | — | Opcional |

### I.2 Figuras nuevas recomendadas (solo si aportan)

1. **Convergencia**: h₁(t) superpuesto para N=10⁴, 4×10⁴ (D5) y Δt, Δt/2;
   separación lisa/fina.
2. **Barrido ε de D3**: pérdida/ε frente a ε.
3. **Isla vs. péndulo**: ancho de isla en J frente a predicción, con barras.
4. **Serie de energía**: ΔE/E(t) para todas las corridas.
5. **η_c frente a g y W₀** (si se hace D.4).

---

## J. Limitaciones y resultados negativos

### J.1 Regímenes donde el método puede fallar

- **a₀ ≥ 1e-2**: la partida silenciosa (mapa fijo) ya no es exacta; hay que
  partir del equilibrio autoconsistente (documentado en PREGUNTAS).
- **Tiempos largos**: el ruido de discreción crece como t^0.9 (a₀=1e-3) o
  t^1.9 (a₀=1e-2) y fija un t_max por resolución.
- **Modos muy cerca del borde**: la cuadratura en J (Nrc=400) tiene una
  resolución finita en frecuencia; modos a distancia menor que el espaciado
  no se resuelven.
- **ε muy pequeño**: la señal cae bajo el piso de ruido; no se puede medir
  el régimen lineal puro indefinidamente.

### J.2 Limitaciones del modelo L=L₀

- No hay dispersión en L: no hay mezcla en L (que en el sistema completo
  domina, según `vlasov_intro.tex` §"Phase mixing en dos direcciones").
- Solo modo l=0 (respiración radial monopolar); no hay asimetrías.
- No hay relajación por dos frecuencias ni resonancias L–r.

### J.3 Conclusiones actuales demasiado amplias

- "η_c≈1.1–1.2" no es universal (B.3).
- "la teoría lineal no predice la pérdida" no está probado (B.1).
- La conexión con Hadžić es motivación, no equivalencia (B.5).

---

## K. Acción priorizada

### Nivel 1 — Esencial (sin esto no es publicable)

1. **Acotar γ lineal del matrix pencil.** *Por qué*: D3 depende de ello.
   *Resultado*: resolución demostrada y cota firme. *Parte del paper*:
   método y D3. *Dificultad*: baja (test sintético). *Depende de*: nada.
2. **Barrido ε de D3** (D.1). *Por qué*: distingue no lineal de lineal.
   *Resultado*: mecanismo probado o descartado. *Paper*: D3.
   *Dificultad*: baja. *Depende de*: 1.
3. **Convergencia Δt y N para D3/D5/D6** (D.2). *Por qué*: sin esto los no
   lineales no valen. *Resultado*: parte fina N⁻¹ᐟ². *Paper*: todos los no
   lineales. *Dificultad*: media. *Depende de*: nada.
4. **Series de energía** (F.1). *Por qué*: control del integrador.
   *Resultado*: deriva acotada. *Paper*: método. *Dificultad*: baja.
   *Depende de*: nada.
5. **Corregir el error del factor en Q del h_k** (B.9). *Paper*: batería.
   *Dificultad*: baja.

### Nivel 2 — Importante (sustancialmente más fuerte)

6. **Ajuste de isla con modelo de péndulo** (F.4). *Resultado*: prueba
   cuantitativa de atrapamiento. *Paper*: D6. *Dificultad*: media.
7. **Barrido η a ε fija** (D.3). *Resultado*: valida ν como parámetro.
   *Paper*: sección no lineal. *Dificultad*: media.
8. **Robustez de η_c** (D.4). *Resultado*: universalidad o ley de escala.
   *Paper*: transición. *Dificultad*: media.
9. **D8 con ε óptimo y sustracción** (B.4). *Resultado*: ventana ampliada.
   *Paper*: dicotomía. *Dificultad*: media.

### Nivel 3 — Futuro

10. **Barrido L₀ en solver lineal** (E.4). *Resultado*: η_c(L₀).
    *Dificultad*: alta.
11. **Continuación analítica** (G.7). *Resultado*: polo vs. transitorio.
    *Dificultad*: alta.
12. **Levantar L fijo** (proyecto separado, ya en PREGUNTAS).

---

## L. Estructura propuesta del artículo final

1. Introducción: problema (autogravedad vs. mezcla), pregunta (η_c y
   régimen no lineal), predicciones falsables.
2. Modelo: reducción a L fijo, isócrono, variables ángulo-acción,
   equilibrio autoconsistente, perturbación.
3. Teoría lineal: banda, mezcla, modo colectivo, η, dicotomía de borde.
4. Método: PIC, validación, diagnóstico h_k, solver lineal, partida
   silenciosa.
5. Validación: PIC vs. lineal (frecuencia, tasa, modos discretos); energía.
6. Resultados lineales: D1, D2, D4, D8; colas y escalas temporales.
7. Resultados no lineales: D5, D6 (atrapamiento), D3 (fragilidad) con
   convergencia y barrido ε.
8. Discusión: η como parámetro, límites de L fijo, relación con Hadžić,
   régimen de O'Neil gravitatorio.
9. Conclusiones separadas en demostrado/interpretado/hipótesis.
10. Apéndice: ruido de discreción, convergencia, reproducción.

---

## M. Evaluación final

**Veredicto: investigación relevante y correctamente diseñada; aparato
metodológico excelente; resultados no lineales interesantes pero no
publicables en su estado actual.**

El punto fuerte es la infraestructura de validación y el diseño
"teoría lineal como referencia + PIC solo donde la lineal falla". El punto
débil es que los tres resultados no lineales (D3, D5/D6, D8) no tienen el
respaldo de convergencia y cota que el estándar exige, y uno de ellos (D3)
depende de una cota γ no demostrada.

Si se completan las cinco acciones esenciales (K, Nivel 1), el trabajo
pasa de "demo convincente" a "contribución sólida" sobre la transición
autogravedad–mezcla en la reducción de L fijo. Si además se hace el ajuste
de isla con el péndulo (K, Nivel 2), el atrapamiento quedaría probado
cuantitativamente, no ilustrado.

**La contribución actual es limitada pero real.** El mensaje honesto para
el autor es: el aparato está listo; lo que falta son cuatro corridas bien
elegidas (D.1, D.2) y una cota γ (K.1), no más análisis de figuras.

**Resultados potencialmente interesantes que el autor puede estar
pasando por alto:**

1. La **fracción de partículas atrapadas** (F.5) como función de ν es un
   observable no reportado y físicamente más directo que la amplitud.
2. La **ley de escala del ancho de isla con ε** (F.4) convierte la figura de
   D6 en un resultado numérico.
3. El **quiebre de pendiente en las colas** (F.7) marca la transición
   mezcla→modo sin depender del ajuste de polos.
4. La **dependencia de η_c con L₀** (E.4) podría ser una ley de escala nueva
   si η no es el único parámetro.