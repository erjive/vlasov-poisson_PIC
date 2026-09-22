# Auditoría de vlasov-poisson_PIC con f = F(r, p_r, t) δ(L − L₀)

**Fecha:** 2026-09-21. **Versión auditada:** rama `clean/comentarios`, commit `a5e7f86` (último
cambio de `src/`: `0652fe8`). **Referencia:** versión original `ff64af7` (2024-05-22) y el código
hermano `VlasovPoisson_PIC_sp` (distribución en L), auditado hoy en
`../VlasovPoisson_PIC_sp/AUDITORIA_FISICA_2026-09-21.md`.

**Reglas.** No se modificó `src/`. Una simulación a la vez. Cada afirmación lleva su prueba y su
cifra. Se compiló una copia aparte con las banderas del Makefile; su desensamblado coincide objeto
por objeto con `objs/` del repositorio (solo `paramfile.o` difiere, por una línea que escribe
`dftype`). Las pruebas unitarias enlazan esos objetos y llaman a las rutinas del código.

**Una corrección durante la auditoría.** Al leer las fuentes filtré las líneas de comentario y con
ellas oculté las directivas `!$OMP`; por eso creí serial el ciclo de `analysish`. Al encontrar una
diferencia entre dos corridas iguales revisé todas las directivas (§C.3). Las conclusiones que
siguen ya lo incorporan.

---

## Resumen

| Categoría | Qué entra |
|---|---|
| **Demostrado matemáticamente** | La reducción por δ(L−L₀): M = 8π²L₀∫∫F dr dp, ρ = (2πL₀/r²)∫F dp, sin factores extra de L, L² ni r²; características ṙ = p, ṗ = L₀²/r³ − ∂Φ/∂r, L̇ = 0; la F reducida cumple Vlasov en (r, p) sin términos adicionales; el código usa esa convención (f = valor de F, masa 8π²L₀Δr_cΔp_c f) |
| **Validado numéricamente** | Fuerza = −∂H/∂r (3·10⁻¹³); órbitas con L₀ > 0 frente a RK4 (órdenes 1.03 / 2.00 / 4.00 / 5.9 y exacto); energía orbital acotada con su orden; L₀ y f constantes; densidad de una F no separable, orden 2 en el interior; colapso frío, orden 2 en Δr; el código actual reproduce al original (8 cifras) con el mismo eps |
| **Consistente, no probado** | Los resultados de mezclado de fase y de Landau con L₀ = 2 (el régimen de producción): los defectos de abajo no los alcanzan o los alcanzan a ~10⁻⁴ relativo en la autogravedad, pero no se repitieron con el esquema corregido |
| **Sospechoso** | Autogravedad cerca del origen o del borde exterior; órbitas con L₀ pequeño a paso fijo |
| **Incorrecto** | L₀ = 0 (masa NaN); `state = gaussian` (NaN); fondos `iso`, `isotrun`, `nfw`, `burkert` (no actúan sobre las partículas); `sphere` con autogravedad (la borra); `null` sin autogravedad (el centrífugo se acumula); reflexión sin cambio de signo de la fuerza; `reduce_arrays` en la evolución (fuerza sin inicializar); nodos no ligados que quedan en la simulación (h_k NaN) |

**Para el artículo** (isócrono, L₀ = 2, `aa`/`aa_quad`, a0 ≤ 10⁻², r_max = 20): las partículas no
bajan de r ≈ 2.6 ni cruzan el origen, así que la reflexión, las imágenes y la interpolación en el
origen no las alcanzan. Lo que sí las alcanza es el error del Runge–Kutta de Poisson, ~(Δr/r)²:
4·10⁻⁴ de la autogravedad a r = 20Δr, es decir ~10⁻⁶ de la fuerza total con a0 = 10⁻².

---

## A. Modelo matemático

Unidades G = M_iso = b = 1; m ≡ 1 en las ecuaciones de movimiento; la masa del gas es `a0`.

    ∂_t F + p ∂_r F + (L₀²/r³ − ∂_rΦ) ∂_p F = 0,        (Vlasov reducida)
    Φ = Φ_ext + Φ_gas,  Φ_ext = −1/(1+√(1+r²)) (isócrono),
    (1/r²) d/dr(r² dΦ_gas/dr) = 4πρ  ⇔  dM/dr = 4πr²ρ,  dΦ_gas/dr = M/r²,  g = −M/r²,
    ρ(r) = (2πL₀/r²) ∫ F dp,   M(r) = ∫₀^r 4πs²ρ ds = 8π²L₀ ∫₀^r ∫ F dp ds,
    ṙ = p,   ṗ = L₀²/r³ − ∂Φ/∂r,   L̇ = 0.

Hamiltoniano por unidad de masa: H = p²/2 + L₀²/(2r²) + Φ(r,t).

**Condiciones iniciales:** nodos regulares en (r, p) (`aa`, `aa_halton`), en (J, Q) (`aa_quad`),
aleatorios (`aa_random`) o de archivo (`checkpoint`); f se normaliza a masa `a0`.
**Frontera en r = 0:** simetría F(r,p) = F(−r,−p), puntos fantasma con ρ par y Φ′ impar,
Φ ≈ Φ₀ + (2π/3)ρ₀r²; reflexión de las partículas que pasan a r < 0.
**Frontera exterior:** Φ + rΦ′ = 0 en el último nodo.

## B. Reducción desde la distribución completa

La medida del espacio de fases en simetría esférica es d³x d³v = 8π² L dr dp dL: el r² de d³x
se cancela con el 1/r² de d²v_t, porque v_t dv_t = L dL/r² (derivación completa en la auditoría
de `_sp`, §A.3). Con f = F(r,p,t) δ(L−L₀):

    M    = 8π² ∫∫∫ F δ(L−L₀) L dr dp dL = 8π² L₀ ∫∫ F dr dp,
    ρ(r) = (2π/r²) ∫∫ F δ(L−L₀) L dp dL = (2π L₀ / r²) ∫ F dp,
    E    = 8π² L₀ ∫∫ (p²/2 + L₀²/2r² + Φ_ext + ½Φ_gas) F dr dp,
    h_k  = 8π² L₀ ∫∫ F Φ̂_k* e^{−ikQ} dr dp.

**Qué sobrevive:** un factor L₀ (de la medida) y el 8π² (de los tres ángulos: dos de posición,
uno del plano tangencial de la velocidad). No aparecen L², ni r², ni ningún jacobiano adicional. La
degeneración angular (la dirección de **L** uniforme en la esfera) es lo que da la simetría
esférica con |**L**| = L₀ fijo.

**La ecuación reducida.** Sustituyendo f = Fδ(L−L₀) en la ecuación en (r, p, L), el término en
∂_L no existe porque L̇ = 0, y la medida L dr dp dL es invariante (divergencia nula, L constante).
F cumple la ecuación de Vlasov bidimensional con L = L₀, sin términos extra.

**La convención del código.** `f(j)` es el valor de F en el nodo, no una densidad respecto de
dr dp: la masa de la partícula es **m_j = 8π² L₀ Δr_c Δp_c f_j** (`initial_data.f90:93` y
equivalentes), ρ usa `factor` = 2πL₀Δr_cΔp_c (`density.f90:52,196`), la energía
8π²L₀Δr_cΔp_c (`energy.f90`), h_k 8π²L₀Δr_cΔp_c (`analysish.f90`). Coincide con el manuscrito:
F₀ ≈ (1/8π²L₀) Σ N_j S(r−r_j) δ(p−p_j), con N_j ↔ 8π²L₀Δr_cΔp_c f_j.

**El límite L₀ → 0.** Con F finita la masa 8π²L₀∫F se anula. Un problema puramente radial con
masa finita exige que F crezca como 1/L₀, es decir, trabajar con la densidad reducida
𝔉 ≡ 8π²L₀F (masa por unidad de dr dp), con ρ = (1/4πr²)∫𝔉 dp, que tiene límite. El código
intenta tratar L₀ = 0 aparte (`factor = 1` en `density.f90:45-48`, `energy.f90`,
`gaussian_fixedL` con otra normalización), pero `initial_data` normaliza dividiendo por L₀ y
`analysish` multiplica por L₀. Resultado medido: masa inicial NaN, h_k NaN (E1).

**Verificación numérica de la reducción (prueba 5).** F(r,p) = g(r)(1 − p²/s(r)²)² con
s(r) = 0.2 + 0.3r (no separable), g elegida para que ρ = ρ₀(1−r²)³, L₀ = 2. Densidad depositada
y campo frente a los analíticos en 0.2 < r < 0.8:

| (N_r, N_p), Δr | error de ρ / ρ₀ | error de g | primer nodo | masa en el campo |
|---|---|---|---|---|
| (50, 20), 0.04 | 2.36·10⁻³ | 9.3·10⁻³ | 29 % | 1.00134 |
| (100, 40), 0.02 | 5.65·10⁻⁴ | 2.2·10⁻³ | 29 % | 1.00036 |
| (200, 80), 0.01 | 1.38·10⁻⁴ | 5.4·10⁻⁴ | 29 % | 1.00009 |
| (400, 160), 0.005 | 3.43·10⁻⁵ | 1.3·10⁻⁴ | 29 % | 1.00002 |

Orden 2.0 en el interior: un factor L₀ duplicado o ausente daría un error de un factor 2. El
primer nodo no converge (E9, sin imágenes, más el sesgo del volumen geométrico en r₁) y la masa
del campo no es exacta (E8).

## C. Correspondencia matemática → código

| Ecuación | Implementación | Verificación |
|---|---|---|
| m_j = 8π²L₀Δr_cΔp_c f_j | `initial_data.f90` (normalización a `a0`) | prueba 5 |
| ρ = (2πL₀/r²)∫F dp, forma de celda | `density.f90:216-246` (`avg_density`): Σ f W_n /(r²Δr + Δr³/12) · 2πL₀Δr_cΔp_c | prueba 5 |
| Poisson | `poisson_rk.f90`: RK2 sobre (Φ, Φ′), arranque (2π/3)ρ₀r₁², Φ(r_N) = −M/r_N | pruebas 6 |
| g(r_j) | `poisson_rk.f90:170-190`: Σ_{j=1}^{N_r} W_n g_j | E10 |
| ṗ ∋ L₀²/r³ | `grav_force.f90:158-163`: `Lfix²·r/(r²+eps²)²` | −∂H/∂r a 3·10⁻¹³ |
| Fondo isócrono | `grav_force.f90:88-110` | −∂H/∂r a 3·10⁻¹³ |
| ṙ = p, kicks | `main.f90` (euler, leapfrog KDK, yoshida4/6 por composición, analytic) | prueba 2 |
| L̇ = 0 | `Lfix` es un parámetro; solo se asigna al leerlo | lectura (`grep`) |
| F(r,p) = F(−r,−p) | `main.f90:311`: (r,p) → (−r,−p) | E3 |
| E | `energy.f90`: 8π²L₀Δr_cΔp_c Σ(p²/2 + pot − ½pot_self) f | §G |
| h_k | `analysish.f90`: 8π²L₀Δr_cΔp_c Σ f B(J) a_k e^{−ikQ} | igual al original ×8π²L₀ a 10⁻⁸ |

**¿Se trata L como coordenada?** No. No hay malla, depósito, interpolación ni integración en L;
`Lfix` es un escalar. No hay operaciones redundantes en L.

**Factores geométricos** (todos rastreados):

| Factor | Dónde | Origen |
|---|---|---|
| 8π² | masa, energía, h_k | ángulos de x (4π) y del plano tangencial de v (2π) |
| L₀ | masa, ρ, energía, h_k | medida L dL integrada contra δ(L−L₀) |
| 1/r² en ρ | `density.f90:126` | d²v_t = 2πL dL/r² |
| r²Δr + Δr³/12 | volumen de celda | 4π∫r² dr sobre la celda (volumen geométrico) |
| 4π en Poisson | `poisson_rk.f90` | ∇²Φ = 4πρ |
| L₀²/r³, L₀²/2r² | `grav_force.f90` | ∂_r de L₀²/2r² |
| Δr_cΔp_c | todas las sumas | cuadratura de la F sobre los nodos |

### C.3 Paralelismo

Los ciclos que mueven partículas, los fondos, el depósito (paralelo por nodo) y la lista de
celdas (conteo en dos pasadas por hilo) son deterministas para un número de hilos dado.
`analysish.f90:140` y `energy.f90:57` usan `REDUCTION(+:...)`, cuyo orden de suma cambia entre
corridas. Medido: dos corridas idénticas dan h_k distintos en 7.9·10⁻²³ (2.6·10⁻¹⁵ relativo) con
r, p, f idénticos bit a bit. Solo afecta a diagnósticos.

## D. Auditoría de optimizaciones

**Método.** Los commits de rendimiento (`fc2dbc5`, `d08aefd`, `f7e79ae`, `bb57363`, `285daa2`,
`1d4c9b1`) son de la época de entrada posicional. En vez de correrlos por pares, se comparó el
original `ff64af7` con el actual sobre el mismo caso: 1800 partículas, estado `aa`, L₀ = 2,
leapfrog, 2000 pasos. Al actual se le impuso el `eps` que el original fija, L₀/(10·pmax) = 0.1.

| Caso | Resultado |
|---|---|
| Sin autogravedad | r, p, f de las 1800 partículas idénticos en las 8 cifras impresas en t = 0 … 50; energía idéntica |
| Con autogravedad (a0 = 10⁻²) | idénticos hasta t = 37.5; en t = 50 difieren 8·10⁻⁸ en el 3.4 % de las filas (orden de suma de la lista de celdas) |
| h_k | actual/original = 8π²L₀ = 157.9137 en los cinco modos, a 10⁻⁸ |
| Energía con autogravedad | original 6.0·10⁻³ de deriva (sin el ½), actual 2.4·10⁻⁶ |

Esto demuestra que, en esos caminos (leapfrog, isócrono, `aa`, autogravedad, `analysish`), las
optimizaciones y limpiezas no cambiaron las ecuaciones. **No quedan cubiertos:** yoshida4/6 y
`analytic` (no existen en el original; se validaron contra RK4), `aa_quad`, `aa_halton`,
`aa_random`, `checkpoint`, la salida `raw`.

| Original | Actual | Equivalentes | Condiciones | Riesgo |
|---|---|---|---|---|
| `-r/sqrt(1+r²)*pot²` | igual en la rama sin autogravedad; `-r/(sq(1+sq)²)` con autogravedad | Sí (ℝ) | todo r | las dos ramas difieren en el último bit |
| `L²r/(r²+eps²)²`, eps = L₀/(10 pmax) | igual, eps = 0 por omisión | **No** con eps ≠ 0 | — | cambio de física (correcto): 1.4·10⁻² en r a t = 50 |
| búsqueda O(N_r N) | lista de celdas | Sí (ℝ) | r ≥ 0 | orden de suma |
| Wn/Sn repetidos | una evaluación por par | Sí, misma secuencia | — | ninguno |
| φ̂_k por partícula | a_k fuera del ciclo | Sí (ℝ) | A independiente de la partícula | último bit |
| suma serial de h_k | `REDUCTION` | Sí (ℝ) | — | no determinismo ~10⁻¹⁵ |
| E sin ½ | E con ½Φ_self | **No** | — | corrige la energía |
| h_k sin 8π²L₀ | con 8π²L₀ | **No** | — | corrige la normalización |

Cerca de r = 0: el término L₀²r/r⁴ es la misma expresión en ambos; con eps = 0 diverge en r → 0,
que es el límite físico para L₀ > 0 (barrera). Nada se reescribió de forma que cambie ese límite.

## E. Errores encontrados

**Estado de las correcciones** (rama `fix/auditoria-L0`):

| Id | Estado | Medición |
|---|---|---|
| E3 | corregido | L₀ = 0 frente a RK4, Δt = 0.02 … 0.0025: leapfrog 1.19·10⁻³ → 1.50·10⁻⁵ (orden 2.00); yoshida4 de orden 2 (1.6·10⁻³ → 2.4·10⁻⁵) a orden 4.00 (1.2·10⁻⁹ → 2.0·10⁻¹³); energía de yoshida4 ×16 por mitad de Δt. `dfstudy__king_quad_5000`, `sg__quad` y `landau__L_a1e-2_n400_e0.1` (L₀ = 2): r, p, f y ρ idénticos bit a bit en todas las instantáneas (751 en la de Landau) |
| E6 | corregido | 70 órbitas ligadas y 30 que escapan (L₀ = 1, isócrono, t = 60, dos reducciones): antes de la corrección las ligadas diferían 3–4·10⁻² en r entre la corrida con reducción y la corrida sin ella (con y sin autogravedad, leapfrog y yoshida4); después, idénticas bit a bit en los cuatro casos. Ninguna corrida de `reproducir/` (106) ni de `paper_runs/` (61) usa `reduceparticles = .true.` |


| Id | Severidad | Archivo | Función | Problema | Evidencia | Corrección |
|---|---|---|---|---|---|---|
| E1 | **HIGH** | `initial_data.f90:93,262,327,429,506,552`; `analysish.f90` | `initial_data`, `analysish` | **L₀ = 0 no funciona:** la normalización divide por L₀ (masa NaN) y h_k multiplica por L₀, mientras `density`/`energy` usan `factor = 1`: tres convenciones distintas | corrida: "Initial total mass = NaN", h_k NaN, salida 0 | Trabajar con 𝔉 = 8π²L₀F (masa por dr dp) en todo el código, o abortar con L₀ = 0 |
| E2 | **HIGH** | `paramfile.f90:370`, `initial_data.f90:56` | `initial_data` | `state = gaussian` (el valor por omisión) no entra en ninguna rama: el código espera `gaussian1`. Partículas en r = 0 con f = 0 | corrida: todo NaN, salida 0 | Unificar el nombre y abortar ante un estado sin rama |
| E3 | **HIGH** (**corregido**) | `main.f90:311` | bucle principal | La reflexión (r,p) → (−r,−p) **no invierte la fuerza**, que se usa en el primer medio kick del paso siguiente | L₀ = 0 frente a RK4: yoshida4 cae a orden 2 (1.6·10⁻³ → 2.4·10⁻⁵), peor que leapfrog; energía ∝ Δt² | `force_part = −force_part` al reflejar (como `_sp`, `a81f4cf`) |
| E4 | **HIGH** | `grav_force.f90:113-131` | `grav_force` | `iso`, `isotrun`, `nfw`, `burkert` escriben solo la malla: **las partículas no sienten el fondo**; con autogravedad además sobrescriben Φ_gas en la malla | partícula en r = 1.5: F = +0.074 (solo centrífugo) en vez de −1.9 … −2.2 | Evaluar en las partículas, en \|r\|, y sumar (como `_sp`, `8f079fd`, `e2de4f2`) |
| E5 | **HIGH** | `grav_force.f90:57-62, 64-76` | `grav_force` | `null` sin autogravedad no reinicia la fuerza (el centrífugo se acumula en cada llamada); `sphere` **asigna** en vez de sumar y borra la autogravedad | F = 0.074, 0.148 en llamadas sucesivas; sphere + autograv.: −0.3704 frente a −0.3770 | Reiniciar y sumar |
| E6 | **HIGH** (**corregido**) | `utils.f90:1061-1138`, `main.f90:445` | `reduce_arrays` | Tras reducir en la evolución, `force_part` queda reasignado **sin inicializar** y no se recalcula: el primer medio kick del paso siguiente usa memoria no inicializada | lectura (`utils.f90:1122`, `main.f90:445`); magnitud no medida | Llamar a `grav_force` después (como `_sp`) |
| E7 | MEDIUM | `initial_data.f90:237-259` (y `aa_halton`) | estado `aa` | Nodos no ligados o NaN van a r = 10⁴ con f = 0, pero `reduce_arrays` solo los quita si queda < 90 %; `analysish` no excluye no ligados | 1 nodo no ligado entre 1800: **todos los h_k NaN** toda la corrida, salida 0 | Quitar siempre esos nodos; excluir E ≥ 0 en `analysish` |
| E8 | MEDIUM | `poisson_rk.f90` | `poisson_rk` | RK2 sobre (Φ, Φ′): pierde la masa de las partículas cerca del origen y la cuenta mal hasta ~5Δr | masa del campo de una partícula: 0 en r₁, 0.60 en Δr, **1.21** en 1.5Δr, 1.007 en 5Δr, 1.0004 en 20Δr; masa total 1.0013 con Δr = 0.04 | Integrar la masa encerrada (`_sp`, `5c0e07f`) |
| E9 | MEDIUM | `density.f90` | `deposit` | Sin imágenes en −r_j: la parte del peso que cae en los fantasmas se pierde | hasta 32 % de la masa de una partícula en r = Δr/4 | Sumar W(r_i + r_j) |
| E10 | MEDIUM | `poisson_rk.f90:170-190` | `poisson_rk` | La interpolación solo usa los nodos 1…N_r: sin espejo en el origen ni solución exterior. g(0⁺) ≠ 0 y más allá de r(N_r) no hay autogravedad | fuerza ~39 veces la exacta en r < 1.5Δr; F = 0 en r = 20.05 (debía ser −2.5·10⁻³) | Nodos espejo y solución exterior (`_sp`, `a81f4cf`) |
| E11 | MEDIUM | `utils.f90:521-569` | `set_timestep` | El paso no resuelve el pericentro de L₀ pequeño | Δt = 0.01: L₀ = 10⁻³ da \|ΔE/E\| = 3.5 (leapfrog), 4·10⁷ (yoshida4); L₀ ≥ 0.1: ≤ 2.5·10⁻⁵ | Cota Δt ≲ η r_p²/L₀ |
| E12 | MEDIUM | `utils.f90:318-359` | `invert_QJ_to_rp` | Newton de Kepler desde η₀ = Q, sin salvaguarda: diverge para e ≳ 0.98 (misma rutina que `_sp`, D3) | auditoría de `_sp`, U6 | Arranque robusto + bisección |
| E13 | LOW | `utils.f90:366-397`, `analysish.f90`, `invert_QJ_to_rp` | mapa AA | Radicandos sin acotar (órbitas circulares) y sin exclusión de no ligados | lectura; igual que `_sp` antes de `e2de4f2` | `max(...,0)`; una sola rutina |
| E14 | LOW | `density.f90:134-146` | `density` | `vlasov_rhomix` suma f·¼π/(r₂² − r₁²): le faltan 8π²L₀Δr_cΔp_c y el volumen es r₂² − r₁² en vez de (4π/3)(r₂³ − r₁³). El diagnóstico no es una densidad | lectura | Como `_sp` (`90ff2fa`) |
| E15 | LOW | `initial_data.f90:82-83, 210-211` | `gaussian1`, `aa` | Nodos en p en el borde derecho de la celda (p_j = p_min + jΔp_c); en `gaussian1` también r desplazado una celda completa ((i+½)Δr_c, el último fuera de la caja) | lectura | Puntos medios |
| E16 | LOW | `utils.f90:56-72` | `set_grid_size` | `int()` trunca N_r; con r_max = 20, Δr = 0.1 el último nodo es 19.95 < r_max | t_borde | Redondeo con tolerancia |
| E17 | INFO | `analysish.f90:140`, `energy.f90:57` | — | Reducciones OpenMP no deterministas (h_k y E, ~10⁻¹⁵) | §C.3 | Sumas por hilo (`_sp`, `e6d2011`) |
| E18 | INFO | `parameters.f90`, `grav_force.f90` | — | `eps` suaviza el centrífugo si el usuario lo fija ≠ 0; la dinámica deja de ser la del isócrono. Por omisión 0 | lectura | Mantener 0 |
| E19 | INFO | `initial_data.f90` | — | Estados muertos (`gaussian2`, `Plummer`, `compact`, `compact2`, `other3`) que `paramfile` no admite | lectura | Quitarlos |
| E20 | INFO | `BUGS_TODO.md` | — | Dice que `eps = Lfix/(10 pmax)` está activo: ya no (eps = 0 por omisión desde `23e03de`) | lectura | Actualizar |

**Autofuerza:** el campo de la malla incluye el de la propia partícula (−m/2r²) y no se resta. Es
lo físicamente correcto (auditoría de `_sp`, D1). Aquí no hay nada que corregir.

## F. Tests analíticos

| Test | Resultado esperado | Resultado obtenido | Error | Estado |
|---|---|---|---|---|
| 1. L₀ = 0: masa y diagnósticos | masa `a0`, h_k finito | masa NaN, h_k NaN | — | **FALLA (E1)** |
| 1. L₀ = 0: órbita radial por el centro | órdenes 2 / 4 | leapfrog 2.0; **yoshida4 2.0** (1.6·10⁻³ → 2.4·10⁻⁵) | — | **FALLA (E3)** |
| 2. L₀ = 1, 30 órbitas del isócrono frente a RK4 | órdenes 1 / 2 / 4 / 6, `analytic` exacto | 1.03 / 2.00 / 4.00 / 5.87; `analytic` 2·10⁻¹³ | — | PASA |
| Hamiltoniano: F = −∂H/∂r | igualdad | 2.8·10⁻¹³ (isócrono), 6.3·10⁻¹³ (sphere), L₀ = 0.7 | — | PASA |
| 3. Energía orbital, potencial fijo | ∝ Δt^k, acotada | ×4, ×16, ×64 por mitad de Δt | 2.4·10⁻¹⁴ (yoshida6) | PASA |
| 3. Energía, L₀ → 0 a Δt fijo | acotada | 3.5 (L₀ = 10⁻³, leapfrog) | — | **FALLA (E11)** |
| 4. L(t) = L₀ | exacta | `Lfix` no se asigna fuera de la lectura | 0 | PASA |
| 5. Densidad de una F no separable | ρ₀(1−r²)³ | orden 2.0 en 0.2 < r < 0.8; primer nodo 29 % | 3.4·10⁻⁵ | PASA en el interior; falla en el origen (E9) |
| 6. Poisson, ρ₀(1−r²)³ cuasi-continuo | Φ, g cerrados | g en la malla **orden 1** (3.0·10⁻² → 3.8·10⁻³), dominado por el origen; orden 2 en el interior | — | PARCIAL (E8) |
| 6. Esfera uniforme, interior | g = −r | orden 1 (4.7·10⁻³ → 5.9·10⁻⁴) | — | PARCIAL (volumen geométrico + E8) |
| 6. Partícula cerca del origen | g → 0 | g ≈ 39 veces la exacta | — | **FALLA (E10)** |
| 6. Partícula fuera de la malla | −m/r² | 0 | 100 % | **FALLA (E10)** |
| Colapso frío (L₀ = 10⁻⁴) | cicloide del continuo | mediana 8·10⁻⁵ (N = 800), sin mejora en N; orden 2 en Δr | — | PASA con piso de malla (90× el de `_sp`) |
| Fondos `iso`, `isotrun`, `nfw`, `burkert` | F = fondo + centrífugo | solo centrífugo, acumulándose | — | **FALLA (E4)** |
| Estado por omisión | partículas válidas | NaN | — | **FALLA (E2)** |

## G. Conservación

| Magnitud | Esperado | Medido |
|---|---|---|
| L₀ | exacta | constante por construcción (parámetro) |
| Masa (Σ f) | exacta | `f` idéntico bit a bit en todas las instantáneas de tres corridas (751 en la de Landau) |
| Masa vista por Poisson | exacta con toda la masa en la malla | error O(Δr²): 1.3·10⁻³ con Δr = 0.04 (E8); nada más allá de r(N_r) (E10) |
| Energía, fondo fijo | acotada | `dfstudy__king_quad_5000`, t = 2000: 9.3·10⁻¹² |
| Energía, autogravedad | error de discretización | Landau, t = 3000, 751 instantáneas: 2.0·10⁻⁷, igual en las dos mitades (sin deriva secular); `sg__quad`, t = 2000: 1.7·10⁻⁷ entre t = 0 y t = 2000 (solo dos instantáneas: no prueba conservación intermedia) |
| Energía, colapso frío | → 0 con Δr | 1.1·10⁻³, 3.0·10⁻⁴, 8.0·10⁻⁵, 2.1·10⁻⁵ para Δr = 0.04 … 0.005 (orden 2): error espacial, no temporal |

Clasificación de las derivas: con fondo fijo, error de integración temporal (escala con Δt^k);
con autogravedad, error espacial del esquema PIC (escala con Δr², no con Δt); el no determinismo
de ~10⁻¹⁵ es punto flotante. No se encontró deriva secular atribuible a un bug en las corridas de
producción.

## H. Convergencia

| Variable | Experimento | Orden observado |
|---|---|---|
| Δt | 30 órbitas, L₀ = 1 | euler 1.03, leapfrog 2.00, yoshida4 4.00, yoshida6 5.87 |
| Δt, L₀ = 0 | órbita por el centro | leapfrog 2.0; yoshida4 **2.0** (E3) |
| Δr_c, Δp (muestreo) | prueba 5, refinamiento conjunto | 2.0 en el interior |
| Δr (Poisson, máximo) | perfil (1−r²)³ | **1.0** (origen); 2.0 en el interior |
| Δr (colapso) | N = 800, dt fijo | 1.9 |
| N (colapso) | Δr = 0.01 | ≈ 0 (piso de malla) |

## Simulaciones del repositorio

Se repitieron `08_verificacion/dfstudy__king_quad_5000`, `09_autogravedad/sg__quad` y
`11_landau/landau__L_a1e-2_n400_e0.1` con el binario actual. Reproducen los h_k guardados en
`exe/` a 2.6·10⁻¹⁵, 8.4·10⁻¹⁶ y 3.2·10⁻¹⁵ relativo. Eso está al nivel del no determinismo de
las reducciones (§C.3): dos corridas del mismo binario ya difieren en esa cantidad. No se buscaron
ni se vieron calentamiento, amortiguamiento o crecimiento artificiales en estas tres corridas más
allá de lo que explican las derivas de §G. No se repitieron las 106 corridas de `reproducir/`.

## I. Conclusión

**Demostrado matemáticamente.** La reducción por δ(L−L₀) y los factores que sobreviven (8π²L₀ y
1/r²); las características con L₀²/r³; la convención del código para f, idéntica a la del
manuscrito.

**Validado numéricamente.** Con L₀ > 0 y fondo isócrono o esfera sin autogravedad, la dinámica es
exactamente la del Hamiltoniano (−∂H/∂r a 10⁻¹³, órdenes 1/2/4/6, energía acotada). La densidad y
el campo de una F conocida convergen con orden 2 lejos del origen. La masa y L₀ se conservan
exactamente. El código actual reproduce al original con el mismo eps: las optimizaciones no
cambiaron las ecuaciones en los caminos comparados.

**Consistente, no probado.** Los resultados de producción con L₀ = 2 (mezclado y Landau): en ese
régimen solo actúa E8, a ~10⁻⁴ de la autogravedad; no se repitieron con el esquema corregido.

**Sospechoso.** Cualquier corrida con autogravedad y masa cerca del origen o del borde, con
L₀ pequeño, con `reduceparticles = .true.`, o con cortes que dejan < 10 % de nodos fuera.

**Incorrecto.** E1 (L₀ = 0), E2 (estado por omisión), E3 (reflexión), E4 y E5 (fondos), E6
(`reduce_arrays`), E7 (nodos no ligados), E10 (interpolación en el origen y fuera de la malla).
Casi todos ya están corregidos en `VlasovPoisson_PIC_sp` y documentados como pendientes de portar
en `BUGS_TODO.md` (sección "Pendiente: consistencia del acoplamiento partícula-malla"); aquí quedan
medidos.

**Recomendación.** Portar de `_sp` el acoplamiento (`a81f4cf`), la integración de la masa
(`5c0e07f`), los fondos (`8f079fd`, `e2de4f2`), la exclusión de no ligados (`732250f`) y las sumas
deterministas (`e6d2011`); resolver L₀ = 0 con 𝔉 = 8π²L₀F; y repetir con eso la corrida de
referencia y la de Landau antes de citarlas en el artículo.
