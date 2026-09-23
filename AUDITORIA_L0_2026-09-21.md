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
| ρ = (2πL₀/r²)∫F dp, forma de celda | `density.f90:216-246` (`avg_density`): Σ f W_n /(r²Δr + Δr³/12) · 2πL₀Δr_cΔp_c (desde E21 la densidad de salida divide por Δr(r² + (n+1)Δr²/12)) | prueba 5 |
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
| r²Δr + Δr³/12 | volumen de celda | 4π∫r² dr sobre la celda (volumen geométrico); desde E21 la salida usa el que cubre W_n, Δr(r² + (n+1)Δr²/12) |
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
| E20 | corregido | Marcadas como desactualizadas en `BUGS_TODO.md` las dos notas que ya no describen el código: `eps` no se deriva de `pmax` desde `23e03de` (es un parámetro de entrada con valor 0 por omisión, y `paramfile` avisa si se pide `eps /= 0` con L₀ /= 0), y `analysish` sí multiplica h_k por 8π²L₀. Se deja el texto original con la corrección fechada, porque el archivo es una bitácora |
| E19 | corregido | Quitados los estados muertos de `initial_data` (`gaussian2`, `Plummer`, `compact`, `compact2`, `other3`), 112 líneas que `paramfile` no admite desde su lista de estados válidos. Los seis estados vivos (`gaussian`, `aa`, `aa_halton`, `aa_quad`, `aa_random` y `checkpoint`), en t = 0: idénticos bit a bit. Compila sin avisos |
| E17 | corregido | `energy` y `analysish` acumulan por hilo y suman las partes en orden de hilo, con reparto `SCHEDULE(STATIC)` para que cada corrida dé a cada hilo las mismas partículas (como `_sp`, `e6d2011`). La reducción de OpenMP combinaba las partes en el orden en que terminaban los hilos. Tres corridas idénticas de `sg__quad` con 4 hilos: antes, 9 y 5 de 11 energías distintas entre corridas y solo 1 de 5 archivos `.tl` idénticos; ahora, las tres idénticas bit a bit (11/11 energías y 5/5 archivos). El cambio mueve la energía 8.6·10⁻¹⁶ relativo respecto del binario anterior. Sigue dependiendo del número de hilos (h_k cambia 3·10⁻¹⁵ entre 1, 2 y 4), que es inevitable al repartir la suma; lo que se gana es reproducibilidad a número de hilos fijo. Costo, ABBA con 4 hilos: 11.98, 11.74, 12.43, 12.68 s, dentro del ruido |
| E22 | corregido | `Sn(1)` vale 1/2 en \|y\| = 1/2 (la media de los dos límites laterales) y `density` ajusta y a ±1/2 cuando está a menos de 8ε·max(1, r_i/Δr_c) de la cara, que es el redondeo de y = (r_i − r_j)/Δr_c. Una partícula sobre una cara de celda, medida como Σ ρ_i 4πr_i²Δr con n = 1 y Δr_c = Δr: antes depositaba **0** de su masa en r_j = 19Δr, 1 en 1Δr, 10Δr y 28Δr, y **2** en 37Δr (y 2 en r_j = 0, contando su imagen); el resultado dependía del último bit. Ahora deposita 1.000 en los cinco radios, a ±2 ulp de la cara, y también en el origen. Solo afecta a la densidad `rho` que se escribe (Poisson usa `avg_rho`, con W_n, que es continua). `sg__quad` (2000 pasos) y `dfstudy__king_quad_5000`: idénticas bit a bit |
| E16 | corregido | `set_grid_size` redondea (rmax − rmin)/Δr en vez de truncarlo, y avisa si no es entero. El cociente casi nunca es exacto en binario, así que la malla tenía una celda de más o de menos según el redondeo: con r_max = 20 y Δr = 0.1 leídos del archivo daba N_r = 201 y el último nodo en 20.05, media celda más allá de r_max, mientras que los mismos números escritos como literales de precisión simple daban 19.95. Con el cambio, N_r = 200 y el último nodo en 19.95, y la malla cubre [0, r_max] exacto. En `sg__quad` (con autogravedad) y `dfstudy__king_quad_5000`: partículas, energías y h_k idénticos bit a bit, y solo los arreglos de malla tienen un punto menos (la celda extra no contenía masa, y la condición de contorno exterior fija el mismo potencial en el interior) |
| E14 | corregido | `vlasov_rhomix` es ahora la densidad media de la cáscara: la masa de las partículas con r₁ ≤ r ≤ r₂, 8π²L₀Δr_cΔp_c Σf, sobre el volumen (4π/3)(r₂³ − r₁³) (como `_sp`, `90ff2fa`). Antes sumaba f/(4π(r₂² − r₁²)), sin la medida y con un área en vez del volumen. `sg__quad` con la cáscara [5, 8] (la de la corrida, [2, 3], está vacía): el valor pasa de 3.534·10⁻⁷ a 9.279·10⁻⁹ en t = 0, un factor constante 0.026258 = 8π²L₀Δr_cΔp_c·3(r₂² − r₁²)/(r₂³ − r₁³), y coincide con el mismo cálculo hecho en Python desde las partículas guardadas a 9·10⁻¹⁰ (las 8 cifras que escribe el archivo), en las 11 instantáneas. El resto de la corrida, idéntico bit a bit (132 datasets y atributos, los otros 4 `.tl`) |
| E15 (`aa`, `aa_halton`) | corregido | Los nodos en p van en el punto medio de la celda, como en r: `aa` los tenía en el borde derecho (p_j = p_min + jΔp_c) y `aa_halton` centraba su desplazamiento en ese borde, media celda corrida (`compact` y `compact2` tienen lo mismo, pero son estados muertos, E19). (a) Caja que corta a F en p (p ∈ [−0.15, 0.15], sin corte en f, N_p = 10 … 320): antes, momento radial neto espurio Σfp/Σf\|p\| = 8.3·10⁻⁴ … 2.4·10⁻⁵, de orden 1 en Δp (la malla incluía p_max y no p_min); ahora 0 a redondeo (≤ 6·10⁻¹⁶). Las cantidades pares (E₀, K₀, \|h_k(0)\|) convergen con orden ≈ 2 en ambos casos, ahora con constantes ≈ 2 veces menores (la regla del extremo coincide con la del trapecio para funciones pares). (b) Las 16 corridas de `paper_runs/` con `aa` y `aa_halton` (caja [1, 10] × [−0.6, 0.6], donde F ≈ 0 en los bordes, así que antes la malla ya era simétrica salvo un nodo con F ≈ 0 y Σfp ≈ 10⁻¹⁶): cambian al nivel del error de cuadratura y del corte al 1 % de f_max. En t = 0, \|h_k\| cambia 1.5·10⁻³, 5·10⁻⁴, ≈ 5·10⁻⁵ y ≈ 7·10⁻⁶ relativo con N ≈ 10², 10³, 10⁴ y 10⁵; E₀ entre 7.6·10⁻⁵ y 2·10⁻⁸; el número de partículas que sobrevive al corte, ±0.1–4 %; sin corte (`cal_aa_nocut`), 2·10⁻⁷. En la evolución (`cal_aa_e2` y `cal_aa_e3`, t ≤ 300) la diferencia no crece: a lo sumo 5·10⁻⁴ y 3·10⁻⁵ de max\|h_k\| |
| E12, E13 | corregido | Una sola función `kepler_eta(Q, e, tol)` en `utils` resuelve Q = η − e sen η para `invert_QJ_to_rp` (integrador `analytic`) y el estado `aa_quad`: primero el Newton de siempre (mismas operaciones); se acepta si salió por tolerancia con un último paso \|g/g′\| ≤ √tol y η ∈ [Q − e, Q + e], y si no, se resuelve con Newton acotado y bisección. Radicandos con `max(·, 0)` y fase de órbita circular (s₁ = s₂) acotada en `invert_QJ_to_rp`, `aa_quad`, los mapas directos de `aa`, `aa_halton` y `aa_random` e `init_action_angle`, como en `analysish`; `init_action_angle` aborta si hay partículas no ligadas. (a) Kepler, 2000 ángulos por e: el Newton sin salvaguarda fallaba desde e = 0.98 (6/2000) hasta ~50/2000 cerca de e = 1, con residuos de hasta 10²³; ahora ninguna falla y residuo ≤ 8.9·10⁻¹⁶ para todo e ≤ 1. Donde el viejo convergía, η idéntico bit a bit salvo 5 de ~21 500 casos (e ≥ 0.98, último paso > √tol, ahora resueltos por la vía segura). (b) Ida y vuelta (Q, J) → (r, p) → (Q, J) con L₀ = 0.25: antes, desde J = 10 (e = 0.992), 29–53 de 2000 fases erradas hasta \|ΔQ\| ≈ π con J correcto; ahora \|ΔQ\| ≤ 6·10⁻¹³ hasta J = 10⁵ (e = 1 − 10⁻¹⁰). (c) Órbita circular exacta (J = 0) con L₀ = 0.25: `invert_QJ_to_rp` daba r = 0; ahora r_c = 0.80410872. Mapa directo a menos de 10⁻⁸ relativo de r_c: Q = NaN en 375, 909 y 783 de 3003 puntos (L₀ = 0.25, 1, 2) con las fórmulas viejas, en ninguno con las protegidas. Partícula no ligada con `analytic`: antes terminaba con código 0, r = 0 y energía infinita; ahora aborta con código 1. Neutralidad: estados iniciales de las 159 configuraciones guardadas (76 `aa_quad`, 83 `aa`, `aa_halton` y `aa_random`, estas últimas con la misma semilla en ambos binarios porque sin `seed` toman una del reloj), idénticos bit a bit; las dos corridas `analytic` de `paper_runs/` completas (6 instantáneas y los 5 `.tl`), idénticas. Costo del integrador `analytic` (ABBA): 14.30 → 14.35 s (una primera versión que volvía a evaluar el residuo costaba +5 %). Observación: `aa_random` no tiene tope de intentos; en una región donde F ≈ 0 (una caja alrededor de la órbita circular, con `dftype = gauss`, F ∝ J² e^(−J/σ)) no termina, con el binario anterior y con el nuevo |
| E11 | corregido | `set_timestep` añade una tercera cota, la condición de Courant con la escala del pericentro: Δt ≤ courant·r_p/v_p = courant·r_p²/L₀, con r_p el menor pericentro de las partículas en el campo de t = 0 (bisección en log r de L₀²/2r² + Φ(r) = E; Φ = fondo cerrado + autogravedad de la malla). `bgpot` y `bgforce` pasan de `grav_force` a `utils` (con el caso del isócrono) para evaluarlo; el traslado es neutral bit a bit (5 fondos y 2 casos con autogravedad). Al arrancar se imprime siempre r_p y Ω_pΔt = L₀Δt/r_p². (a) Ley medida (una órbita en el isócrono): el pico de energía en cada pericentro es max\|ΔE/E\| ≈ 0.03 (Ω_pΔt)^p, con coeficiente 0.031 (leapfrog, p = 2) y 0.032 (yoshida4, p = 4). (b) Esa órbita con courant = 0.5, T = 50, error fuera del pericentro (r > 0.1): L₀ = 10⁻³: leapfrog 3.5 → 1.9·10⁻⁴, yoshida4 4.0·10⁷ → 1.4·10⁻⁴; L₀ = 10⁻⁴: 0.38 → 2.8·10⁻⁴ y 56 → 1.9·10⁻⁴. El Δt elegido es 0.5·r_p²/L₀ con el r_p analítico (1.746·10⁻³). Con courant = 0.25 y 0.125 (L₀ = 10⁻³): 1.7·10⁻⁸ y 3.6·10⁻⁹ (leapfrog), 1.6·10⁻⁹ y 2.0·10⁻¹⁴ (yoshida4). Con L₀ ≥ 10⁻² la cota no se activa y nada cambia. (c) Las 109 corridas de `reproducir/` y `paper_runs/`: Δt idéntico, la cota nunca se activa. `sg__quad`, `dfstudy__king_quad_5000` y la órbita con L₀ = 10⁻³ (cota activa): idénticas bit a bit a las corridas previas a la línea informativa. (d) L₀ = 0.25, el mínimo que se usará: la cota no se activa (órbita: 0.48 frente a Δt = 0.01; `sg__quad`: ≥ 0.75 frente a 0.1) y la cota instantánea queda constante en toda la corrida (el potencial no se hunde: la barrera detiene el colapso). La precisión la fija courant: en `sg__quad` con L₀ = 0.25, Ω_pΔt ≈ 0.26 y max\|ΔE/E\| = 2.6·10⁻⁵ (1.7·10⁻⁷ con L₀ = 2). **Resultado negativo:** con autogravedad el Δt fijo de t = 0 no es una cota si el potencial se hunde. Colapso frío con L₀ = 10⁻⁴: Δt 2.5·10⁻³ → 5.0·10⁻⁵ (hasta t = 0.8 el resultado no cambia, 1.78·10⁻⁵ → 1.80·10⁻⁵ con N = 200, 1.04·10⁻⁶ → 9.3·10⁻⁷ con N = 800, a 45 veces el costo); atravesando el centro (t = 1.5) fallan los dos, max\|ΔE/E₀\| = 32 antes y 1.8·10¹⁷ ahora (85 de 200 capas expulsadas). Entre t = 0 y el rebote (t = 1.108) Φ_min pasa de −1.5 a −43.7 y la cota instantánea de 5.0·10⁻⁵ a 5.9·10⁻⁷ (≈ 68 000 pasos hasta el rebote con Δt por paso). Costo de recalcular Δt en cada paso (ABBA, 1 hilo): con la bisección +1000 % (`sg__quad`, 10⁴ partículas), +3400 % (N5e4, 50 384) y +1000 % (colapso, 800); con la cota cerrada r_p ≥ L₀/√(2(E_max − Φ_min)), +3.3 %, +2.3 % y dentro del ruido. No se implementó (Δt variable deja de ser simpléctico y con L₀ ≥ 0.25 no hace falta) |
| E21 | corregido | Por decisión del usuario, la densidad de salida (`density`) divide la masa depositada por el volumen que cubre W_n, V_i = 4πΔr(r_i² + (n+1)Δr²/12), en vez del de celda. `avg_density`, la entrada de Poisson, escribe ahora en un arreglo propio con las mismas operaciones: la dinámica no cambia y en t = 0 ya no pisa la densidad de salida. (a) Esfera uniforme, capas cuasi-continuas: error en los nodos 1–4 de +0.250/+0.036/+0.013/+0.007 (n = 1), +0.50/… (n = 2) y +0.75/… (n = 3), igual con Δr = 0.04 y 0.01, a 5·10⁻⁵ (n = 1, muestreo) y 10⁻¹⁵ (n = 2, 3). (b) ρ₀(1−r²)³: primer nodo +0.25 → −0.0031 (Δr = 0.04) y −0.00015 (Δr = 0.01), ahora converge con orden 2. **Resultado negativo:** en el interior (0.2 < r < 0.8) el error es 1.7–2 veces mayor (n = 1: 2.0·10⁻³ → 3.5·10⁻³ con Δr = 0.04, 1.25·10⁻⁴ → 2.2·10⁻⁴ con 0.01), sigue de orden 2. El volumen de celda añade al error de suavizado σ²Δr²(2ρ′/r + ρ″/2) el término positivo (σ² − 1/12)Δr²ρ/r², que en un perfil decreciente lo compensa en parte; con ρ ∝ r² (creciente) pasa lo contrario, y V_W da un error interior 10–13 % menor (2.29·10⁻³ → 2.08·10⁻³ con n = 1). (c) Una partícula: Σ ρ_i V_i = m exacto con el volumen nuevo. Neutralidad (1 hilo): `sg__quad` (4000 pasos, 21 instantáneas, malla incluida), `landau__L_a1e-2_n400_e0.1` completa, `dfstudy__king_quad_5000` y colapso frío: todo idéntico bit a bit salvo `avg_rho`, que cambia exactamente en (r² + Δr²/12)/(r² + (n+1)Δr²/12) (a 4·10⁻¹⁵), también en t = 0 con autogravedad |
| E1 | mitigado | Por decisión del usuario, `validate` (`paramfile.f90`) detiene toda corrida con L₀ = 0, con el plan anotado en el código: trabajar en todo el código con 𝔉 = 8π²L₀F (masa por dr dp), que tiene límite finito en L₀ → 0; las ramas `Lfix == 0` de `initial_data`, `density` y `energy` quedan marcadas como inalcanzables. Antes, `dfstudy__king_quad_5000` con L₀ = 0: "Initial total mass = NaN" y código de salida 0; ahora aborta con el mensaje y código 1. Colapso frío con L₀ = 10⁻⁴ y `dfstudy__king_quad_5000` (L₀ = 2): idénticas bit a bit. Ninguna corrida guardada usa L₀ = 0. Las órbitas radiales con L₀ = 0 de las pruebas de E3 y E4 se hicieron antes de este cambio (vía `checkpoint`) y ya no se pueden repetir con el binario actual |
| E3 | corregido | L₀ = 0 frente a RK4, Δt = 0.02 … 0.0025: leapfrog 1.19·10⁻³ → 1.50·10⁻⁵ (orden 2.00); yoshida4 de orden 2 (1.6·10⁻³ → 2.4·10⁻⁵) a orden 4.00 (1.2·10⁻⁹ → 2.0·10⁻¹³); energía de yoshida4 ×16 por mitad de Δt. `dfstudy__king_quad_5000`, `sg__quad` y `landau__L_a1e-2_n400_e0.1` (L₀ = 2): r, p, f y ρ idénticos bit a bit en todas las instantáneas (751 en la de Landau) |
| E6 | corregido | 70 órbitas ligadas y 30 que escapan (L₀ = 1, isócrono, t = 60, dos reducciones): antes de la corrección las ligadas diferían 3–4·10⁻² en r entre la corrida con reducción y la corrida sin ella (con y sin autogravedad, leapfrog y yoshida4); después, idénticas bit a bit en los cuatro casos. Ninguna corrida de `reproducir/` (106) ni de `paper_runs/` (61) usa `reduceparticles = .true.` |
| E2, E15 (`gaussian`) | corregido | La rama se llamaba `gaussian1` y el lector solo admite `gaussian`: antes, todas las partículas en r = 0 con f = 0 y NaN; ahora, nodos en los puntos medios ([1.05, 4.95] × [−0.95, 0.95] en una caja [1, 5] × [−1, 1]), masa = a0, energía 6.8·10⁻¹⁴ y \|h₀\| constante en una caja ligada. Densidad frente a ρ = (2πL₀/r²)∫F dp: orden 1.88, 1.98, 2.00 con dos cáscaras por celda (con media cáscara por celda da orden 1 por muestreo, no por el código). Un estado sin rama ahora aborta. `aa_quad`, `aa`, `aa_halton`, `checkpoint` y `aa_random`: idénticos bit a bit. Nota: con nodos no ligados en la caja, h_k sigue en NaN (E7) |
| E7 (y E13 en `analysish`) | corregido | (a) Un nodo no ligado que quedaba en r = 10⁴ (menos del 10 % descartado): antes h_k NaN; ahora se elimina en t = 0 y h_k es finito. (b) `checkpoint` con 30 órbitas ligadas y una no ligada: antes NaN; ahora h_k es 30/31 del de las 30 solas en las 15 cifras, en todos los instantes y modos. (c) 50 órbitas casi circulares, 6 con radicando negativo por redondeo: antes NaN; ahora finito, cociente 30/80 exacto. Neutralidad, con 1 hilo: las 16 corridas de `paper_runs/` con `aa` y `aa_halton`, idénticas bit a bit en t = 0 (partículas y h_k); evoluciones de `N1e3`, `cal_aa_halton_e1`, `N5e4_selfgrav`, `dfstudy__king_quad_5000` y Landau, idénticas en partículas y en los cuatro archivos de h_k |
| E5 | corregido | (a) `null` sin autogravedad: F(1.5) = 0.07407 en dos llamadas seguidas (antes 0.074, 0.148). Partícula libre con L₀ = 1 frente a r(t)² = (r₀+p₀t)² + (L₀t/r₀)², Δt = 0.02, 0.01, 0.005: antes, errores de 7.95·10², 1.28·10³ y 2.07·10³ (el centrífugo acumulado domina); ahora leapfrog 1.25·10⁻⁴, 3.13·10⁻⁵, 7.82·10⁻⁶ (orden 2.00) y yoshida4 1.13·10⁻⁷, 7.06·10⁻⁹, 4.41·10⁻¹⁰ (orden 4.00). (b) `sphere` con autogravedad: F(1.5) = −0.377030, la suma esperada (antes −0.3704, solo el fondo); la malla lleva también el campo total, como con el isócrono. La esfera se evalúa en \|r\| con fuerza impar (antes r < −1 tomaba la rama interior). Neutralidad, con 1 hilo: `sphere` sin autogravedad (30 órbitas con r₀ ∈ [0.3, 3], L₀ = 1, leapfrog y yoshida4, 41 instantáneas) y `dfstudy__king_quad_5000` (isócrono): idénticas bit a bit, h_k incluido |
| E4 | corregido | Los cinco fondos analíticos (`sphere`, `iso`, `isotrun`, `nfw`, `burkert`) pasan por `add_background`: actúan sobre las partículas, se suman a la autogravedad y se evalúan en \|r\| con fuerza impar (mismas formas que `_sp`). (a) r = ±1.5, L₀ = 0.5: antes F = ±0.07407 y ±0.14815 en dos llamadas (solo el centrífugo, acumulándose); ahora ∓1.92593, ∓0.69216, ∓2.17510 y ∓1.47504, lo esperado, iguales en las dos llamadas. F = −dΦ/dr por diferencia central (h = 10⁻⁴) hasta 4·10⁻⁸, el error de truncamiento. Con autogravedad, F(1.5) coincide con autogravedad + fondo + centrífugo y la malla lleva Φ_gas + Φ_bg exacto (antes le faltaba Φ_gas, 2.7·10⁻²). (b) 30 órbitas con L₀ = 1, r₀ ∈ [0.3, 3], T = 20: antes los cuatro fondos daban lo mismo, las partículas escapaban (r ≈ 3·10⁴) y la energía divergía; ahora quedan ligadas y max\|ΔE/E\| converge con orden 2.00 (leapfrog) y 4.00 (yoshida4) en los cuatro. (c) Órbitas radiales (L₀ = 0) que cruzan el origen, frente a una referencia RK4 por regiones (h = 10⁻⁴, error ≲ 10⁻¹²): `isotrun` y `burkert`, orden 2.00 y 4.00. **Resultado negativo:** `nfw` baja a orden ≈ 1 y `sphere` da órdenes erráticos en las órbitas que cruzan r = 1 (las que no, orden 2.00 y 4.00). El código reproduce exactamente (diferencia 0) un leapfrog sobre la recta con la fuerza impar, así que viene del método de paso fijo: la fuerza de `nfw` salta en r = 0 (−8 → +8) y la derivada de la de `sphere` salta en r = 1. No afecta a L₀ > 0, que no llega al origen. Neutralidad (1 hilo): `sphere` sin autogravedad (las corridas de E5), `sphere` e isócrono con autogravedad (`sg__quad`, 2000 pasos, malla incluida), `dfstudy__king_quad_5000` y `null` (partícula libre, 6 corridas): idénticas bit a bit, energías y h_k incluidos. Ninguna corrida guardada usa estos cuatro fondos |
| E9 | corregido | `density` y `avg_density` suman la imagen de cada partícula en (−r_j, −p_j), W((r_i + r_j)/Δr), con signo opuesto en la corriente (como `_sp`). Una partícula sola, Δr = 0.1, masa en la malla Σ ρ_i 4π(r_i²Δr + Δr³/12) / m: antes 0.50 en r_j = 0 (n = 1, 2, 3), 0.75 / 0.72 / 0.68 en Δr/4 y 0.25 / 0.28 / 0.32 en r_j = −Δr/4 (a mitad de paso); ahora 1.000000000 en todas las posiciones y órdenes. El depósito de (r_j, p_j) y de (−r_j, −p_j) es ahora idéntico (diferencia 0 en ρ, ρ promedio y corriente; antes hasta el 100 %). Capas cuasi-continuas (64 por Δr), esfera uniforme: primer nodo +0.234 / +0.450 / +0.656 → +0.250 / +0.500 / +0.750 (n = 1, 2, 3), independiente de Δr; el resto, igual (orden 2 en 0.2 < r < 0.8). Las imágenes solo aportan la masa cercana al origen, que va como r²; el error del primer nodo es el sesgo del volumen geométrico (E21). Neutralidad (1 hilo): `sg__quad` (2000 pasos), `dfstudy__king_quad_5000` y `landau__L_a1e-2_n400_e0.1` completa (751 instantáneas): idénticas bit a bit, malla, energías y `.tl` incluidos (sus partículas no bajan de r ≈ 3.5). Costo, ABBA con `sg__quad`: 8.26, 7.98, 8.03, 8.22 s, sin aumento |
| E8 | corregido | `poisson_rk` integra la masa encerrada, dM/dr = 4πr²ρ y dΦ/dr = M/r², en forma cerrada sobre cada intervalo con ρ lineal entre nodos (esquema de `_sp`, `5c0e07f`). Para que el campo vea exactamente la masa depositada, la recta pasa por ρ̃_k = m_k / (4πΔr(r_k² + Δr²/6)) (el segundo momento de la recta), calculada a partir de `avg_rho`; la densidad que se escribe sigue con el volumen de celda, así que E21 ya no entra al campo con n = 1. (a) Una partícula, masa vista −F(r_N)r_N²: antes 0 en r_j ≤ Δr/2, 0.60 en Δr, 1.21 en 1.5Δr, 1.0074 en 5Δr y 1.0004 en 20Δr; ahora 1.00000000 en todas las posiciones, n = 1, 2, 3. (b) ρ₀(1−r²)³ cuasi-continuo, n = 1: error máximo de g en r < 0.9 de orden 1 (3.2·10⁻² → 4.1·10⁻³, el primer nodo +25 %) a orden 2.00 (8.97·10⁻³ → 1.41·10⁻⁴, primer nodo → 0); masa vista 1.0013 → 1.0000000000; en 0.2 < r < 0.8, igual que antes (±8 %). Esfera uniforme, n = 1: error de g 5.0·10⁻³ → 2.8·10⁻⁶. Con n = 2 y 3 el primer nodo queda +20 % y +40 % (la recta supone el peso lineal); la masa es exacta. (c) Colapso frío (L₀ = 10⁻⁴) frente a la cicloide, error mediano: N = 200, 400, 800 con Δr = 0.01: 9.7·10⁻⁵, 8.3·10⁻⁵, 8.0·10⁻⁵ (piso de malla) → 1.8·10⁻⁵, 4.5·10⁻⁶, 1.0·10⁻⁶ (orden 2 en 1/N, el de cáscaras con su autogravedad); N = 800 con Δr = 0.04 … 0.005: 1.2·10⁻³ … 2.1·10⁻⁵ → 1.0·10⁻⁶ en todos (sin piso de malla; el nivel de `_sp`). Error de energía a la mitad (8.0·10⁻⁵ → 4.0·10⁻⁵ con Δr = 0.01), sigue de orden 2 en Δr. (d) Corridas de referencia, 1 hilo: `landau__L_a1e-2_n400_e0.1` cambia h_k en 2.3·10⁻⁶ y 4.3·10⁻⁶ relativo al máximo y E₀ en 3.2·10⁻⁷; max\|ΔE/E₀\| 2.0·10⁻⁷ → 1.7·10⁻⁷. `sg__quad` (t = 2000): h_k 1.4·10⁻⁵ y 3.8·10⁻⁵, E₀ 5.9·10⁻⁷. **Resultado negativo:** en `sg__quad` max\|ΔE/E₀\| pasa de 1.7·10⁻⁷ a 5.3·10⁻⁷, por un salto en los primeros 20 unidades de tiempo (5.8·10⁻⁸ → 3.6·10⁻⁷); después la fluctuación (σ = 2.3·10⁻⁸ → 3.1·10⁻⁸) y la deriva son del mismo orden. En `_sp` el mismo cambio dio 1.94·10⁻⁷ → 2.49·10⁻⁷. Revisado con E10 (que deja `sg__quad` idéntica): el salto es error espacial, no temporal. En t = 20 vale 3.580, 3.583 y 3.583·10⁻⁷ con Δt, Δt/2 y Δt/4, y 3.6·10⁻⁷, 1.2·10⁻⁷ y 3.6·10⁻⁸ con Δr, Δr/2 y Δr/4 (orden ≈ 1.6). Con el RK2 era 5.8·10⁻⁸, −1.2·10⁻⁸ y −1.0·10⁻⁸: más chico en Δr = 0.1, pero sin converger. Sin autogravedad `poisson_rk` no se llama: `dfstudy__king_quad_5000`, idéntica bit a bit |
| E10 | corregido | La interpolación a las partículas usa todos los nodos del soporte de W_n: los j ≤ 0 como espejo del nodo 1−j (potencial par, fuerza impar) y los que pasan de r(N_r) con la solución exterior Φ = Φ(r_N)r_N/r, F = F(r_N)(r_N/r)² (esquema de `_sp`, `a81f4cf`). Para una partícula cuyo soporte cae en 1…N_r las operaciones son las mismas. (a) Partícula de prueba fuera de la malla (r(N_r) = 19.95, masa 1 en r = 1): F en r = 20.0, 20.05, 20.1, 25: antes −1.26·10⁻³, −7.5·10⁻⁹, 0, 0; ahora −2.5000·10⁻³, −2.4875·10⁻³, −2.4752·10⁻³, −1.6000·10⁻³ = −m/r². (b) Fuerza sobre partículas con r < 1.5Δr, ρ₀(1−r²)³ cuasi-continuo: error relativo con n = 1 de 31 veces la exacta a 6.3·10⁻³ … 4.1·10⁻⁵ (Δr = 0.04 … 0.005); con n = 2 y 3, de 38 y 47 veces a 0.20 y 0.32 (el sesgo del primer nodo de E8). El error máximo sobre las partículas en r < 0.9 con n = 1 pasa de orden 1 (6.4·10⁻² → 8.1·10⁻³) a orden 2.00 (1.26·10⁻² → 1.98·10⁻⁴). En el interior, igual. (c) `sg__quad` con r_max = 8 (hasta 53 % de la masa fuera de la malla, t = 400): max\|ΔE/E₀\| 6.0·10⁻⁴ → 4.1·10⁻⁴; la masa fuera de la malla sigue sin gravitar entre sí, así que la malla debe contener el sistema. (d) Colapso frío: sin cambio (sus capas en 0.1 < r₀ < 0.9 no bajan de r ≈ 0.064 en t = 0.8). Neutralidad (1 hilo): `sg__quad` completa (101 instantáneas) y `landau__L_a1e-2_n400_e0.1` completa (751): idénticas bit a bit, energías y `.tl` incluidos. Costo, ABBA con `sg__quad` (2000 pasos): 7.89, 7.26, 7.10, 7.90 s, ≈ 9 % menos (el bucle recorre menos nodos) |


| Id | Severidad | Archivo | Función | Problema | Evidencia | Corrección |
|---|---|---|---|---|---|---|
| E1 | **HIGH** (**mitigado**: aborta; la convención 𝔉 queda pendiente) | `initial_data.f90:93,262,327,429,506,552`; `analysish.f90` | `initial_data`, `analysish` | **L₀ = 0 no funciona:** la normalización divide por L₀ (masa NaN) y h_k multiplica por L₀, mientras `density`/`energy` usan `factor = 1`: tres convenciones distintas | corrida: "Initial total mass = NaN", h_k NaN, salida 0 | Trabajar con 𝔉 = 8π²L₀F (masa por dr dp) en todo el código, o abortar con L₀ = 0 |
| E2 | **HIGH** (**corregido**) | `paramfile.f90:370`, `initial_data.f90:56` | `initial_data` | `state = gaussian` (el valor por omisión) no entra en ninguna rama: el código espera `gaussian1`. Partículas en r = 0 con f = 0 | corrida: todo NaN, salida 0 | Unificar el nombre y abortar ante un estado sin rama |
| E3 | **HIGH** (**corregido**) | `main.f90:311` | bucle principal | La reflexión (r,p) → (−r,−p) **no invierte la fuerza**, que se usa en el primer medio kick del paso siguiente | L₀ = 0 frente a RK4: yoshida4 cae a orden 2 (1.6·10⁻³ → 2.4·10⁻⁵), peor que leapfrog; energía ∝ Δt² | `force_part = −force_part` al reflejar (como `_sp`, `a81f4cf`) |
| E4 | **HIGH** (**corregido**) | `grav_force.f90:113-131` | `grav_force` | `iso`, `isotrun`, `nfw`, `burkert` escriben solo la malla: **las partículas no sienten el fondo**; con autogravedad además sobrescriben Φ_gas en la malla | partícula en r = 1.5: F = +0.074 (solo centrífugo) en vez de −1.9 … −2.2 | Evaluar en las partículas, en \|r\|, y sumar (como `_sp`, `8f079fd`, `e2de4f2`) |
| E5 | **HIGH** (**corregido**) | `grav_force.f90:57-62, 64-76` | `grav_force` | `null` sin autogravedad no reinicia la fuerza (el centrífugo se acumula en cada llamada); `sphere` **asigna** en vez de sumar y borra la autogravedad | F = 0.074, 0.148 en llamadas sucesivas; sphere + autograv.: −0.3704 frente a −0.3770 | Reiniciar y sumar |
| E6 | **HIGH** (**corregido**) | `utils.f90:1061-1138`, `main.f90:445` | `reduce_arrays` | Tras reducir en la evolución, `force_part` queda reasignado **sin inicializar** y no se recalcula: el primer medio kick del paso siguiente usa memoria no inicializada | lectura (`utils.f90:1122`, `main.f90:445`); magnitud no medida | Llamar a `grav_force` después (como `_sp`) |
| E7 | MEDIUM (**corregido**) | `initial_data.f90:237-259` (y `aa_halton`) | estado `aa` | Nodos no ligados o NaN van a r = 10⁴ con f = 0, pero `reduce_arrays` solo los quita si queda < 90 %; `analysish` no excluye no ligados | 1 nodo no ligado entre 1800: **todos los h_k NaN** toda la corrida, salida 0 | Quitar siempre esos nodos; excluir E ≥ 0 en `analysish` |
| E8 | MEDIUM (**corregido**) | `poisson_rk.f90` | `poisson_rk` | RK2 sobre (Φ, Φ′): pierde la masa de las partículas cerca del origen y la cuenta mal hasta ~5Δr | masa del campo de una partícula: 0 en r₁, 0.60 en Δr, **1.21** en 1.5Δr, 1.007 en 5Δr, 1.0004 en 20Δr; masa total 1.0013 con Δr = 0.04 | Integrar la masa encerrada (`_sp`, `5c0e07f`) |
| E9 | MEDIUM (**corregido**) | `density.f90` | `deposit` | Sin imágenes en −r_j: la parte del peso que cae en los fantasmas se pierde | hasta 32 % de la masa de una partícula en r = Δr/4 | Sumar W(r_i + r_j) |
| E10 | MEDIUM (**corregido**) | `poisson_rk.f90:170-190` | `poisson_rk` | La interpolación solo usa los nodos 1…N_r: sin espejo en el origen ni solución exterior. g(0⁺) ≠ 0 y más allá de r(N_r) no hay autogravedad | fuerza ~39 veces la exacta en r < 1.5Δr; F = 0 en r = 20.05 (debía ser −2.5·10⁻³) | Nodos espejo y solución exterior (`_sp`, `a81f4cf`) |
| E11 | MEDIUM (**corregido**) | `utils.f90:521-569` | `set_timestep` | El paso no resuelve el pericentro de L₀ pequeño | Δt = 0.01: L₀ = 10⁻³ da \|ΔE/E\| = 3.5 (leapfrog), 4·10⁷ (yoshida4); L₀ ≥ 0.1: ≤ 2.5·10⁻⁵ | Cota Δt ≲ η r_p²/L₀ |
| E12 | MEDIUM (**corregido**) | `utils.f90:318-359` | `invert_QJ_to_rp` | Newton de Kepler desde η₀ = Q, sin salvaguarda: diverge para e ≳ 0.98 (misma rutina que `_sp`, D3) | auditoría de `_sp`, U6 | Arranque robusto + bisección |
| E13 | LOW (**corregido**) | `utils.f90:366-397`, `analysish.f90`, `invert_QJ_to_rp` | mapa AA | Radicandos sin acotar (órbitas circulares) y sin exclusión de no ligados | lectura; igual que `_sp` antes de `e2de4f2` | `max(...,0)`; una sola rutina |
| E14 | LOW (**corregido**) | `density.f90:134-146` | `density` | `vlasov_rhomix` suma f·¼π/(r₂² − r₁²): le faltan 8π²L₀Δr_cΔp_c y el volumen es r₂² − r₁² en vez de (4π/3)(r₂³ − r₁³). El diagnóstico no es una densidad | lectura | Como `_sp` (`90ff2fa`) |
| E15 | LOW (**corregido**) | `initial_data.f90:82-83, 210-211` | `gaussian1`, `aa` | Nodos en p en el borde derecho de la celda (p_j = p_min + jΔp_c); en `gaussian1` también r desplazado una celda completa ((i+½)Δr_c, el último fuera de la caja) | lectura | Puntos medios |
| E16 | LOW (**corregido**) | `utils.f90:56-72` | `set_grid_size` | `int()` trunca N_r; con r_max = 20, Δr = 0.1 el último nodo es 19.95 < r_max | t_borde | Redondeo con tolerancia |
| E17 | INFO (**corregido**) | `analysish.f90:140`, `energy.f90:57` | — | Reducciones OpenMP no deterministas (h_k y E, ~10⁻¹⁵) | §C.3 | Sumas por hilo (`_sp`, `e6d2011`) |
| E18 | INFO | `parameters.f90`, `grav_force.f90` | — | `eps` suaviza el centrífugo si el usuario lo fija ≠ 0; la dinámica deja de ser la del isócrono. Por omisión 0 | lectura | Mantener 0 |
| E19 | INFO (**corregido**) | `initial_data.f90` | — | Estados muertos (`gaussian2`, `Plummer`, `compact`, `compact2`, `other3`) que `paramfile` no admite | lectura | Quitarlos |
| E20 | INFO (**corregido**) | `BUGS_TODO.md` | — | Dice que `eps = Lfix/(10 pmax)` está activo: ya no (eps = 0 por omisión desde `23e03de`) | lectura | Actualizar |
| E21 | MEDIUM (**corregido**) | `density.f90` | `density`, `avg_density` | (Hallado al corregir E9.) El volumen r_i²Δr + Δr³/12 es el de la celda, no el que cubre W_n, cuyo segundo momento es (n+1)Δr²/12: una densidad uniforme sale (r_i² + (n+1)Δr²/12)/(r_i² + Δr²/12) veces la verdadera. En el primer nodo, +25 / +50 / +75 % (n = 1, 2, 3) sin converger con Δr; en el nodo i cae como Δr²/r_i² | t_E9: esfera uniforme, 0.250 / 0.500 / 0.750 exactos | V_i = 4πΔr(r_i² + (n+1)Δr²/12), como `_sp`. Decisión pendiente: el manuscrito usa el volumen geométrico. Desde E8 solo afecta a la densidad que se escribe: Poisson reconstruye la suya a partir de las masas de los nodos |
| E22 | LOW (**corregido**) | `functions.f90:23-29` | `Sn(1)` | (Hallado al corregir E9.) El top-hat vale 1 en \|y\| ≤ 1/2, cerrado en los dos extremos: una partícula justo en una cara de celda, o en r = 0 con su imagen, cuenta doble en `rho` (solo diagnóstico de salida; Poisson usa `avg_rho`) | t_E9: masa 2 en r_j = 2Δr y en r_j = 0 | Peso 1/2 en \|y\| = 1/2 |

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
