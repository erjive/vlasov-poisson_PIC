# Auditoría del proyecto `vlasov-poisson_PIC`

Fecha de la auditoría: 2026-09-25
Alcance: código fuente (`src/`, `Makefile`), documentación interna
(`BUGS_TODO.md`, `AUDITORIA_L0_2026-09-21.md`) y los dos documentos del
experimento eta (`docs/experimento_eta/bateria_eta.tex`,
`docs/demo_eta/demo_eta.tex`).

---

## 1. El código: qué es y cómo funciona

Un PIC (partícula-en-celda) en Fortran 90 con OpenMP para la **ecuación de
Vlasov–Poisson en simetría esférica a momento angular fijo**.

### 1.1 Modelo reducido

La densidad en el espacio fase se escribe

```
f(r, p_r, L, t) = F(r, p_r, t) · δ(L − L₀)
```

La medida de fase es `d³x d³v = 8π² L dr dp dL`, de la que salen todos los
factores del código:

- Masa: `M = 8π² L₀ ∫∫ F dr dp`
- Densidad: `ρ(r) = (2πL₀/r²) ∫ F dp`
- Peso de partícula nodal: `m_j = 8π² L₀ Δr_c Δp_c f_j`

La dinámica se reduce a un problema 2D en `(r, p_r)`:

```
dr/dt = p
dp/dt = L₀²/r³ − ∂Φ/∂r
```

con `H = p²/2 + L₀²/(2r²) + Φ_ext + Φ_gas`.

### 1.2 Estructura de `src/`

```
main.f90          bucle principal e integradores
parameters.f90    parámetros globales
paramfile.f90     lector "nombre = valor" + línea de comandos + validación
arrays.f90        arrays de partículas y de malla
utils.f90         grilla, memoria, RNG, pericentro, kepler_eta,
                  reduce_arrays, guardado
functions.f90     núcleos B-spline W_n y Sn
distribution.f90  F0(Q,J) en variables ángulo-acción
initial_data.f90  estados iniciales: gaussian, aa, aa_halton, aa_quad,
                  aa_random, checkpoint
density.f90       depósito a la malla (rho, avg_rho, curr, rhomix)
grav_force.f90    fondos analíticos + centrífugo + autogravedad
poisson_rk.f90    Poisson por masa encerrada + interpolación
energy.f90        energías
analysish.f90     proyección h_k sobre funciones de prueba
hdf5_io.f90       salida HDF5
raw_io.f90        salida binaria cruda
```

**Detalle del Makefile:** compila TODOS los `.f90` de `src/` automáticamente
(`wildcard`), con los módulos en orden fijo. En consecuencia, los archivos
sin extensión `src/poisson`, `src/poisson_ps` y `src/reduce_arrays`
**no compilan**: son versiones viejas reemplazadas por `poisson_rk.f90` y por la
rutina en `utils.f90`. Candidatas a limpieza del repo.

### 1.3 Ciclo de vida de una corrida

1. Lee parámetros (archivo + overrides de línea de comandos) y **aborta si
   `Lfix = 0`** (la convención `f = F·δ(L−L0)` se rompe en ese límite; el plan
   es migrar a `𝔉 = 8π²L0·F`).
2. Arma la grilla radial (staggered para `rmin = 0`).
3. Genera partículas según `state`.
4. Si `integrator = "analytic"`, guarda `(Q3, J3)` iniciales.
5. Inicialización: `density → grav_force → energy → analysish → set_timestep`.
6. Bucle temporal con `euler`, `leapfrog` (kick-drift-kick, simpléctico),
   `yoshida4/6` (composiciones por tabla) o `analytic`
   (avance exacto `Q3(t) = Q3(0) + ω(J3)t`). `rk4` aborta a propósito
   (no simpléctico, reintroduce deriva secular).
7. Tras cada paso: reflexión en el origen **con cambio de signo de la fuerza**
   (E3), depósito/Poisson según `forcetype`, salidas, `analysish`,
   `reduce_arrays`.

### 1.4 La parte PIC

- **Depósito:** peso B-spline `W_n` con **imagen en (−r_j, −p_j)** para
  conservar masa cerca del origen; lista de celdas (counting sort paralelo)
  para evitar barridos O(N_r·N_part).
- **Poisson:** integra la **masa encerrada** con forma cerrada sobre cada tramo
  (ρ lineal entre nodos); el campo ve exactamente la masa depositada.
- **Interpolación:** nodos espejo para `j ≤ 0` y solución exterior
  `Φ ∝ 1/r` más allá de la malla.

### 1.5 Diagnósticos `h_k`

`analysish.f90` proyecta el ensamble sobre gaussianas de prueba y extrae
modos `k = 0…4`:

```
h_k = 8π² L₀ · Δr_c Δp_c · Σ F(Q,J) Φ̂_k(J) e^{−ikQ}
```

La cuadratura en Q se factoriza exactamente (`A(Q)·B(J)`) y se calcula una
vez por corrida (~10× más rápido). Guarda magnitud y **complejo** (fase).

### 1.6 Estado del proyecto

- Rama auditada: `fix/auditoria-L0`.
- Auditoría exhaustiva en `AUDITORIA_L0_2026-09-21.md` con hallazgos E1–E22,
  casi todos corregidos y verificados bit a bit.
- Pendientes en `PREGUNTAS_ABIERTAS.md`: mecanismo del ruido de discreción
  creciente, extensión de las mediciones de Landau, régimen `a0 ≥ 1e-2`.

---

## 2. Verificación de consistencia de los documentos eta

Los números clave cierran entre sí y con cuentas independientes.

### 2.1 Estructura del isócrono

- `c = (L₀+√(L₀²+4))/2 = 2.4142` con `L₀=2`; `Ω(0) = c⁻³ = 0.07106`,
  que coincide con el `0.0711` citado.
- La identidad epicíclica enmarcada `κ² = 4πρ(r_c) + GM(<r)/r_c³` es
  correcta (Φ″+3Φ′/r con Poisson y Φ′=GM/r²).
- Verificación numérica con el isócrono: `ρ(r_c)=1.09e-4`,
  `M(<r_c)=0.6966`, `κ=0.0711`. Todo cuadra.

### 2.2 Decaimiento libre gaussiano

- De `|h₁|/|h₁(0)| = e^{-σ_Ω²t²/2}`: `σ_Ω` inferido 2.73, 2.81, 2.88, 3.00,
  3.07×10⁻³ para t=400–1400, con deriva del 12%.
- `ΔΩ/σ_Ω ≈ 6.2–7.0` contra `6.1` esperado: convincente.

### 2.3 Parámetro de O'Neil

- Para A4 (η=1): `ν≈11√ε` da `3.01` para ε=0.075 (reportan 3.0) y `10.08`
  para ε=0.84 (reportan 10): consistente.

### 2.4 Resonancia fuera de la banda

- L5: x=−0.01 → distancia 1.75×10⁻⁴, `J_r≈J_t+0.002`.
- L6: x=−0.19 → 4.3×10⁻³, `J_r≈0.164`.
- Ambas cuadran con lo citado. La física "la amplitud finita puebla la
  resonancia de un modo discreto pegado al borde" sale de números
  internamente consistentes.

### 2.5 Escalas temporales

- τ₁ = 2π/ΔΩ da 277–556, compatible con "280–560".
- Periodo de rebote de D5: 5700 ≈ 14τ₁ (τ₁ ≈ 395): cierra.

**Conclusión parcial: no se encontraron inconsistencias aritméticas internas.
Los números se reproducen entre sí.**

---

## 3. Juicio de árbitro sobre el estudio eta

### 3.1 Relevancia

**Sí, es relevante y oportuna.** La pregunta —si la autogravedad de una
componente de masa pequeña puede sostener una oscilación sin amortiguar, y
qué significa η_c— toca tres líneas vivas:

1. Amortiguamiento de Landau en sistemas esféricos (Antonov, Kalnajs).
2. La dicotomía del exponente del borde de
   Hadžić–Rein–Schrecker–Straub (ARMA 2025).
3. El régimen de atrapamiento de O'Neil cerca de una transición con γ→0.

La reducción a L fijo es extrema pero legítima: fondo exacto, variables
ángulo-acción limpias y un solver lineal como referencia precisa. El diseño
"teoría lineal como referencia + PIC solo donde la lineal falla" es la
práctica correcta.

### 3.2 Fortalezas

1. **Separación de escalas bien construida.** η como cociente entre
   corrimiento y ancho de banda, ligado al espectro esencial del operador de
   Antonov, es sensato.
2. **Uso de ν en vez de ε.** D1 (ε=1, ν=0.17, sigue lineal al 0.1%)
   ilustra muy bien que la linealidad la decide ν, no ε.
3. **Colas de mezcla de fases calculadas en serio.** La "regla de las colas"
   con potencias t⁻⁽α⁺¹⁾, los puntos de retorno, la componente regenerada en
   Ω_max y su verificación contra una corrida libre es análisis de alta
   calidad.
4. **Partida silenciosa en (Q,J).** Elimina ruido de muestreo de raíz:
   elección correcta para medir señales a 10⁻³–10⁻⁴.
5. **Validación del régimen lineal.** La PIC da frecuencia y tasa al 1%
   (D1, D2) y frecuencia de modos discretos a 10⁻⁴ (D3, D4): es lo que un
   árbitro quiere ver antes de creer lo no lineal.

### 3.3 Críticas mayores

1. **El resultado central sobre modos discretos no tiene γ lineal acotada
   con firmeza.** Para D3: γ_PIC≈5×10⁻⁵, γ_lineal "<10⁻⁵". La afirmación
   "la teoría lineal no predice la pérdida del 19%" depende de que "<10⁻⁵"
   sea una cota real y no el límite de resolución del matrix pencil. Si la
   γ lineal verdadera es 4–5×10⁻⁵, la pérdida de D3 es lineal y no prueba
   nada nuevo. **Es el punto que más debilita el trabajo.** Hay que demostrar
   la resolución de γ del solver lineal y hacer el barrido en ε
   (D3 a ε→0.003); si la pérdida desaparece cuando ω_b<Ω_min−ω, el
   mecanismo queda probado.

2. **Sin tests de convergencia, los resultados no lineales no son
   publicables aún.** D3 (19%), D5 (saturación a 44τ₁) y D6 (isla)
   descansan en una sola resolución. La "parte fina" de D5 (hasta 2×10⁻²)
   puede ser ruido de discretización amplificado por 1/ε. Hay que mostrar
   que esa parte fina baja como N⁻¹/² y que la parte lisa no cambia al
   partir Δt y cuadruplicar N. Lo que el documento lista como "lo que
   sigue" debe ir en el paper, no después.

3. **La dicotomía del borde es frágil por rango dinámico.** La diferencia
   de D8 (g=1) queda en 2.6×10⁻³, apenas ~13× el piso de 2×10⁻⁴, y se
   pierde a ~15τ₁. Ventana estrecha. Con ε≈0.3–0.5 y sustracción de la
   parte estática podría extenderse; hoy es sugestivo, no concluyente.

4. **η como ordenador de la transición es provisional.** η_c≈1.1–1.2 es
   una medición empírica del solver lineal, no una predicción. La
   robustez de η_c es una conjetura hasta completar el bloque V.

5. **El modelo de L fijo y la comparación con Hadžić et al.** El documento
   es honesto al marcar el modelo como "extremadamente anisótropo", pero
   conviene ser más explícito: solo el borde Ω_min es representativo de un
   sistema isótropo. La dicotomía del exponente tal como la prueban Hadžić
   y colaboradores es para un problema distinto; la analogía es motivación,
   no evidencia de equivalencia.

### 3.4 Críticas menores

1. **Error conceptual sobre el h_k del Fortran.** La batería dice que el
   factor `e^{-sin²(Q/2)/s_Q²}` "solo multiplica cada armónico por una
   constante": falso. Es una función periódica de Q que convoluciona los
   armónicos; solo es constante por armónico en el límite s_Q→∞. No afecta
   al estudio pero debe corregirse.
2. **Separación "lisa/fina" ad hoc.** El suavizado sobre 0.5 en r (5 celdas)
   es una elección no justificada; falta robustez frente a otras escalas.
3. **Ajustes de polo sobre piso.** Extraer ω de D5 en 30–44τ₁
   (Δω=5.7×10⁻⁴ sobre señal 10⁻²) está al límite; faltan incertidumbres
   de ajuste.
4. **Conservación de energía por corrida no reportada.** Con islas y
   atrapamiento a tiempos largos, mostrar la deriva de energía de cada
   corrida como control del integrador.

### 3.5 Veredicto

**Investigación relevante y correctamente diseñada; demo convincente como
prueba de concepto y como validación del régimen lineal; los resultados no
lineales nuevos son todavía insuficientes para publicación.**

Aprobaría como "revisiones mayores" con estas condiciones:

1. Acotar firmemente la γ del solver lineal y hacer el barrido en ε de D3.
2. Tests de convergencia en Δt y N para D3/D5/D6.
3. Extender D8 con ε óptima para salir del piso.
4. Corregir el error del factor en Q del h_k del Fortran.

Con eso, el trabajo sería un aporte sólido a la cuestión
autogravedad-vs-mezcla de fases cerca de la transición.

---

## 4. Observaciones adicionales sobre el código

1. **Tres archivos muertos en `src/`**: `poisson`, `poisson_ps` y
   `reduce_arrays` (sin extensión) no compilan. Son residuos del desarrollo.
2. **`output_format="raw"`** existe como alternativa a HDF5: ~20× más
   rápido en aislamiento, pero solo ~6–13% en corridas realistas; no es
   autodescriptivo y necesita un lector propio.
3. **Optimización mayor ya aplicada**: `analysish` factorizado
   (cuadratura en Q fuera del bucle de partículas): 10.4×; barridos
   seriales fusionados: la corrida completa pasó de 333 s a 14 s (23×).
4. **El piso de h_k estuvo dominado históricamente por `eps`** (suavizado
   del centrífugo); corregido (eps=0 por omisión). Con eps=0 el límite
   pasó a ser la fase del leapfrog; `yoshida4` resuelve los modos altos.
5. **`integrator="analytic"`** es una herramienta de validación excelente
   (avance exacto, sin error de integración); requiere setup integrable y
   aborta con mensaje claro si no se cumple.
6. **`aa_random` sin tope de intentos** si F≈0 en la región muestreada:
   pendiente.
7. **La densidad de salida usa el volumen que cubre W_n** (E21), mientras
   que Poisson reconstruye su propia `avg_rho`; la dinámica es consistente
   pero conviene documentarlo bien para evitar confusiones de usuario.