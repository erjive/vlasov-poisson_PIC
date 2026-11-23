# The passive μ: definition, sensitivity to the cut and convergence

22 November 2026. Follows mu_conventions.md. No file of the manuscript or of the code was changed. The tests were run
with a script in the session scratchpad that imports `reproducir/scripts/mu_convenciones.py`, `reducido.py` and
`bulto.py` without modifying them; the script and its output are not versioned.

## Summary

- The passive amplitude is half the change of the potential **inside** the cloud, G M_gas ⟨1/r⟩, between the initial
  pericenter and the first apocenter. It does not isolate a coherent oscillation.
- The cut does not change A. It changes ΔJ ΔΩ, and with it μ: +18 % at 5 percent and +32 % at 10 percent for the
  clouds of Section V. The ratios between these clouds change by less than 2 percent, so the comparisons between
  families of Section V do not depend on the cut. Comparisons between shapes do.
- A is converged to about 10⁻⁵ in the five families tested. The radial window [3, 15] cuts the plateau in the Kepler
  potential and lowers A by 1.8 percent there.
- Precision: μ is numerically exact to three digits within a convention. Its absolute value carries a convention
  dependence of about 30 percent. Between families with the same shape, the change agrees to about 0.05 in μ.

## 1. How A is computed for the passive cloud

**Datum and particles.** The initial DF is Eq. (42), F₀(Q, J) = α exp[−sin²(Q/2)/σ_Q²] f(J), with f(J) = J² exp(−J²/σ_J²)
or one of the shapes of VII.E. It is sampled on the lattice of Eq. (41), J_a = (a − 1/2) 6σ_J/N_J and
Q_b = (b − 1/2) 2π/N_Q, with N_J × N_Q = 400 × 100. Points where either factor is below 10⁻⁷ of its maximum are
dropped. Since the cells have equal volume, the mass fractions are

$$ w_{ab} = \frac{F_0(Q_b, J_a)}{\sum_{a'b'} F_0(Q_{b'}, J_{a'})}, \qquad \sum_{ab} w_{ab} = 1 . $$

**Motion.** The gas moves in the external potential alone, with L = L₀ = 2. The exact solution is

$$ J_{ab}(t) = J_a, \qquad Q_{ab}(t) = Q_b + \Omega(J_a)\,t, \qquad \Omega(J) = \left[J + c(L_0)\right]^{-3}, $$

Eq. (23), with c = (L₀ + (L₀² + 4)^(1/2))/2 in the isochrone potential and c = L₀ in the Kepler potential.

**Radii.** Each lattice point is a thin shell of radius r_ab(t) = r(Q_ab(t), J_a). In the isochrone potential
(Binney and Tremaine, section 3.5.2), with E = −1/[2(J + c)²], C = −1/(2E) − 1 and e = [1 − (L₀²/C)(1 + 1/C)]^(1/2),

$$ \eta - \frac{eC}{C+1}\sin\eta = Q, \qquad s = 2 + C\,(1 - e\cos\eta), \qquad r = \sqrt{s(s-2)} . $$

In the Kepler potential, Eq. (43): η − e sin η = Q and r = a(1 − e cos η), with a = (J + L₀)² and
e = (1 − L₀²/a)^(1/2). Kepler's equation is solved with 60 Newton steps.

**Density.** No grid and no smoothing are used. The density is the sum of the shells,

$$ \rho(t, r) = M_{gas}\sum_{ab} w_{ab}\,\frac{\delta\!\left(r - r_{ab}(t)\right)}{4\pi r^2} . $$

**Potential.** Eq. (21) for this density gives, exactly,

$$ \Phi_{gas}(t, r) = -G M_{gas}\left[\frac{1}{r}\sum_{r_{ab}(t) < r} w_{ab} + \sum_{r_{ab}(t) \ge r} \frac{w_{ab}}{r_{ab}(t)}\right] . $$

It is evaluated on 241 radii r_i equally spaced in [3, 15] (Δr = 0.05), at the times t_k = k Δt with Δt = 0.25.

**Observable and interval.** With T = 2π/Ω_min the period of the outermost orbit of the band,

$$ A = \max_{r_i \in [3,15]}\ \frac{1}{2}\left[\max_{0 \le t_k \le T}\Phi_{gas}(t_k, r_i) - \min_{0 \le t_k \le T}\Phi_{gas}(t_k, r_i)\right], \qquad \mu^2 = \frac{A}{\Delta\Omega\,\Delta J} . $$

Since A ∝ M_gas, μ = c (M_gas/M)^(1/2), with c fixed by the shape, σ_Q, σ_J and the potential.

**What the observable is, in fact.** In every family the maximum over r falls on the edge r = 3 of the window. In the
isochrone potential no shell comes below r = 3: the smallest pericenter is 3.4 (σ_J = 0.2) to 4.1 (σ_J = 0.05). The
potential there is the constant value inside all the shells, and A reduces to

$$ A = \frac{G M_{gas}}{2}\left[\max_{0\le t\le T}\left\langle r^{-1}\right\rangle(t) - \min_{0\le t\le T}\left\langle r^{-1}\right\rangle(t)\right], \qquad \left\langle r^{-1}\right\rangle(t) = \sum_{ab}\frac{w_{ab}}{r_{ab}(t)} . $$

With the window [1, 15] the result is the same to all digits in these families. In the Kepler potential, 40 percent
of the mass has its pericenter below r = 3, the window misses the plateau, and A is 1.8 percent lower than with the
window [1, 15] (0.07498 against 0.07637 per unit mass).

**Time interval.** The minimum of Φ_gas occurs at t = 0, when the cloud is at its pericenter. The maximum occurs near
the first passage through the apocenter. A reaches its final value within [0, T/2] and does not change for windows up
to 3T. With [0, T/4] it is 5 to 24 percent smaller, depending on the family.

**Coherent oscillation and other components.** The procedure does not separate them. It takes the whole
peak-to-peak change at fixed r. Compared with the potential of the mixed state (the same actions with the angles
uniform), the excursions are not symmetric:

| Family | Φ_gas(0, 3)/M_gas (pericenter) | Mixed state | Largest value in [0, T] | Down / up excursion |
| --- | --- | --- | --- | --- |
| Reference | −0.2233 | −0.1590 | −0.1239 | −0.064 / +0.035 |
| σ_J = 0.05 | −0.2100 | −0.1663 | −0.1379 | −0.044 / +0.028 |
| σ_J = 0.2 | −0.2393 | −0.1460 | −0.1061 | −0.093 / +0.040 |
| σ_Q = 0.4 | −0.2037 | −0.1590 | −0.1275 | −0.045 / +0.031 |
| Kepler (r = 3) | −0.3223 | −0.2242 | −0.1723 | −0.098 / +0.052 |

The downward excursion is 1.4 to 2.3 times the upward one: the concentration of the cloud at its pericenter raises
⟨1/r⟩ more than its first apocenter lowers it. During the first half period the two ends of the band also drift apart
by ΔΩ T/2 ≈ 1.2 radians in the reference cloud, but these edges carry little mass: the mass-weighted dipole loses only 1.5 percent by then (A_interpretation.md, section 3), so the shear is not the cause of the asymmetry. A
is therefore the largest change of the potential inside the cloud during its first half period, not the amplitude of
a sinusoidal oscillation. For a single compact shell this is the quantity of Eq. (58), which gives 0.0509 per unit mass
for the reference cloud, against 0.04967 here.

## 2. Sensitivity to the cut

The cut p is the fraction of the maximum of f(J) that defines [J_lo, J_hi]. The time window T = 2π/Ω_min changes with
p, but A is reached within T/2 for every p tested, so **A does not depend on the cut** (identical to all digits). The
change of μ comes from ΔΩ ΔJ alone.

| Family | Cut | J_lo – J_hi | Mass inside | ΔJ | ΔΩ | τ₁ | c | μ / μ(1 %) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Reference, σ_Q = 0.2, σ_Q = 0.4 | 1 % | 0.0061 – 0.2764 | 99.82 % | 0.2703 | 0.01919 | 327 | 3.094 / 3.007 / 2.710 | 1 |
| | 5 % | 0.0137 – 0.2397 | 98.87 % | 0.2260 | 0.01637 | 384 | 3.664 / 3.561 / 3.209 | 1.184 |
| | 10 % | 0.0196 – 0.2211 | 97.40 % | 0.2016 | 0.01473 | 427 | 4.090 / 3.975 / 3.582 | 1.322 |
| σ_J = 0.05 | 1 % | 0.0030 – 0.1382 | 99.82 % | 0.1351 | 0.01066 | 589 | 5.001 | 1 |
| | 5 % | 0.0068 – 0.1198 | 98.87 % | 0.1130 | 0.00901 | 697 | 5.949 | 1.190 |
| | 10 % | 0.0098 – 0.1106 | 97.40 % | 0.1008 | 0.00808 | 778 | 6.653 | 1.330 |
| σ_J = 0.2 | 1 % | 0.0122 – 0.5527 | 99.82 % | 0.5406 | 0.03172 | 198 | 1.972 | 1 |
| | 5 % | 0.0274 – 0.4793 | 98.87 % | 0.4519 | 0.02743 | 229 | 2.319 | 1.176 |
| | 10 % | 0.0391 – 0.4423 | 97.40 % | 0.4031 | 0.02482 | 253 | 2.581 | 1.309 |
| Kepler | 1 % | 0.0061 – 0.2764 | 99.82 % | 0.2703 | 0.03909 | 161 | 2.664 | 1 |
| | 5 % | 0.0137 – 0.2397 | 98.87 % | 0.2260 | 0.03345 | 188 | 3.149 | 1.182 |
| | 10 % | 0.0196 – 0.2211 | 97.40 % | 0.2016 | 0.03014 | 208 | 3.513 | 1.319 |
| J⁴ | 10 % | 0.0515 – 0.2594 | 97.15 % | 0.2079 | 0.01438 | 437 | 4.373 | 1.359 |
| (J−0.05)₊² | 10 % | 0.0633 – 0.2390 | 97.32 % | 0.1757 | 0.01222 | 514 | 5.111 | 1.332 |
| Ring | 10 % | 0.1045 – 0.1955 | 96.81 % | 0.0910 | 0.00632 | 994 | 10.03 | 1.415 |
| Water bag | any | 0.05 – 0.25 | 100 % | 0.2000 | 0.01395 | 451 | 4.502 | 1 |
| J⁰ | 10 % | 0 – 0.1517 | 96.81 % | 0.1517 | 0.01188 | 529 | 4.328 | 1.381 |

(At 5 percent the factors of the other shapes are 1.206, 1.191, 1.240, 1 and 1.221.)

**Effect on the change.** The runs that bracket the change (Table III) take these values of μ:

| Family | 1 % | 5 % | 10 % |
| --- | --- | --- | --- |
| Reference | 0.438 – 0.489 | 0.518 – 0.579 | 0.578 – 0.647 |
| σ_J = 0.05 | 0.354 – 0.500 | 0.421 – 0.595 | 0.470 – 0.665 |
| σ_J = 0.2 | 0.394 – 0.558 | 0.464 – 0.656 | 0.516 – 0.730 |
| σ_Q = 0.2 | 0.425 – 0.601 | 0.504 – 0.712 | 0.562 – 0.795 |
| σ_Q = 0.4 | 0.542 – 0.766 | 0.642 – 0.908 | 0.716 – 1.013 |
| Kepler | 0.377 – 0.533 | 0.445 – 0.630 | 0.497 – 0.703 |
| Common part, σ_Q ≤ 0.2 | 0.438 – 0.489 | 0.518 – 0.579 | 0.578 – 0.647 |

In every case the common part is the interval of the reference family. The deviation of m_c from the reference
family at equal μ (σ_J = 0.05, 0.2 and Kepler, five runs with a core) has a mean of 8.6, 9.5 and 9.9 percent for the
three cuts.

**Is the 1 percent cut physical?** No. It is an operational convention. The band at 1 percent holds 99.8 percent of the
mass, at 10 percent 97.4 percent: removing 2.4 percent of the mass from the tails changes ΔJ ΔΩ by a factor 0.57 and μ
and τ₁ by about 30 percent. Nothing in the dynamics singles out 1 percent; the band is set by the tails of the DF. The
consequences are of two kinds.

- For clouds of one shape (the J² clouds of Section V, in either potential and for every σ_J and σ_Q tested) the
  factor μ(p)/μ(1 %) is the same within 1.6 percent. The relative statements of Section V (same change for all
  families, collapse of m_c) do not depend on the cut.
- Absolute values do: "μ ≳ 0.4", "the change lies between 0.44 and 0.49" and τ₁ hold only for 1 percent. Between
  shapes the factors differ by up to 1.42 (ring) against 1 (water bag), so the comparison of shapes at equal μ in
  VII.E depends on the cut. I have not recomputed the crossings of VII.E with other cuts.

A definition based on the mass distribution (a quantile or a moment of Ω) would not depend on the tails. I have not
computed it.

## 3. Convergence

Relative change of a = A/M_gas with respect to the base computation (400 × 100 points, Δt = 0.25, 241 radii in [3, 15]):

| Change | Reference | σ_J = 0.05 | σ_J = 0.2 | σ_Q = 0.4 | Kepler |
| --- | --- | --- | --- | --- | --- |
| Base a | 0.049673 | 0.036036 | 0.066640 | 0.038095 | 0.074982 |
| N_J × 2 | 6 × 10⁻⁹ | 6 × 10⁻⁹ | 6 × 10⁻⁹ | 6 × 10⁻⁹ | 1.3 × 10⁻⁷ |
| N_Q × 2 | 1 × 10⁻⁹ | 1 × 10⁻⁹ | 2 × 10⁻⁹ | < 10⁻¹⁴ | 5.8 × 10⁻⁶ |
| Δt = 0.1 | 9 × 10⁻⁷ | 6 × 10⁻⁶ | 7 × 10⁻⁷ | 3 × 10⁻⁷ | 1.2 × 10⁻⁵ |
| Radii × 2 | 0 | 0 | 0 | 0 | 0 |
| Snapshots every 10 (as the PIC output) | −6 × 10⁻⁵ | −3.3 × 10⁻³ | −2.5 × 10⁻³ | −1.4 × 10⁻⁴ | −3.3 × 10⁻⁴ |
| 400 × 25 and snapshots every 10 | −1.3 × 10⁻³ | −4.2 × 10⁻³ | −5.2 × 10⁻³ | −1.4 × 10⁻⁴ | −2.3 × 10⁻⁴ |
| Radial window [1, 15] | 0 | 0 | 0 | 0 | +1.8 × 10⁻² |

All the tests requested are available, and convergence holds in the five families. The only effect above 10⁻³ that is
not numerical is the radial window in the Kepler potential: a convention of the definition, not a discretization error.
The computation for the σ_Q = 0.2 family and for the shapes of VII.E was not repeated with finer parameters; their
values come from the base computation.

## Conclusion: what precision μ has

| Source | Size in μ | Kind |
| --- | --- | --- |
| Lattice, time step and radii of the passive computation | about 10⁻⁵ or less | Numerical |
| PIC-like sampling (snapshots every 10, N_Q = 25) | up to 0.3 % | Numerical, only if A is measured that way |
| Radial window [3, 15] in the Kepler potential | 0.9 % | Convention |
| Run instead of passive amplitude | 1 to 3 % for M_gas/M ≤ 0.03; up to 47 % at 0.3 | Convention (mixes in the outcome) |
| Leading order instead of passive | −15 to +4 % | Convention (approximation) |
| Cut at 5 or 10 % instead of 1 % | +18 and +32 % (J² clouds); 0 to +42 % between shapes | Convention |
| Location of the change between J² families at equal μ | brackets 0.14 – 0.16 wide, all containing 0.44 – 0.49; collapse of m_c within ±19 %, or about 0.06 in μ | Physical (finite width of the band) and limited resolution in mass |
| μ²/Γ between J² families | 0.77 to 0.96 (σ_J), 0.64 to 0.85 (σ_Q, Kepler) | Physical |

In consequence:

- Within one convention, μ is a precise label: three digits are exact, and steps of 0.01 between masses of one family
  are meaningful.
- Between families of the same shape, a common location of the change can be claimed to about ±0.05 in μ, not better.
- The absolute value of μ is a convention. Statements such as "μ ≳ 0.4" or "between 0.44 and 0.49" must state the
  definition, including the cut, and do not have a meaning independent of it.
- Comparisons of μ between shapes depend on the cut by up to 40 percent.

Not checked: a cut-free (mass-based) definition of the band; the crossings of VII.E with other cuts; whether A
measured with the potential inside the cloud is the quantity that controls the frequency shift in the heuristic of
V.A.
