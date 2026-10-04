# Oscillation threshold of isotropic polytropes and King models from the Mathur operator

4 October 2026. Feasibility test of direction P1 of `docs/literature/spherical_vp_review.md`.

## 1. Result

For the isotropic polytropes f₀ = (E₀ − E)ᵏ of the self-gravitating Vlasov–Poisson system,
the gain of the self-consistent loop at the lower edge of the band of radial frequencies,
λ_edge, which is the norm of the Mathur operator at the edge of the principal gap, crosses
one at

    k* = 1.24253 ± 0.00002.

For k < k* there is a discrete radial mode below the band, and for k > k* there is none.
For the King models f₀ = e^{E₀−E} − 1 with y(0) = κ, the threshold is κ* = 2.0495: there is
a mode for κ < κ*.

Straub (2024) placed the threshold of the polytropes between 1.2 and 1.3 by inspection of
time series, could not classify k = 1.25, and fitted by hand a boundary that gives
12/π² = 1.2159 for isotropic models. For the King models he could state only that
κ ≤ 3/2 oscillates and κ ≥ 5/2 damps. Both thresholds are now fixed by an eigenvalue
computation that takes seconds.

## 2. Method

**Operator.** For a radial perturbation of f₀(E, L), the angular momentum of each particle
is conserved and the linearized equation is the same for each L,

    (kΩ − ω) a_k = k (∂f₀/∂J) φ_k,   ∂f₀/∂J = Ω ∂f₀/∂E,   Ω = Ω(E, L),

with J the radial action and φ_k the harmonic of δΦ along the orbit. The even part
b_k = a_k + a_{−k} obeys b_k = R_k φ_k, R_k = 2k²Ω ∂_J f₀/(k²Ω² − ω²). The mass element is
8π² L dL dJ dQ, so the Poisson equation gives

    φ_k(J, L) = Σ_k′ ∫ L′dL′ dJ′ M_kk′ b_k′,
    M_kk′ = 4π ∬ cos kQ cos k′Q′ G(r, r′) dQ dQ′,   G = −1/max(r, r′).

Compared with the fixed-|L| case only the measure changes, L₀ dJ → L dL dJ. With
∂f₀/∂E < 0 and ω below the band the operator is similar to |D|^{1/2}(−M)|D|^{1/2}, where
D = L dL dE 2k²Ω ∂_E f₀/(k²Ω² − ω²) is the weight of each orbit node. Writing
−M = 4π P(−G)Pᵀ, with P the deposit of cos kQ along each orbit on a radial grid, the
non-zero eigenvalues are those of

    4π Cᵀ (Pᵀ|D|P) C,    −G = CCᵀ,

a matrix of the size of the radial grid whatever the number of orbits. This is the Mathur
operator in configuration space (Hadžić, Rein and Straub 2022). λ(ω) is its largest
eigenvalue; λ_edge = λ(Ω_b), with Ω_b the smallest radial frequency of the support. A mode
exists if and only if λ_edge > 1, and its frequency solves λ(ω_d) = 1.

**Orbits.** With s = r², dt = ds/(2√g), g(s) = 2(E − Φ)s − L², which has simple zeros at
both turning points also when L → 0. The substitution s = s_m + s_a sin θ gives
dt = dθ/(2√h) with h = g/((s − s₋)(s₊ − s)) smooth, and the midpoint rule in θ converges
fast. The turning points come from bisection on each side of the circular orbit.

**Steady states.** y = E₀ − U₀ solves y″ + 2y′/r = −4πρ(y), y(0) = κ, with
ρ(y) = 4π√2 ∫₀^y Φ(η)(y − η)^{1/2} dη; for the polytropes this is the Lane–Emden equation of
index k + 3/2, and κ = 1 by scale invariance. The lowest frequency Ω_b is that of the radial
orbit of energy E₀, which goes from the centre to the surface.

**Nodes.** E₀ − E = κ s^q with s on Gauss–Legendre nodes in (0, 1), which accumulate at the
edge; L = L_c(E) x with x on Gauss–Legendre nodes and L_c the angular momentum of the
circular orbit. Default: 48 × 24 orbits, 128 points in θ, 400 radial nodes, 16 harmonics.

The code is `reproducir/scripts/lambda_L.py`.

## 3. Validation

**Fixed-|L| limit.** With f₀ = F(E) g(L), g narrow around L₀ and ∫ g L dL = L₀, in the total
potential of three fixed-|L| equilibria of the point-mass setting, the new operator must
give the λ of `lambda_borde.py`. The orbit integration, the nodes and the radial grid are
different in the two codes.

| equilibrium | ω/Ω_min | fixed L | σ_L/L₀ = 10⁻³ | σ_L/L₀ = 10⁻² |
|---|---|---|---|---|
| k = 2, a₀ = 1 | 0.5 | 0.497970 | 0.497963 | 0.498431 |
| | 0.9 | 0.674400 | 0.674389 | 0.674874 |
| | 1 | 0.819933 | 0.820117 | 0.822143 |
| k = 3, a₀ = 1 | 0.5 | 0.532715 | 0.532719 | 0.534172 |
| | 1 | 0.738042 | 0.738046 | 0.739771 |
| k = 2, a₀ = 0.1 | 0.5 | 0.097866 | 0.097864 | 0.097960 |
| | 1 | 0.144571 | 0.144568 | 0.144632 |

The agreement is 10⁻⁵ away from the edge. At the edge of the k = 2 case it is 2 × 10⁻⁴,
because with a spread in L the edge of the band moves slightly.

**Mass.** 8π² Σ L dL dE T f₀ over the nodes equals the mass of the Lane–Emden solution to
10⁻⁸ or better for the polytropes and the King models, except k = 3 (8 × 10⁻⁸).

**Antonov stability.** λ(0) is between 0.57 and 0.61 in all models, below one.

**An independent measurement.** Ramming and Rein (2018, Table 3) measured the period of the
oscillation of the isotropic polytrope with k = 1 in nonlinear PIC runs: T = 2.761 for
y(0) = 1. The mode found here has ω_d = 2.27769, that is T = 2π/ω_d = 2.7586. The
difference is 0.09 %.

**Convergence** of k*:

| orbits | θ | radial nodes | harmonics | q | λ(k = 1.2) | λ(k = 1.25) | k* |
|---|---|---|---|---|---|---|---|
| 48 × 24 | 128 | 400 | 16 | 2 | 1.048476 | 0.992052 | 1.242530 |
| 96 × 48 | 128 | 400 | 16 | 2 | 1.048476 | 0.992052 | 1.242530 |
| 48 × 24 | 256 | 400 | 16 | 2 | 1.048495 | 0.992070 | 1.242547 |
| 48 × 24 | 128 | 800 | 16 | 2 | 1.048463 | 0.992039 | 1.242517 |
| 48 × 24 | 128 | 1600 | 16 | 2 | 1.048460 | 0.992036 | 1.242514 |
| 48 × 24 | 128 | 400 | 6 | 2 | 1.048069 | 0.991555 | 1.242082 |
| 48 × 24 | 128 | 400 | 10 | 2 | 1.048467 | 0.992040 | 1.242519 |
| 48 × 24 | 128 | 400 | 24 | 2 | 1.048476 | 0.992052 | 1.242530 |
| 48 × 24 | 128 | 400 | 16 | 3 | 1.048477 | 0.992053 | 1.242530 |
| 96 × 48 | 256 | 1600 | 24 | 2 | 1.048478 | 0.992053 | 1.242531 |

Only the number of harmonics matters below 10; the other changes are at most 2 × 10⁻⁵.

## 4. Results

**Polytropes.** x_d = 1 − ω_d/Ω_b is the distance of the mode below the edge.

| k | R | Ω_b | Ω_max/Ω_b | λ(0) | λ_edge | ω_d | x_d |
|---|---|---|---|---|---|---|---|
| 0.25 | 0.3763 | 5.39886 | 2.26 | 0.613 | 13.12 | 3.92663 | 0.273 |
| 0.5 | 0.4648 | 4.13085 | 2.62 | 0.614 | 4.358 | 3.23026 | 0.218 |
| 0.75 | 0.5689 | 3.17843 | 3.07 | 0.608 | 2.199 | 2.71331 | 0.146 |
| 1 | 0.6940 | 2.44355 | 3.64 | 0.603 | 1.3736 | 2.27769 | 0.0679 |
| 1.1 | 0.7515 | 2.19629 | 3.92 | 0.601 | 1.1878 | 2.11550 | 0.0368 |
| 1.15 | 0.7821 | 2.08118 | 4.07 | 0.600 | 1.1133 | 2.03503 | 0.0222 |
| 1.2 | 0.8141 | 1.97132 | 4.23 | 0.599 | 1.0485 | 1.95372 | 0.0089 |
| 1.22 | 0.8274 | 1.92878 | 4.30 | 0.599 | 1.0250 | 1.92052 | 0.0043 |
| 1.24 | 0.8408 | 1.88703 | 4.37 | 0.598 | 1.0027 | 1.88631 | 0.00038 |
| 1.25 | 0.8477 | 1.86643 | 4.40 | 0.598 | 0.9921 | — | — |
| 1.3 | 0.8828 | 1.76628 | 4.58 | 0.597 | 0.9428 | — | — |
| 1.5 | 1.0415 | 1.40829 | 5.43 | 0.594 | 0.7999 | — | — |
| 2 | 1.6346 | 0.75093 | 8.96 | 0.588 | 0.6442 | — | — |
| 3 | 6.6840 | 0.09841 | 55.8 | 0.578 | 0.5793 | — | — |

Units: G = 1, E₀ − U₀(0) = 1 and f₀ = (E₀ − E)ᵏ without a prefactor. For k < 1/2 the
integrand is unbounded at the corner of the support and the default nodes are not enough:
λ_edge at k = 0.25 is 11.76 with them and 13.12 with 192 × 96 orbits; the value in the
table is the converged one. From k = 0.75 on, the default nodes are converged to 2 × 10⁻⁴.

**King models.**

| κ | R | Ω_b | λ_edge | ω_d | x_d |
|---|---|---|---|---|---|
| 0.5 | 1.1353 | 1.03202 | 1.2774 | 0.98210 | 0.0484 |
| 1 | 0.6576 | 2.45449 | 1.1839 | 2.38270 | 0.0293 |
| 1.5 | 0.4736 | 4.05105 | 1.0938 | 4.00124 | 0.0123 |
| 1.75 | 0.4173 | 4.88501 | 1.0504 | 4.85828 | 0.0055 |
| 2 | 0.3738 | 5.72923 | 1.0082 | 5.72602 | 0.00056 |
| 2.25 | 0.3392 | 6.57358 | 0.9673 | — | — |
| 2.5 | 0.3110 | 7.40764 | 0.9279 | — | — |
| 3 | 0.2685 | 9.00055 | 0.8537 | — | — |
| 4 | 0.2174 | 11.55025 | 0.7274 | — | — |

## 5. Comparison with the literature

| statement | source | here |
|---|---|---|
| polytropes oscillate for 0 < k ≤ 1.2 and damp for k ≥ 1.3 | Straub 2024, Observation 4.1 | k* = 1.2425 |
| k = 1.25: "probably fully damped", "rather slow" damping, "there might also be an undamped part" | Straub 2024 | λ_edge = 0.992: no mode, 0.8 % below the threshold |
| the boundary fitted by hand gives k = 12/π² = 1.2159 | Straub 2024, Observation 4.4 | λ_edge(1.2159) = 1.0297; the threshold is 0.027 higher |
| k* ≈ 1.25 | Straub, thesis, outlook | 1.2425 |
| near k = 1.2 in the limit of small redshift of Einstein–Vlasov; "not feasible to numerically determine the threshold" | Wolfschmidt 2023 | 1.2425 is the Newtonian value |
| King models oscillate for 1/2 ≤ κ ≤ 3/2 and damp for 5/2 ≤ κ ≤ 4 | Straub 2024, Observation 4.3 | κ* = 2.0495 |
| for κ > 1 the frequency is too close to the edge to decide | Straub 2024 | x_d = 1.2 % at κ = 1.5, 0.06 % at κ = 2 |
| the sufficient criterion of Kunze holds only for k below about 0.03 | Straub 2024, Section 4.3 | λ_edge > 1 up to k = 1.2425 |
| period of the k = 1 oscillation, 2.761 | Ramming and Rein 2018 | 2.7586 |

## 6. Why time integration could not decide

The mode separates from the response of the edge of the band after a time of order
1/(Ω_b − ω_d), that is 1/(2π x_d) periods of the slowest orbit: 2 periods for k = 1, 18
for k = 1.2, 37 for k = 1.22 and 420 for k = 1.24. Above the threshold the damping is
correspondingly slow. Within 0.05 of k* a time series of some tens of periods cannot tell
the two cases apart.

Near the threshold x_d is close to linear in k* − k, with a coefficient that decreases
slowly, from 0.26 at k = 1.1 to 0.15 at k = 1.24.

**The edge law with a distribution in L.** λ_edge − λ(Ω_b(1 − δ)) ∝ δ^p, with the measured p:

| k | 0.25 | 0.5 | 0.75 | 1 | 1.25 | 1.5 | 2 |
|---|---|---|---|---|---|---|---|
| p at δ = 10⁻⁴ | 0.25 | 0.47 | 0.68 | 0.84 | 0.94 | 0.98 | 1.00 |

The exponent tends to min(k, 1), against min(k − 1, 1) with fixed |L|. The lowest frequency
is at the corner E = E₀, L = 0 of the support, where Ω − Ω_b ≈ a(E₀ − E) + bL², and the
integral over L removes one power of the singularity. Two consequences: λ_edge is finite
for every k > 0 (it is 13.1 at k = 0.25), so the dichotomy at k = 1 of the fixed-|L|
theorem of Hadžić et al. (2023) becomes a threshold at finite k; and the binding near the
threshold is linear, not a high power.

## 7. What remains before a paper

1. **Literature.** Read Mathur (1990) and the book of Kunze (2021) in full, to confirm that
   the operator was not evaluated numerically for these models.
2. **A check in time.** A linearized time-domain solver with L, or runs of the code with a
   distribution in L, for k = 1.1, 1.2 and 1.3, long enough to see the mode at k = 1.2
   (about 200 periods of the slowest orbit).
3. **Other families.** Shells with L₀ > 0, where Straub finds oscillations always; the
   anisotropic polytropes with ℓ > 0, where the band reaches zero frequency and the modes
   are embedded, so that this criterion does not apply directly.
4. **The relativistic value.** The Newtonian limit of the Einstein–Vlasov threshold found
   near 1.2 should be 1.2425.
5. **The value itself.** Whether 1.24253 has a closed form is not known.

## 8. Reproduction

From `reproducir/scripts/`, with four threads:

| command | output | time |
|---|---|---|
| `python3 lambda_L.py validar` | the fixed-\|L\| limit | 2 min |
| `python3 lambda_L.py politropos` | `exe/lambda_L/politropos.txt`, k* | 2 min |
| `python3 lambda_L.py king` | `exe/lambda_L/king.txt`, κ* | 1 min |
| `python3 lambda_L.py convergencia` | `exe/lambda_L/convergencia.txt` | 7 min |

## References

- Hadžić, Rein, Straub, *On the existence of linearly oscillating galaxies*, ARMA 243 (2022), [arXiv:2102.11672](https://arxiv.org/abs/2102.11672)
- Hadžić, Rein, Schrecker, Straub, *Damping versus oscillations for a gravitational Vlasov–Poisson system*, ARMA 249:45 (2025), [arXiv:2301.07662](https://arxiv.org/abs/2301.07662)
- Kunze, *A Birman–Schwinger principle in galactic dynamics*, [arXiv:2202.04555](https://arxiv.org/abs/2202.04555)
- Mathur, *Existence of oscillation modes in collisionless gravitating systems*, MNRAS 243, 529 (1990)
- Ramming, Rein, *Oscillating solutions of the Vlasov–Poisson system — a numerical investigation*, Physica D 365 (2018), [arXiv:1602.07989](https://arxiv.org/abs/1602.07989)
- Straub, *Numerical experiments on stationary, oscillating, and damped spherical galaxy models*, Physica D (2024), [arXiv:2405.01235](https://arxiv.org/abs/2405.01235)
- Straub, *Pulsating Galaxies*, thesis, Bayreuth 2024, [epub 7639](https://epub.uni-bayreuth.de/id/eprint/7639)
- Wolfschmidt, *Stability and Oscillations of Star Clusters in General Relativity*, thesis, Bayreuth 2023, [epub 7337](https://epub.uni-bayreuth.de/id/eprint/7337)
