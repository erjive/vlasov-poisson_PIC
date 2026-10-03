# Hadžić and collaborators on the gravitational Vlasov–Poisson system: what to study numerically

This note collects what Hadžić, Rein, Schrecker, Straub and Moreno have proved about the
linearized gravitational Vlasov–Poisson system, and the numerical studies of the same group.
It compares their setting with that of our PIC code and lists the statements that can be
tested numerically. The PDFs are stored locally in `exe/referencias/hadzic/` (not versioned);
the arXiv links are below.

## 1. The papers

### Theory

| ref. | paper | setting | main result |
|---|---|---|---|
| [HRS21] | Hadžić, Rein, Straub, *On the existence of linearly oscillating galaxies*, ARMA 243 (2022), [arXiv:2102.11672](https://arxiv.org/abs/2102.11672) | spherical Antonov-stable steady states (polytropes, King) and plane-symmetric analogues | essential spectrum of the Antonov operator; principal gap; Birman–Schwinger (Mathur) criterion for an eigenvalue in the gap; existence of oscillating modes (rigorous in the planar case, under a monotonicity assumption on the period function in the spherical case) |
| [HRSS23] | Hadžić, Rein, Schrecker, Straub, *Damping versus oscillations for a gravitational Vlasov–Poisson system*, ARMA 249:45 (2025), [arXiv:2301.07662](https://arxiv.org/abs/2301.07662) | **our reduction**: all particles with the same \|L\|, point mass at the centre, small steady states ε(E₀ − E)ᵏ | sharp dichotomy in the edge exponent: k > 1 damps, 1/2 < k ≤ 1 does not (a discrete eigenvalue in the gap) |
| [HRSS24] | Hadžić, Rein, Schrecker, Straub, *Quantitative phase mixing for Hamiltonians with trapping*, [arXiv:2405.17153](https://arxiv.org/abs/2405.17153) | pure transport (no self-gravity in the perturbation) in the potential of a steady state; includes the fixed-\|L\| (1+1) case | decay rates of ∂ₜU, the force and ∂ₜρ, limited by the elliptic point (the circular orbit) |
| [HS25] | Hadžić, Schrecker, *On quantitative linear gravitational relaxation*, [arXiv:2505.14856](https://arxiv.org/abs/2505.14856) | full linearized system, point mass, distributions in (E, L) with ε(E₀−E)^μ (L−L₀)^ν, μ > 2, ν > 1 | first quantitative decay for the linearized system; the rate improves with the regularity of the steady state and of the data |
| [HM24] | Hadžić, Moreno, *On absence of embedded eigenvalues and stability of BGK waves*, [arXiv:2412.07025](https://arxiv.org/abs/2412.07025) | electrostatic plasma, periodic BGK waves with trapped electrons | no embedded eigenvalues, non-quantitative damping; the method uses the same action–angle tools |

### Numerical studies of the same group

| ref. | paper | main observations |
|---|---|---|
| [RR18] | Ramming, Rein, *Oscillating solutions of the Vlasov–Poisson system — A numerical investigation*, Physica D 365 (2018), [arXiv:1602.07989](https://arxiv.org/abs/1602.07989) | spherical perturbations of stable polytropes oscillate, periodically or with damping; an Eddington–Ritter relation between the period and the central density |
| [S24] | Straub, *Numerical experiments on stationary, oscillating, and damped spherical galaxy models*, [arXiv:2405.01235](https://arxiv.org/abs/2405.01235) | systematic study, linear and nonlinear, of spherical models with distributed L (see Section 3.4) |

## 2. The setting of [HRSS23] and ours

The radial system with fixed squared angular momentum L and a point mass M is

∂ₜf + w∂ᵣf − (U′ + M/r² − L/r³)∂_w f = 0,  ρ = (π/r²)∫f dw,

with steady states f = ε(E₀ − E)ᵏ₊ and E = w²/2 + Ψ(r), Ψ = U − M/r + L/(2r²). Their
linearized operator is the second-order Antonov operator 𝓐 = −T² − R acting on the part of
the perturbation that is odd in w; an eigenvalue λ of 𝓐 is a frequency squared, λ = ω².

| | [HRSS23] | our PIC runs (demo η, battery η) |
|---|---|---|
| background | point mass M | isochrone of unit mass and scale |
| angular momentum | squared modulus L | modulus L₀ = 2 (their L = 4) |
| steady state | ε(E₀ − E)ᵏ, a polytrope in the energy | lowered Maxwellian A·E_γ(g, (E_t − E)/T); also J²(J_t−J)² and the Gaussian in J |
| edge exponent | k | g |
| mass | ε → 0 (theorems valid for ε < ε₀(k), not estimated) | finite, a₀ = 0.003–0.5 |
| band | wide: the single-gap condition requires T_max ≥ 2T_min, so the bands k = 1, 2, … overlap and the only gap is the principal one, (0, Ω_min²) | narrow: Ω_max/Ω_min ≈ 1.2, so there are also gaps between consecutive bands |
| criterion for a discrete mode | Mathur operator M_λ (Birman–Schwinger): an eigenvalue in the gap exists iff M_λ ≥ 1 for some λ in the gap (Prop. 5.3) | λ(ω), the largest eigenvalue of the loop gain, and λ_edge = λ(Ω_min) (`lambda_borde.py`) |

**Our λ(ω) is their Mathur operator.** Their kernel is
K_λ(r,s) = (c/rs) Σ_j ∫ |ϕ′(E)| sin(2πjθ(r,E)) sin(2πjθ(s,E)) / [T(E)(4π²j²/T(E)² − λ)] dE,
with c a numerical constant and θ ∈ [0, 1) the angle,
which is the radial-kernel form of our |R|^{1/2}(−M)|R|^{1/2} with λ = ω². In both, the
operator norm grows towards the edge of the gap. For k ≤ 1 it diverges there, because
∫|ϕ′|/(Ω² − Ω_min²) dE ∝ ∫(E₀ − E)^{k−2} dE. Their proof of the dichotomy (Theorems 5.4 and
5.5) is the statement λ_edge < 1 for k > 1 and λ_edge = ∞ for k ≤ 1 at small ε. What our
numerics add is the value of λ_edge at finite mass, where it can exceed 1 also for k > 1.

**The code can run their exact setting.** `BGtype = "sphere"` is a uniform ball of mass 1
and radius 1, so it is a point mass for r > 1. With M = 1 and L₀ = 2 the circular orbit is at
r = L₀² = 4, with energy E_min = −1/8, and the single-gap condition is
−2^{−2/3}/8 = −0.079 < E₀ < 0. With E₀ = −0.05, for example, the shell extends from 2.25 to
17.7 in the limit ε → 0, so the particles never enter r < 1, but the radial grid must
extend beyond r = 20. What is missing is in the tools, which assume the isochrone:
`aa_numerico.py`, `equilibrio.py`, `lineal.py` and `lambda_borde.py` need a point-mass
background, and `equilibrio.py` needs the form ε(E₀ − E)ᵏ (the limit W₀ → 0 of the lowered
Maxwellian).

## 3. What they report that can be tested numerically

### 3.1 The dichotomy [HRSS23, Thm. 1.2]

For ε < ε₀(k): if 1/2 < k ≤ 1 there is at least one eigenvalue of 𝓐 in the principal gap, so
the perturbation does not damp; if k > 1 the point spectrum is empty and the field damps in
the time-averaged sense lim (1/T)∫₀ᵀ‖∇U‖² dt = 0 (RAGE theorem, no rate). The King model
does not damp.

- *Already seen with the isochrone*: D1 (g = 2) and D8 (g = 1) at η = 0.1; λ_edge = ∞ for
  g ≤ 1; the discrete modes of the King edge at weak coupling are bound by 10⁻⁸–10⁻²⁵ band
  widths, so in a finite run they look like a response at Ω_min that decays ever more slowly.
- *To do*: the same in their setting, point mass and ε(E₀ − E)ᵏ, scanning k = 0.75, 1, 1.5, 2, 3
  at small ε, with λ_edge, the linear solver and PIC runs.

### 3.2 The threshold mass ε₀(k), which the theorems do not estimate

For k > 1, damping is proved only below ε₀(k). λ_edge gives a numerical estimate of where it
stops: the mass at which λ_edge = 1. In the isochrone, for g = 2, it is a₀ ≈ 0.091. In their
setting it can be computed as a function of k, together with the way λ_edge diverges as
k → 1⁺.

### 3.3 No embedded eigenvalues for k > 1 [HRSS23, Thm. 4.5]

There are no eigenvalues inside the essential spectrum at small ε. In our runs all the
discrete modes lie below the band, and [S24, Obs. 4.2] reports the same for isotropic
polytropes. A numerical test is to look, at moderate mass, for undamped components at
frequencies inside the band. In our narrow-band setting the gaps between consecutive bands,
which their single-gap condition excludes, are also worth a look.

### 3.4 Pure transport [HRSS24, Thm. 1.3]

In the fixed-|L| case, with data f₀ = |ϕ′(E)|g₀, the solution of the pure transport equation
satisfies

| quantity | bound |
|---|---|
| \|∂ₜU(t, r)\|, r ≠ r* | t^{−min(2, k)}, with a constant that grows as log\|r − r*\| |
| \|∂ₜU(t, r)\|, r > r*, k ≥ 5/2 | t^{−5/2}, with a constant ∝ (r − r*)⁻¹ |
| ‖∂ₜ∂ᵣU‖∞ (the force) | t^{−min(3/2, k)}, likely optimal |
| ‖∂ₜρ‖∞ | t^{−min(1/2, k − 1/2)} |

In the 1+1 case these rates do not improve with smoother data, because of the elliptic point
(the circular orbit). The decay is faster outside r* than inside. These bounds can be
compared directly with PIC runs without self-gravity (`autointeraction = .false.`) and with
`colas_libres.py`, which already gave tails t^{−5/2} from the turning points and t^{−2} from
Ω_max in the battery.

### 3.5 Quantitative relaxation of the full linearized system [HS25, Thm. 1.1]

With a point mass and steady states ε(E₀−E)^μ(L−L₀)^ν, μ > 2, ν > 1, the force decays as
(1+t)^{−b} for b ≤ K − 1, K = min(μ − 1, ν, k), with k the regularity of the data. The rate
improves with regularity. The authors expect this improvement not to hold in the fixed-|L|
case. Our runs can measure the decay rate of the linear response in the fixed-|L| case for
k > 1, and how it degrades as λ_edge → 1.

### 3.6 Spectrum, gap and Eddington–Ritter relation [HRS21]

The essential spectrum is {4π²n²/T(E,L)²}, the principal gap is (0, 4π²/sup T²), and an
eigenvalue λ in the gap is an oscillation of period P = 2π/√λ, longer than every radial
period. For galaxy pulsations an Eddington–Ritter relation, P ρ(0)^{1/2} ≈ const, holds
along families of steady states ([HRS21, Sec. 3.4], [RR18]). In our shells the analogue is
the frequency of the discrete mode against the density of the shell; ω(a₀) is already known
for the reference family.

### 3.7 Straub's numerical observations [S24]

- Isotropic polytropes: undamped for 0 < k ≤ 1.2, fully damped for 1.3 ≤ k ≤ 3. At finite
  mass the threshold is not at k = 1. Our λ_edge shows the same effect: at a₀ ≈ 0.076,
  g = 1.5 already has a discrete mode.
- No eigenvalue embedded in the essential spectrum.
- King models: undamped for small κ, fully damped for large κ. In our model a King edge always
  has λ_edge = ∞. A possible reading, to be tested, is that the "damped" King models have a
  mode bound so weakly that no simulation can separate it from the edge, as for our G1a.
- Anisotropic polytropes (k, ℓ): undamped oscillations if and only if (π²/12)k − (π/3)ℓ < 1.
- Shells with L₀ not too small always keep an undamped part.
- The nonlinear and linearized systems behave alike for weak perturbations; pure transport
  behaves differently. In our runs this holds where ν ≪ 1; near the transition (D3, D5, D6)
  the nonlinear runs depart from linear theory.

## 4. Proposed numerical programme

1. **Their setting in our tools.** Point-mass background in `aa_numerico.py`, `equilibrio.py`,
   `lineal.py` and `lambda_borde.py`; the form ε(E₀ − E)ᵏ in `equilibrio.py`; a radial grid
   beyond r = 20 in the PIC runs. Check the single-gap condition and the monotonicity of the
   period function T(E), which their proofs need.
2. **Linear map.** λ_edge(k, ε) for k = 0.75–3, the threshold ε₀(k) for k > 1, the binding of
   the mode for k ≤ 1, and the comparison with the time-domain linear solver.
3. **PIC runs for intuition.** At small ε, k = 0.75, 1, 1.5 and 2: damping against
   oscillation in δΦ, h₁ and the phase space; the same runs without self-gravity to compare
   the decay with [HRSS24].
4. **Beyond the theorems.** Masses above ε₀(k), where the dichotomy is not proved; the gaps
   between bands in a narrow-band setting; the nonlinear regime near the threshold.
