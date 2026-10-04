# Literature review: oscillation versus damping in gravitational Vlasov–Poisson systems

4 October 2026

## 1. Question and scope

The η battery, the η demo and the report on the setting of Hadžić, Rein, Schrecker and
Straub contain six results (Section 4). This review asks which of them are known, which
are new, and where an extension would answer a question that the literature leaves open.

Sources:

- the reference lists of seven papers of the Hadžić–Rein group (about 190 distinct works);
- the works that cite them, from OpenAlex and Semantic Scholar on 3–4 October 2026
  (HRS22: 18 citing works, HRSS23: 12, Kunze's notes: 9, S24: 3, HS25: 3, HM24: 2, HRSS24: 1);
- targeted searches in mathematics, astrophysics, statistical physics and spectral theory;
- two doctoral theses from Bayreuth (Straub 2024, Wolfschmidt 2023).

Depth of reading. Sections of the full text were read for HRS22, HRSS23, HRSS24, HS25, RR18,
S24, the two theses and Kunze's lecture notes. All other works were read in abstract only.
Not accessed: the full text of Mathur (1990), Weinberg (1991), Louis (1992), Klaus and Simon
(1980), Sellwood and Pryor (1998), and Kunze's book.

## 2. The literature

### 2.1 The Bayreuth–London line

| work | content | relation to our results |
|---|---|---|
| Mathur 1990 | Radial perturbations of spherical systems and the 1D system: a reduction of the eigenvalue problem that shows that oscillation modes exist | origin of the operator whose norm is our λ |
| Ramming, Rein 2018 [RR18] | Nonlinear PIC runs near polytropes: oscillations, periodic or damped; period–density relation | first numerical evidence of both behaviours |
| Hadžić, Rein, Straub 2022 [HRS22] | Essential spectrum of the linearized operator from the radial periods; rigorous Birman–Schwinger–Mathur principle; eigenvalue in the principal gap for some steady states | the criterion λ_edge > 1 is their principle at the edge of the gap |
| Kunze 2021 (book), 2022 (notes) | Birman–Schwinger principle for the Antonov operator; sufficient criteria for an eigenvalue | the notes end with open questions; the first is "Do some numerics" |
| Günther, Rein, Straub 2022/2025 | Birman–Schwinger principle for Einstein–Vlasov; stable shells around a black hole | same method in general relativity |
| Hadžić, Rein, Schrecker, Straub 2023/2025 [HRSS23] | Point mass, fixed \|L\|: damping for k > 1 and oscillation for 1/2 < k ≤ 1, below a mass ε₀(k) that is not estimated | the setting of our report |
| Moreno, Rioseco, Van Den Bosch 2023 | Absolutely continuous spectrum; criteria for the existence and number of oscillating modes, with an attractive external potential | analytic criteria; no numerical evaluation (0 occurrences of "numeric") |
| Straub 2024 [S24] and thesis | Time-domain study, linearized and nonlinear, of polytropes, King models and shells; thresholds found by inspection of the time series | states the open problems of Section 3 |
| Wolfschmidt 2023 (thesis) | The same for Einstein–Vlasov; threshold near k = 1.2 at small redshift | same ambiguity near the threshold |
| Straub, Wolfschmidt 2024 | Neural network that predicts the stability of Einstein–Vlasov steady states from 10⁴ models | a data-driven classifier already exists in this group |
| Hadžić, Rein, Schrecker, Straub 2024 [HRSS24] | Pure transport: decay bounds, t^{−min(2,k)} for ∂ₜU in the fixed-\|L\| case | upper bounds; no self-gravity |
| Hadžić, Schrecker 2025 [HS25] | Linearized system with a point mass, distributed L: quantitative decay, which improves with regularity | first decay rates with self-gravity |
| Hadžić, Moreno 2024 [HM24] | No embedded eigenvalues for BGK waves; damping without a rate | same method in plasmas |
| Rein 2023; Pausader 2026 | Reviews | context |

After 2024 Straub's papers are on neural operators and physics-informed networks. S24 has
three citing works (Moreno et al., Chaturvedi and Luk, Kunze and Ortega); none is numerical.

### 2.2 Other mathematical work

| work | content |
|---|---|
| Rioseco, Sarbach 2020; Moreno, Rioseco, Van Den Bosch 2022 | Phase mixing in a central potential and in an anharmonic well, with rates |
| Chaturvedi, Luk 2021 | 1D transport in a confining potential: decay t⁻² of ∂ₜφ |
| Chaturvedi, Luk 2024 | Vlasov–Poisson with an external Kepler potential: linear phase mixing in 3D and nonlinear phase mixing in spherical symmetry on long time scales |
| Weder 2025 | Scattering theory for the plane-symmetric Antonov operator around polytropes and King states |
| Nguyen 2023 | Homogeneous plasmas with compact velocity support: a survival threshold κ₀² = 4π∫u²μ/(Υ² − u²)du; oscillations below it; the decay law at the threshold depends on how μ vanishes at the maximal speed Υ |
| Bedrossian, Masmoudi, Mouhot 2022; Han-Kwan, Nguyen, Rousset 2021 | Linearized Vlasov–Poisson on ℝ³: Klein–Gordon-type oscillations at low wave number |
| Faou, Horsin, Rousset 2021; Faou, Rousset, Seetohul 2026 | Inhomogeneous states of the HMF model: algebraic linear damping; nonlinear dynamics for long times |

Nguyen's threshold is the homogeneous analogue of λ_edge = 1: an integral over the
equilibrium with a singularity at the edge of the band, finite or not according to how the
equilibrium vanishes there.

### 2.3 Statistical physics

Barré, Olivetti and Yamaguchi (2010, 2011) study perturbations of inhomogeneous states of
the HMF model. After a transient Landau damping, the perturbation decays algebraically with
exponent −2 and a definite frequency. They test this with a weighted-particle code that
reduces the finite-size noise. Barré and Yamaguchi (2013) classify the branch singularities in
two dimensions, which covers spherical self-gravitating systems, and test the result on an
advection equation for the isochrone model.

### 2.4 Astrophysics

**Radial pulsations and modes of spheres.** Louis and Gerhard (1988) build a self-consistent
nonlinear oscillating mode of the isochrone sphere. Sellwood and Pryor (1998) study the
pulsation modes of spherical systems. Vandervoort (1983, 2003, 2004) formulates matrix methods
for radial oscillations and for stationary oscillations of the van Kampen type. Weinberg (1994) finds weakly damped l = 1 modes in King models;
Heggie, Breen and Varri (2020) see them in N-body runs; Fouvry and Prunet (2022) and Petersen
et al. (2024, LinearResponse.jl) compute damped modes by analytic continuation of the
response matrix. Lau and Binney (2021) and Polyachenko, Shukhman and Borodina (2021) discuss
van Kampen modes against Landau-damped waves. These works look for weakly damped modes by
continuation of the response matrix, mostly for l = 1. As far as the abstracts show, none
evaluates the self-adjoint criterion below the band for radial modes, or relates the result
to the exponent of the edge.

**The one-dimensional slab.** This is the same class as our fixed-|L| model: one degree of
freedom, a band of orbital frequencies, an external potential and self-gravity. Mathur (1990)
treats it together with the radial case. Weinberg (1991) shows by linear analysis and N-body
runs that a 1D disc may have undamped modes. Widrow and Bonner (2015) find normal modes or
Landau-damped oscillations depending on the distribution function. Darling and Widrow (2019)
apply dynamic mode decomposition to simulations of the isothermal plane and find nearly
undamped modes when the external potential is about four times the self-gravity. The topic
is active because of the Gaia phase spiral: Widrow (2023), Asano and Antoja (2025), who measure
a delay of the phase mixing caused by self-gravity, and the workshop summary of Frankel et
al. (2026).

### 2.5 Spectral theory

Simon (1976) shows that a weak attractive potential always binds in one and two dimensions,
with binding energy of order λ² and exp(−c/λ). Klaus and Simon (1980) study the eigenvalue
e(λ) that is absorbed at a threshold coupling λ₀ and find a leading power that depends on the
dimension. The standard results are e ∝ (λ − λ₀)² in three dimensions, (λ − λ₀)/|ln(λ − λ₀)|
in four, and linear above; these exponents are quoted from memory and must be checked
against the paper. Fassari and Klaus (1998) treat thresholds at the band edges of periodic
operators.

### 2.6 Numerical methods

Weighted particles for linear response (Barré, Olivetti, Yamaguchi 2011), perturbation
particles (Leeuwin, Combes, Binney 1993) and δf particle-in-cell methods carry the perturbation
in the weights, as our quiet start does. Dynamic mode decomposition has been used on
collisionless slabs (Darling, Widrow 2019). Steady states as fixed points of a
mass-preserving iteration are analysed by Andréasson, Kunze and Rein (2024).

## 3. Open problems stated in the literature

1. **Numerical evaluation of the Mathur operator.** Kunze's notes (2022) list "Do some
   numerics" as the first open question. Straub's thesis (2024, Outlook VII, "Numerical
   Analysis of the Mathur Operator") says it "would certainly be helpful to be able to
   numerically analyse the Mathur operator M_λ and the number M".
2. **The threshold exponent of isotropic polytropes.** S24 finds undamped oscillations for
   0 < k ≤ 1.2 and damping for k ≥ 1.3, and does not classify k = 1.25: "the damping is
   rather slow and there might also be an undamped part". The thesis gives k* ≈ 1.25 and
   notes that proving monotonicity "could be easier than actually determining whether M is
   > 1 or ≤ 1 for a fixed steady state". For Einstein–Vlasov, Wolfschmidt writes that "it is
   not feasible to numerically determine the threshold value of k", and that they "do not
   know whether k = 1.2 is of any deeper analytical meaning or if it is an artifact of
   numerical inaccuracy".
3. **The linear relation between exponents.** For anisotropic polytropes
   f ∝ (E₀ − E)ᵏ Lˡ, S24 (Observation 4.4) gives the empirical boundary
   (π²/12)k − (π/3)ℓ < 1, "fitted by manual trial and error". For ℓ > 0 the essential
   spectrum is [0, ∞) and the modes are embedded.
4. **The rigorous criteria are far from sharp.** The sufficient criterion of Kunze holds for
   isotropic polytropes only for k below about 0.03, and for shells with L₀ not too small
   (S24, Section 4.3). The observed threshold is k ≈ 1.2.
5. **Small shells around a point mass** with k + ℓ ≤ 0: existence of pulsations not proved
   (thesis, Outlook VIII). The threshold ε₀(k) of HRSS23 is not estimated.
6. **Decay rates.** HRSS24 call their rate for the force "likely optimal". HS25 state that
   the improvement of the decay with regularity is unlikely in 1+1 dimensions, and that it is
   unclear whether the rate t⁻² of the radial case can be improved without the point mass.

## 4. Where our results stand

| our result | closest prior work | assessment |
|---|---|---|
| R1. λ_edge, the norm of the Mathur operator at the edge of the gap, computed numerically for 61 equilibria (isochrone and point mass, fixed \|L\|) and checked against the linear solver and PIC runs | principle: Mathur 1990, HRS22, Kunze 2021; numerics: none found | apparently new; it is open problem 1. The fixed-\|L\| model limits its reach |
| R2. Threshold masses of the HRSS23 setting: a₀ ≈ 0.53 (k = 1.25), 0.93 (k = 1.5), above 1 for k ≥ 2 | HRSS23 (ε₀ not estimated) | apparently new; of interest mainly to that group |
| R3. Edge law λ_edge − λ(Ω_min − δ) ∝ δ^{min(k−1,1)} and the binding laws, with c the coupling (proportional to the mass): δ_d ∝ c^{1/(1−k)} for k < 1, exp(−const/c) for k = 1, (λ_edge − 1)^{1/(k−1)} for 1 < k < 2 | the divergence for k ≤ 1 is in the proof of HRSS23; the laws are those of Simon and Klaus–Simon with d = 2k; Nguyen 2023 in homogeneous plasmas | the analogy and the laws do not appear in HRS22, HRSS23, Kunze's notes, Straub's thesis or Moreno et al. (no occurrence of "weak coupling", "coupling constant" or Klaus). Apparently new in this setting; the mathematics is classical |
| R4. Observability: the mode separates from the edge after t ≈ 1/(Ω_min − ω_d); at small mass the modes of k ≤ 1 lie less than 10⁻⁶ band widths below the edge, and every k decays as t^{−(k+1)} | none found | apparently new, and it explains open problem 2: near the threshold the binding is a high power of the distance to it, so no time integration can decide |
| R5. Algebraic tails t^{−(k+1)} with self-gravity | Barré et al. 2011 (exponent −2), Barré, Yamaguchi 2013, HRSS24 bounds, HS25 | partly known; a numerical complement to open problem 6 |
| R6. PIC results: quiet start in angle–action variables with the D − Z subtraction; grid and particle noise separated; finite-amplitude loss of the modes next to the edge | weighted particles (Barré et al. 2011), perturbation particles (Leeuwin et al. 1993), δf methods | known ideas; supporting material |
| R7. η is not universal | Darling, Widrow 2019 use a similar ratio (external potential to self-gravity, about 4:1) | a negative result; useful as a remark: such a ratio cannot be the criterion |

Dynamic mode decomposition and a neural-network classifier have both been applied to this
problem (Darling and Widrow 2019; Straub and Wolfschmidt 2024), so neither would be new as a
method.

## 5. Possible papers

**A. Numerical analysis of the Mathur operator.** Compute λ_edge for the steady states of
S24: isotropic polytropes, King models and shells. The expected results are the threshold
exponent k* as the root of λ_edge(k) = 1, a test of the value 12/π², the classification of
k = 1.25, and, with R3 and R4, the reason why the time-domain studies cannot resolve it. R1
and R2 give the validation in the fixed-|L| case, where theory, linear solver and PIC runs
agree. Work needed: extend `lambda_borde.py` from one action to two, (J_r, L), for radial
perturbations of isotropic states. The anisotropic relation of open problem 3 involves
embedded eigenvalues and is outside the criterion below the band; it would be an outlook.
Natural journals: Physica D (RR18 and S24 appeared there) or Classical and Quantum Gravity.

**B. Persistent modes of the vertical slab.** The slab is a one-action problem, so the
present code applies after changing the Green function and the background. The question,
when a self-gravitating slab in an external potential keeps an undamped mode instead of
winding into a phase spiral, is of current interest (Gaia). Weinberg (1991), Widrow and Bonner
(2015) and Darling and Widrow (2019) found modes case by case; a criterion with the laws R3
was not found. This option needs its own review of the phase-spiral literature (Banik,
Weinberg and van den Bosch; Widrow 2023) before any work.

**C. A note on decay exponents** (R5) is possible but small.

Option A answers a question that the authors of the field have written down, and most of
its validation is done. A feasibility test on isotropic polytropes decides it.

## 6. Limits of this review

- Mathur (1990) was read only in abstract. If that paper already evaluates the criterion for
  polytropes, R1 loses part of its novelty. HRS22 and S24, which build on it, do not mention
  such an evaluation.
- Kunze's book was not read; the lecture notes, written a year later, still ask for
  numerics.
- The citing lists come from two indexes that disagree (OpenAlex finds 6 citing works for
  HRSS23, Semantic Scholar 12). Works from 2026 may be missing.
- The astrophysical literature on radial modes before 2000 was covered through the reference
  lists of HRS22 and S24, not searched independently.

## References

Group of Hadžić and Rein, and Bayreuth:

- [RR18] Ramming, Rein, Physica D 365 (2018), [arXiv:1602.07989](https://arxiv.org/abs/1602.07989)
- [HRS22] Hadžić, Rein, Straub, ARMA 243 (2022), [arXiv:2102.11672](https://arxiv.org/abs/2102.11672)
- [HRSS23] Hadžić, Rein, Schrecker, Straub, ARMA 249:45 (2025), [arXiv:2301.07662](https://arxiv.org/abs/2301.07662)
- [HRSS24] Hadžić, Rein, Schrecker, Straub, [arXiv:2405.17153](https://arxiv.org/abs/2405.17153)
- [HM24] Hadžić, Moreno, Ann. Henri Poincaré (2026), [arXiv:2412.07025](https://arxiv.org/abs/2412.07025)
- [HS25] Hadžić, Schrecker, [arXiv:2505.14856](https://arxiv.org/abs/2505.14856)
- [S24] Straub, Physica D (2024), [arXiv:2405.01235](https://arxiv.org/abs/2405.01235)
- Straub, *Pulsating Galaxies*, thesis, Bayreuth 2024, [epub 7639](https://epub.uni-bayreuth.de/id/eprint/7639)
- Wolfschmidt, *Stability and Oscillations of Star Clusters in General Relativity*, thesis, Bayreuth 2023, [epub 7337](https://epub.uni-bayreuth.de/id/eprint/7337)
- Kunze, *A Birman–Schwinger Principle in Galactic Dynamics*, Birkhäuser 2021; lecture notes, Class. Quantum Grav. (2022), [arXiv:2202.04555](https://arxiv.org/abs/2202.04555)
- Günther, Rein, Straub, ARMA (2025), [arXiv:2204.10620](https://arxiv.org/abs/2204.10620)
- Günther et al., Class. Quantum Grav. (2020), [arXiv:2009.08163](https://arxiv.org/abs/2009.08163)
- Günther, Straub, Rein, ApJ (2021), [arXiv:2105.05556](https://arxiv.org/abs/2105.05556)
- Straub, Wolfschmidt, Class. Quantum Grav. 41 (2024), [arXiv:2310.08253](https://arxiv.org/abs/2310.08253)
- Rein, Class. Quantum Grav. (2023), [arXiv:2305.02098](https://arxiv.org/abs/2305.02098)
- Andréasson, Kunze, Rein, [arXiv:2412.01544](https://arxiv.org/abs/2412.01544)
- Kunze, Ortega, [arXiv:2512.07746](https://arxiv.org/abs/2512.07746)
- Mathur, MNRAS 243, 529 (1990), [record](https://www.osti.gov/etdeweb/biblio/6932499)

Other mathematics:

- Moreno, Rioseco, Van Den Bosch, [arXiv:2201.07019](https://arxiv.org/abs/2201.07019) and [arXiv:2305.05749](https://arxiv.org/abs/2305.05749)
- Chaturvedi, Luk, [arXiv:2109.12402](https://arxiv.org/abs/2109.12402) and [arXiv:2409.14626](https://arxiv.org/abs/2409.14626)
- Weder, [arXiv:2501.04175](https://arxiv.org/abs/2501.04175)
- Nguyen, J. Funct. Anal. (2026), [arXiv:2305.08672](https://arxiv.org/abs/2305.08672)
- Faou, Horsin, Rousset, [arXiv:2105.02484](https://arxiv.org/abs/2105.02484); Faou, Rousset, Seetohul, [arXiv:2609.16998](https://arxiv.org/abs/2609.16998)
- Rioseco, Sarbach, Class. Quantum Grav. 37, 195027 (2020)
- Bedrossian, Masmoudi, Mouhot, SIAM J. Math. Anal. 54 (2022); Han-Kwan, Nguyen, Rousset (2021)

Statistical physics:

- Barré, Olivetti, Yamaguchi, J. Stat. Mech. P08002 (2010); J. Phys. A 44, 405502 (2011), [arXiv:1104.1890](https://arxiv.org/abs/1104.1890)
- Barré, Yamaguchi, J. Phys. A 46, 225501 (2013), [arXiv:1210.8040](https://arxiv.org/abs/1210.8040)

Astrophysics:

- Weinberg, ApJ 373, 391 (1991), [ADS](https://ui.adsabs.harvard.edu/abs/1991ApJ...373..391W/abstract); ApJ 421, 481 (1994), [arXiv:astro-ph/9306020](https://arxiv.org/abs/astro-ph/9306020)
- Louis, Gerhard, MNRAS 233, 337 (1988), [journal](https://academic.oup.com/mnras/article/233/2/337/969454); Louis, MNRAS 258, 552 (1992), [journal](https://academic.oup.com/mnras/article/258/3/552/1081867)
- Sellwood, Pryor, Highlights of Astronomy 11 (1998), [chapter](https://link.springer.com/chapter/10.1007/978-94-011-4778-1_10)
- Vandervoort, MNRAS 339, 537 (2003), [arXiv:astro-ph/0207259](https://arxiv.org/abs/astro-ph/0207259); MNRAS (2004), [doi](https://doi.org/10.1111/j.1365-2966.2004.07361.x)
- Heggie, Breen, Varri, MNRAS 492 (2020), [arXiv:2002.00783](https://arxiv.org/abs/2002.00783)
- Fouvry, Prunet, MNRAS 509 (2022), [arXiv:2105.01371](https://arxiv.org/abs/2105.01371)
- Petersen, Roule, Fouvry, Pichon, Tep, [arXiv:2311.10630](https://arxiv.org/abs/2311.10630)
- Lau, Binney, [arXiv:2104.07044](https://arxiv.org/abs/2104.07044) and [arXiv:2106.10297](https://arxiv.org/abs/2106.10297)
- Polyachenko, Shukhman, Borodina, [arXiv:2101.08287](https://arxiv.org/abs/2101.08287)
- Hamilton, Fouvry, [arXiv:2402.13322](https://arxiv.org/abs/2402.13322)
- Widrow, Bonner, MNRAS 450 (2015), [arXiv:1503.05741](https://arxiv.org/abs/1503.05741)
- Darling, Widrow, MNRAS (2019), [doi](https://doi.org/10.1093/mnras/stz2539)
- Widrow, MNRAS 522 (2023), [arXiv:2302.14524](https://arxiv.org/abs/2302.14524)
- Asano, Antoja, A&A (2026), [arXiv:2510.11801](https://arxiv.org/abs/2510.11801)
- Frankel et al., [arXiv:2603.09015](https://arxiv.org/abs/2603.09015)

Spectral theory:

- Simon, Ann. Phys. 97, 279 (1976)
- Klaus, Simon, Ann. Phys. 130, 251 (1980), [record](https://authors.library.caltech.edu/records/95ns3-zwc43)
- Fassari, Klaus, J. Math. Phys. 39, 4369 (1998), [journal](https://pubs.aip.org/aip/jmp/article-abstract/39/9/4369/441228/Coupling-constant-thresholds-of-perturbed-periodic)
