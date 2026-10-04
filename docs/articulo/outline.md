# Outline of an article that joins the fixed-|L| and the dispersed-L codes

4 October 2026. A proposal to discuss with the coauthors. It builds on the draft
`VlasovPoisson_PIC_sp/Vlasov_Poisson_evolutions/main.md`, keeps its Sections 2–4, and
replaces its results by what was learned since.

## 1. The article in one paragraph

Working title: *Phase mixing against self-gravity for a spherical collisionless gas in a
central potential: fixed and dispersed angular momentum.*

Without self-gravity a gas in a central potential mixes in phase space (Rioseco and Sarbach
2020). With self-gravity the response to a perturbation of an equilibrium is Landau damped
below a threshold in mass and oscillates without damping above it. With all particles at
the same angular momentum the threshold is where the gain of the self-consistent loop
reaches one. A spread in angular momentum removes the oscillation that fixed |L| allows.
The article shows this with a particle-in-cell code, a linearized solver and the loop gain,
which agree with each other.

## 2. What the contribution is, and what it is not

Small and honest:

1. A quantitative account of where self-gravity stops the mixing of a shell in an isochrone
   or point-mass potential, with three methods that agree: particles, a linearized solver in
   time, and an eigenvalue criterion.
2. The effect of a spread in L: it suppresses the unwinding of a phase-space spiral and
   removes the divergence of the gain at the edge of the band, so that the mode of the
   fixed-|L| case disappears.
3. The numerical conditions needed to measure these effects with particles.

Not new, and to be said so: phase mixing and Landau damping; the Birman–Schwinger–Mathur
criterion (Mathur 1990; Hadžić, Rein and Straub 2022); perturbations carried by the weights
(Barré et al. 2011; Leeuwin et al. 1993).

## 3. Sections

| # | section | content | source that already exists |
|---|---|---|---|
| 1 | Introduction | mixing and its role; what self-gravity changes (damping or oscillation: Hadžić et al., Straub); the two reductions, fixed \|L\| and dispersed L; what is done here | draft §1 (comments of Olivier); `docs/hadzic/literature_review.md` |
| 2 | Model | Vlasov–Poisson with an external potential in spherical symmetry; fixed \|L\| (one degree of freedom) and a distribution in L; units; angle–action variables of the total potential; the observables h_k and δΦ | draft §2–3, shortened |
| 3 | Linear theory | linearized equation; free mixing and the exponent of the edge; the loop gain λ(ω) and the criterion λ_edge > 1; the same with a distribution in L (the measure becomes L dL dJ) | demo §1.5; `docs/mathur/threshold_isotropic.md` §2 |
| 4 | Numerical methods | particles as shells, grid, deposit, integrator; quiet start on a lattice in the angle–action variables of the self-consistent equilibrium; perturbation in the weights and reference run; construction of the equilibria; linearized solver; evaluation of the operator | draft §4; demo §2.1; report of the Hadžić setting §2 |
| 5 | Validation | (a) no self-gravity: exact solution, both codes; (b) stationarity of the equilibria, error of order Δr²; (c) particles against the linearized solver, 1 %; (d) the operator with dispersion in the limit of fixed L, 10⁻⁵; (e) how long the quiet start lasts against N | draft §4.4; demo §5; `lambda_L.py validar`; `relajacion_N.py` |
| 6 | Results with fixed \|L\| | 6.1 weak self-gravity: Landau damping, ω and γ; 6.2 the threshold: λ_edge against the mass, and the discrete mode; 6.3 the edge: exponent, algebraic tails, how weakly the mode is bound; 6.4 finite amplitude, briefly | intro document of the fixed-L code; battery η §3.9; demo §3; report of the Hadžić setting §3 |
| 7 | Results with a spread in L | 7.1 free mixing: the unwinding of a spiral is suppressed; 7.2 the gain with dispersion: finite at the edge, threshold raised; 7.3 particles with dispersion: the mode of fixed L becomes damped | `_sp` intro document, section on the spiral; `lambda_L.py dispersion`; 7.3 is to be done |
| 8 | Discussion | limits (only radial perturbations, noise, range of L); relation with the theorems (fixed L: Hadžić et al. 2023; distributed L: Hadžić and Schrecker 2025); outlook | reviews in `docs/` |
| 9 | Conclusions | the three points of Section 2 | — |
| A | Appendices | the numerical angle–action map; normalization of the weights; reproduction | demo §2.1; scripts |

Length: about 20 pages.

## 4. Numbers that carry each section

**Section 5.**

- Operator with dispersion against fixed L: 0.497963 against 0.497970, 0.738046 against
  0.738042.
- Quiet start at a₀ = 1: usable to t ≈ 1500 with 10⁴ particles, beyond 4000 with 4 × 10⁴.

**Section 6.**

- Landau damping at a₀ = 0.01: ω = 0.05705, γ = 5.08 × 10⁻³ with particles; γ = 5.01–5.07 × 10⁻³
  with the linearized solver.
- Threshold masses. Isochrone, Wilson model: 0.091 (reference width), 0.042 (narrow support),
  above 0.22 (wide). Point mass: 0.53 (k = 1.25), 0.93 (k = 1.5), above 1 (k ≥ 2).
- Discrete modes, particles against λ(ω) = 1: 0.14156 against 0.14169 (k = 0.75) and 0.14766
  against 0.14785 (k = 1); damping of k = 2, 6.8 × 10⁻³ against 6.5 × 10⁻³.
- Small mass: every k decays as t^{−(k+1)}; the mode of k ≤ 1 lies less than 10⁻⁶ band widths
  below the edge and is not observable.
- Finite amplitude: at t = 500 the amplitude is 0.97, 0.84 and 0.38 of the linear one for
  ε = 0.03, 0.1 and 0.3.

**Section 7.**

- Unwinding of a spiral without self-gravity: |h₁| grows by 2.2 with fixed L and falls to
  2.2 × 10⁻³ with σ_L = 0.2.
- Gain at the edge, k = 0.75: with fixed L it diverges (14.4 at 10⁻⁵ of the edge); with a
  spread of 3 % it is 3.07 at a₀ = 1 and 0.77 at a₀ = 0.1, where the mode is lost.

## 5. Figures and tables

| | content | state |
|---|---|---|
| F1 | phase space of a pulse, free mixing | draft |
| F2 | h_k, exact against particles, both codes | draft; `_sp` |
| F3 | Landau damping: \|h₁\| of particles and of the linearized solver | exists |
| F4 | λ_edge against the mass for the families, with the threshold | exists (`lambda_mapa.pdf` and battery) |
| F5 | discrete mode: \|h₁\| does not decay; damped case beside it | exists (`pic_a1.pdf`) |
| F6 | finite amplitude: envelopes for three ε, and with Δr/2 and 4N | exists (`pic_barrido.pdf`) |
| F7 | unwinding with and without a spread in L | exists in the `_sp` document |
| F8 | λ near the edge for several spreads | to plot from `exe/lambda_L/dispersion.txt` |
| F9 | particles with a spread in L: mode against damping | to be done |
| T1 | numerical parameters of the runs | to assemble |
| T2 | validation errors | to assemble |
| T3 | thresholds and frequencies of the modes, three methods | exists |
| T4 | gain against the spread | exists |

## 6. What is missing

| | task | why | effort |
|---|---|---|---|
| 1 | equilibria with a compact edge in energy and a spread in L, in `equilibrio_L.py` (the forms exist in the fixed-L code) | Section 7.3 and a self-consistent Section 7.2 | 3–5 days |
| 2 | four to six runs of the code with L: one mass, spreads from 0 to 30 % | F9 | 1–2 days of computing |
| 3 | the gain with dispersion with the potential recomputed | the present table uses the potential of fixed L | 1–2 days |
| 4 | repeat the validation and the Landau run of the draft with the audited codes, in common units | the draft predates the fixes, and uses the isochrone map, which leaves a static residue | 2–3 days |
| 5 | convergence in N, Δr and Δt of each headline number | most exist; put them in one table | 1–2 days |
| 6 | read Mathur (1990) and Kunze's book | to state correctly what is known about the criterion | — |

If tasks 1–2 turn out harder than expected, Section 7.3 is dropped and Section 7 rests on
the exact free solution and on the operator. The article is then weaker but still complete.

## 7. What stays out

- The parameter η: one sentence, as a quantity that orders the transition within a family
  but does not locate it.
- The design of the battery of 30 runs, the videos (supplementary material at most), the
  relaxation against N.
- The non-equilibrium blob of the draft (its Sections 5.1–5.2). It is a different question,
  the self-binding of a clump, and its classification into mixed and not mixed used the
  isochrone map. To decide with the coauthors whether it becomes a short section or stays out.
- The thresholds of isotropic polytropes and King models (k* = 1.2425). They do not need the
  particle code and address another readership; a second, short article. Here, one sentence
  in the outlook.

## 8. Statements to avoid

- That the runs reproduce the dichotomy at k = 1: at small mass it is not observable.
- Anything at late times below the noise of the reference run.
- That the method is new.
- That the result of Section 7.2 is self-consistent, until task 3 is done.

## 9. Journals

Classical and Quantum Gravity, where the article of Rioseco and Sarbach on mixing appeared,
as its continuation with self-gravity; or Physica D, where the numerical studies of this
problem appeared (Ramming and Rein 2018; Straub 2024).
