# Literature review: oscillation versus damping in gravitational Vlasov–Poisson systems

4 October 2026; revised 8 October 2026

## 1. Question and scope

The η battery, the η demo and the report on the setting of Hadžić, Rein, Schrecker and
Straub contain the results listed in Section 4. This review asks which of them are known, which
are new, and where an extension would answer a question that the literature leaves open.

Sources:

- the reference lists of seven papers of the Hadžić–Rein group (about 190 distinct works);
- the works that cite them, from OpenAlex and Semantic Scholar on 3–4 October 2026
  (HRS22: 18 citing works, HRSS23: 12, Kunze's notes: 9, S24: 3, HS25: 3, HM24: 2, HRSS24: 1);
- targeted searches in mathematics, astrophysics, statistical physics and spectral theory;
- two doctoral theses from Bayreuth (Straub 2024, Wolfschmidt 2023);
- in the revision of 8 October: about one hundred queries to the arXiv interface on 7 and 8
  October 2026, the reference lists of S24, Widrow and Bonner (2015) and Karpov et al.
  (2021), and searches in plasma physics, accelerator physics, vortex dynamics and cold
  atoms.

Revision of 8 October 2026. Sections 2.7 to 2.11 are new, Sections 2.1 to 2.6 and 3 were
extended, and the table of Section 4 was extended and corrected for two results that changed
on 7 October:
the threshold masses, after the quadrature at the edge of the support was corrected, and the
exponent of the algebraic tail with self-gravity, which is −k and not −(k + 1). The reading of
Mathur (1990), Weinberg (1991) and Louis (1992) changed the assessment of several results:
for the one-dimensional slab, the numerical location of the modes in the gaps, the models at
which they enter the continuum and the limit of their amplitude were already there
(Section 2.4).

Depth of reading. Sections of the full text were read for HRS22, HRSS23, HRSS24, HS25, RR18,
S24, the two theses and Kunze's lecture notes, and, in the revision, for Moreno et al.
(2023), Nguyen (2023), Chaturvedi and Luk (2024), Barré et al. (2011), Barré and Yamaguchi
(2013), Widrow and Bonner (2015), Fouvry and Prunet (2022), Petersen et al. (2024),
Polyachenko and Shukhman (2026), Stucchi and Lauber (2025), Burov (2012), Karpov et al.
(2021), the three papers on which the criterion rests in one dimension, Mathur (1990),
Weinberg (1991) and Louis (1992), and Hénon (1973) and Louis and Gerhard (1988). All other
works were read in abstract only. The works marked with † in the list of
references were seen only through the summary of the abstract returned by a search engine,
and must be checked before they are cited. Not accessed: the full text of Klaus and Simon
(1980) and of Sellwood and Pryor (1998), and Kunze's book.

## 2. The literature

### 2.1 The Bayreuth–London line

| work | content | relation to our results |
|---|---|---|
| Mathur 1990 | Radial perturbations of spherical systems and the 1D system, with a fixed external field: a mode with real frequency in a gap exists if and only if a symmetric Hilbert–Schmidt kernel has the eigenvalue one. Example with a distribution that vanishes linearly at the edge (Section 2.4) | origin of the operator whose norm is our λ. The example is the case k = 1 of our edge law |
| Ramming, Rein 2018 [RR18] | Nonlinear PIC runs near polytropes: oscillations, periodic or damped; period–density relation | first numerical evidence of both behaviours |
| Hadžić, Rein, Straub 2022 [HRS22] | Essential spectrum of the linearized operator from the radial periods; rigorous Birman–Schwinger–Mathur principle; eigenvalue in the principal gap for some steady states | the criterion λ_edge > 1 is their principle at the edge of the gap |
| Kunze 2021 (book), 2022 (notes) | Birman–Schwinger principle for the Antonov operator; sufficient criteria for an eigenvalue | the notes end with open questions; the first is "Do some numerics" |
| Günther, Rein, Straub 2022/2025 | Birman–Schwinger principle for Einstein–Vlasov; stable shells around a black hole | same method in general relativity |
| Hadžić, Rein, Schrecker, Straub 2023/2025 [HRSS23] | Point mass, fixed \|L\|: damping for k > 1 and oscillation for 1/2 < k ≤ 1, below a mass ε₀(k) that is not estimated | the setting of our report |
| Moreno, Rioseco, Van Den Bosch 2023 | Absolutely continuous spectrum; criteria for the existence and number of oscillating modes, with an attractive external potential | analytic criteria; no numerical evaluation (0 occurrences of "numeric") |
| Straub 2024 [S24] and thesis | Time-domain study, linearized and nonlinear, of polytropes, King models and shells; thresholds found by inspection of the time series | states the open problems of Section 3 |
| Wolfschmidt 2023 (thesis) | The same for Einstein–Vlasov; threshold near k = 1.2 at small redshift | same ambiguity near the threshold |
| Straub, Wolfschmidt 2024 | Neural network that predicts the stability of Einstein–Vlasov steady states from 10⁴ models | a data-driven classifier already exists in this group |
| Hadžić, Rein, Schrecker, Straub 2024 [HRSS24] | Pure transport with initial data \|φ′(E)\| g₀: decay bounds, t^{−min(2,k)} for ∂ₜU in the fixed-\|L\| case (Theorem 1.3) | upper bounds; no self-gravity. Data that vanish at the edge like F′ already decay as t^{−k} |
| Hadžić, Schrecker 2025 [HS25] | Linearized system with a point mass, distributed L: quantitative decay, which improves with regularity | first decay rates with self-gravity |
| Hadžić, Moreno 2024 [HM24] | No embedded eigenvalues for BGK waves; damping without a rate | same method in plasmas |
| Rein 2023; Pausader 2026 | Reviews | context |

After 2024 Straub's papers are on neural operators and physics-informed networks. S24 has
three citing works (Moreno et al., Chaturvedi and Luk, Kunze and Ortega); none is numerical.

The quantity that we call λ_edge is defined in this line as a supremum. In HRS22 (Theorem
8.11 and Remark 8.12) the norm M_λ of the Mathur operator increases with λ in the principal
gap G = ]0, min σ_ess[, an eigenvalue exists in G if and only if M_λ ≥ 1 for some λ in G, and
the number M := sup_{λ∈G} M_λ decides: M > 1 gives at least one eigenvalue and M < 1 none;
for M = 1 the answer depends on whether the supremum is attained. Kunze's notes (Lemma 9.11
and Theorem 10.3) write μ₁(λ) for the largest eigenvalue of the Birman–Schwinger operator,
show that it is increasing and convex, define μ_* as its limit when λ tends from below to
δ₁² = min σ_ess, which is also its supremum, and prove that μ_* > 1 if and only if the best
constant λ_* of the Antonov inequality lies below δ₁²; in that case λ_* is the eigenvalue. HRSS23 (Proposition 5.3) call
it the Birman–Schwinger–Mathur criterion. In these works the spectral parameter λ is the
square of the frequency. None of them gives a numerical value of M or μ_*.

### 2.2 Other mathematical work

| work | content |
|---|---|
| Rioseco, Sarbach 2020; Moreno, Rioseco, Van Den Bosch 2022 | Phase mixing in a central potential and in an anharmonic well, with rates |
| Chaturvedi, Luk 2021 | 1D transport in a confining potential: decay t⁻² of ∂ₜφ |
| Chaturvedi, Luk 2024 | Vlasov–Poisson with an external Kepler potential: linear phase mixing in 3D and nonlinear phase mixing in spherical symmetry on long time scales. The data are supported on bounded orbits with energies away from zero |
| Weder 2025 | Scattering theory for the plane-symmetric Antonov operator around polytropes and King states |
| Nguyen 2023 | Homogeneous plasmas with compact velocity support: a survival threshold κ₀² = 4π∫u²μ/(Υ² − u²)du; oscillations below it; the decay law at the threshold depends on how μ vanishes at the maximal speed Υ |
| Bedrossian, Masmoudi, Mouhot 2022; Han-Kwan, Nguyen, Rousset 2021 | Linearized Vlasov–Poisson on ℝ³: Klein–Gordon-type oscillations at low wave number |
| Faou, Horsin, Rousset 2021; Faou, Rousset, Seetohul 2026 | Inhomogeneous states of the HMF model: algebraic linear damping; nonlinear dynamics for long times |
| Perez, Aly 1996 | Linear stability of a spherical collisionless system whose particles also feel an external potential: if ∂f₀/∂E < 0 it is stable with respect to all spherical perturbations |
| Pausader, Widmayer 2020 | Small radial perturbations of a point charge, repulsive case: global solutions with modified scattering, by action-angle coordinates |
| Velozo Ruiz, Velozo Ruiz 2023; Bigorgne et al. 2023 | Small data with the external potential −\|x\|²/2, which has unstable trapping: sharp decay and modified scattering |
| Eo 2026 | Einstein–Vlasov near Schwarzschild with matter on bounded geodesics: quantitative linear phase mixing, and nonlinear phase mixing up to the time scale of the echoes |
| Bedrossian 2016 | On the periodic box, the theorem of Mouhot and Villani does not extend in general to Sobolev regularity in the gravitational case, because of echoes |
| Lin, Zeng 2011 | BGK waves exist arbitrarily close to homogeneous states in W^{s,p} with s < 1 + 1/p; under the Penrose condition there are none for larger s |
| Stucchi, Lauber 2025 | Roots of the dispersion relation for Maxwellian, kappa and cut-off distributions (Section 2.11) |
| Kunze, Ortega 2025; Sridhar 1989† | Exact time-periodic solutions: around the Kurth solution, and with uniform density |

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

Criteria and thresholds. Campa and Chavanis (2010) and Ogawa (2013) give necessary and
sufficient conditions of linear stability for inhomogeneous states of the HMF model, as
explicit inequalities. Bachelard et al. (2011) treat models in which the inhomogeneity is
created by an external field: they write the Vlasov equation in action-angle variables,
derive a dispersion relation, obtain the growth rate and the stability threshold, and compare
with N-body simulations. Chavanis (2013) shows that, in the linear regime, the water-bag
distribution responds to a step perturbation with permanent oscillations, as a barotropic gas
does, whereas other distributions are Landau damped. Chavanis (2012) applies the Nyquist
method to the Jeans problem and gives the damping rate and the period for isothermal and
polytropic distributions. Barré, Métivier and Yamaguchi (2016) study non-oscillating
bifurcations of inhomogeneous states: the resonances are strongly suppressed, and the
saturation generalizes the trapping scaling of plasmas.

Coupled oscillators. In the Kuramoto model below the critical coupling, with a frequency
distribution of compact support, the linearized relaxation is exponential at intermediate
times and slower than exponential at long times (Strogatz, Mirollo and Matthews 1992†).
Fernandez, Gérard-Varet and Giacomin (2014) and Dietert (2016) prove the damping of the order
parameter; for sufficiently regular distributions Dietert finds exponential decay in the
stable case and finitely many eigenmodes otherwise.

Trapped particles with long-range forces. Olivetti et al. (2009) compute the breathing mode
of trapped particles with power-law interactions from a dynamical ansatz and support it with
simulations. Chalony et al. (2012) realize a one-dimensional attraction of gravitational type
in a cold strontium gas and measure the density profile, the size of the cloud and its
breathing oscillations.

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

Petersen et al. (2024, Appendix C2) find that the frequency of the l = 1 damped mode of the
isotropic isochrone reported by Fouvry and Prunet, 0.0143 − 0.00142 i in units of the
frequency scale, drifts towards the origin as the number of resonances increases. They
conclude that it seems to be the neutral translation mode, not converged, and find no robust
damped mode with l = 1 to 4 in the isotropic isochrone and Plummer models. That number
should not be quoted as a physical mode.

**N-body studies of radial oscillations.** Hénon (1973), read in the scan, uses 1000 concentric
shells to test polytropes against spherical perturbations and finds them stable almost down
to the index n = 1/2; at n = 1/2 he sees collective oscillations and judges the case
slightly unstable. He states that in his
earlier papers (Hénon 1967, 1968) all the shells had the same angular momentum, so that the
plane of radius and radial velocity described the system completely, and that the stratified
sampling of the initial data, which he had already used in 1968, is similar to the quiet
start of plasma simulations. According to Louis and Gerhard (1988), Hénon (1968) found
virtually undamped oscillations of such a system of shells under certain circumstances. The
papers of 1967 and 1968 were not accessed. Miller and Smith (1994)† report normal-mode oscillations that continue
undamped long after the initial transients, with the kinetic energy oscillating by up to
10 % of its mean; Wachlin and Muzzio (1997)† confirm them with a perturbation-particle
method. David and Theuns (1989)† find long-lived radial pulsations whose lifetime depends
strongly on the number of particles. S24 cites Sweatman (1993) for damped oscillations near
the Plummer sphere and Namboodiri (2000) for polytropes; these two were not accessed.

**The model with fixed |L|.** The model goes back to Hénon (1967, 1968), as described above.
Klinko and Miller (2002) simulate N concentric mass shells with
a fixed magnitude of the angular momentum, the particle version of the fixed-|L| model; the
abstract mentions no external potential. Their subject is a phase transition between
quasi-uniform and core–halo states. They report that the equilibration follows a power law
and that there is strong evidence of long-lived collective oscillations in the supercritical
region. Gargar (2011, thesis) studies systems of few and of many shells, the latter in the
Vlasov limit. Destri (2014) integrates the spherically symmetric Vlasov–Poisson system for
the collapse of warm dark matter halos and attributes their hollow cores to the conservation
of the squared angular momentum.

**Matrix method for radial perturbations.** Polyachenko and Shukhman (2026) build two
families of potential–density pairs for l = 0 in systems of infinite extent, in which the
condition of zero perturbed mass is part of the basis. Their introduction states that
current studies of spherical models focus on the damping of perturbations in stable systems
and on the possible existence of Landau quasi-modes. They cite Weinberg (1994) and
Polyachenko et al. (2021), and no work of the Bayreuth–London line. Vandervoort (2004)
formulates a matrix method in a Lagrangian representation and computes the frequencies of the
lowest radial modes of a family of spherical models.

**Linear response.** Nelson and Tremaine (1997) derive the response operator of a general
inhomogeneous system. The damping of a collective mode is given by the anti-Hermitian part of
a polarization operator, and without the self-gravity of the response the expressions reduce
to the formulae of perturbation theory in action-angle variables. Weinberg (1997) computes
the amplification of the Poisson noise by self-gravity: a factor of about six in the dipole
power for a King model and fifteen for a Hernquist model. Weinberg (2022) finds that small
wiggles in the distribution function can destabilize the weakly damped dipole modes. Ng and
Bhattacharjee (2021) show that with a weak collision operator the Landau modes become true
eigenmodes. Jalali and Tremaine (2011) compute the normal modes of near-Keplerian discs with
finite elements: slow modes exist for an arbitrarily small disc mass as long as the
self-gravity of the disc dominates the apsidal precession.

**Finite amplitude and resonant trapping.** Louis and Gerhard (1988), read in the text,
construct a self-consistent distribution function for a radial oscillation of the isochrone
sphere with an amplitude of 10 % in the central density. Many orbits are trapped by the
resonant families, in which the frequency of the perturbation is a multiple of the radial
frequency. The stable closed resonant orbits oscillate against the perturbation and damp it,
and the non-resonant orbits support it. The self-consistent solution therefore has holes in
the energy distribution around the resonances, which a distribution with df/dE < 0 does not
have, and they conclude that such oscillations are difficult to set up. Vandervoort (2003) builds stationary oscillations of small finite amplitude by canonical
perturbation theory. Chiba and Schönrich (2022) describe the orbits trapped at the
resonances of a slowly evolving bar: the libration frequency falls towards the separatrix,
the trapped density winds up into a spiral, and the torque oscillates and damps. Hamilton et
al. (2022) add diffusion to the kinetic equation near a resonance; without it they recover
the result of Tremaine and Weinberg that the friction vanishes by phase mixing of the
resonant orbits.

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

Mathur (1990), read in the text. The system has a fixed external field, which he takes
harmonic (a halo much larger than the system), and finite extent. For the radial problem the
system is a set of shells with angular momentum L, and L has a positive minimum. A mode with
real frequency ω in a gap exists if and only if a symmetric Hilbert–Schmidt kernel K_ω has the
eigenvalue one. He argues that, among the K_ω of a gap, the largest eigenvalue is attained at
the right end of the gap. The example is f₀ = μ(1 − 2E) for 0 ≤ E ≤ 1/2 with μ infinitesimal,
a distribution that vanishes linearly at the edge (k = 1 in our notation). The quadratic form
of the kernel at the band edge diverges logarithmically (his Eq. A13), and a mode follows by
continuity for small μ. For the radial case he takes a narrow Gaussian in L and nearly
circular orbits, which reduces the kernel to the one-dimensional one; our fixed-|L| model is
the same limit without the restriction to nearly circular orbits. He gives no numerical
evaluation. He expects the modes to disappear when the distribution has a long tail, and
notes that a larger spread in L shortens the gaps.

Mathur's paper is not consistent about where the mode lies. Result (i) of its Section 4.1 and
the discussion state that no mode can lie below the lowest orbital frequency, and that the
modes lie in the other gaps. The example, however, is placed "in the gap between S₁ and S₋₁",
which is the principal gap, at its right end 2π/T(E_max); Appendix A announces the gap between
S₁ and S₂; and Eq. (4.1) has the opposite sign of the sum of the terms n and −n of Eq. (3.16).
HRS22 prove that the operator is nonnegative in the principal gap and place the eigenvalue
there. In our reading the principal gap is the correct one for the example: in a harmonic
external field the rigid oscillation of the whole system is an exact mode at the frequency of
the field, which Mathur mentions in a footnote, and that frequency lies below all the orbital
frequencies, because the self-gravity of the system only adds to the restoring force.

Weinberg (1991), read in the text. He writes the dispersion relation with the matrix method
of Kalnajs and locates the modes with real frequency as its zeros in the gaps. The models are
truncated isothermal slabs of depth W₀, a variant with f′(E_max) = 0, and polytropes, with a
harmonic and an anharmonic halo. His Table 1 gives the ends of the gaps and the modes, and
marks the gaps that contain no mode. Except in the principal gap, the mode lies near the
right end of its gap. A sharp edge favours a mode, because the integrand depends on df/dE; a
distribution that goes smoothly to zero at the edge "will have more difficulty sustaining a
self-consistent response", which he finds in the variant with f′(E_max) = 0 and in the
polytropes of large index. A halo widens the gaps. In N-body runs with 64 000 particles a mode
imposed on the initial data keeps its frequency for about 100 dynamical times. The amplitude
is limited: the symmetric breathing mode reaches about 3 % of the background force, and the
next mode, which is closer to the continuum (8.146 against 8.16, where the first has 5.38
against 5.44), reaches 0.5 %. He attributes the limit to the nonlinear broadening of the
mode, whose wings overlap the continuum and damp, and observes a slow beat at the difference
between the frequency of the mode and the edge of the continuum. He adds that a small measure
of orbits at higher energy, which closes a gap, only turns the undamped mode into a weakly
damped one, and that a system without gaps may still have ranges of frequency in which the
damping is slow.

Louis (1992), read in the text. He closes the moment equations, a fluid model in two
approximations, which reduces the eigenvalue problem to an ordinary differential equation, and
compares with N-body simulations. The external potential is bx². His Figure 1 is a mode
diagram: the frequencies of the modes against the polytropic index n and against the depth W₀
of King models, together with the smallest and largest orbital frequencies and their
harmonics. The modes lie close to the lower edge of a continuum segment of the same parity;
in a symmetric slab a mode couples only to the harmonics of its parity. With increasing
central concentration the modes cross into the continuum, the higher ones first. The
fundamental breathing mode disappears at n ≈ 2 for polytropes and at W₀ ≈ 2 for King models,
partly in disagreement with the table of Weinberg. An error of 1 % in the frequency is
already significant, because the modes are so close to the continuum. With an external
potential the bands become narrower and the modes move closer to the lower edge, while the
critical values of n and W₀ "remain almost unchanged". Inside the continuum the modes
continue as slightly damped density waves, whose damping rate is very small close to the
lower edge, because the phase-space density of the least bound particles is small. At finite
amplitude the trapping of resonant particles leads to oscillations of the wave amplitude,
and "it is not clear whether the amplitude will necessarily decay to zero". He calls the
examples of Mathur "somewhat artificial" and announces an extension to spherical geometry.

Widrow and Bonner (2015), read in the text, state the role of the gap. In the untruncated
Spitzer sheet the orbital frequencies range from Ω_c down to zero, there are no gaps, and
"all coherent oscillations of the Spitzer sheet are damped"; they compute these oscillations
by analytic continuation of the response matrix. In the lowered sheet the frequencies have a
lower bound Ω_min, and for W = 2 there are no resonant particles for 0 < ω ≲ 0.72 Ω_c, the
principal gap. The modes found with the method of Mathur and Weinberg, who restrict the
search to real frequencies in the gaps, lie slightly below nΩ_min.

One difference with our model follows from these papers. In the slabs of Weinberg and Louis,
which are symmetric about the midplane, the principal gap contains only the translation or
seiche mode, and the breathing mode lies below the segment of the second harmonic, because a
mode couples only to the harmonics of its parity. The radial problem has no such symmetry, and the mode lies in
the principal gap, below Ω_min, as in HRS22.

Banik, Weinberg and van den Bosch (2022, 2023) give a perturbative formalism for the phase
mixing of the response of a slab and of a disc, without the self-gravity of the response.
Chiba et al. (2025) derive the galactic analogue of the plasma echo in angle-action variables
for the vertical motion. Sellwood (1996) finds that razor-thin discs have a discrete spectrum
of neutral bending modes, which decay by Landau damping when the disc has random motion
normal to its plane. Miller et al. (2023)† review the simulations of one-dimensional gravity;
Colombi and Touma (2014) and Joyce and Worrakitpoonpon (2010) study the relaxation to
quasi-stationary states in one dimension.

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

Sellwood (1983)† develops quiet starts for disc simulations: the mode frequencies agree with
linear theory within 2 % for cold discs and about 10 % for warm discs, and starts that are
not quiet give errors above 50 % with the same number of particles. Wachlin, Rybicki and
Muzzio (1993) describe another perturbation-particle method (known here from its title). Among the direct solvers,
Yoshikawa, Yoshida and Umemura (2012) integrate the Vlasov–Poisson system in six dimensions
and test the stability of King spheres and the Landau damping against linear theory; Halle,
Colombi and Peirani (2019) compare a spherical Vlasov solver, a shell code and an N-body
code; Colombi and Touma (2014) use the water-bag method in one dimension. Macridin et al.
(2015)† analyze particle simulations of bunched beams with dynamic mode decomposition and
obtain the shapes, the frequency shifts and the damping rates of the modes. Boine-Frankenheim
and Egenolf (2024) use the conservation of the energy, the entropy and the phase-space area
to assess particle tracking with space charge, and stress the choice of the cut-off harmonic
of a grid-based solver for the long-term accuracy.

### 2.7 Plasmas

The standard numerical test is the linear Landau damping of a small sinusoidal perturbation
of a Maxwellian. The frequency and the damping rate are fitted to the field energy and
compared with the root of the dispersion relation. For k = 0.5 the values quoted in the
literature are ω = 1.4156 and γ = −0.1533; they were not checked against a specific paper.
Finn et al. (2023) note that approximations to the dispersion relation are not accurate
enough for this comparison when the parameters are varied.

Finite amplitude. Manfredi (1997)† finds numerically that, after some initial damping, the
field oscillates around an approximately constant value and a vortex survives in phase
space; Isichenko (1997) had predicted an algebraic decay, E ∝ 1/t. Lancellotti and Dorning
(1998)† show that critical initial states separate the initial conditions that damp to zero
from those that evolve to a nonzero state. Brunetti, Califano and Pegoraro (2000)† determine
the parameter that marks the transition between the regimes of Landau and of O'Neil; the
long-time state is a superposition of two BGK waves. Ivanov, Cairns and Robinson (2004)† find
a critical initial amplitude ε*: below it the damping is the linear one, at ε* the field
decays as t^{−3.26}, and above it the field regrows and oscillates. For a Maxwellian, ε*
corresponds to a bounce frequency equal to the Landau rate. Klimas, Viñas and Araneda
(2017)† compute critical exponents of this transition. Danielson, Anderegg and Driscoll
(2004)† measure both regimes in a pure electron plasma: at low amplitude the rate agrees with
linear theory, and at larger amplitude the wave damps, regrows and approaches a steady state,
as predicted by O'Neil (1965). Ouyang, Zhu and Ng (2025) identify the asymptotic state as a
BGK structure with several waves.

Confined plasmas. Anderegg et al. (2016)† measure Landau damping at the harmonics of the
bounce frequency in a trapped plasma whose potential varies along the axis.

### 2.8 Accelerator physics: loss of Landau damping

After the one-dimensional slab of Section 2.4, the longitudinal dynamics of a bunch in a
synchrotron is the closest analogue of our reduced model that was found. The particles move in the potential well of the radio-frequency
system. The well is nonlinear, so that the synchrotron frequencies fill a band. The
self-force comes from the impedance of the machine, for instance space charge or the
inductance of the chamber. The stationary distribution has compact support, and the family in
common use is binomial, proportional to (1 − E/E_max)^μ in the energy E of the synchrotron
oscillation. The linearized Vlasov equation is written in action-angle variables.

In this field "Landau damping is considered to be lost when the frequency of the coherent
bunch oscillations moves outside the incoherent frequency band" (Karpov, Argyropoulos and
Shaposhnikova 2021). According to that paper, Lebedev (1968) wrote the first self-consistent
matrix equation, Sacherer (1973) and Hofmann and Pedersen (1979) gave approximate criteria,
Chin, Satoh and Yokoya (1983) introduced the van Kampen modes, with the threshold reached
when a discrete mode emerges from the band, Oide and Yokoya (1990) gave a numerical
eigenvalue method, and Burov (2010) combined the two to compute thresholds.

The results that correspond to ours are the following.

- Criterion. The threshold is reached when the largest eigenfrequency equals the edge of the
  band. Karpov et al. obtain an analytic threshold from the first-order expansion
  det(1 + εX) ≈ 1 + ε tr X, that is, with the trace of the operator in the place of its
  largest eigenvalue.
- Dependence on the edge. For the binomial family and a constant inductive impedance the
  threshold is ζ_th = πφ_max⁵/[32 μ(μ + 1) χ], with φ_max the largest phase amplitude in the bunch and χ a function of
  μ and of the cut-off of the impedance. It vanishes when the impedance has no cut-off at high
  frequency, and for the water-bag distribution (μ = 0) there is no threshold. Burov (2012)
  finds the thresholds "extremely sensitive to the small-argument behaviour of the bunch
  distribution function", that is, to its form at the end of the band where the mode emerges.
  Burov (2021) shows that for a repulsive inductive impedance the discrete spectrum consists
  of infinitely many modes with real frequencies that accumulate at the edge of the band.
- Compact support. In Burov (2012), discrete modes of this type "may only appear if the
  distribution function is of a finite width", and "even a tiny tail covering the coherent
  frequency yields Landau damping, killing that discrete mode".
- Reported quantities. The threshold intensity; below it, the damping time of the response
  to a rigid kick, which diverges at the threshold; above it, the amplitude of the residual
  oscillations. Karpov et al. compute them by expanding the kick in van Kampen modes, and
  compare with particle simulations with tens of millions of macro-particles and with
  measurements in the LHC.
- Simulations. Boine-Frankenheim and Shukla (2005) confirm the simple criteria by particle
  tracking. Boine-Frankenheim and Egenolf (2024) obtain the spectrum of the bunch
  oscillations and locate the point where the frequency of the dipole mode emerges from the
  incoherent spectrum. Intelisano, Damerau and Karpov (2025) extend the thresholds to
  double-harmonic systems and compare with measurements in two synchrotrons.

Two differences, in our reading. With a constant inductive impedance the force is
proportional to the derivative of the line density, and the sum over harmonics diverges
without a cut-off; this is why the threshold vanishes. The gravitational force of our model is
proportional to the enclosed mass, and the threshold is finite for k > 1. For the repulsive
case treated by Karpov et al. the mode emerges above the band, at the frequency of the centre
of the bunch; in our attractive case it emerges below the band, at the frequency of the outer
edge.

### 2.9 Two-dimensional vortices

The linearized Euler equation around a circular vortex has the same structure, with the
rotation frequency of the fluid in the place of the orbital frequency. Schecter et al.
(2000)† study the damping of an elliptical perturbation in theory and in an experiment with
an electron plasma. An impulse excites a quasi-mode that decays exponentially at a rate
proportional to the gradient of the vorticity at the critical radius, where the fluid rotates
in resonance with the wave; the quasi-mode is not an eigenmode. Balmforth, Llewellyn Smith
and Young (2001)† find a critical forcing amplitude: below it the quasi-mode decays, and
above it nonlinear effects stop the decay and cat's eyes form. Bedrossian, Coti Zelati and
Vicol (2017) prove the inviscid damping of the velocity field around a monotone vortex, with
the rates t⁻¹ and t⁻². Isichenko (1997) predicts t^{−5/2} for the stream function of a
perturbed shear flow.

### 2.10 What the studies report

Numerical studies report measured quantities. RR18 and S24 plot the kinetic and potential
energies against time and classify the solutions by inspection. S24 measures the fundamental
period p and compares (2π/p)² with the essential spectrum, and does not fit damping rates.
Plasma studies fit ω and γ and compare them with the dispersion relation; at finite amplitude
they give a critical amplitude. Studies of clusters give the complex frequency of a mode and
compare its period with the power spectrum of the density centre of N-body runs (Weinberg
1994; Heggie, Breen and Varri 2020). Barré et al. (2011) give the exponent and the frequency
of the algebraic tail. For the slab, Weinberg (1991) tabulates the ends of the gaps and the
frequencies of the modes, and Louis (1992) draws the frequencies of the modes against a
parameter of the model together with the band edges. In accelerator physics the reported
quantities are the threshold intensity, the damping time below it and the residual amplitude
above it.

Analytical studies define a quantity and compare it with one: the number M of HRS22 and μ_*
of Kunze (Section 2.1); the norm that bounds the number of modes in Moreno et al. (2023); the
determinant det[I − M(ω)] of the response matrix in stellar dynamics, where Fouvry and Prunet
(2022) also plot the largest eigenvalue of [I − M(ω)]⁻¹ against the real frequency; the
survival threshold of Nguyen (2023); and the trace criterion of Karpov et al. (2021).

### 2.11 Supports that are not compact

The criterion needs a gap below the band. In the isochrone and Kepler potentials the orbital
frequency tends to zero as the action grows, so that a support that reaches the escape energy
leaves no gap.

- All the mathematical works on the gravitational case assume compact support: HRS22,
  HRSS23, HRSS24, Kunze, Moreno et al., and Chaturvedi and Luk (2024). S24 names the steady
  states with unbounded support as an extension to be done.
- In the untruncated Spitzer sheet all coherent oscillations are damped (Widrow and Bonner
  2015, Section 2.4). Mathur (1990) expects the modes to disappear when a long tail is added
  to the distribution. Weinberg (1991) argues that a small measure of orbits added at high
  energy turns an undamped mode into a weakly damped one.
- Nguyen (2023): for equilibria that are positive for all velocities the survival threshold
  is zero and there are no pure oscillations. At small wave number the Landau rate is of order
  |k|^{2N₀−3} for equilibria that decay as the power −N₀ of the energy, and of order
  |k|⁻³exp(−c/|k|²) for Gaussians: the faster the equilibrium decays, the weaker the damping.
- Stucchi and Lauber (2025): with a cut-off at the velocity v_c, a wave with phase velocity
  above v_c has no resonant particles and is not damped. Replacing the step by a smooth
  function changes the set of roots, and differently for different smooth functions.
- Burov (2012): a tail of the distribution that covers the coherent frequency restores the
  damping (Section 2.8).
- Barré and Yamaguchi (2013): for the self-consistent isochrone model without the
  self-gravity of the perturbation, the weakly bound orbits give a "singularity at infinity"
  and the slowest tail, t^{−2/3} with zero frequency. In two action dimensions, a frequency
  that decays as |J|^{−a} and an integrand that decays as |J|^{−b} give t^{−(b−2)/a}.
- Petersen et al. (2024) find no robust damped modes with l = 1 to 4 in the isotropic
  isochrone and Plummer models (Section 2.4). For radial oscillations near the Plummer sphere
  only the N-body result cited by S24 was found.
- Our own cases are the two Gaussian equilibria of the η battery, F ∝ J²exp(−J²/σ_J²) cut at
  6σ_J. The reference run (σ_J = 0.10, a₀ = 0.01) is damped, with a frequency inside the
  band; without self-gravity its decay is Gaussian in time over four orders of magnitude. Case
  E2 (σ_J = 0.05, a₀ = 0.069) has the frequency 0.06695, below the 1 % band
  [0.06842, 0.07871], and is still damped, with γ = 2.0 × 10⁻⁴. A formula for this rate in
  terms of F′ at the resonant action, the analogue of Landau's formula and of the one of
  Nelson and Tremaine (1997), was not derived.

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
7. **Decay rates with self-gravity.** S24 (Outlook) writes that for the full linearized
   system "no quantitative damping results are available", and that the numerics should
   determine the optimal decay rates and their dependence on the initial data and on the
   steady state. HS25, which appeared later, gives the first analytical rates.
8. **Unbounded supports and external potentials** are named by S24 as extensions of the
   numerical methods.
9. **Number of modes and value of the supremum.** Moreno et al. (2023): when the integral
   that bounds the number of modes diverges, eigenvalues occur for arbitrarily weak coupling,
   and "whether they are finite in number or not, remains an open question". Kunze's notes
   ask to "determine the limit μ_*" in terms of the steady state, and whether λ_* = δ₁² can
   happen.

## 4. Where our results stand

| our result | closest prior work | assessment |
|---|---|---|
| R1. λ_edge, the norm of the Mathur operator at the edge of the gap (the number M of HRS22 and μ_* of Kunze), computed numerically for 61 equilibria (isochrone and point mass, fixed \|L\|) and checked against the linear solver and PIC runs | principle: Mathur 1990, HRS22, Kunze 2021. Numerics for the gravitational radial problem: none found. For the one-dimensional slab, Weinberg (1991) locates the modes numerically as zeros of the dispersion relation in the gaps and tabulates them for three families, with and without a halo, and Louis (1992) gives a mode diagram with the values of the polytropic index and of W₀ at which the modes enter the continuum. In accelerator physics the same criterion is evaluated routinely (Section 2.8) | new only for the radial problem with a central potential, as far as this search goes; it is open problem 1. In one dimension the numerical evaluation exists since 1991. The fixed-\|L\| model limits its reach. The slab papers and the beam literature must be cited: criterion, numerical method and reported quantities coincide |
| R2. Threshold masses of the HRSS23 setting: a₀ ≈ 0.514 (k = 1.25), 0.920 (k = 1.5), 1.53 (k = 2), 2.38 (k = 3). These are the values of 7 October; they replace 0.53, 0.93 and "above 1 for k ≥ 2" | HRSS23 (ε₀ not estimated) | apparently new; of interest mainly to that group |
| R3. Edge law λ_edge − λ(Ω_min − δ) ∝ δ^{min(k−1,1)} and the binding laws, with c the coupling (proportional to the mass): δ_d ∝ c^{1/(1−k)} for k < 1, exp(−const/c) for k = 1, (λ_edge − 1)^{1/(k−1)} for 1 < k < 2 | the divergence for k ≤ 1 is in the proof of HRSS23; the laws are those of Simon and Klaus–Simon with d = 2k; Nguyen 2023 in homogeneous plasmas. Mathur (1990): logarithmic divergence at the edge for a distribution that vanishes linearly, the case k = 1. Weinberg (1991): a sharp edge favours a mode and a smooth one hinders it. In beam physics the threshold depends on the exponent μ of the binomial distribution (Karpov et al. 2021) and on the form of the distribution where the mode emerges (Burov 2012) | the analogy and the laws do not appear in HRS22, HRSS23, Kunze's notes, Straub's thesis or Moreno et al. (no occurrence of "weak coupling", "coupling constant" or Klaus). The qualitative role of the edge is known since 1990 for the slab, and in beam physics. The laws themselves, with the exponent min(k − 1, 1) and the binding laws, were not found; the mathematics is classical |
| R4. Observability: the mode separates from the edge after t ≈ 1/(Ω_min − ω_d); at small mass the modes of k ≤ 1 lie less than 10⁻⁶ band widths below the edge | Louis 1992: the modes are so close to the continuum that an error of 1 % in the frequency prevents locating where they cross. Karpov et al. 2021: the damping time diverges at the threshold, so that the damping "can be effectively lost even below the threshold" | apparently new for the gravitational thresholds, and it explains open problem 2: near the threshold the binding is a high power of the distance to it, so no time integration can decide |
| R5. Algebraic tails: t^{−(k+1)} for the free mixing of a perturbation proportional to F_eq, and t^{−k} with self-gravity at late times. Corrected on 7 October: the earlier statement, t^{−(k+1)} with self-gravity, holds only as a transient at small mass | HRSS24, Theorem 1.3: t^{−min(2,k)} for pure transport with data \|φ′(E)\| g₀, which vanish at the edge like F_eq′. Barré et al. 2011 and Barré, Yamaguchi 2013: the tail is set by the strongest branch singularity of the response function or of the initial datum. HS25 | the exponent −k is not new as a rate. What appears to be new is that self-gravity turns the tail of a perturbation proportional to F_eq from t^{−(k+1)} into t^{−k}, with the numerical evidence. It bears on open problems 6 and 7. The derivation is formal |
| R6. PIC results: quiet start in angle–action variables with the D − Z subtraction; grid and particle noise separated; finite-amplitude loss of the modes next to the edge, described by a pendulum at the resonance | quiet starts (Hénon 1968, as described in Hénon 1973; Sellwood 1983), weighted particles (Barré et al. 2011), perturbation particles (Leeuwin et al. 1993), δf methods. Weinberg 1991: the amplitude of a mode close to the continuum is limited by its nonlinear broadening (3 % and 0.5 % in two cases), with a beat at the difference frequency. Louis 1992: trapping leads to oscillations of the amplitude. Louis and Gerhard 1988: the resonant orbits oscillate against the perturbation and damp it. Also O'Neil 1965; Ivanov et al. 2004 (critical amplitude); Balmforth et al. 2001 (vortices); Chiba and Schönrich 2022 (pendulum, oscillating torque); Karpov et al. 2021 (residual amplitude) | known ideas; supporting material. That the closeness to the continuum limits the amplitude is in Weinberg (1991), as a conjecture supported by two runs. A formula for the threshold amplitude was not found in this search |
| R7. η is not universal | Darling, Widrow 2019 use a similar ratio (external potential to self-gravity, about 4:1) | a negative result; useful as a remark: such a ratio cannot be the criterion |
| R8. The reduced model: shells with fixed \|L\| in an external central potential | Hénon 1967, 1968: shells that all have the same angular momentum, with virtually undamped oscillations under certain circumstances (through Hénon 1973 and Louis and Gerhard 1988). Mathur 1990: shells with a narrow distribution of L in an external field, as the limit that reduces the radial problem to one dimension. Klinko and Miller 2002 simulate concentric shells with fixed \|L\| and report long-lived collective oscillations; Hénon 1973 (shell model); HRSS23 and HRSS24 (analysis) | the model is that of Hénon (1967, 1968), and it is the limit used by Mathur (1990). No numerical linear analysis or threshold was found for it |

Dynamic mode decomposition and a neural-network classifier have both been applied to this
problem (Darling and Widrow 2019; Straub and Wolfschmidt 2024), so neither would be new as a
method. Dynamic mode decomposition is also used on beam simulations (Macridin et al. 2015).

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
Such a paper should cite the beam literature of Section 2.8, whose way of presenting the
results (threshold, damping time below it, residual amplitude above it) is a model.

**B. Persistent modes of the vertical slab.** The slab is a one-action problem, so the
present code applies after changing the Green function and the background. The question,
when a self-gravitating slab in an external potential keeps an undamped mode instead of
winding into a phase spiral, is of current interest (Gaia). Weinberg (1991) and Louis (1992)
located the modes and the models at which they enter the continuum, and Widrow and Bonner
(2015) the damped oscillations; the laws R3 were not found. This option needs its own review of the phase-spiral literature (Banik,
Weinberg and van den Bosch; Widrow 2023) before any work.

**C. A note on decay exponents** (R5) is possible but small.

Option A answers a question that the authors of the field have written down, and most of
its validation is done. A feasibility test on isotropic polytropes decides it.

## 6. Limits of this review

- Mathur (1990), Weinberg (1991) and Louis (1992) were read in the revision. Mathur's text is
  not consistent about the gap that contains the mode (Section 2.4); what is said here about
  it is our reading and should be checked by a second reader. None of the three evaluates the
  criterion for spherical systems; Louis announces it as future work, and no such work by him
  was found.
- The papers of Hénon of 1967 and 1968, in which the shells with the same angular momentum
  appear, were not accessed; what is said about them comes from Hénon (1973) and from Louis
  and Gerhard (1988).
- Kunze's book was not read; the lecture notes, written a year later, still ask for
  numerics.
- The citing lists come from two indexes that disagree (OpenAlex finds 6 citing works for
  HRSS23, Semantic Scholar 12). Works from 2026 may be missing.
- The astrophysical literature on radial modes before 2000 was covered through the reference
  lists of HRS22 and S24. In the revision, the N-body studies of the 1980s and 1990s were
  added from summaries of their abstracts.
- The works marked † were seen only through the summary returned by a search engine.
- The accelerator literature was covered only for the longitudinal loss of Landau damping,
  through four papers and their reference lists. The transverse problem with space charge and
  the microwave instability were not searched.
- The identification of our λ(ω) with the norm M_λ of HRS22 rests on the fact that the
  operators AB and BA have the same nonzero spectrum. It has not been written out.
- Cosmological simulations and Vlasov–Maxwell systems were not searched.

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
- Mathur, MNRAS 243, 529 (1990), [record](https://www.osti.gov/etdeweb/biblio/6932499), [scan](https://articles.adsabs.harvard.edu/pdf/1990MNRAS.243..529M)

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

- Weinberg, ApJ 373, 391 (1991), [ADS](https://ui.adsabs.harvard.edu/abs/1991ApJ...373..391W/abstract), [scan](https://articles.adsabs.harvard.edu/pdf/1991ApJ...373..391W); ApJ 421, 481 (1994), [arXiv:astro-ph/9306020](https://arxiv.org/abs/astro-ph/9306020)
- Louis, Gerhard, MNRAS 233, 337 (1988), [journal](https://academic.oup.com/mnras/article/233/2/337/969454), [scan](https://articles.adsabs.harvard.edu/pdf/1988MNRAS.233..337L)
- Louis, MNRAS 258, 552 (1992), [journal](https://academic.oup.com/mnras/article/258/3/552/1081867), [scan](https://articles.adsabs.harvard.edu/pdf/1992MNRAS.258..552L)
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

Added in the revision of 8 October 2026. A dagger marks the works seen only through a
search-engine summary of the abstract.

Mathematics:

- Perez, Aly, MNRAS 280, 689 (1996), [arXiv:astro-ph/9511103](https://arxiv.org/abs/astro-ph/9511103); Perez, Alimi, Aly, Scholl, MNRAS 280, 700 (1996), [arXiv:astro-ph/9511090](https://arxiv.org/abs/astro-ph/9511090)
- Pausader, Widmayer, [arXiv:2008.08013](https://arxiv.org/abs/2008.08013)
- Velozo Ruiz, Velozo Ruiz, [arXiv:2304.12017](https://arxiv.org/abs/2304.12017); Bigorgne, Velozo Ruiz, Velozo Ruiz, [arXiv:2310.17424](https://arxiv.org/abs/2310.17424)
- Eo, [arXiv:2609.13394](https://arxiv.org/abs/2609.13394)
- Bedrossian, [arXiv:1605.06841](https://arxiv.org/abs/1605.06841)
- Lin, Zeng, Commun. Math. Phys. (2011), [arXiv:1003.3005](https://arxiv.org/abs/1003.3005)
- Stucchi, Lauber, J. Plasma Phys. 91, E44 (2025), [arXiv:2411.06769](https://arxiv.org/abs/2411.06769)
- † Sridhar, MNRAS 238, 1159 (1989), [journal](https://academic.oup.com/mnras/article/238/4/1159/1037490)

Statistical physics, oscillators and trapped gases:

- Campa, Chavanis, J. Stat. Mech. P06001 (2010), [arXiv:1003.2378](https://arxiv.org/abs/1003.2378)
- Ogawa, [arXiv:1301.1130](https://arxiv.org/abs/1301.1130)
- Bachelard, Staniscia, Dauxois, De Ninno, Ruffo, J. Stat. Mech. P03022 (2011), [arXiv:1010.4647](https://arxiv.org/abs/1010.4647)
- Chavanis, Eur. Phys. J. B 85, 229 (2012), [arXiv:1002.0291](https://arxiv.org/abs/1002.0291); Eur. Phys. J. Plus 128, 38 (2013), [arXiv:1209.5987](https://arxiv.org/abs/1209.5987)
- Barré, Métivier, Yamaguchi, Phys. Rev. E 93, 042207 (2016), [arXiv:1511.07645](https://arxiv.org/abs/1511.07645)
- † Strogatz, Mirollo, Matthews, Phys. Rev. Lett. 68, 2730 (1992), [journal](https://link.aps.org/doi/10.1103/PhysRevLett.68.2730)
- Fernandez, Gérard-Varet, Giacomin, [arXiv:1410.6006](https://arxiv.org/abs/1410.6006); Dietert, J. Math. Pures Appl. 105, 451 (2016), [arXiv:1411.3752](https://arxiv.org/abs/1411.3752)
- Olivetti, Barré, Marcos, Bouchet, Kaiser, [arXiv:0907.4423](https://arxiv.org/abs/0907.4423)
- Chalony, Barré, Marcos, Olivetti, Wilkowski, [arXiv:1202.1258](https://arxiv.org/abs/1202.1258)

Astrophysics:

- Hénon, Astron. Astrophys. 24, 229 (1973), [record](https://www.osti.gov/biblio/4381418), [scan](https://articles.adsabs.harvard.edu/pdf/1973A%26A....24..229H)
- Hénon, Mém. Soc. R. Sci. Liège (5) 15, 243 (1967); Bull. Astron. Paris (3) 3, 241 (1968): taken from the reference list of Hénon (1973), not accessed
- † Miller, Smith, Celest. Mech. Dyn. Astron. 59, 161 (1994), [journal](https://link.springer.com/article/10.1007/BF00692131)
- † Wachlin, Muzzio, Celest. Mech. Dyn. Astron. 67, 225 (1997), [record](http://sedici.unlp.edu.ar/handle/10915/137950?show=full)
- † David, Theuns, MNRAS (1989), [record](https://www.osti.gov/etdeweb/biblio/7242590)
- Sweatman, MNRAS 261, 497 (1993); Namboodiri, Celest. Mech. Dyn. Astron. 76, 69 (2000); Vandervoort, ApJ 273, 511 (1983); Wachlin, Rybicki, Muzzio, MNRAS 262, 1007 (1993): taken from the reference list of S24, not accessed
- Klinko, Miller, [arXiv:astro-ph/0201488](https://arxiv.org/abs/astro-ph/0201488)
- Gargar, thesis, [arXiv:1101.4877](https://arxiv.org/abs/1101.4877)
- Destri, Phys. Rev. D 90, 123531 (2014), [arXiv:1409.6244](https://arxiv.org/abs/1409.6244)
- Polyachenko, Shukhman, [arXiv:2609.04012](https://arxiv.org/abs/2609.04012)
- Nelson, Tremaine, [arXiv:astro-ph/9707161](https://arxiv.org/abs/astro-ph/9707161)
- Weinberg, [arXiv:astro-ph/9707206](https://arxiv.org/abs/astro-ph/9707206); [arXiv:2209.06846](https://arxiv.org/abs/2209.06846)
- Ng, Bhattacharjee, ApJ 923, 271 (2021), [arXiv:2109.07806](https://arxiv.org/abs/2109.07806)
- Jalali, Tremaine, [arXiv:1110.4551](https://arxiv.org/abs/1110.4551)
- Sellwood, [arXiv:astro-ph/9604123](https://arxiv.org/abs/astro-ph/9604123)
- Banik, Weinberg, van den Bosch, [arXiv:2208.05038](https://arxiv.org/abs/2208.05038); Banik, van den Bosch, Weinberg, [arXiv:2303.00034](https://arxiv.org/abs/2303.00034)
- Chiba, Schönrich, MNRAS 513, 768 (2022), [arXiv:2109.10910](https://arxiv.org/abs/2109.10910)
- Hamilton, Tolman, Arzamasskiy, Duarte, [arXiv:2208.03855](https://arxiv.org/abs/2208.03855)
- Chiba, Ding, Hamilton, Kunz, Tremaine, MNRAS (2025), [arXiv:2506.16512](https://arxiv.org/abs/2506.16512)
- † Miller, Manfredi, Pirjol, Rouet, Class. Quantum Grav. 40, 073001 (2023), [journal](https://iopscience.iop.org/article/10.1088/1361-6382/acb8fb)
- Colombi, Touma, [arXiv:1404.5175](https://arxiv.org/abs/1404.5175); Joyce, Worrakitpoonpon, [arXiv:1012.5042](https://arxiv.org/abs/1012.5042)
- García-Perciante, Guerrero, Núñez, Sarbach, [arXiv:2610.06913](https://arxiv.org/abs/2610.06913): phase space mixing in an integrable axisymmetric potential, with N-particle simulations and without self-gravity

Plasmas:

- † Manfredi, Phys. Rev. Lett. 79, 2815 (1997), [journal](https://link.aps.org/doi/10.1103/PhysRevLett.79.2815)
- Isichenko, Phys. Rev. Lett. 78, 2369 (1997), [arXiv:chao-dyn/9612021](https://arxiv.org/abs/chao-dyn/9612021)
- † Lancellotti, Dorning, Phys. Rev. Lett. 81, 5137 (1998), [journal](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.81.5137)
- † Brunetti, Califano, Pegoraro, Phys. Rev. E 62, 4109 (2000), [journal](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.62.4109)
- † Ivanov, Cairns, Robinson, Phys. Plasmas 11, 4649 (2004), [journal](https://pubs.aip.org/aip/pop/article-abstract/11/10/4649/261102/Wave-damping-as-a-critical-phenomenon)
- † Danielson, Anderegg, Driscoll, Phys. Rev. Lett. 92, 245003 (2004), [doi](https://doi.org/10.1103/PhysRevLett.92.245003)
- † Klimas, Viñas, Araneda, J. Plasma Phys. 83 (2017), [journal](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/simulation-study-of-landau-damping-near-the-persisting-to-arrested-transition/5229878E27140B3B9339FE9CE97EBA05)
- Ouyang, Zhu, Ng, [arXiv:2512.17269](https://arxiv.org/abs/2512.17269)
- † Anderegg, Affolter, Kabantsev, Dubin, Ashourvan, Driscoll, Phys. Plasmas 23, 055706 (2016), [journal](https://pubs.aip.org/aip/pop/article/23/5/055706/964519/Bounce-harmonic-Landau-damping-of-plasma-waves)
- Finn, Knepley, Pusztay, Adams, [arXiv:2303.12620](https://arxiv.org/abs/2303.12620)
- O'Neil, Phys. Fluids 8, 2255 (1965)

Accelerator physics:

- Karpov, Argyropoulos, Shaposhnikova, Phys. Rev. Accel. Beams 24, 011002 (2021), [arXiv:2011.07985](https://arxiv.org/abs/2011.07985)
- Burov, [arXiv:1207.5826](https://arxiv.org/abs/1207.5826), [arXiv:1208.4338](https://arxiv.org/abs/1208.4338), [arXiv:1209.2715](https://arxiv.org/abs/1209.2715); Phys. Rev. Accel. Beams 24, 064401 (2021), [arXiv:2103.07523](https://arxiv.org/abs/2103.07523)
- Boine-Frankenheim, Egenolf, [arXiv:2410.13698](https://arxiv.org/abs/2410.13698)
- Intelisano, Damerau, Karpov, [arXiv:2502.14548](https://arxiv.org/abs/2502.14548)
- Karpov, [arXiv:2210.00080](https://arxiv.org/abs/2210.00080); Karpov, Shaposhnikova, [arXiv:2309.06638](https://arxiv.org/abs/2309.06638)
- † Macridin, Burov, Stern, Amundson, Spentzouris, Phys. Rev. ST Accel. Beams 18, 074401 (2015), [record](https://www.osti.gov/pages/biblio/1201004-simulation-transverse-modes-intrinsic-landau-damping-bunched-beams-presence-space-charge)
- Taken from the reference list of Karpov et al. (2021), not accessed: Lebedev, Atomic Energy 25, 851 (1968); Sacherer, IEEE Trans. Nucl. Sci. 20, 825 (1973); Hofmann, Pedersen, IEEE Trans. Nucl. Sci. 26, 3526 (1979); Chin, Satoh, Yokoya, Part. Accel. 13, 45 (1983); Oide, Yokoya, KEK Preprint 90-10 (1990); Boine-Frankenheim, Shukla, Phys. Rev. ST Accel. Beams 8, 034201 (2005); Boine-Frankenheim, Chorniy, Phys. Rev. ST Accel. Beams 10, 104202 (2007)

Fluids:

- † Schecter, Dubin, Cass, Driscoll, Lansky, O'Neil, Phys. Fluids 12, 2397 (2000), [text](https://www.cora.nwra.com/~schecter/pubs/schecter00_pf.pdf)
- † Balmforth, Llewellyn Smith, Young, J. Fluid Mech. 426, 95 (2001), [journal](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/disturbing-vortices/30D4DB569B343F2C793AB26C4A1533EC)
- Bedrossian, Coti Zelati, Vicol, [arXiv:1711.03668](https://arxiv.org/abs/1711.03668)

Numerical methods:

- † Sellwood, J. Comput. Phys. 50, 337 (1983), [journal](https://www.sciencedirect.com/science/article/abs/pii/002199918390102X)
- Yoshikawa, Yoshida, Umemura, [arXiv:1206.6152](https://arxiv.org/abs/1206.6152)
- Halle, Colombi, Peirani, Astron. Astrophys. 621, A8 (2019), [arXiv:1701.01384](https://arxiv.org/abs/1701.01384)
