# Spherically symmetric Vlasov–Poisson with particle methods: literature and research directions

4 October 2026

The code is treated here as a tool. The aim is to find physics questions for which a
spherical particle-in-cell (PIC) code gives new results, not to describe the code.

Two earlier documents are part of this review and are not repeated: the search of
20 September on echoes, anisotropy and ℓ = 0 modes
(`VlasovPoisson_PIC_sp/docs/literatura_2026.md`) and the review of oscillation against
damping (`docs/hadzic/literature_review.md`).

## 1. Scope, sources and level of verification

Scope: the gravitational Vlasov–Poisson system in spherical symmetry, solved with PIC,
shell or other particle methods, and the phase-space (Eulerian or semi-Lagrangian) and
N-body work on the same problems. Electrostatic and relativistic work is included where the
method or the mathematics carries over.

Sources: keyword searches on OpenAlex and arXiv with the combinations of the request and
others (shell code, Hénon sphere, splitting scheme, self-similar infall, kinetic blocking);
reference lists and citing works of the papers found; the two earlier reviews.

Verification. Every work in the reference list was found in OpenAlex, arXiv or a publisher
page; its title, authors and year come from that record. One recollection was wrong and was
corrected by the check (the third author of the rotation-induced transition paper is
Prokhorenkov). Unless stated otherwise, what is said about a paper comes from its abstract.
Sections of the full text were read for Halle, Colombi and Peirani (2019), Sellwood (2015),
and the works of the earlier review. Statements marked *inference* are mine and are not in
the papers.

## 2. The codes

| | `vlasov-poisson_PIC` | `VlasovPoisson_PIC_sp` |
|---|---|---|
| phase space | (r, p_r), all particles with the same \|L\| = L₀ | (r, p_r, L), a distribution in L |
| field | monopole on a radial grid, M(<r); self-gravity and analytic backgrounds (isochrone, uniform ball or point mass, NFW, Burkert) | the same |
| particles | spherical shells; quiet start on a regular lattice in the angle–action variables of the equilibrium, perturbation carried by the weights | the same, with a quadrature in L |
| companion tools | self-consistent equilibria, angle–action map, linearized time-domain solver, loop gain λ(ω) (the norm of the Mathur operator) | equilibria F(J, L), exact solution without self-gravity, true angle–action map |
| limits | L₀ ≥ 0.25 in use; only ℓ = 0 | only ℓ = 0; linear solver with L not written |

Both codes see only spherically symmetric perturbations. Non-spherical instabilities are
outside them (Section 3.6).

## 3. The literature by scientific question

### 3.1 Equilibria and their stability

*Studied.* Existence and stability of steady states f(E, L). The classical results are
Antonov's criterion and the theorem of Doremus, Feix and Baumann (1973): a distribution that
decreases with energy is stable against radial perturbations whatever its dependence on L.
Nonlinear orbital stability is proved by Guo and Rein (2007) and by Lemou, Méhats and
Raphaël (2011).

*Numerical.* Hénon (1973) evolved generalized polytropes with an N-body code with enforced
spherical symmetry and found the first spherical equilibria that are unstable on a
dynamical time. Barnes, Goodman and Hut (1986) confirmed this radial instability with a code
that kept forces to quadrupole order, and found two non-radial ones. Dattathri et al. (2025)
show that double-power-law spheres with an inflection of the distribution function,
df/dE > 0, grow an ℓ = 1 mode that saturates into a long-lived soliton, which traps
particles and erodes the bump; they compare it with the bump-on-tail instability.

*Open.* The nonlinear saturation of the purely radial instability, and whether its end state
is a pulsating solution of the kind built by Louis and Gerhard (1988), was not found
studied.

*Distance to the codes.* Close. The radial instability is the one instability that a
spherical code sees, and λ(ω) of `lambda_borde.py` already handles non-monotonic F.

### 3.2 Linear response: phase mixing, damping and discrete modes

Covered in `docs/hadzic/literature_review.md`. In short: the group of Hadžić and Rein
proves when radial perturbations damp or oscillate, with a Birman–Schwinger principle due
to Mathur; Ramming and Rein (2018) and Straub (2024) study it with a spherical PIC code; the
threshold exponent of isotropic polytropes (k ≈ 1.2–1.3) is not determined, and the
numerical evaluation of the Mathur operator is listed as open by Kunze and by Straub.
Decay rates for pure transport and for the linearized system are in Hadžić et al. (2024),
Hadžić and Schrecker (2025), Chaturvedi and Luk (2021, 2024), Rioseco and Sarbach (2020).
Rozier and Errani (2024) add a result of a different kind (Section 3.9).

*Distance.* This is where the codes have been used. The fixed-|L| results are in the
report `docs/hadzic/hadzic_informe.pdf`.

### 3.3 Cold collapse and violent relaxation

*Studied.* The collapse of a cold or cool sphere and the state it reaches. Hénon (1964)
followed a spherical cluster with a shell model. Van Albada (1982) showed with N-body runs
that cold irregular collapses give the r^{1/4} law when the collapse factor is large.
Aarseth, Lin and Papaloizou (1988) derived that the collapse factor of a cold uniform sphere
of N particles scales as N^{1/3}; Boily, Athanassoula and Kroupa (2002) verified it in
spherical symmetry and found N^{1/6} for axisymmetric configurations. Joyce, Marcos and
Sylos Labini (2009) found that the ejected energy grows as N^{1/3} and the ejected mass
about logarithmically. Sylos Labini (2012) separates mild from violent relaxation by the
initial virial ratio: the violent case gives a density tail r⁻⁴. Roy and Perez (2004) and
Trenti, Bertin and van Albada (2005) give systematic N-body studies.

*In phase space and in spherical symmetry.* Fujiwara (1983) integrated the collisionless
Boltzmann equation for spherical systems with a splitting scheme and found that after the
collapse of a uniform sphere the core is partially degenerate, not Maxwellian. Hozumi,
Fujiwara and Kan-ya (1996) and Hozumi, Burkert and Fujiwara (2000) used the same method:
the core shrinks as the initial radial anisotropy grows, so anisotropy controls the cusp.
Merrall and Henriksen (2003) find Gaussian velocity distributions in the centre with a
Vlasov integrator and a tree code. Colombi et al. (2015) compare a spherical Vlasov solver
with Gadget for Hénon spheres of virial ratio 0.5 and 0.1: the agreement is very good, but
for the colder case the N-body result is not converged even with 10⁸ particles. Halle,
Colombi and Peirani (2019) use a Vlasov solver in (r, v_r, j), a shell code with 10⁷ shells
and N-body runs for power-law spheres with Gaussian velocities. They separate three phases:
violent relaxation to a quasi-steady state whose phase-space density is a smooth spiral
consistent with self-similar predictions; a phase in which small-scale radial instabilities
destroy the spiral without changing coarse-grained properties; and, for the cool and steep
cases, the radial-orbit instability.

*Open.* In Halle et al. the time at which the small-scale radial instabilities appear
depends on the grid in the Vlasov code and on N in the particle codes, so they are seeded
numerically. Whether the spiral is physically unstable, and with what growth rate, is not
settled. They also note that Henriksen and Widrow (1997), with few shells, may have missed
the intermediate phase.

*Distance.* Medium. Cold collapse needs the centre and small pericentres, which the
fixed-|L| code avoids only for L₀ not small. A grid code is not the best tool for cold
caustics (Section 5).

### 3.4 The relaxed state: statistical theories

Lynden-Bell (1967) proposed a statistical mechanics of violent relaxation. Tremaine, Hénon
and Lynden-Bell (1986) introduced H-functions as mixing measures. Arad and Lynden-Bell (2005)
showed an inconsistency in these theories, and Arad and Johansson (2005) tested the
Lynden-Bell and Nakamura theories with test particles in N-body runs. Levin, Pakter and
Rizzato (2008) find that initial conditions that satisfy the virial condition relax to a
Lynden-Bell distribution with a cutoff, and that the others oscillate, evaporate mass and
form a core and a halo whose mass they predict without free parameters; the mechanism is a
parametric resonance with the bulk oscillation. Teles, Levin and Pakter (2011) do the same
for sheets. Williams and Hjorth (2010) and Pontzen and Governato (2013) derive the relaxed
distribution from maximum entropy in energy or in action space; the former compare it with
a shell code that relaxes only in energy and leaves the angular momentum of each shell
unchanged. Sylos Labini and
Capuzzo-Dolcetta (2020) conclude that the quasi-stationary states are not universal and
depend on the energy exchanged during the collapse.

A disagreement. Beraldo e Silva et al. (2017, 2019) estimate the entropy in N-body collapses
and in orbit ensembles in fixed potentials, find a fast increase that converges with N, and
conclude that the Vlasov–Poisson equation does not describe violent relaxation, the cause
being discreteness. Colombi et al. (2015) find close agreement between a Vlasov solver and
N-body runs for the same kind of collapse. The two statements are not about the same
quantity (an entropy estimate against the phase-space density), and they have not been
confronted in one controlled setting.

*Distance.* Medium. The parametric-resonance picture and the entropy question can be posed
in spherical symmetry with shells; they need non-equilibrium initial data.

### 3.5 Halo formation in spherical symmetry

Fillmore and Goldreich (1984) and Bertschinger (1985) gave the self-similar solutions of
secondary infall with radial orbits. Angular momentum was added by White and Zaritsky
(1992), Sikivie, Tkachev and Wang (1997), Nusser (2001), Le Delliou and Henriksen (2003) and
Zukin and Bertschinger (2010); the inner slope depends on how angular momentum is assigned
at turnaround. Łokas and Hoffman (2000) obtain inner slopes between r⁻² and r⁻²·³ with
radial orbits and state that angular momentum makes the profile shallower. Henriksen and
Widrow (1997, 1999) follow the collapse with a shell code: a self-similar phase, then an
instability of the similarity solution that relaxes the system. Lu, Mo and Katz (2006) use
one-dimensional simulations and show that isotropization during fast accretion gives
ρ ∝ r⁻¹. Vogelsberger, Mohayaee and White (2011) run the same initial conditions in three
dimensions: the haloes are triaxial, and the structure remembers the initial conditions.

*Distance.* Far. These problems need an expanding background or continuous infall, and
cold initial data. They are well covered by shell codes.

### 3.6 The radial-orbit instability and what spherical symmetry removes

Merritt and Aguilar (1985) and Barnes, Goodman and Hut (1986) established the instability
of radially anisotropic spheres numerically; Polyachenko and Shukhman (2015) and Maréchal
and Perez (2011) discuss its nature. Barnes, Lanzel and Williams (2009) find that the
initial anisotropy together with the virial ratio decide whether it occurs, and Bellovary et
al. (2008) that an isotropic velocity dispersion suppresses it. Huss, Jain and Steinmetz
(1999) and MacMillan, Widrow and Henriksen (2006) ran collapses with and without the
non-radial forces: without them the density is a pure power law, the energy distribution is
close to Boltzmann and the orbits are radially biased; with them the profile becomes of the
NFW form and the centre isotropic.

Consequence for a spherical code. A spherical run is a controlled experiment in which this
instability is switched off. That is useful as a comparison, as in MacMillan et al., but a
result obtained for a radially anisotropic model in spherical symmetry may not survive in
three dimensions. Results should be claimed for models that are stable to non-spherical
perturbations, or stated as constrained to spherical symmetry.

### 3.7 Fixed angular momentum: the shell model and its phase transition

The model of the fixed-|L| code has a literature in statistical mechanics. Miller and
Youngkins (1998) and Youngkins and Miller (2000) study concentric spherical mass shells and
find two phases, quasi-uniform and centrally concentrated, in mean-field theory and in
dynamical simulation. Klinko and Miller (2000) add the conserved sum of squared angular
momenta. Klinko, Miller and Prokhorenkov (2001) take every particle with specific angular
momentum of the same magnitude l, prove that the entropy is bounded, and find a phase
transition when l falls below a critical value, in the microcanonical and canonical
ensembles. Klinko and Miller (2002, 2004) simulate N rotating shells of fixed l and see the
transition between quasi-uniform and core–halo states, with power-law relaxation. Di Cintio
and Ciotti (2011) use shells with a force 1/r^α and find that virialization and phase-mixing
times depend on α.

*Open.* These works study thermal equilibrium, reached by collisions among few shells in a
confined system. The collisionless evolution of the same model, the quasi-stationary state
reached after violent relaxation as a function of l, was not found.

*Distance.* Very close: this is the model of the code. A confining wall is needed for
the thermodynamic comparison.

### 3.8 Time-dependent potentials: adiabatic and impulsive response

Young (1980) computed the adiabatic growth of a central black hole in a star cluster;
Quinlan, Hernquist and Sigurdsson (1995) confirm the cusps with N-body runs and find them
insensitive to the initial anisotropy. Gondolo and Silk (1999) and Ullio, Zhao and
Kamionkowski (2001) apply it to dark-matter spikes. Sellwood and McGaugh (2005) implement
Young's algorithm for the compression of haloes by baryons and find that haloes with random
motion resist compression. Pontzen and Governato (2012) explain cores by repeated fast
changes of the central potential, with an impulsive approximation; Burger and Zavala (2019)
show that radial actions are conserved when the core forms adiabatically and change when it
forms impulsively. Errani et al. (2025) study the response of tracer populations to
impulsive and adiabatic changes of the potential. Kandrup, Vass and Sideris (2003) find that
an oscillating potential produces transient chaos by resonance between orbital and driving
frequencies.

*Distance.* Close after a small change: a background that depends on time. The response in
action space is the natural diagnostic of the codes.

### 3.9 Return to equilibrium after a perturbation

Rozier and Errani (2024) show with the matrix method that the phase-space distribution of
the final virialized state follows from the initial disequilibrium without computing the
evolution. For energy-truncated Hernquist spheres, a model of tidal stripping, the linear
prediction agrees with N-body runs at the per cent level where the response is linear, and
departs where it is not. Errani and Navarro (2021) give the N-body phenomenology of tidal
remnants.

*Distance.* Close. The fixed-|L| linear solver already gives the time evolution; the static
part of the response is the prediction to compare.

### 3.10 Echoes and nonlinear effects

Chiba et al. (2025) derive echo theory in angle–action variables and apply it to a
one-dimensional model of the vertical motion of the Galaxy (earlier notes, Section 1.5).
O'Neil (1965) gives the nonlinear evolution of Landau damping by trapping. For stellar
systems, trapping is studied for bars and for the dipole mode (Dattathri et al. 2025); no
study of trapping for radial modes of spheres was found, which may reflect the search terms.

### 3.11 Finite N: relaxation, noise and kinetic theory

Weinberg (1998) predicts with the dressed-particle formalism that self-gravity amplifies
Poisson noise, by about six for the dipole of a King model. Lau and Binney (2019) measure in
10⁴ realizations of a cluster an amplification of more than ten, dominated by the dipole
mode, and conclude that local-scattering theory is qualitatively wrong. Sellwood (2015)
measures relaxation in spherical N-body models with four field methods: diffusion rates scale
as N⁻¹ in inhomogeneous models, but as N^{−1/2} in the uniform sphere, where the modes
excited by shot noise are almost undamped and dominate the relaxation. Hamilton and
Heinemann (2020, 2023) unify noise and waves, and prove for a homogeneous system that the
divergence of the wake at marginal stability is cancelled by the mode, so that calculations
with the Balescu–Lenard equation that ignore modes may need revision.

One-dimensional systems. Joyce and Worrakitpoonpon (2010) find that sheets relax in a time
proportional to N, faster for colder initial states. Roule, Fouvry and Pichon (2022)
validate the Balescu–Lenard equation for one-dimensional self-gravitating systems and find
that collective effects reduce the diffusion tenfold. Fouvry (2022) and Fouvry and Roule
(2023) show that in one-dimensional inhomogeneous systems with a monotonic frequency profile
and only 1:1 resonances the Balescu–Lenard flux vanishes exactly: the relaxation is then
driven by 1/N² effects (kinetic blocking).

From the full text. The exact blocking of Fouvry (2022) holds for pairwise interactions
that depend on the difference of the angles, for which only 1:1 resonances exist; the
systems simulated are particles on a sphere and classical spins. For self-gravitating
sheets, Roule et al. (2022) find a quasi-blocking, a flux 10⁵ times smaller than the
diffusion coefficients, and give four reasons: with a monotonic Ω(J) the resonances k = k′
are local and give no flux; a symmetry forbids k and k′ of different parity; the finite
range of frequencies requires k/k′ ≤ Ω(0)/Ω(J); and the couplings fall as 1/k². What is
blocked is the flux, the change of F(J). The diffusion of each particle is not: the local
resonances contribute to it, and it scales as 1/N.

*Inference.* The interaction between two shells, −1/max(r, r′), is not a function of the
difference of the angles, so resonances with k ≠ k′ are allowed, and there is no parity
rule because pericentre and apocentre are not equivalent. The third reason of Roule et al.
applies: a resonance k:k′ needs Ω_max/Ω_min ≥ k′/k, so the 1:2 resonance exists only when
the band is wider than an octave, which is the single-gap condition of Hadžić et al.
(2023). The flux of this system should therefore be strongly suppressed for narrow bands
and of order 1/N for wide ones. The band width as a control parameter is the new element;
the mechanism is theirs.

*Distance.* Very close. Each particle of the code is a shell, so the code is a realization
of this system.

### 3.12 The relativistic counterpart

Shapiro and Teukolsky (1985) and Rasio, Shapiro and Teukolsky (1989) solve the spherical
Einstein–Vlasov system with particles and in phase space. Rein, Rendall and Schaeffer
(1998), with a particle code, and Olabarrieta and Choptuik (2002), with a particle-mesh
method in which each particle is a shell with its own angular momentum, find type I critical
behaviour at the threshold of black-hole formation; the critical solutions are unstable
static solutions. Andréasson and Rein (2006) relate the threshold to the stability of
steady states. Akbarian and Choptuik (2014) integrate in phase space, find no universality
of the critical solution, and clarify the role of angular momentum. Rein and Rodewis (2003)
prove the convergence of a PIC scheme for this system, and Schaeffer (1987) of a particle
scheme for Vlasov–Poisson.

### 3.13 Electrostatic problems with spherical symmetry

Boella, Coppa and D'Angola (2017) present gridless particle techniques with shells or rings
for collisionless plasmas with spherical or axial symmetry, suited to regions that change
fast, such as Coulomb explosions. The method is the plasma version of the shell code.

## 4. Canonical papers

| paper | question | method | result | relevance |
|---|---|---|---|---|
| Lynden-Bell 1967 | why collisionless systems relax | statistical theory | violent relaxation, coarse-grained equilibrium | background |
| Hénon 1964, 1973 | collapse and stability of spheres | shells, enforced spherical symmetry | first unstable spherical equilibria (radial) | direct: same symmetry |
| Doremus, Feix, Baumann 1973 | radial stability of f(E, L) | energy principle | df/dE < 0 is enough | direct: defines what can be unstable |
| van Albada 1982 | end state of cold collapse | N-body | r^{1/4} law for large collapse factor | background |
| Fujiwara 1983 | phase space of collapse | splitting scheme in (r, v_r, L) | degenerate core | direct: same reduced system |
| Barnes, Goodman, Hut 1986 | instabilities of spheres | N-body to quadrupole order | radial instability confirmed; two non-radial ones | direct |
| Mathur 1990 | existence of oscillations | reduction of the eigenvalue problem | criterion | direct |
| Henriksen, Widrow 1997 | cold spherical collapse | shell code | self-similar phase, then instability | medium |
| Weinberg 1994, 1998 | weakly damped modes, noise | matrix method | modes persist; noise amplified | direct for finite N |
| Levin, Pakter, Rizzato 2008 | end state of relaxation | theory and N-body | core–halo from parametric resonance | medium |
| Klinko, Miller, Prokhorenkov 2001 | thermodynamics with fixed l | mean field | transition at a critical l | direct: same model |
| MacMillan, Widrow, Henriksen 2006 | role of non-radial forces | runs with and without them | power law without, NFW with | direct: meaning of a spherical run |
| Sellwood 2015 | numerical relaxation | N-body, four field methods | N⁻¹, but N^{−1/2} with undamped modes | direct |
| Colombi et al. 2015; Halle et al. 2019 | Vlasov against N-body | spherical Vlasov solver, shells, tree code | agreement; three dynamical phases | direct: closest numerical work |
| Hadžić et al. 2022, 2023 | oscillation or damping | spectral theory | dichotomy in the edge exponent | direct |
| Straub 2024 | which models oscillate | spherical PIC | thresholds by inspection | direct: closest physical work |
| Fouvry 2022; Fouvry, Roule 2023 | relaxation in 1D | kinetic theory, N-body | kinetic blocking | direct |
| Rozier, Errani 2024 | final state after a perturbation | matrix method | linear prediction at per cent level | direct |
| Chiba et al. 2025 | echoes in galaxies | angle–action theory, 1D model | echoes and their damping | direct for the code with L |

## 5. Numerical methods compared

| method | field | strong points | weak points | used by |
|---|---|---|---|---|
| shell code | exact, M(<r) by sorting | no grid, any dynamic range, natural for cold collapse | force discontinuous at shell crossings; each shell is a strong scatterer | Hénon; Henriksen and Widrow; MacMillan et al.; Halle et al. |
| spherical PIC | monopole on a radial grid | cheap, low diffusion, weights can carry the perturbation, smooth force | field error of order Δr²; grid limits cold caustics; particle noise | Ramming and Rein; Straub; this work |
| Eulerian or semi-Lagrangian in (r, v_r, L) | grid | no particle noise, fine phase-space detail | diffusion and aliasing of filaments; cost N_r N_v N_L (2048 × 2048 × 128 in Halle et al.) | Fujiwara; Hozumi; Rasio et al.; Colombi et al.; Akbarian and Choptuik |
| N-body in 3D, or expansions with ℓ > 0 | tree, grid or basis | sees non-spherical instabilities | noise; cost; convergence not reached for cold cases with 10⁸ particles | van Albada; Barnes et al.; Sellwood |
| matrix (linear response) | basis functions | exact linear answer, no noise | linear only; slow convergence for ℓ = 0 (Polyachenko and Shukhman 2026) | Weinberg; Fouvry and Prunet; Rozier and Errani |

Where a spherical PIC code with a quiet start is the better tool:

- small perturbations of an equilibrium followed for tens of dynamical times, where the
  weights carry the signal and a reference run removes the common error;
- questions in which the number of shells is the physical parameter (Section 3.11);
- parameter scans: a run costs minutes, so hundreds of runs are possible;
- diagnostics in the true angle–action variables, particle by particle.

Where it is worse: cold collapse and caustics (shell or phase-space-sheet methods);
fine-grained phase-space structure at late times (limited by recurrence of the lattice and
by noise); anything non-spherical; small pericentres, which force a small time step.

## 6. Angular momentum, anisotropy and perturbations under numerical control

**Angular momentum.**

| choice | what it gives | control |
|---|---|---|
| one value L₀ | an exact reduction to one degree of freedom; the centre is empty | time step set by the pericentre, Ω_p Δt small; L₀ ≥ 0.25 in use |
| a quadrature in L (N_L nodes) | a distribution F(J, L), anisotropy by the width σ_L | recurrence in L at T = 2π/(k \|∂Ω/∂L\| ΔL), measured at t = 1500 with N_L = 16; Halle et al. use 128–512 slices |
| random or low-discrepancy sampling in L | no recurrence | shot noise in L |
| L assigned by a rule, L ∝ √(GMr) at turnaround | the secondary-infall models of Section 3.5 | belongs to collapse problems |

**Anisotropy.** Families with known equilibria: F(J, L) with a Gaussian in L; f = L^{−2β}g(E);
Osipkov–Merritt; generalized polytropes (E₀ − E)ᵏ Lˡ. By the theorem of Doremus et al. the
radial stability depends only on the sign of ∂f/∂E, so anisotropy changes frequencies and
damping rates but does not create a radial instability. Radially anisotropic models are
candidates for the radial-orbit instability, which the code does not see (Section 3.6);
tangential anisotropy is safe in that respect. Hozumi et al. (1996) warn that the initial
anisotropy is a poor indicator of that instability in a collapse.

**Perturbations.**

| way | property | check |
|---|---|---|
| in the weights, δF = ε F s(J) cos kQ on the lattice | no shot noise; mass and mean profile unchanged | s(J) must vanish at the circular orbit as J^{k/2} or faster; linearity by halving ε |
| difference with a reference run of ε = 0 | removes the common discretization error | same particles, grid and time step in both runs |
| an impulsive external potential δΦ(r)δ(t) | a perturbation that a physical event would produce; two kicks give an echo | amplitude in the linear range |
| a dilation r → (1 + ε)r | excites the scale-invariant solution of Polyachenko and Shukhman (2023), which is exact | use it as a test of the code |
| a virial ratio different from one | non-equilibrium initial data | no quiet start in angle–action variables; a lattice in (r, v_r, L) instead |

The amplitude must be chosen against O'Neil's parameter ν = ω_b/γ: near a threshold the
damping rate tends to zero, so the evolution is nonlinear for any finite ε at late times
(battery η, Section 1; measured in the demo).

## 7. Parameters

| parameter | what the literature suggests | regime that is little explored |
|---|---|---|
| L₀ (one value) | a phase transition in thermal equilibrium at a critical l (Klinko et al. 2001); sets the width of the band, Ω_max/Ω_min | collisionless evolution as a function of l |
| width in L, anisotropy | changes cusps after collapse (Hozumi et al. 2000); cusp around a black hole insensitive to it (Quinlan et al. 1995); suppresses the unwinding of a spiral (measured, earlier notes) | damping of ℓ = 0 modes against anisotropy |
| mass of the component, a₀ | threshold λ_edge = 1 for a discrete mode | relaxation and noise across the threshold |
| edge exponent k | dichotomy at small mass (Hadžić et al. 2023); threshold near 1.2 for isotropic polytropes (Straub 2024) | the exact threshold |
| virial ratio | mild or violent relaxation (Sylos Labini 2012); core–halo against Lynden-Bell (Levin et al. 2008) | with a fixed L₀ |
| perturbation amplitude ε | linear below ν ≈ 1; trapping above | radial modes of spheres |
| perturbation scale (harmonic k in Q) | recurrence time and decay exponent depend on k | higher harmonics |
| number of particles N | N⁻¹ or N^{−1/2} (Sellwood 2015); N⁻² under kinetic blocking (Fouvry and Roule 2023) | shells with fixed L |
| Δr, Δt, N_L | field error of order Δr² (measured); recurrence | — |

Varying a parameter gives a paper only where the literature predicts a threshold or a
scaling to test. Of the list, that holds for a₀ and k (the threshold), N (the scaling of the
relaxation), L₀ through the width of the band, and σ_L (echo suppression).

## 8. Analytical predictions that can be turned into tests

| prediction | source | test |
|---|---|---|
| a discrete mode exists if and only if λ_edge > 1; its frequency solves λ(ω) = 1 | Mathur; Hadžić et al. | done for fixed L; to do with a distribution in L |
| binding laws near the edge, by class of k | battery η; Simon; Klaus and Simon | measure δ_d against the coupling |
| decay t^{−(k+1)} of the edge response; bounds t^{−min(2,k)} | Hadžić et al. 2024 | slopes of the envelope, done at small mass |
| the flux in J is suppressed when no resonance k ≠ k′ fits in the band; the diffusion of each particle is not | Fouvry 2022; Roule et al. 2022 | change of F(J) against N for narrow and wide bands, with equal-mass random sampling |
| amplification of noise by self-gravity | Weinberg 1998; Lau and Binney 2019 | power spectrum of δΦ in unperturbed runs |
| final state after a perturbation from linear theory | Rozier and Errani 2024 | static part of the response against ε |
| echo amplitude, and its loss with a spread in L, exp[−½k²(∂Ω/∂L)²σ_L²t*²] | Chiba et al. 2025; earlier notes | two kicks, against σ_L |
| collapse factor ∝ N^{1/3} for a cold uniform sphere | Aarseth et al. 1988 | a test of the code with L₀ → small |
| core–halo masses from the resonance with the bulk oscillation | Levin et al. 2008 | collapse with a virial ratio different from one |
| adiabatic invariance of J under slow changes; Young's algorithm | Young 1980; Sellwood and McGaugh 2005 | grow a central mass slowly and fast |
| critical l of the fixed-l transition | Klinko et al. 2001 | confined system, long runs with few shells |
| the scale-invariant ℓ = 0 solution | Polyachenko and Shukhman 2023 | exact control for any equilibrium with one length |
| virial theorem; conservation of energy and of F along orbits | — | standard checks |

## 9. Observables

Used in the literature and available in the runs:

- phase-space density in (r, p_r) and in (Q, J); harmonics f_k(J, t), now stored for the
  34 runs of the Hadžić setting;
- the moment h_k with a test function, and its poles;
- δΦ(r, t) and its power spectrum;
- the diffusion of the action of each particle, ⟨ΔJ²⟩(t), and the flux in J;
- energy distribution N(E), and its change (Pontzen and Governato; MacMillan et al.);
- velocity moments and anisotropy β(r); pseudo-phase-space density ρ/σ³;
- virial ratio and its oscillation (bulk mode);
- entropy estimates and H-functions (Tremaine et al. 1986; Beraldo e Silva et al.);
- recurrence times; energy error.

Observables that may hold unused information in the existing runs:

- the displacement of J in the reference runs (ε = 0) against N: read in Section 13; it
  measures the breakdown of the quiet start, not the relaxation law of Section 3.11;
- the profile in J of f₁ for the discrete modes, to compare with the eigenfunction;
- the static part f₀(J, ∞) − f₀(J, 0) against ε: the final state of Section 3.9;
- the second harmonic and the mean shift of J at large ε: the trapping regime.

## 10. Physics against numerical effects

| phenomenon | artefact that imitates it | test |
|---|---|---|
| a mode that does not damp | static residue from the wrong angle–action map; recurrence of the lattice | true map of the self-consistent potential; T_rec beyond the run |
| slow damping of a mode | dip of the amplitude at finite ε; noise of the reference run | scan in ε; four times more particles |
| non-stationary equilibrium | field error of order Δr²; growth from discreteness | Δr/2 at fixed Δt; 4N |
| instability of a phase-space spiral | seeded by the grid or by Poisson noise (Halle et al.) | growth rate against N and grid; seed a known perturbation |
| relaxation, entropy growth | discreteness (Beraldo e Silva et al.); kinetic blocking changes the N law | ensembles at several N; random and quiet starts |
| growth of noise at late times | breakdown of the ordered start, as the multibeam instability of plasma codes (Gitomer and Adam 1976) | growth rate against the number of rows in J; random start for comparison |
| cusp or core after collapse | softening, time step at pericentre, N^{1/3} collapse factor | L₀ and Δt scans; compare with a shell code |
| enhanced response near a threshold | nonlinearity (ν large) | ε → 0 extrapolation |
| loss of mass or energy | outer boundary of the grid | larger r_max |

## 11. Gaps

For each gap: why it is one, and how sure that is.

1. **The exact oscillation threshold.** Stated as open by Kunze (2022), Straub (2024) and
   Wolfschmidt (2023). A genuine gap. Risk: Mathur (1990) was read only in abstract.
2. **Relaxation of shells with fixed L, and the role of the band width.** The papers were
   read. The mechanism, a finite frequency range that excludes resonances, is described by
   Roule et al. (2022) for sheets. A system in which the band ratio is a control parameter
   that crosses 2 was not found. A gap of moderate size: an application of a known
   mechanism, with a clear test.
3. **Fluctuations and relaxation across the threshold of a discrete mode.** Sellwood's
   N^{−1/2} and the cancellation of Hamilton and Heinemann point to open theory. No test in
   a system where the threshold can be crossed continuously was found. Probably genuine.
4. **Nonlinear evolution of radial modes near the threshold.** Straub reports the same
   qualitative behaviour for the linearized and the nonlinear system; the demo measured
   trapping and a loss that depends on ε. Absent from the papers found; the search for
   trapping in stellar dynamics was not exhaustive.
5. **Saturation of the radial instability.** Hénon and Barnes et al. found the instability;
   the saturation of the ℓ = 1 analogue was described in 2025. The ℓ = 0 case was not found.
   May be simply unsearched.
6. **Collisionless dynamics of the fixed-l model.** The equilibrium thermodynamics is done.
   A gap, of interest mainly to statistical physics.
7. **Whether the phase-space spiral of a quasi-steady state is physically unstable.** Posed
   by Halle et al., who call it rather academic because the radial-orbit instability
   dominates in three dimensions.
8. **Entropy growth in violent relaxation: discreteness or coarse graining.** Two groups
   disagree. A real disagreement, but partly about definitions.
9. **ℓ = 0 damping against anisotropy, and echoes with a spread in L.** From the earlier
   notes; the linear theory for ℓ = 0 is being developed now (Polyachenko and Shukhman 2026).

## 12. Research directions

Scores from 1 to 5; cost: 5 is cheap.

| | direction | relevance | originality | feasibility | cost | analytic support | clear result | publishable | sum |
|---|---|---|---|---|---|---|---|---|---|
| P1 | exact oscillation threshold of polytropes and King models | 4 | 4 | 3 | 3 | 5 | 5 | 4 | 28 |
| P2 | relaxation of fixed-L shells: kinetic blocking, band width and the discrete mode | 4 | 3 | 4 | 4 | 4 | 4 | 3 | 26 |
| P3 | echoes and unwinding with a spread in L | 3 | 4 | 5 | 5 | 4 | 4 | 3 | 28 |
| P4 | nonlinear evolution of modes near the threshold | 3 | 3 | 5 | 5 | 3 | 3 | 3 | 25 |
| P5 | final state after a perturbation: linear prediction and its failure | 3 | 2 | 5 | 5 | 5 | 4 | 3 | 27 |
| P6 | radial instability of a non-monotonic F and its saturation | 3 | 4 | 3 | 4 | 3 | 3 | 3 | 23 |
| P7 | response to a time-dependent potential in action space | 4 | 2 | 4 | 4 | 4 | 3 | 3 | 24 |
| P8 | collisionless dynamics of the fixed-l model against mean field | 2 | 4 | 3 | 3 | 4 | 3 | 2 | 21 |
| P9 | violent relaxation and entropy with a quiet start | 3 | 2 | 2 | 3 | 2 | 2 | 2 | 16 |

**P1. Exact oscillation threshold.** Question: for which steady states does a radial
perturbation oscillate for ever. Literature: Hadžić et al., Kunze, Straub. Hypothesis:
λ_edge(k) crosses one near k = 1.2–1.3 for isotropic polytropes, and the mode is bound by a
high power of the distance to the threshold, which is why time integration cannot decide.
Initial data: the equilibria of Straub. Control parameters: k, and κ for King. Observables:
λ_edge, ω_d, the binding δ_d; in PIC runs, f₁(J, L, t) and h₁. Numerical tests: resolution of
the operator in J, L and harmonics; PIC runs with ε → 0 and 4N. Interesting: a sharp k* and a
test of 12/π². Inconclusive: a threshold that moves with the resolution. Novelty: physical
and numerical; high, subject to gap 1.

**P2. Relaxation of fixed-L shells.** Question: how does a finite number of shells relax,
and does the law change with the width of the band and with the presence of a discrete
mode. Literature: Fouvry (2022), Fouvry and Roule (2023), Roule et al. (2022), Sellwood
(2015), Hamilton and Heinemann (2023). Hypothesis: the flux in J, the rate of change of F(J), is of order
1/N when Ω_max/Ω_min exceeds 2 and much smaller below, while the diffusion of each particle
scales as 1/N in both cases; near λ_edge = 1 the dressing enhances both. Initial data:
equilibria sampled at random, with equal masses; with unequal weights the exchange of
energy between masses adds a flux of its own (Sellwood 2015). Control parameters: N, the
band ratio through J_t/L₀, a₀. Observables: the change of F(J) (a Wasserstein distance
between cumulative distributions), ⟨ΔJ²⟩(t), spectrum of δΦ. Expected: a flux that drops
when the band closes the 1:2 resonance, softened by higher resonances (k:k+1 needs a ratio
above (k+1)/k). Numerical tests: Δr, Δt, ensembles, random against
quiet start. Interesting: a flux that depends on the band ratio as predicted, or a peak at the
threshold. Uninteresting: no dependence; even then the runs calibrate the numerical
relaxation of every other project. Novelty: physical but moderate (gap 2).

**P3. Echoes with a spread in L.** As T1 of the earlier notes. Hypothesis: a spread in L
suppresses the echo without collisions, as exp[−½k²(∂Ω/∂L)²σ_L²t*²]. Novelty: physical,
moderate; the mechanism is dephasing in a second action and may be regarded as expected.

**P4. Nonlinear evolution near the threshold.** Uses the demo and series 4 of the Hadžić
setting. Hypothesis: the deficit of amplitude scales between ε and ε², and trapping sets in
at ν ≈ 1. Weakness: no analytic law to test beyond O'Neil's estimate. It is better as a
section of P1 than as a paper.

**P5. Final state after a perturbation.** The linear prediction exists (Rozier and Errani)
and was tested with N-body runs. The contribution would be the range of validity in ε with
a noise-free code. Mostly a methodological addition.

**P6. Radial instability and saturation.** Build a non-monotonic F(J), compute the growth
rate with λ(ω), follow the saturation. Interesting only if the end state is a persistent
pulsation that can be related to the oscillating solutions of Section 3.2.

**P7. Time-dependent potential.** The physics is known from N-body work (Burger and Zavala);
the contribution would be the resolved response in action space and its comparison with
linear response. Parameter exploration unless a sharp prediction is found.

**P8. Fixed-l model.** A statistical-physics paper at most.

**P9. Violent relaxation and entropy.** Not recommended: the best tools are shell and
phase-space codes, and the question is partly one of definitions.

Ideas that sound interesting and are unlikely to give a strong paper: dark-matter spikes
and adiabatic contraction (done with N-body and Young's algorithm); secondary infall with
angular momentum (done with shell codes, and needs cosmology); the radial-orbit instability
(not visible in spherical symmetry).

## 13. Minimal campaigns for the first three

**P1.**
- Baseline: isotropic polytrope f ∝ (E₀ − E)ᵏ with k = 1, self-gravitating, radial
  perturbation.
- Operator: extend `lambda_borde.py` to nodes in (J, L); about 60 × 30 nodes and 6
  harmonics; largest eigenvalue by iteration.
- Sweep: k = 0.5, 0.75, 1, 1.1, 1.2, 1.25, 1.3, 1.5, 2; King models over the range of
  Straub.
- Validation: the fixed-L limit must reproduce the 61 equilibria already computed; a linear
  time-domain solver with L for three values of k; PIC runs with the code with L for k = 1,
  1.25 and 1.5 at ε = 0.03 and 0.1.
- Convergence: nodes in J and L doubled; harmonics 4, 6, 8; for PIC, 4N and Δr/2.
- Figures: λ_edge(k) with the crossing; δ_d(k) on a logarithmic scale; time series for
  three k with the linear solution; the observability number (Ω_min − ω_d)T.
- Finding of interest: k* to three digits and a binding law that explains the undecided
  window 1.2–1.3.

**P2.**
- Baseline: point-mass background, L₀ = 2, polytrope in energy with k = 2 (no discrete mode
  up to a₀ = 1), a₀ = 0.1, random sampling.
- Sweep: Ω_max/Ω_min = 1.2, 1.5, 2.0, 2.5 (J_t = 0.126, 0.289, 0.52, 0.714);
  N = 10³, 2 × 10³, 4 × 10³, 8 × 10³, 1.6 × 10⁴; 16 to 64 realizations each. Then a₀ across
  the threshold for k = 1.25 (a₀ from 0.3 to 0.8).
- Length: until ⟨ΔJ²⟩ is measurable, 10⁴ to 10⁵ time units. A run of 10⁴ particles and
  10⁵ steps takes about ten minutes; the ensembles suit the two cluster nodes.
- Diagnostics: the change of F(J) and ⟨ΔJ²⟩(t) in the true map; their N dependence;
  spectrum of δΦ.
- Convergence: Δr/2, Δt/2; quiet start for contrast.
- Figures: flux and diffusion coefficient against N for each band ratio; flux against the
  band ratio; both against a₀ across the threshold.
- Finding of interest: a flux that rises when the band opens the 1:2 resonance, with a
  diffusion coefficient that does not change.

**First reading of the existing reference runs** (`reproducir/scripts/relajacion_N.py`,
output in `exe/relajacion/`). The runs have a quiet start and unequal weights, so they do
not test the hypothesis. They measure the noise of the scheme. d₂ is the rms change of the
action since t = 0, weighted by mass, in units of 10⁻³ J_t; W₁ is the Wasserstein distance
between the cumulative mass distributions in J, the net change of F(J).

| run | N | d₂ at t_fin/4, t_fin/2, t_fin | growth after the first quarter |
|---|---|---|---|
| point mass, k = 1.25, a₀ = 1 (t_fin = 4000) | 10⁴ | 1.36, 6.27, 37.2 | t^{2.6} |
| the same, 4N | 4 × 10⁴ | 0.50, 0.58, 2.49 | exponential, rate 6.8 × 10⁻⁴ |
| the same, Δr/2 | 10⁴ | 1.71, 9.78, 42.5 | t^{2.4} |
| point mass, k = 2, a₀ = 1 (no discrete mode) | 10⁴ | 1.06, 4.36, 29.9 | t^{2.8} |
| point mass, k = 2, a₀ = 0.01 (t_fin = 8000) | 10⁴ | 0.013, 0.023, 0.036 | t^{0.74} |
| isochrone A4, a₀ = 0.075 (t_fin = 17200) | 10⁴ | 0.25, 0.41, 0.50 | t^{0.5} |
| the same, 4N | 4 × 10⁴ | 0.12, 0.12, 0.12 | none |
| isochrone M5, a₀ = 0.5, Δr/2 (t_fin = 5000) | 10⁴ | 1.33, 2.15, 7.19 | exponential, rate 4.4 × 10⁻⁴ |
| the same, 4N | 4 × 10⁴ | 0.21, 0.31, 0.40 | slow, about t^{0.5} |

- At a₀ = 1 the displacement grows faster than ballistically, reaches 3.7 % of J_t at
  t = 4000, and is the same for k = 1.25, 1.5 and 2 and with Δr/2. It does not depend on the
  presence of a discrete mode, nor on the grid.
- With 4N it is 12 to 15 times smaller at the end, and the run with 10⁴ particles reaches
  the final value of the 4N run at t = 1500. A diffusion by Poisson noise would be 2 times
  smaller. With a₀ = 0.01 it is a thousand times smaller and grows as t^{0.74}.
- The net change of F(J) is much smaller than the displacement: W₁ = 0.86 against
  d₂ = 37 at t = 4000 (k = 1.25, 10⁴ particles), and W₁ stays at 0.29 with 4N. The
  particles are reshuffled. This is what diffusion gives (a change of F of second order in
  the displacement), so it says nothing for or against blocking.
- The growth is the breakdown of the ordered start, strong where self-gravity is strong.
  Its dependence on the number of rows and its fast growth resemble the multibeam
  instability of ordered loadings in plasma simulations (Gitomer and Adam 1976); this was
  not tested. It sets the usable time of a quiet start: about t = 1500 with 10⁴ particles
  and beyond 4000 with 4 × 10⁴ at a₀ = 1.

**P3.**
- Baseline: isochrone background, no self-gravity, spiral or two kicks, exact solution
  available.
- Sweep: σ_L from 0 to 0.2; then a₀ = 10⁻³ and 10⁻².
- Resolution: N_L chosen so that the recurrence in L is beyond the echo time.
- Diagnostics: |h₁| at the echo time against σ_L; comparison with the Gaussian law.
- Finding of interest: the law holds without free parameters, and self-gravity changes it
  in a measurable way.

## 14. Roadmap

**A. Well understood.** Stability of f(E, L) with df/dE < 0 against radial perturbations.
The qualitative course of cold collapse and violent relaxation, and its dependence on the
virial ratio. Self-similar infall and the effect of angular momentum on the inner slope.
That the radial-orbit instability shapes haloes and is absent in spherical symmetry.
Adiabatic growth of a central mass. Existence of oscillating modes and of damping in the
linear theory.

**B. Still interesting.** The exact oscillation threshold. Relaxation and noise in systems
with weakly damped or undamped modes. Kinetic blocking in a gravitational system with a
tunable band. The nonlinear fate of radial modes. The saturation of the radial instability.

**C. Suited to spherical PIC.** Small perturbations of equilibria over tens of dynamical
times; questions in which the number of shells is the parameter; scans of hundreds of runs;
diagnostics in angle–action variables.

**D. In the runs already made.** The breakdown time of the quiet start against N (read: Section 13); the
profiles f₁(J) of the discrete modes; the static response against ε; the trapping regime of
the demo; the suppression of unwinding by σ_L.

**E. First parameter to vary.** N, with random equal-mass starts, in one equilibrium with a
narrow band and one with a wide band. It is cheap and decides P2. The existing runs cannot
do it.

**F. Strongest direction.** P1, the exact threshold, with P4 as its nonlinear section and
P2 as the finite-N counterpart.

**G. First experiments.** (1) Done: the ε = 0 runs read for the N law (Section 13). (2) Four
random-start ensembles, two band ratios and two values of N. (3) Done: the feasibility test of λ_edge with a
distribution in L gives k* = 1.2425 for isotropic polytropes (`docs/mathur/threshold_isotropic.md`).

**H. Analytical prediction to test.** That the flux in J is suppressed when no resonance
with k ≠ k′ fits in the band, and λ_edge(k*) = 1.

**I. Essential validation.** The true angle–action map; a reference run with ε = 0; scans
in ε, N and Δr at fixed Δt; recurrence times in J and in L beyond the run; energy
conservation; agreement with the linear solver before any nonlinear claim.

**J. Central question of a paper.** When does a spherical collisionless system keep
oscillating, and what ends the oscillation: phase mixing, trapping or particle noise. P1
answers the first part, P4 and P2 the second.

## 15. Limits of this review

- Most statements about papers rest on abstracts.
- The search on trapping for radial modes was not exhaustive; gaps 4 and 5 need the papers
  read before any claim. Fouvry (2022) and Roule et al. (2022) were read in full; Fouvry and
  Roule (2023) only in abstract.
- The literature before 1990 was reached through reference lists.
- OpenAlex refused requests during the last part of the search; the remaining checks were
  made on arXiv and on publisher pages.
- The inference of Section 3.11 on the band width is mine. It has not been checked against
  the kinetic equations for this system.

## References

All entries were verified as described in Section 1. The works on oscillation and damping
are listed in `docs/hadzic/literature_review.md`.

- Aarseth, Lin, Papaloizou (1988), *On the collapse and violent relaxation of protoglobular clusters*, The Astrophysical Journal. [doi](https://doi.org/10.1086/165895)
- Akbarian and Choptuik (2014), *Critical collapse in the spherically symmetric Einstein-Vlasov model*, Physical Review D. [doi](https://doi.org/10.1103/physrevd.90.104023)
- Andréasson and Rein (2006), *A numerical investigation of the stability of steady states and critical phenomena for the spherically symmetric Einstein-Vlasov system*, Classical and Quantum Gravity.
- Arad and Lynden-Bell (2005), *Inconsistency in theories of violent relaxation*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2005.09133.x)
- Arad and Johansson (2005), *A numerical comparison of theories of violent relaxation*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2005.09293.x)
- Barnes, Hut, Goodman (1986), *Dynamical instabilities in spherical stellar systems*, The Astrophysical Journal. [doi](https://doi.org/10.1086/163786)
- Barnes, Lanzel, Williams (2009), *The radial orbit instability in collisionless N-body simulations*, The Astrophysical Journal. [doi](https://doi.org/10.1088/0004-637x/704/1/372)
- Bellovary, Dalcanton, Babul et al. (2008), *The Role of the Radial Orbit Instability in Dark Matter Halo Formation and Structure*, The Astrophysical Journal. [doi](https://doi.org/10.1086/591120)
- Beraldo e Silva, de Siqueira Pedra, Junior et al. (2017), *The Arrow of Time in the Collapse of Collisionless Self-gravitating Systems: Non-validity of the Vlasov–Poisson Equation during Violent Relaxation*, The Astrophysical Journal. [doi](https://doi.org/10.3847/1538-4357/aa876e)
- Beraldo e Silva, de Siqueira Pedra, Valluri et al. (2019), *The Discreteness-driven Relaxation of Collisionless Gravitating Systems: Entropy Evolution in External Potentials, N-dependence, and the Role of Chaos*, The Astrophysical Journal. [doi](https://doi.org/10.3847/1538-4357/aaf397)
- Bertschinger (1985), *Self-similar secondary infall and accretion in an Einstein-de Sitter universe*, The Astrophysical Journal Supplement Series. [doi](https://doi.org/10.1086/191028)
- Boella, Coppa, D’Angola et al. (2017), *Gridless particle technique for the Vlasov–Poisson system in problems with high degree of symmetry*, Computer Physics Communications. [doi](https://doi.org/10.1016/j.cpc.2017.11.004)
- Boily, Athanassoula, Kroupa (2002), *Scaling up tides in numerical models of galaxy and halo formation*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1046/j.1365-8711.2002.05372.x)
- Burger and Zavala (2018), *The nature of core formation in dark matter haloes: adiabatic or impulsive?*, [arXiv:1810.10024](https://arxiv.org/abs/1810.10024)
- Chiba, Ding, Hamilton et al. (2025), *Galactic echoes*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/staf1463)
- Di Cintio and Ciotti (2011), *Relaxation of spherical systems with long-range interactions: a numerical investigation*, [arXiv:1103.5436](https://arxiv.org/abs/1103.5436)
- Colombi, Sousbie, Peirani et al. (2015), *Vlasov versus N-body: the Hénon sphere*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/stv819)
- Dattathri, van den Bosch, Weinberg et al. (2025), *The dipole instability in gravitational N-body systems: a natural explanation for lopsidedness and off-centered nuclei in galaxies*, [arXiv:2505.23905](https://arxiv.org/abs/2505.23905)
- Delliou and Henriksen (2003), *Non-radial motion and the NFW profile*, [arXiv:astro-ph/0307046](https://arxiv.org/abs/astro-ph/0307046)
- Doremus, Feix, Baumann (1973), *Stability of a self gravitating system with phase space density function of energy and angular momentum*, Astronomy and Astrophysics 29, 401.
- Errani and Navarro (2020), *The asymptotic tidal remnants of cold dark matter subhalos*, [arXiv:2011.07077](https://arxiv.org/abs/2011.07077)
- Errani, Walker, Rozier et al. (2025), *Impulsive mixing of stellar populations in dwarf spheroidal galaxies*, [arXiv:2502.19475](https://arxiv.org/abs/2502.19475)
- Fillmore and Goldreich (1984), *Self-similar gravitational collapse in an expanding universe*, The Astrophysical Journal. [doi](https://doi.org/10.1086/162070)
- Fouvry (2022), *Kinetic theory of one-dimensional inhomogeneous long-range interacting N-body systems at order 1/N² without collective effects*, [arXiv:2207.05349](https://arxiv.org/abs/2207.05349)
- Fouvry and Roule (2023), *Kinetic blockings in long-range interacting inhomogeneous systems*, [arXiv:2306.04613](https://arxiv.org/abs/2306.04613)
- Fujiwara (1984), *Integration of the Collisionless Boltzmann Equation for Spherical Stellar Systems*, Publications of the Astronomical Society of Japan. [doi](https://doi.org/10.1093/pasj/35.4.547)
- Gitomer and Adam (1976), *Multibeam instability in a Maxwellian simulation plasma*, Physics of Fluids 19, 719. [journal](https://pubs.aip.org/aip/pfl/article-abstract/19/5/719/837379/Multibeam-instability-in-a-Maxwellian-simulation)
- Gondolo and Silk (1999), *Dark Matter Annihilation at the Galactic Center*, Physical Review Letters. [doi](https://doi.org/10.1103/physrevlett.83.1719)
- Guo and Rein (2007), *A Non-Variational Approach to Nonlinear Stability in Stellar Dynamics Applied to the King Model*, Communications in Mathematical Physics. [doi](https://doi.org/10.1007/s00220-007-0212-8)
- Halle, Colombi, Peirani (2019), *Phase-space structure analysis of self-gravitating collisionless spherical systems*, Astronomy and Astrophysics 621, A8, [arXiv:1701.01384](https://arxiv.org/abs/1701.01384)
- Hamilton and Heinemann (2020), *Noise and waves: a unified kinetic theory for stellar systems*, [arXiv:2011.14812](https://arxiv.org/abs/2011.14812)
- Hamilton and Heinemann (2023), *The linear response of stellar systems does not diverge at marginal stability*, [arXiv:2304.07275](https://arxiv.org/abs/2304.07275)
- Henriksen and Widrow (1997), *Self-Similar Relaxation of Self-Gravitating Collisionless Particles*, Physical Review Letters. [doi](https://doi.org/10.1103/physrevlett.78.3426)
- Henriksen and Widrow (1999), *Relaxing and virializing a dark matter halo*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1046/j.1365-8711.1999.02124.x)
- Hozumi, Fujiwara, Kan‐ya (1996), *Growth of Velocity Dispersions for Collapsing Spherical Stellar Systems*, Publications of the Astronomical Society of Japan. [doi](https://doi.org/10.1093/pasj/48.3.503)
- Hozumi, Burkert, Fujiwara (2000), *The origin and formation of cuspy density profiles through violent relaxation of stellar systems*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1046/j.1365-8711.2000.03058.x)
- Huss, Jain, Steinmetz (1999), *How universal are the density profiles of dark halos?*, The Astrophysical Journal.
- Hénon (1964), *L'évolution initiale d'un amas sphérique*, Annales d'Astrophysique.
- Hénon (1973), *Numerical experiments on the stability of spherical stellar systems*, Astronomy and Astrophysics 24, 229.
- Joyce, Marcos, Sylos Labini (2009), *Energy ejection in the collapse of a cold spherical self-gravitating cloud*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2009.14922.x)
- Joyce and Worrakitpoonpon (2010), *Relaxation to thermal equilibrium in the self-gravitating sheet model*, Journal of Statistical Mechanics Theory and Experiment. [doi](https://doi.org/10.1088/1742-5468/2010/10/p10012)
- Kandrup, Vass, Sideris (2002), *Transient chaos and resonant phase mixing in violent relaxation*, [arXiv:astro-ph/0211056](https://arxiv.org/abs/astro-ph/0211056)
- Klinko and Miller (2000), *Mean field theory of spherical gravitating systems*, Physical Review E. [doi](https://doi.org/10.1103/physreve.62.5783)
- Klinko, Miller, Prokhorenkov (2001), *Rotation-induced phase transition in a spherical gravitating system*, Physical Review E. [doi](https://doi.org/10.1103/physreve.63.066131)
- Klinko and Miller (2002), *Angular momentum induced phase transition in spherical gravitational systems: N-body simulations*, Physical Review E 65, 056127, [arXiv:astro-ph/0201488](https://arxiv.org/abs/astro-ph/0201488)
- Klinko and Miller (2004), *Dynamical study of a first order gravitational phase transition*, Physics Letters A 333, 187, [record](https://www.osti.gov/etdeweb/biblio/20616541)
- Lau and Binney (2019), *Relaxation of spherical stellar systems*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/stz2567)
- Lemou, Méhats, Raphaël (2011), *Orbital stability of spherical galactic models*, Inventiones mathematicae. [doi](https://doi.org/10.1007/s00222-011-0332-9)
- Levin, Pakter, Rizzato (2008), *Collisionless relaxation in gravitational systems: From violent relaxation to gravothermal collapse*, Physical Review E. [doi](https://doi.org/10.1103/physreve.78.021130)
- Lu, Mo, Katz et al. (2006), *On the origin of cold dark matter halo density profiles*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2006.10270.x)
- Lynden-Bell (1967), *Statistical Mechanics of Violent Relaxation in Stellar Systems*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/136.1.101)
- MacMillan, Widrow, Henriksen (2006), *On Universal Halos and the Radial Orbit Instability*, The Astrophysical Journal. [doi](https://doi.org/10.1086/508602)
- Maréchal and Perez (2011), *Radial Orbit Instability: Review and Perspectives*, Transport Theory and Statistical Physics. [doi](https://doi.org/10.1080/00411450.2011.654750)
- Merrall and Henriksen (2003), *Relaxation of a Collisionless System and the Transition to a New Equilibrium Velocity Distribution*, The Astrophysical Journal. [doi](https://doi.org/10.1086/377249)
- Merritt and Aguilar (1985), *A numerical study of the stability of spherical galaxies*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/217.4.787)
- Miller and Youngkins (1998), *Phase Transition in a Model Gravitating System*, Physical Review Letters. [doi](https://doi.org/10.1103/physrevlett.81.4794)
- Nusser (2001), *Self-similar spherical collapse with non-radial motions*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1046/j.1365-8711.2001.04527.x)
- O'Neil (1965), *Collisionless damping of nonlinear plasma oscillations*, Physics of Fluids 8, 2255. [doi](https://doi.org/10.1063/1.1761193)
- Olabarrieta and Choptuik (2001), *Critical phenomena at the threshold of black hole formation for collisionless matter in spherical symmetry*, Physical Review D. [doi](https://doi.org/10.1103/physrevd.65.024007)
- Polyachenko and Shukhman (2015), *On the nature of the radial orbit instability in spherically symmetric collisionless stellar systems*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/stv844)
- Polyachenko and Shukhman (2023), *Scale-invariant mode in collisionless spherical stellar systems*, [arXiv:2311.05551](https://arxiv.org/abs/2311.05551)
- Polyachenko and Shukhman (2026), *Two sets of potential-density basis pairs for the study of radial perturbations in collisionless spherical stellar systems*, [arXiv:2609.04012](https://arxiv.org/abs/2609.04012)
- Pontzen and Governato (2012), *How supernova feedback turns dark matter cusps into cores*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2012.20571.x)
- Pontzen and Governato (2013), *Conserved actions, maximum entropy and dark matter haloes*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/sts529)
- Quinlan, Hernquist, Sigurðsson (1995), *Models of Galaxies with Central Black Holes: Adiabatic Growth in Spherical Galaxies*, The Astrophysical Journal. [doi](https://doi.org/10.1086/175295)
- Rasio, Shapiro, Teukolsky (1989), *Solving the Vlasov equation in general relativity*, The Astrophysical Journal. [doi](https://doi.org/10.1086/167785)
- Rein, Rendall, Schaeffer (1998), *Critical collapse of collisionless matter: A numerical investigation*, Physical Review D. [doi](https://doi.org/10.1103/physrevd.58.044007)
- Rein and Rodewis (2003), *Convergence of a particle-in-cell scheme for the spherically symmetric Vlasov-Einstein system*, Indiana University Mathematics Journal. [doi](https://doi.org/10.1512/iumj.2003.52.2363)
- Rioseco and Sarbach (2020), *Phase space mixing in an external gravitational central potential*, Classical and Quantum Gravity. [doi](https://doi.org/10.1088/1361-6382/ababb3)
- Roule, Fouvry, Pichon et al. (2022), *Long-term relaxation of one-dimensional self-gravitating systems*, Physical Review E. [doi](https://doi.org/10.1103/physreve.106.044118)
- Roy and Perez (2003), *Dissipationless collapse of a set of N massive particles*, [arXiv:astro-ph/0310871](https://arxiv.org/abs/astro-ph/0310871)
- Rozier and Errani (2024), *Collisionless Relaxation from Near-equilibrium Configurations: Linear Theory and Application to Tidal Stripping*, The Astrophysical Journal. [doi](https://doi.org/10.3847/1538-4357/ad4c6e)
- Schaeffer (1987), *Discrete approximation of the Poisson-Vlasov system*, Quarterly of Applied Mathematics.
- Sellwood and McGaugh (2005), *The Compression of Dark Matter Halos by Baryonic Infall*, [arXiv:astro-ph/0507589](https://arxiv.org/abs/astro-ph/0507589)
- Sellwood (2015), *Relaxation in N-body simulations of spherical systems*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/stv1846)
- Shapiro and Teukolsky (1985), *Relativistic stellar dynamics on the computer. I - Motivation and numerical method*, The Astrophysical Journal. [doi](https://doi.org/10.1086/163587)
- Sikivie, Tkachev, Wang (1997), *Secondary infall model of galactic halo formation and the spectrum of cold dark matter particles on Earth*, Physical Review D. [doi](https://doi.org/10.1103/physrevd.56.1863)
- Sylos Labini (2012), *Violent and mild relaxation of an isolated self-gravitating uniform and spherical cloud of particles*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2012.21019.x)
- Sylos Labini and Capuzzo-Dolcetta (2020), *Properties of self-gravitating quasi-stationary states*, [arXiv:2009.11624](https://arxiv.org/abs/2009.11624)
- Teles, Levin, Pakter (2011), *Statistical Mechanics of 1d Self-Gravitating Systems: The Core-Halo Distribution*, [arXiv:1109.3780](https://arxiv.org/abs/1109.3780)
- Tremaine, Henon, Lynden-Bell (1986), *H-functions and mixing in violent relaxation*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1093/mnras/219.2.285)
- Trenti, Bertin, van Albada (2004), *A family of models of partially relaxed stellar systems II. Comparison with the products of collisionless collapse*, University of Groningen research database (University of Groningen / Centre for Information Technology). [doi](https://doi.org/10.48550/arxiv.astro-ph/0411541)
- Ullio, Zhao, Kamionkowski (2001), *A Dark-Matter Spike at the Galactic Center?*, [arXiv:astro-ph/0101481](https://arxiv.org/abs/astro-ph/0101481)
- van Albada (1982), *Dissipationless galaxy formation and the r^{1/4} law*, Monthly Notices of the Royal Astronomical Society 201, 939. [doi](https://doi.org/10.1093/mnras/201.4.939)
- Vogelsberger, Mohayaee, White (2011), *Non-spherical similarity solutions for dark halo formation*, Monthly Notices of the Royal Astronomical Society. [doi](https://doi.org/10.1111/j.1365-2966.2011.18605.x)
- Weinberg (1997), *Fluctuations in finite N equilibrium stellar systems*, [arXiv:astro-ph/9707206](https://arxiv.org/abs/astro-ph/9707206)
- White and Zaritsky (1992), *Models for Galaxy halos in an open universe*, The Astrophysical Journal. [doi](https://doi.org/10.1086/171552)
- Williams and Hjorth (2010), *Statistical mechanics of collisionless orbits. ii. structure of halos*, The Astrophysical Journal. [doi](https://doi.org/10.1088/0004-637x/722/1/856)
- Young (1980), *Numerical models of star clusters with a central black hole. I - Adiabatic models*, The Astrophysical Journal. [doi](https://doi.org/10.1086/158553)
- Youngkins and Miller (2000), *Gravitational phase transitions in a one-dimensional spherical system*, Physical Review E. [doi](https://doi.org/10.1103/physreve.62.4583)
- Zukin and Bertschinger (2010), *Self-similar spherical collapse with tidal torque*, Physical Review D. [doi](https://doi.org/10.1103/physrevd.82.104044)
- Łokas and Hoffman (2000), *Formation of Cuspy Density Profiles: A Generic Feature of Collisionless Gravitational Collapse*, The Astrophysical Journal. [doi](https://doi.org/10.1086/312928)
