# Literature check for the cloud article

Date: 2026-11-15. Purpose: place the epicyclic model (Section 7), the absence of a
threshold in mu (Section 6) and the persistent core before the article claims anything
as new. All references below were checked in Crossref or arXiv on this date. Statements
marked "abstract only" were checked against the abstract, not the full text.

## 1. Averaging over a common fast phase and the point-vortex form

- **Tep, Fouvry and Pichon 2026**, A&A 714, A177, arXiv:2603.11238 ("Very long-term
  relaxation of harmonic 1D self-gravitating systems"). N-body relaxation of 1D systems;
  harmonic cores relax on N^2 t_dyn. In the perspectives of v2 they propose to
  time-average the Hamiltonian of a harmonic core over the common orbital period,
  "following an approach similar to Zelnikov & Kuskov (2016)", and state in words that
  after a canonical transformation "the resulting Hamiltonian closely resembles that of
  two-dimensional point vortices". No equations are given. Read in the text (v2, through
  a fetch tool).
- **Zelnikov and Kuskov 2016**, MNRAS 455, 3597, doi:10.1093/mnras/stv2389
  ("Super-Chandrasekhar dynamical friction in a constant-density spherical star system").
  Time averaging in a harmonic core. Not read.
- **Touma, Tremaine and Kazandjian 2009**, MNRAS 394, 1085,
  doi:10.1111/j.1365-2966.2009.14409.x (Gauss's method, softened; N-ring code).
- **Sridhar and Touma 2016**, MNRAS 458, 4129, doi:10.1093/mnras/stw542 (orbit-averaged
  collisionless Boltzmann equation over the fast Kepler phase). Abstract only.

Reading for the article: averaging over a phase which is common to all orbits is the
basis of secular dynamics, and the point-vortex form has been noted in words for harmonic
1D cores. We found no precedent for (i) the explicit averaged equations for shells near a
circular orbit of an anharmonic potential, with the differential rotation Omega(J) - kappa_0,
(ii) the single parameter Gamma proportional to mu^2, (iii) the closed-form omega_c, J_res and
omega_b, (iv) a parameter-free comparison with full Vlasov-Poisson runs, and (v) the
extension to distributed L. The article should claim these, and not the averaging or
the vortex analogy as such.

## 2. Slow modes of near-Keplerian discs (the closest analogue of the core)

- **Tremaine 2001**, AJ 121, 1776, doi:10.1086/319398. Slow m = 1 modes of low-mass discs
  around a massive body; all stable; g-modes only in narrow rings or when an external
  potential dominates the precession; a finite number (sometimes zero) of discrete modes
  and a continuum of singular modes. Abstract only.
- **Jalali and Tremaine 2012**, MNRAS 421, 2368, doi:10.1111/j.1365-2966.2012.20469.x.
  Slow modes exist for arbitrarily small disc mass as long as the self-gravity of the
  disc is the dominant source of apsidal precession; hybrid modes from resonant
  interaction with singular van Kampen modes. Abstract only.
- **Bannikova et al. 2026**, arXiv:2604.11528. N-body: a global slow m = 1 mode forms
  spontaneously in a self-gravitating torus around a dominant mass and persists over
  secular times. Abstract only; probably not needed.

Analogy: in the zeta-plane the core is a concentration which rotates rigidly at
omega_c - kappa_0, slowly compared with kappa_0. Omega(J) - kappa_0 plays the role of the
apsidal precession due to the external potential, and mu^2 compares the rotation induced
by self-gravity with it, which is the comparison that decides whether slow modes exist.

## 3. Threshold of nonlinear Landau damping

- **Lancellotti and Dorning 1998**, PRL 81, 5137, doi:10.1103/PhysRevLett.81.5137
  ("Critical initial states in collisionless plasmas"). Not read.
- **Ivanov and Cairns 2006** (already cited as aIiC06).
- **Klimas, Vinas and Araneda 2017**, J. Plasma Phys. 83, 905830405,
  doi:10.1017/S002237781700054X. Transition from persisting to arrested damping;
  near the transition strength trapping plays no role; no phase-space evidence of a
  second-order transition. Abstract only.

## 4. Persistent clusters

- **Barre and Yamaguchi 2009**, PRE 79, 036208, doi:10.1103/PhysRevE.79.036208. Small
  traveling clusters in HMF models persist (quasi-asymptotic periodic states); the
  parameter region is much smaller in the attractive case. Abstract only.
- Already cited: Levin et al. 2014 (core-halo), Joyce and Worrakitpoonpon 2011 (sheet
  model), Darling and Widrow 2019, Widrow 2023, Asano and Antoja 2026 (phase spirals),
  Schecter et al. 2000, Schecter and Dubin 2001 (vortices).

## Proposed changes (after Erik's review of Sections 6 and 7)

1. Section 7, first paragraph, after "We refer to the resulting model as the epicyclic
   model": "Averaging over a phase which is common to all the orbits is the basis of the
   secular dynamics of nearly Keplerian systems~\cite{jTsTmK09,sSjT16}, and for
   one-dimensional systems in a harmonic potential it has been noted that the averaged
   Hamiltonian resembles the one of point vortices~\cite{kTjFcP26}."
2. Section 7, paragraph on two-dimensional fluids: add the analogy with the slow modes of
   near-Keplerian discs~\cite{sT01,mJsT12}, which exist when the self-gravity of the disc
   dominates the precession of the orbits.
3. Section 6, line 116: cite Lancellotti and Dorning 1998 together with Ivanov and Cairns
   2006 for the threshold of Langmuir waves.
4. Conclusions: state what is new as in Section 1 above.
