# Design of the test of the second regime (loss of the core for massive clouds)

Written on 2026-11-16, after the outcome recorded in `prediction_L_scaling.md` and before any of
the runs proposed here. Nothing of this has been run, and the term of stage 1 has not been
implemented in the N-shell model. The expectations of the section "Expectations" are recorded
with the commit of this file.

## What has to be explained

For clouds with Gamma <= 0.6 half of the core is lost at chi = 0.20 to 0.30. For the two largest
masses the PIC simulations give chi = 0.29 (M_gas/M = 0.07) and 0.36 (0.10), whereas the
epicyclic model with the same particles gives 0.21 and 0.22: the difference between PIC and model
at the half is 12, 27 and 39 percent for M_gas/M = 0.05, 0.07 and 0.10.

In the PIC simulations the ratio B2 between the spread |dOmega/dL| sigma_L and the rotation
g m_c M_gas/rho of the core is 0.21, 0.18 and 0.20 at the half for these three masses. Since
rho/sqrt(J_med) = 0.303, B2 = 0.303 chi/m_c, so "B2 = 0.2" is the rule chi_half = 0.66 m_c, with
m_c the mass of the core for equal angular momenta.

## Three hypotheses

- A. There is no second regime: chi_half = 0.27 +- 0.05 and the value 0.36 is a fluctuation.
- B. Empirical rule: chi_half = 0.66 m_c for massive cores. Since m_c saturates (0.50, 0.54 and
  0.47 for M_gas/M = 0.07, 0.10 and 0.30; its bound is 0.73), chi_half saturates at 0.33 to 0.36.
- C. Cause: the frequency shift by the enclosed mass. A shell with a larger L has a larger
  circular orbit and encloses more gas on average, which raises its frequency by
  kappa_u times that mass (Eq. (UniformShift); kappa_u = 0.1454 in the isochrone potential and
  2/L_0^3 = 0.25 in the Kepler potential, L_0 = 2) and counteracts the decrease of Omega with L.
  The model omits this term (section 8A). Its effect grows with the radial density of the gas at
  the circular orbits, hence with the mass, without saturating.

The existing runs already speak against B as a mechanism: the model contains the rotation
induced by the core, and its halves stay at chi = 0.21 and 0.22 where rule B, applied to the
cores of the model (0.42 and 0.49), would give 0.28 and 0.32. For two shells the term of C has
been checked against the exact equations of motion (`reducido_pruebas.py`, P7): the error of the
relative phase drops from 2-6 percent to 0.4-1.5 percent. A direct measurement of the rotation of
each layer in the existing runs (phase of h_1 of each layer) was tried and is too noisy to decide.

## Stage 1: the model with the enclosed mass (no new PIC simulations)

Add to the N-shell model the term

    d zeta_i/dt += i kappa_u zeta_i sum_j m_j [1/2 + asin(clip(d_ij/|zeta_i - zeta_j|, -1, 1))/pi],
    d_ij = sqrt(kappa_0) [R(L_i) - R(L_j)],

which is the average mass enclosed by shell i. For equal angular momenta it is a common rotation
and changes nothing. Checks before any run: equal L reproduces the present model; the pair test
P7 with the N-shell code; P is conserved; convergence with the time step.

Runs: the 28 runs with distributed L of the seven series in which the excess of the model can be
compared with the PIC up to the half (M_gas/M = 0.03, 0.05, 0.07, 0.10; J05M010, J20M080, KM040),
with the same particles. Cost: about twice the batch `ley` of the model, of the order of 12 hours
of jobs with 32 CPUs.

## Stage 2: a larger mass, predicted in advance

Isochrone, sigma_J = 0.1, M_gas/M = 0.15 (Gamma = 1.72): equal angular momenta and three widths.

| chi | s_L | sigma_L/L_0 | N_L | particles | PIC run |
|---|---|---|---|---|---|
| 0.25 | 0.077 | 2.72 % | 65 | 136 500 | 1.1 h |
| 0.40 | 0.124 | 4.38 % | 115 | 241 500 | 1.9 h |
| 0.55 | 0.170 | 6.01 % | 175 | 367 500 | 2.9 h |

The corrected model of stage 1 is solved first for this series (same particles for chi = 0.25
and 0.40, a lattice with half the rows for 0.55; about 2, 7 and 4 hours), and its half is
recorded before the PIC runs are submitted. The analysis of the largest run needs about 64 GB.

## Expectations

Stage 1, if C is the cause: with the term, the half of the model lies within 12 percent of the
one of the PIC in each of the seven series. In particular it rises from 0.22 to at least 0.31 for
M_gas/M = 0.10 (PIC 0.357) and from 0.21 to at least 0.26 for 0.07 (PIC 0.291), and it changes by
less than 0.02 for J05M010 (0.211; PIC 0.218). If the half of the corrected model stays below 0.26
for M_gas/M = 0.10, C is refuted as the cause.

Stage 2, half of the core for M_gas/M = 0.15:

| hypothesis | chi at the half | sigma_L/L_0 |
|---|---|---|
| A | 0.22 - 0.32 | 2.4 - 3.5 % |
| B | 0.33 - 0.36 (0.66 m_c, m_c = 0.50 - 0.54) | 3.6 - 3.9 % |
| C | the value of the corrected model; rough estimate 0.45 - 0.55 | 4.9 - 6.0 % |

The rough estimate for C uses the difference between PIC and model at the half, which is
3.9 M_gas/M for the two largest masses (0.27 and 0.39), extrapolated to 0.58 at M_gas/M = 0.15:
0.22/(1 - 0.58) = 0.52. It is an extrapolation and not a result of the model.

## What each outcome means for the article

- C holds in both stages: the model with one further term, derived from the same expansion,
  reproduces the loss of the core at all the masses, and the two regimes of section 8D become one
  criterion with an effective spread of the frequencies. This is the strongest outcome.
- C holds in stage 1 but the PIC of stage 2 follows B: the term explains the trend up to
  M_gas/M = 0.10, and something else limits the resistance of more massive cores.
- C fails in stage 1: the growth of chi_half with the mass has another origin (the leading order
  in the eccentricity, or the self-gravity being too strong for the average). Section 8D then
  states the two regimes as an observation of the PIC, with rule B as its description.

## Not included

- M_gas/M = 0.20 (Gamma = 2.29): widths up to 11 percent and more than 500 000 particles; only if
  stage 2 confirms C. Its number of layers has to be recomputed with a finer reference.
- The Kepler potential at large Gamma: B and C differ there by about 15 percent, which is less
  than the scatter among the families.
