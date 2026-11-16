# Prediction for the scaling of the loss of the core with the spread of the angular momenta

Recorded on 2026-11-16, before the runs of this batch were submitted (the commit of this file
precedes them). Scripts: `VlasovPoisson_PIC_sp/reproducir/scripts/nube_L.py` (batches `ley` and
`ley_grande`, steps `frontera` and `figura_chi`) and `reproducir/scripts/reducido.py` (table `LEY`,
series `ley_*`).

## The law under test

Section 8 of the article shows that a spread sigma_L of the angular momenta removes the dense core
which the cloud of Eq. (Pulse) leaves at the circular orbit, and that the spread which is needed
grows with the mass of the gas. The estimate given there compares the spread
|dOmega/dL| sigma_L of the frequencies of the layers with different L with the rotation
g m_c M_gas/rho which the core induces at its edge. At the spread at which the excess
m_c - m_mix of the mass of the core has fallen to one half, this ratio is 0.60, 0.37 and 0.21 for
M_gas/M = 0.02, 0.03 and 0.05: it does not give the dependence on the mass.

The same spread, compared with the rotation which the whole cloud induces at its distance
sqrt(J_med) from the circular orbit, defines

    chi := |dOmega/dL| sigma_L sqrt(J_med) / (g M_gas) = c'(L_0) sigma_L / (Gamma J_med),

with |dOmega/dL| evaluated at the circular orbit (J = 0, L = L_0), sigma_L = s_L/sqrt(2) the root
mean square width of the factor exp(-(L - L_0)^2/s_L^2), J_med the median of the action variable of
the cloud (0.1087 for sigma_J = 0.1, on the lattice of 400 x 100 points), g = G/(pi R_0^2
sqrt(kappa_0)) and Gamma = g M_gas/(|Omega_0'| J_med^(3/2)) as in Eq. (Gamma). Since
Omega = (J + c(L))^-3, dOmega/dL = c'(L) dOmega/dJ, which gives the second form. chi is the time
sqrt(J_med)/(g M_gas) in which the gas is transported to the circular orbit, divided by the time in
which the layers dephase, in the same way as 1/Gamma is the same transport time divided by the
time 1/(|Omega_0'| J_med) of the differential rotation in J.

With the existing runs (`nube_L.py frontera`, isochrone, sigma_J = 0.1, L_0 = 2) the excess of the
core mass falls to one half at

| M_gas/M | Gamma | sigma_L/L_0 at one half | ratio with the rotation of the core | chi (PIC) | chi (model) |
|---|---|---|---|---|---|
| 0.02 | 0.229 | 0.30 % | 0.60 | 0.205 | 0.181 |
| 0.03 | 0.343 | 0.58 % | 0.37 | 0.268 | 0.225 |
| 0.05 | 0.572 | 1.00 % | 0.21 | 0.276 | 0.242 |

In these three series only the mass varies. The dependences of chi on J_med and on the constants
of the potential (g, |dOmega/dL|) have not been tested.

## Hypotheses

- (S) Similarity. For clouds of the shape of Eq. (Pulse) with sigma_Q = 0.1, the mass of the core
  relative to its value for equal angular momenta depends on the mass of the gas, on sigma_J, on
  sigma_L and on the external potential only through Gamma and chi.
- (L) Linearity. For Gamma >= 0.3 this relative excess is nearly a function of chi alone: it falls
  to one half at chi = 0.27 +- 0.05, such that the spread which removes the core grows linearly
  with the mass of the gas.

## New runs

Five families, each with equal angular momenta and four spreads with chi = 0.15, 0.25, 0.35 and
0.5 (s_L rounded). Lattice of 200 x 50 points in (Q, J) as in Section 8 (2100 particles per layer),
N_L the smallest number which satisfies the bound 7e-4 of `nube_L.py reticulas`, runs up to
30 tau_1. The epicyclic model is solved with the same particles (series `ley_*`); its results
are not known when this file is written.

| family | potential | sigma_J | M_gas/M | Gamma | tau_1 | s_L (chi; N_L) |
|---|---|---|---|---|---|---|
| M070 | isochrone | 0.1 | 0.07 | 0.800 | 327 | 0.022 (0.153; 20), 0.036 (0.250; 28), 0.050 (0.347; 40), 0.072 (0.499; 60) |
| M100 | isochrone | 0.1 | 0.10 | 1.143 | 327 | 0.031 (0.150; 24), 0.052 (0.252; 40), 0.072 (0.349; 60), 0.103 (0.500; 90) |
| J05M010 | isochrone | 0.05 | 0.01 | 0.323 | 589 | 0.0044 (0.151; 10), 0.0073 (0.251; 14), 0.0102 (0.350; 16), 0.0146 (0.501; 20) |
| J20M080 | isochrone | 0.2 | 0.08 | 0.323 | 198 | 0.0175 (0.150; 12), 0.029 (0.249; 16), 0.041 (0.352; 20), 0.058 (0.498; 28) |
| KM040 | Kepler | 0.1 | 0.04 | 0.335 | 161 | 0.0077 (0.150; 12), 0.0129 (0.251; 14), 0.018 (0.350; 18), 0.0257 (0.499; 26) |

The last three families have the value of Gamma of the series with M_gas/M = 0.03 (0.343), but
another J_med or another potential. With equal angular momenta they are the clouds J05_0.01,
J20_0.08 and K_0.04 of Section 6 (m_c = 0.186, 0.309 and 0.234 on their lattices of 400 x 25 and
800 x 25 points). Before the batch, the code with distributed L was checked against those three
runs with equal L up to t = 2000 (`nube_L.py consistencia CJ05 CJ20 CK04`): |h_1|(t) agrees to
within 9e-5 of its maximum in each interval of time, and the amplitude of the potential and mu coincide to all the digits
printed.

## Predictions for the PIC runs

Primary. In each family the excess m_c - m_mix falls to one half of its value for equal angular
momenta (linear interpolation in s_L, as in `nube_L.py frontera`) at chi between 0.22 and 0.32,
that is, at

| family | sigma_L/L_0 at one half, predicted | if it fell at 0.58 % as for M_gas/M = 0.03 |
|---|---|---|
| M070 | 1.12 - 1.63 % (1.38 % for chi = 0.27) | chi = 0.11 |
| M100 | 1.60 - 2.33 % (1.97 %) | chi = 0.08 |
| J05M010 | 0.23 - 0.33 % (0.28 %) | chi = 0.56 |
| J20M080 | 0.91 - 1.32 % (1.11 %) | chi = 0.14 |
| KM040 | 0.40 - 0.58 % (0.49 %) | chi = 0.32 |

Secondary. At the four spreads the relative excess is 0.76, 0.55, 0.35 and 0.17, each within
+- 0.15. These are the values of the series with M_gas/M = 0.03 and 0.05 at chi = 0.15, 0.25, 0.35
and 0.5 (0.80 and 0.72; 0.55 and 0.54; 0.33 and 0.37; 0.16 and 0.18).

For the model the expectation is the one of the existing series: the same dependence on chi, with
the half 10 to 20 percent below the one of the PIC.

## What would refute it

- (S) is refuted if the half of J05M010 or J20M080 lies outside 0.22 <= chi <= 0.32. These two
  families are the decisive ones: if the spread which removes the core did not scale with J_med,
  they would give chi = 0.56 and 0.14. KM040 tests the constants of the potential, but it does
  not discriminate against that alternative (0.32).
- (L) is refuted if the half of M070 or M100 lies outside the same interval. The estimate of the
  article, with its ratio fixed at the value 0.21 of M_gas/M = 0.05 and the masses of the core
  0.479 and 0.516 of the runs with equal L, gives chi = 0.33 and 0.36 for these two families.

## Caveats

- For M_gas/M = 0.07 and 0.10 the core oscillates with an amplitude which is larger than its size.
  m_c is the mass inside a single disc for the whole gas, so it decreases both when the layers
  lose their cores and when they keep them but are out of phase. From this batch on the mass of
  the core of each layer is measured as well (`nube_L.py medir`, nucleo_capas.npz).
- chi neglects the softening of the interaction by the spread of the radii of the circular orbits.
  Its parameter sqrt(kappa_0) R'(L_0) sigma_L/sqrt(2 J_med) reaches 0.17 in this batch (M100,
  chi = 0.5); for M_gas/M = 0.03 it was shown to be irrelevant (Section 8C).
- The model is of leading order in the eccentricity of the orbits, and for sigma_J = 0.2 the ratio
  between the amplitude of the radial oscillation and R_0 is 0.41.
- The coherence between the layers depends on N_L and is not part of the prediction.

## Outcome

Added on 2026-11-16, after the runs (PIC: jobs 228533 to 228557; model: 228558 to 228562; all
finished, largest relative change of the energy 3.7e-5). Numbers from `nube_L.py frontera` and
`figura_chi`. The text above this section was not changed.

| family | Gamma | sigma_L/L_0 at one half | chi (PIC) | chi (model) | primary prediction (0.22 - 0.32) |
|---|---|---|---|---|---|
| M070 | 0.800 | 1.48 % | 0.291 | 0.211 | inside |
| M100 | 1.143 | 2.60 % | 0.357 | 0.217 | outside |
| J05M010 | 0.323 | 0.22 % | 0.218 | 0.211 | 0.002 below the lower end |
| J20M080 | 0.323 | 1.25 % | 0.303 | 0.239 | inside |
| KM040 | 0.335 | 0.49 % | 0.270 | 0.216 | inside |

- Primary: three of five inside, one at the lower end, one outside.
- Secondary (relative excess at chi = 0.15, 0.25, 0.35, 0.5; predicted 0.76, 0.55, 0.35, 0.17, each
  +- 0.15): M070 0.79, 0.57, 0.40, 0.20 (all inside); M100 0.86, 0.68, 0.51, 0.23 (one outside);
  J05M010 0.70, 0.40, 0.19, 0.01 (two outside); J20M080 0.84, 0.62, 0.39, 0.19 (all inside); KM040
  0.82, 0.54, 0.33, 0.10 (all inside).
- (S), the scaling with J_med and with the potential, holds against the alternative stated above:
  the families with sigma_J = 0.05 and 0.2 give chi = 0.218 and 0.303, not 0.56 and 0.14. It is not
  exact: at Gamma of about 0.33 the half lies at chi = 0.218, 0.268 and 0.303 for sigma_J = 0.05, 0.1
  and 0.2, and at 0.270 in the Kepler potential.
- (L), the linear growth with the mass up to Gamma = 1.1, is refuted by M100. The estimate of the
  article with its ratio fixed at 0.21, which was given above as the alternative (chi = 0.33 and
  0.36), is the one which holds for the largest mass (0.291 and 0.357 measured): the ratio with
  the rotation of the core at the half is 0.21, 0.18 and 0.20 for M_gas/M = 0.05, 0.07 and 0.10.
- Model: its halves lie at chi = 0.21 to 0.24 in all the families, 3, 21, 20, 27 and 39 percent
  below those of the PIC (J05M010, J20M080, KM040, M070, M100); the expectation was 10 to 20 percent.
- The mass of the core of each layer differs from m_c by -2 to +6 percent for chi <= 0.25, by -1 to
  +23 percent for chi = 0.35 and by +20 to +51 percent for chi = 0.5.
- With equal angular momenta, on the lattice of 200 x 50 points, m_c = 0.180, 0.274, 0.215, 0.501
  and 0.541 for J05M010, J20M080, KM040, M070 and M100 (0.186, 0.309, 0.234, 0.479 and 0.516 in the
  runs of Section 6 on their lattices).
