# Prediction for the second batch of initial profiles

Recorded on 2026-11-15, before the runs of this batch were submitted (the commit of this file
precedes them). Scripts: `reproducir/scripts/formas.py` and `reducido.py` (series `formas2_*`),
and `VlasovPoisson_PIC_sp/reproducir/scripts/nube_L.py` (batch `formas2`).

## The law under test

The first batch (`formas.py colapso`) compared the profiles J^2, J^4 and (J - 0.05)_+^2 times
exp(-J^2/s_J^2). With

    Gamma = g M_gas / (|Omega_0'| J_med^(3/2)),

J_med the median action of the cloud, the epicyclic model gives m_c = 0.1 at Gamma = 0.24-0.28,
m_c = 0.3 at Gamma = 0.52-0.54 and |h_1|/h_0 = 0.1 at Gamma = 0.31-0.33 for the three profiles,
whereas the values of mu at these points differ by 48-62 percent. The nine PIC runs of the
first batch lie on the same curve (m_c = 0.104 at Gamma = 0.229, 0.105 at 0.258, 0.200 at
0.334, 0.215 at 0.343).

Claim: for clouds which are compact in angle (s_Q = 0.1) and whose DF vanishes at the circular
orbit, the core and its oscillation are determined by Gamma, not by mu.

## New profiles

- ring: exp(-(J - 0.15)^2/0.03^2), narrow band (Delta J = 0.129) far from the circular orbit
  (J_med = 0.150). Its mu is large and its Gamma small, so the two parameters predict opposite
  outcomes.
- water bag: constant for 0.05 <= J <= 0.25 (J_med = 0.149), sharp edges.
- the hollow profile in the Kepler potential (model only): Gamma contains g and Omega_0' of the
  potential.
- J^0 at smaller masses (model, and one PIC run): its DF does not vanish at the circular orbit,
  so the transport to the circular orbit, on which the law rests, is not needed.

Criteria: no core if m_c < 1.5 m_mix and |h_1|/h_0 < 0.03; core if m_c > 2 m_mix; persistent
oscillation if |h_1|/h_0 > 0.1 (averages as in `formas.py tabla`). m_mix = 0.094 (ring), 0.025
(water bag), 0.056 (J^0).

## Predictions for the PIC runs

| run | Gamma | mu | prediction from Gamma | expectation from mu (article cloud at the same mu) |
|---|---|---|---|---|
| G_anillo_150 | 0.15 | 1.04 | no core, no oscillation | core with m_c of about 0.5 |
| G_anillo_300 | 0.30 | 1.46 | core with m_c 0.12-0.25 | m_c of about 0.55 |
| G_anillo_450 | 0.45 | 1.80 | core with m_c 0.22-0.38, oscillation | m_c of about 0.55 |
| G_meseta_150 | 0.15 | 0.65 | no core, no oscillation | core with m_c of about 0.4 |
| G_meseta_450 | 0.45 | 1.13 | core with m_c 0.22-0.38, oscillation | m_c of about 0.5 |
| G_p0_150 | 0.15 | 0.19 | outside the law; expected: core (m_c > 0.11) | no core |

For the model series the prediction is that m_c and |h_1| of the ring, the water bag and the
Kepler run follow, as functions of Gamma, the curve of the first batch within about 20 percent.
