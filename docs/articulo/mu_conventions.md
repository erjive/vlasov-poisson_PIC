# The definition of μ: three conventions compared

22 November 2026. Follows observations 3.1 to 3.3 and 6.2 of review_section_V.md. No file of the manuscript was
changed. Numbers from `reproducir/scripts/mu_convenciones.py` (steps `convergencia`, `libre`, `tabla`, `colapso`;
raw output in exe/demo_eta/mu/mu_convenciones.txt, not versioned).

## Recommendation

Use the **passive** convention everywhere: A is the amplitude of the potential of the cloud moving in the external
potential alone, ΔΩ and ΔJ the band of Eq. (45). It is fixed by the initial datum, it applies to every shape, width
and potential, it is what Section VII.E and Table VIII already use, and it collapses the families of Section V better
than the other two. With it, μ = c (M_gas/M)^(1/2) exactly, with c = 3.094 for the reference cloud.

**Update of 23 November 2026.** The scripts now use the passive μ with A = (G M_gas/2) Δ⟨1/r⟩, without the window 3 ≤ r ≤ 15.
Nothing changes for the isochrone clouds, whose smallest pericenter is 3.4. In the Kepler potential c goes from 2.664 to
2.688 and μ grows by 0.9 percent: M_gas/M = 0.01, 0.02, 0.04 and 0.08 give 0.269, 0.380, 0.538 and 0.760. The Kepler values
in the tables below keep the earlier window.

## The three conventions

All three use μ² = A/(ΔΩ ΔJ) with the band at one percent of the factor in J. They differ in A and ΔΩ.

| Convention | A | ΔΩ | Where it is used now |
| --- | --- | --- | --- |
| Run | Measured in the self-gravitating run, first radial period, PIC snapshots | Band of the external potential | V (Tables I to III, Fig. 5), VI, VIII, Table VII |
| Passive | Cloud moving in the external potential alone, exact mixing and exact shell potential | Band of the external potential | VII.E, Table VIII, Fig. 11 |
| Leading order | π g M_gas (2σ_J)^(1/2), Eq. (58) | \|Ω₀′\| ΔJ | Eq. (59): μ² = 0.69 Γ |

| Criterion | Run | Passive | Leading order |
| --- | --- | --- | --- |
| Fixed by the initial datum | No: includes the response of the cloud in its first period | Yes | Yes |
| Defined for any shape, σ_Q and potential | Yes | Yes | No: only J² exp(−J²/σ_J²), and it ignores σ_Q |
| μ² proportional to M_gas | Only for M_gas/M ≲ 0.025 | Exactly | Exactly |
| Families of Section V with a core: deviation of m_c from the reference family at equal μ (5 runs) | Mean 11.9 %, from −12 to +28 % | Mean 8.6 %, from −12 to +19 % | Mean 17.6 %, from −14 to +45 % |
| Change of the published values | None | Reference cloud: +1 to +3 % for M_gas/M ≤ 0.03, +10 % at 0.1, +47 % at 0.3 | −14 to +33 % over all the runs (−8 % at M_gas/M = 0.02) |

The leading-order convention is proportional to Γ^(1/2), so its row on the collapse of the families also holds for Γ.

## Evidence

**The passive amplitude is converged.** For the reference cloud, a = A/M_gas = 0.049673 with 400 × 100 points and
time step 0.25. Doubling N_J, doubling N_Q, the step 0.1 or twice as many radii change the fifth digit at most. With
snapshots every 10 time units, as in the PIC output, a = 0.04967.

**The measured amplitude falls with the mass because of self-gravity.** The ratio of the amplitude measured in the
runs to the passive one, and the corresponding ratio of μ:

| M_gas/M | 0.001 | 0.01 | 0.02 | 0.025 | 0.03 | 0.05 | 0.1 | 0.3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_run/A_passive | 0.997 | 0.982 | 0.964 | 0.956 | 0.947 | 0.910 | 0.820 | 0.465 |
| μ_run/μ_passive | 0.999 | 0.991 | 0.982 | 0.978 | 0.973 | 0.954 | 0.906 | 0.682 |

At M_gas/M = 10⁻³ the run reproduces the passive amplitude within 0.3 percent; 0.12 percent of it comes from the
lattice with N_Q = 25. The rest of the deficit grows with the mass, so the run convention mixes part of the outcome
into μ.

**μ² and Γ.** In the passive convention μ²/Γ is a constant of each family:

| Family | Reference | Kepler | σ_J = 0.05 | σ_J = 0.2 | σ_Q = 0.2 | σ_Q = 0.4 | J⁴ | (J−0.05)₊² | ring | water bag | J⁰ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| μ²/Γ | 0.838 | 0.848 | 0.773 | 0.962 | 0.791 | 0.642 | 1.430 | 1.843 | 7.125 | 2.871 | 0.249 |

Eq. (59) gives 0.69 for the J² shape at leading order. The spread between σ_J = 0.05 and 0.2 (0.77 to 0.96) is the
finite width of the band, which the passive convention keeps and Γ does not.

## Recalculated values

Passive c (μ = c (M_gas/M)^(1/2)): reference 3.094; σ_Q = 0.2: 3.007; σ_Q = 0.4: 2.710; σ_J = 0.05: 5.001;
σ_J = 0.2: 1.972; Kepler: 2.664; J⁴: 3.218; (J−0.05)₊²: 3.836; ring: 7.092; water bag: 4.502; J⁰: 3.134.

Section V, reference family (σ_Q = σ_J = 0.1, isochrone):

| M_gas/M | Γ | μ run (now) | μ passive | μ leading order |
| --- | --- | --- | --- | --- |
| 0.001 | 0.011 | 0.098 | 0.098 | 0.089 |
| 0.01 | 0.114 | 0.307 | 0.309 | 0.281 |
| 0.015 | 0.171 | 0.374 | 0.379 | 0.344 |
| 0.02 | 0.229 | 0.430 | 0.438 | 0.397 |
| 0.021 | 0.240 | 0.440 | 0.448 | 0.407 |
| 0.022 | 0.252 | 0.450 | 0.459 | 0.417 |
| 0.023 | 0.263 | 0.460 | 0.469 | 0.426 |
| 0.024 | 0.274 | 0.470 | 0.479 | 0.435 |
| 0.025 | 0.286 | 0.478 | 0.489 | 0.444 |
| 0.03 | 0.343 | 0.522 | 0.536 | 0.486 |
| 0.04 | 0.457 | 0.597 | 0.619 | 0.562 |
| 0.05 | 0.572 | 0.660 | 0.692 | 0.628 |
| 0.07 | 0.800 | 0.764 | 0.819 | 0.743 |
| 0.1 | 1.143 | 0.886 | 0.979 | 0.888 |
| 0.3 | 3.430 | 1.156 | 1.695 | 1.538 |

The runs on other lattices (q, j, h, r, t) have the same passive μ; their measured μ differs by less than 0.002.

Section V, other families (Table III and Table VII):

| Family | M_gas/M | Γ | μ run (now) | μ passive | μ leading order |
| --- | --- | --- | --- | --- | --- |
| σ_Q = 0.2 | 0.01 | 0.114 | 0.298 | 0.301 | 0.281 |
| | 0.02 | 0.229 | 0.418 | 0.425 | 0.397 |
| | 0.04 | 0.457 | 0.578 | 0.601 | 0.562 |
| | 0.08 | 0.915 | 0.781 | 0.851 | 0.794 |
| σ_Q = 0.4 | 0.01 | 0.114 | 0.268 | 0.271 | 0.281 |
| | 0.02 | 0.229 | 0.375 | 0.383 | 0.397 |
| | 0.04 | 0.457 | 0.518 | 0.542 | 0.562 |
| | 0.08 | 0.915 | 0.693 | 0.766 | 0.794 |
| σ_J = 0.05 | 0.005 | 0.162 | 0.351 | 0.354 | 0.334 |
| | 0.01 | 0.323 | 0.493 | 0.500 | 0.472 |
| | 0.02 | 0.647 | 0.687 | 0.707 | 0.668 |
| σ_J = 0.2 | 0.02 | 0.081 | 0.274 | 0.279 | 0.236 |
| | 0.04 | 0.162 | 0.384 | 0.394 | 0.334 |
| | 0.08 | 0.323 | 0.529 | 0.558 | 0.472 |
| Kepler | 0.01 | 0.084 | 0.262 | 0.266 | 0.240 |
| | 0.02 | 0.167 | 0.367 | 0.377 | 0.340 |
| | 0.04 | 0.335 | 0.510 | 0.533 | 0.480 |
| | 0.08 | 0.669 | 0.691 | 0.753 | 0.679 |

Section VI (3200 × 100 points; the leading order is defined only for the J² shape):

| Run | M_gas/M | Γ | μ run | μ passive |
| --- | --- | --- | --- | --- |
| J², Γ = 0.025 / 0.05 / 0.1 / 0.15 | 0.00219 / 0.00437 / 0.00875 / 0.01312 | 0.025 / 0.050 / 0.100 / 0.150 | 0.145 / 0.204 / 0.287 / 0.351 | 0.145 / 0.205 / 0.289 / 0.354 |
| J⁴, Γ = 0.05 / 0.1 | 0.00692 / 0.01385 | 0.050 / 0.100 | 0.266 / 0.375 | 0.268 / 0.379 |
| (J−0.05)₊², Γ = 0.05 / 0.1 | 0.0063 / 0.0126 | 0.050 / 0.101 | 0.303 / 0.426 | 0.305 / 0.431 |
| ring, Γ = 0.05 / 0.1 | 0.00714 / 0.01427 | 0.050 / 0.101 | 0.596 / 0.837 | 0.599 / 0.847 |
| water bag, Γ = 0.05 / 0.1 | 0.00703 / 0.01406 | 0.050 / 0.099 | 0.375 / 0.528 | 0.377 / 0.534 |

Section VII.E (Table VIII already uses the passive values, reproduced here to three digits):

| Shape | M_gas/M | Γ | μ run | μ passive |
| --- | --- | --- | --- | --- |
| J² | 0.02 / 0.03 / 0.05 | 0.229 / 0.343 / 0.572 | 0.430 / 0.522 / 0.661 | 0.438 / 0.536 / 0.692 |
| J⁴ | 0.01849 / 0.02774 / 0.04623 | 0.134 / 0.201 / 0.335 | 0.431 / 0.525 / 0.667 | 0.438 / 0.536 / 0.692 |
| (J−0.05)₊² | 0.01301 / 0.01952 / 0.03253 | 0.104 / 0.156 / 0.260 | 0.433 / 0.528 / 0.674 | 0.438 / 0.536 / 0.692 |
| ring | 0.02117 / 0.04234 / 0.06351 | 0.149 / 0.299 / 0.448 | 1.015 / 1.410 / 1.694 | 1.032 / 1.459 / 1.787 |
| water bag | 0.02109 / 0.06327 | 0.149 / 0.447 | 0.643 / 1.076 | 0.654 / 1.132 |
| J⁰ | 0.00376 | 0.148 | 0.191 | 0.192 |

Section VIII (equal angular momenta; the runs with distributed L keep the label of their mass): M_gas/M = 0.01,
0.02, 0.03, 0.05, 0.07 and 0.1 give μ = 0.309, 0.438, 0.536, 0.692, 0.819 and 0.979 (now 0.307, 0.430, 0.522,
0.661, 0.765 and 0.887). The families of the χ test give 0.500 (σ_J = 0.05, M_gas/M = 0.01), 0.558 (σ_J = 0.2,
0.08) and 0.533 (Kepler, 0.04).

## What would change in the text (not applied)

- V.A: A computed for the passive cloud; "μ = 3.05 (M_gas/M)^(1/2) to within 2 percent" becomes μ = 3.09 (M_gas/M)^(1/2)
  exactly. The measured amplitude stays as a remark: it agrees within 5 percent for M_gas/M ≤ 0.025 (4.4 percent at
  0.025).
- V.A: "μ = 0.10, 0.31, 0.89 and 1.16" becomes 0.10, 0.31, 0.98 and 1.70; the masses 0.02 to 0.025 give 0.44 to 0.49.
- V.A: μ²σ_J^(3/2)M/M_gas lies between 0.28 and 0.35 (now 0.26 and 0.34). A cloud four times wider in J needs a mass
  6.4 times larger for the same μ, not about eight.
- V.B to V.E and Section IX: the change lies between μ = 0.44 and 0.49 (now 0.43 and 0.48). All runs with μ ≤ 0.45
  have |h₁|/h₀ ≤ 2.1 × 10⁻²; all with μ ≥ 0.49 and σ_Q ≤ 0.2 have at least 8 × 10⁻² (lattice of Table II). No core
  (m_c ≤ 0.048) for μ ≤ 0.39 in every family.
- Table III brackets: reference 0.44–0.49; σ_J = 0.05: 0.35–0.50; σ_J = 0.2: 0.39–0.56; σ_Q = 0.2: 0.43–0.60;
  σ_Q = 0.4: 0.54–0.77; Kepler: 0.38–0.53. The four brackets with σ_Q ≤ 0.2 still contain 0.44–0.49.
- V.E: m_c of the families with σ_J = 0.05, 0.2 and Kepler differs from the reference family at equal μ by −12 to
  +19 percent (now −12 to +28); the reference family reaches these m_c at values of μ that differ by less than 0.07
  (unchanged). With σ_Q = 0.2, m_c is smaller by 37 to 41 percent; with σ_Q = 0.4 by 66 to 72 percent.
- VI.A: μ between 0.15 and 0.35 (unchanged at two digits); VI.C: Γ = 0.17 is μ = 0.38 (now 0.37).
- VII.B: Eq. (59) stays as the leading order; the values it gives are about 9 percent below the passive ones (now
  "about 8 percent below those of Section V"). The exact passive factor is μ² = 0.84 Γ for the reference cloud.
- VII.D, Table VII: the μ column takes the passive values above.
- VII.E: the sentence that distinguishes the μ of this subsection from that of V.A is no longer needed.
- VIII: M_gas/M = 0.01, 0.02, 0.03 and 0.05 are μ = 0.31, 0.44, 0.54 and 0.69.
- Figures with μ on an axis or in a label: 2, 3, 4, 5, 9, 10, and those of Section VIII.

## Open points

- The collapse comparison rests on five runs with a core in three families; it favors the passive convention but is
  not a strong test.
- Whether Section V should state its results in Γ instead (observation 6.1) is a separate decision. If it does, μ
  keeps the passive definition for the remark that relates it to τ₁ and to trapping.
