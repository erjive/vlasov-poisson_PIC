# What the amplitude A in μ² = A/(ΔΩ ΔJ) means

22 November 2026. Follows mu_conventions.md and mu_precision.md. No file of the manuscript or of the code was changed.
The computations are in scripts of the session scratchpad, which only read the existing modules and runs. They are
not versioned. No PIC simulation was run. G = 1 and amplitudes are per unit mass of the gas, unless stated otherwise.

## Conclusion

**A is a reasonable heuristic for the strength of the perturbation, but it is not directly a trapping scale.** As a
criterion for the core it works only as an indicator within one shape of the DF.

- **Not a direct scale.** A is evaluated inside the cloud, where the potential is uniform in r, exerts no force and
  holds no gas. The part of the perturbation that acts on the gas is the secular Hamiltonian of section 1.
- **A reasonable proxy.** In the 11 clouds tested, the first harmonic of the secular Hamiltonian is 0.19 to 0.24 of A.
  Up to that nearly constant factor, A measures the strength of the perturbation.
- **The form √(|Ω′|A) is right in its scaling and wrong in its prefactor and location.** For a compact shell the
  libration frequency is 0.34 √(|Ω′|A), and the island lies at the phase opposite to the cloud, not on it.
- **Pendulum trapping does not explain the core.** At M_gas/M = 0.02 a libration takes 6.5 τ₁, while the core is visible
  from 3 τ₁. The criterion "island as wide as the band" also fails against the PIC runs of the ring and the water bag.
- **What decides is the denominator.** Replacing A by the secular amplitude changes almost nothing. Replacing the band
  ΔΩ ΔJ by the distance |Ω₀′| J_med² collapses the shapes.

| Statement | Status |
| --- | --- |
| μ = √(s̄ A)/ΔΩ with s̄ = ΔΩ/ΔJ | Identity |
| A = (G M_gas/2) Δ⟨1/r⟩ in the isochrone clouds | Exact (the smallest pericenter, 3.4, lies above r = 3) |
| A uniform Φ → Φ + f(t) does not change the motion | Exact |
| Secular kernel tends to g\|ζ − ζ′\| as J → 0 | Derived; approached numerically |
| Libration about ψ = π at the frequency of the secular Hamiltonian | Verified with test particles, restricted problem, five masses (0.1 to 1.7 %) |
| Saddles at distance ∝ M from a compact shell | Verified (exponent 1.003) |
| Merger of the fixed points at Γ* = 0.544 for a point cloud | Derived; verified as J → 0 (0.540). It is 0.47 at J₀ = 0.1 |
| \|V₁\|/A = 0.19 to 0.24 | Computed at t = 0, first order in M, cloud taken as rigid |
| The island-overlap criterion fails for the ring and the water bag | Tested against PIC runs |
| The core does not form by pendulum trapping | Strong inference from time scales; not demonstrated |
| The merger at Γ* sets the mass of the core | Conjecture, not tested |

## 1. Relation with trapping

**Identity.** With the mean slope s̄ := ΔΩ/ΔJ of the band,

$$ \mu = \frac{\omega_A}{\Delta\Omega}, \qquad \omega_A := \sqrt{\bar{s}\,A} . $$

This says nothing yet about whether ω_A is a frequency of the system.

**Slow Hamiltonian.** A test shell with variables (Q, J) of the external potential has
H = E(J) + Φ_gas(t, r(Q, J)). The source shells have phases Q′ = Q′₀ + Ω(J′)t. When |Ω − Ω′| ≪ κ₀, only the terms that
depend on ψ = Q − Q′ are slow. To first order in M_gas, in the frame that rotates with a cloud of frequency Ω_c,

$$ K(\psi, J) = E(J) - \Omega_c J + M_{gas} V(\psi, J), \qquad V(\psi, J) = \sum_b w_b\, \bar{K}(\psi - Q_b; J, J_b), $$

$$ \bar{K}(\psi; J, J') = -G \left\langle \frac{1}{\max\{ r(\theta + \psi, J),\, r(\theta, J') \}} \right\rangle_\theta . $$

To leading order in the eccentricity K̄ = const + uniform term + g|ζ − ζ′|, as in Eq. (54). For two shells with the same
action, with ρ = (2J)^(1/2),

$$ g|\zeta - \zeta'| = 2g\rho\left|\sin\frac{\psi}{2}\right| = 2g\rho\left[\frac{2}{\pi} - \frac{4}{\pi}\sum_{n\ge 1}\frac{\cos n\psi}{4n^2 - 1}\right] . $$

Since A₅₈ = π g M ρ (Eq. 58), the slow potential is much smaller than A:

| Quantity | Leading order | Exact, shell with J₀ = 0.1 |
| --- | --- | --- |
| Half depth of V in ψ | g M ρ = 0.318 A | 0.283 A |
| First harmonic \|V₁\| | 8 g M ρ/(3π) = 0.270 A | 0.233 A |

The fast change of the potential at fixed r averages out for the gas. Only this part acts over many orbits.

**Pendulum.** Near J_c, with only the first harmonic,

$$ K \simeq -\tfrac{1}{2}|\Omega'|\,\delta J^2 - |V_1|\cos\psi, \qquad \omega_1 = \sqrt{|\Omega'||V_1|} . $$

Since Ω′ < 0, the elliptic point is at **ψ = π**, opposite to the cloud, and ψ = 0 is hyperbolic. This is the same
structure as co-orbital motion in a Keplerian potential. Near ψ = π the exact leading-order potential is
2gMρ cos(x/2), with x = ψ − π, so ω_π² = |Ω′| g M ρ/2 and ω_π = 0.40 √(|Ω′|A₅₈).

**Test with particles.** A shell of mass m and J₀ = 0.1 moves on its unperturbed orbit in the isochrone potential. Test
particles start near ψ = π and are integrated in its exact, time-dependent field. ω_H is the frequency predicted by
the Hessian of K with the exact kernel.

| m | ω_H | Measured | Measured/ω_H | ω_H/√(\|Ω′\|A) | Range of ψ (x₀ = 0.3) |
| --- | --- | --- | --- | --- | --- |
| 0.0025 | 1.031 × 10⁻³ | 1.033 × 10⁻³ | 1.002 | 0.339 | 2.84 – 3.45 |
| 0.005 | 1.463 × 10⁻³ | 1.467 × 10⁻³ | 1.003 | 0.340 | 2.83 – 3.45 |
| 0.01 | 2.081 × 10⁻³ | 2.089 × 10⁻³ | 1.004 | 0.342 | 2.83 – 3.45 |
| 0.02 | 2.961 × 10⁻³ | 2.988 × 10⁻³ | 1.009 | 0.344 | 2.81 – 3.47 |
| 0.04 | 4.194 × 10⁻³ | 4.265 × 10⁻³ | 1.017 | 0.344 | 2.77 – 3.51 |

The secular Hamiltonian is therefore correct for this problem. The form ω ∝ (|Ω′|M)^(1/2) is justified, with a
prefactor of 0.34 instead of 1, **for a prescribed compact source and for gas near the opposite phase**.

**The region around the cloud is not a pendulum.** The kernel has a cone at the source, so the source is an elliptic
point of K with a rotation speed g M independent of the distance. The closed region around it ends at two saddles at
ψ = 0, at

$$ |\delta J_s| = \frac{g M}{\rho |\Omega'|} , $$

**linear in M**, not proportional to M^(1/2). With the exact kernel and J₀ = 0.1, the half width is 1.14 times this value and
its exponent in M is 1.003. There is no single trapping frequency here. The averaging also fails at distances
d ≲ g M/κ₀ from the source.

**Conditions for ω_trap ~ √(|Ω′|A) to be a real scale:**

1. first order in M_gas, and slow averaging (|Ω − Ω_c| ≪ κ₀);
2. a dominant smooth harmonic, far from the cone of the source;
3. a width of the island small compared with the scale on which Ω′ varies;
4. a source that stays coherent during a libration, ω τ_mix ≫ 1, as in O'Neil's theory.

**Condition 4 fails.** At M_gas/M = 0.02, ω_H τ₁ = 0.97 rad: the libration period is 2120 = 6.5 τ₁. The core is visible
from 3 τ₁ (V.B). By contrast, the transport time across the distance of the cloud at the speed g M is
ρ/(g M) ≈ 640 ≈ 2 τ₁, compatible with that time. This comparison of time scales is an estimate.

The manuscript contains a case where the pendulum does apply: Eq. (64), for halo gas far from a core. There
|ζ − ζ_c| ≈ |ζ| − |ζ_c| cos(ϑ − ϑ_c) is a pure first harmonic. That is the late-time resonance of the core with the halo,
not the formation of the core. I have not tested it here.

## 2. Which potential

The quantities below are not equivalent. Each is computed for the passive cloud in its first period T, in units of A.

| Quantity | What it measures | Reference | Range (11 clouds) |
| --- | --- | --- | --- |
| A, potential inside the cloud | Change of G M⟨1/r⟩; no force, no gas there | 1 | 1 |
| Half range at the radius of highest mixed density | Field where the gas is densest (inner turning point) | 0.922 | 0.71 – 0.98 |
| Half range at r = R₀ | Field at the circular orbit, where the core forms | 0.503 | 0.50 – 0.59 |
| Largest \|Φ − Φ_mixed\| | Deviation from the mixed state | 1.294 | 1.17 – 1.40 |
| rms of Φ − Φ_mixed at the densest radius | Typical deviation | 0.645 | 0.53 – 0.70 |
| Along an orbit at J_med, phase 0 | Potential felt by gas moving with the cloud | 1.017 | 1.01 – 1.08 |
| Along an orbit, phase π/2 | Same, a quarter of a turn away | 0.533 | 0.39 – 0.66 |
| Along an orbit, phase π | Same, opposite phase | 0.214 | 0.17 – 0.28 |
| Half depth of V(ψ, J_med) | Secular Hamiltonian | 0.245 | 0.22 – 0.26 |
| First harmonic \|V₁(J_med)\| | Secular Hamiltonian, the pendulum amplitude | 0.221 | 0.19 – 0.24 |
| \|V₁\| averaged over the mass | Same, over all actions | 0.211 | 0.18 – 0.23 |

- **Inside the cloud.** A uniform potential exerts no force, and adding any f(t) to Φ leaves the motion unchanged. A is
  fixed only by the condition Φ → 0 at infinity. It carries a dynamical meaning only through its proportionality to
  the other quantities.
- **Region of highest density.** It is close to A, except for the J⁰ cloud (0.71). It is still a change at fixed r and
  includes fast terms.
- **Deviation from the mixed state.** It is the natural size of the departure from equilibrium, but it is mostly fast.
- **Along an orbit.** It depends on the phase by a factor of about five. A gas element moving with the cloud sees mainly
  its own radial motion in the field. This variation does not by itself change the action.
- **Hamiltonian in action-angle variables.** This is the quantity that enters the equations for the slow motion, and it
  is the appropriate one for trapping. It is what ω_H measures in section 1.

## 3. Temporal asymmetry

Inside the cloud the potential is exactly

$$ \left\langle \frac{1}{r} \right\rangle(t) = \sum_J W_J\left[ b_0(J) + 2\sum_{n\ge 1} b_n(J)\,\hat{f}_n\cos\left(n\Omega(J)t\right) \right], $$

with b_n the Fourier coefficients of 1/r(Q, J) in Q and f̂_n those of the cloud in angle. The term b₀ is the mixed state.
Between the pericenter and the first apocenter, odd harmonics change sign and even ones do not. Hence **the peak to
peak difference cancels the even harmonics**, which carry the asymmetry.

| Cloud | Down/up excursion | Dipole at apocenter/at t = 0 | Share of A: n = 1 | n = 2 | n ≥ 3 | A / dipole at t = 0 |
| --- | --- | --- | --- | --- | --- | --- |
| Reference | 1.83 | −0.985 | 90.0 % | 0.8 % | 9.2 % | 1.10 |
| σ_Q = 0.2 | 1.72 | −0.985 | 92.3 % | 0.8 % | 6.9 % | 1.08 |
| σ_Q = 0.4 | 1.42 | −0.985 | 97.1 % | 0.5 % | 2.4 % | 1.02 |
| σ_J = 0.05 | 1.54 | −0.996 | 94.8 % | 0.2 % | 5.0 % | 1.05 |
| σ_J = 0.2 | 2.33 | −0.947 | 81.2 % | 3.4 % | 15.4 % | 1.20 |
| Kepler (plateau) | 1.93 | −0.979 | 88.2 % | 1.2 % | 10.6 % | 1.12 |
| J⁴ | 1.99 | −0.985 | 87.5 % | 0.9 % | 11.6 % | 1.13 |
| (J−0.05)₊² | 1.95 | −0.988 | 88.3 % | 0.7 % | 11.0 % | 1.13 |
| Ring | 1.97 | −0.997 | 88.3 % | 0.2 % | 11.6 % | 1.13 |
| Water bag | 2.00 | −0.978 | 87.2 % | 1.3 % | 11.6 % | 1.13 |
| J⁰ | 1.58 | −0.986 | 93.7 % | 0.7 % | 5.6 % | 1.06 |

- The asymmetry comes from the shape of 1/r, through n = 2. It does **not** come from phase mixing: the mass-weighted
  dipole loses only 0.3 to 5 percent by the first apocenter.
- **Correction to mu_precision.md.** I had attributed part of the asymmetry to the shear of the cloud. The band edges do
  drift 1.2 rad, but they carry little mass. I have corrected that sentence of the note.
- The peak to peak value is thus a good measure of the **radial dipole**: A is 1.02 to 1.20 times its amplitude at t = 0.
- It remains a poor measure of trapping, for the reasons of sections 1 and 2. Its robustness concerns the dipole, not the
  scale of the slow motion. A cleaner definition of the same quantity is the amplitude of the n = 1 term at t = 0.

## 4. Relation with Eq. (58) and with Γ

Eq. (58), A₅₈ = π g M (2σ_J)^(1/2), follows from A through these hypotheses:

| Hypothesis | What happens for a finite cloud |
| --- | --- |
| All the gas at one action, J = σ_J | The mean of ρ over the mass replaces (2σ_J)^(1/2) |
| Small eccentricity: 1/r linear in ξ | e = 0.28 at J = 0.1; the exact dipole is 14 % below π g⟨ρ⟩ f̂₁ in the reference cloud |
| Pure dipole | Harmonics n ≥ 3 add 2 to 15 % to A |
| No shear in half a period | Holds: the dipole loses 0.3 to 5 % |
| Compact in angle | With σ_Q = 0.4, A falls to 0.77 of its value at σ_Q = 0.1 |
| Window reaches the plateau | Holds in the isochrone potential; A is 1.8 % low in the Kepler potential |

These errors partly compensate: A/A₅₈(σ_J) = 0.976 for the reference cloud, and A/A₅₈(J_med) lies between 0.89 and 0.96 for
all the clouds with σ_Q = 0.1.

**Relation with Γ.** With A₅₈(J) = π g M (2J)^(1/2), the definition of Γ, Eq. (57), is identical to

$$ \Gamma = \frac{A_{58}(J_{med})}{\sqrt{2}\,\pi\,|\Omega_0'|\,J_{med}^2}, \qquad \mu^2 = \frac{A}{\bar{s}\,\Delta J^2} . $$

The factor 0.69 of Eq. (59) becomes 0.838 for the passive cloud because of two corrections: A/A₅₈(σ_J) = 0.976, and
|Ω₀′|/s̄ = 1.244, since the mean slope of the band is smaller than the slope at the circular orbit.

**What Γ means for a compact source.** Γ/√2 is the ratio between the shear time 1/(|Ω₀′|J) and the transport time ρ/(gM).
It also gives the size of the closed region: |δJ_s|/J = Γ/√2. In the frame of a point cloud the two fixed points between
the circular orbit and the cloud merge at

$$ \Gamma^* = \frac{2\sqrt{2}}{3\sqrt{3}} = 0.544 . $$

With the exact kernel, and removing the common frequency shift, Γ* = 0.540, 0.536, 0.528, 0.504, 0.468 and 0.407 for
J₀ = 0.005, 0.01, 0.02, 0.05, 0.1 and 0.2. These results are for a point source. The real cloud extends to the circular orbit,
and its width is of the order of its distance. That m_c reaches 0.3 at Γ ≈ 0.52 to 0.55 may be a coincidence.

**Which part of μ fails.** For ten PIC runs with a core, q is the value of each parameter divided by its value in the
reference cloud with the same m_c. A perfect parameter gives q = 1. The table shows the rms of ln q.

| Parameter | All 10 | 5 runs: σ_J and Kepler | 5 runs: other shapes |
| --- | --- | --- | --- |
| μ² = A/(ΔΩ ΔJ) | 1.09 | 0.13 | 1.53 |
| \|V₁\|/(ΔΩ ΔJ) | 1.08 | 0.19 | 1.51 |
| Γ | 0.17 | 0.19 | 0.14 |
| A/(\|Ω₀′\| J_med²) | 0.18 | 0.21 | 0.14 |
| \|V₁\|/(\|Ω₀′\| J_med²) | 0.22 | 0.28 | 0.14 |
| Half depth of V/(\|Ω₀′\| J_med²) | 0.22 | 0.28 | 0.14 |

Changing the numerator from A to the secular amplitude gives no gain. Changing the denominator from the band to the
distance reduces the error between shapes by a factor of about ten. The sample is small and uses one measure, m_c.

**Island overlap.** The usual criterion is that the pendulum island is as wide as the band, 4(|V₁|M/s̄)^(1/2) = ΔJ. It gives the
mass M_W = ΔΩ ΔJ/(16|V₁|).

| Shape | M_W with \|V₁\| | M_W with A | PIC: largest mass without a core | PIC: smallest mass with a core |
| --- | --- | --- | --- | --- |
| J² | 0.0295 | 0.0065 | 0.015 | 0.02 |
| J⁴ | 0.0283 | 0.0060 | 0.0277 | 0.0462 |
| (J−0.05)₊² | 0.0196 | 0.0042 | 0.0195 | 0.0325 |
| Ring | 0.0057 | 0.0012 | 0.0212 | 0.0423 |
| Water bag | 0.0146 | 0.0031 | 0.0211 | 0.0633 |

With |V₁|, the criterion is off for J² (too high), and for the ring and the water bag (too low). The water bag does not depend on
the cut, since its band has sharp edges. With A, the criterion fails for all the shapes.

## 5. Determination

| Option | Answer |
| --- | --- |
| Direct physical interpretation as a scale of the perturbation | **No.** It is a change of a uniform potential, without force. Its direct meaning is the radial dipole of the cloud |
| Reasonable heuristic approximation | **Yes, for the strength of the perturbation**: \|V₁\| = (0.19 to 0.24) A in the 11 clouds. **No, as a trapping frequency**: wrong prefactor, phase and time scale |
| Mainly an empirical indicator | **Yes, for the core**: μ orders the runs within one shape, and it does not order different shapes |

Consequences for earlier recommendations:

- Observation 3.3 of review_section_V.md proposed to restore the sentence that motivates μ through ω_b = (|Ω′|A)^(1/2). **I withdraw
  it in that form.** That sentence would state a mechanism which these tests do not support. If the root is kept, μ should
  be presented as a heuristic ratio.
- The passive convention remains the better definition of A. These results do not change that comparison.

**What should be computed to test the mechanism** (in order of cost; none is done):

1. **Trapped mass in the frozen field.** From the initial datum, compute K(ψ, J), its separatrices and the mass of the cloud inside
   them, and compare with m_c. No runs are needed.
2. **Origin of the core gas.** Follow back the particles that end in the core, in (ψ, J), using the existing snapshots. If they start inside
   a definite region of K, the mechanism is identified. If not, the frozen picture fails.
3. **Measured secular harmonic.** From the PIC snapshots, measure V₁(J, t) and the accumulated phase ∫ω₁(t) dt, to test whether it
   reaches order one for the masses with a core.
4. **Model with the interaction truncated.** Solve the epicyclic model keeping only the harmonic n = 1 of the kernel, and then with the full
   cone. If the core disappears with n = 1 alone, the pendulum is excluded as the mechanism. It requires new solutions of the model.

**Limits of this analysis.** The secular results are to first order in the mass and treat the source as prescribed or rigid. The test
particles validate the secular Hamiltonian for one compact shell, not for the self-consistent cloud. The test of the parameters uses
ten runs.
