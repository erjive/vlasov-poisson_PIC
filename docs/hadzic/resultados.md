# First results in the setting of Hadžić, Rein, Schrecker and Straub

These are the notes of the first series (a₀ = 0.01). The full results, with the runs at
a₀ = 1 and ε = 0.5, the η table and the figures, are in `hadzic_informe.pdf`.

Setting of [HRSS23] (see `hadzic_base.md`): all particles with |L| = L₀ = 2, a point mass of
mass 1 at the centre, and steady states that are polytropes in the energy,
F = A (E_t − E)ᵏ, with the edge at the action J_t = 0.7. Without self-gravity this is
E_t = −0.0686, inside their single-gap range (−0.079 < E_t < 0), with
Ω_max/Ω_min = 2.46. Everything below is produced by `reproducir/scripts/hadzic.py`
(steps `lineal`, `preparar`, `correr`, `analizar`); the outputs are in `exe/hadzic/`.

## 1. The loop gain λ_edge(k, a₀)

λ_edge is the gain of the self-consistent loop at the lower edge of the band, that is,
the norm of their Mathur operator at the edge of the principal gap. A discrete mode below
the band exists if and only if λ_edge > 1. δ_d = (Ω_min − ω_d)/ΔΩ is the distance of the
mode below the edge.

| k \ a₀ | 0.01 | 0.03 | 0.1 | 0.3 | 0.6 | 1.0 |
|---|---|---|---|---|---|---|
| 0.75 | ∞ | ∞ (δ_d ≈ 2×10⁻⁶) | ∞ (δ_d ≈ 10⁻⁴) | ∞ (δ_d = 0.004) | ∞ (δ_d = 0.027) | ∞ (δ_d = 0.075) |
| 1 | ∞ | ∞ | ∞ | ∞ (δ_d = 0.0002) | ∞ (δ_d = 0.0085) | ∞ (δ_d = 0.044) |
| 1.25 | 0.025 | 0.075 | 0.240 | 0.643 | 1.112 (δ_d = 0.0001) | 1.558 (δ_d = 0.019) |
| 1.5 | 0.018 | 0.054 | 0.170 | 0.454 | 0.759 | 1.050 (δ_d = 0.0015) |
| 2 | 0.016 | 0.046 | 0.145 | 0.371 | 0.609 | 0.820 |
| 3 | 0.016 | 0.047 | 0.145 | 0.358 | 0.567 | 0.738 |

- **The dichotomy holds at small mass.** For k > 1, λ_edge ∝ a₀ (0.016–0.025 at a₀ = 0.01),
  as their bound M_λ ≤ Cε says, so there is no discrete mode. For k ≤ 1, λ diverges at the
  edge: as ln(1/δ) for k = 1 and as a power for k = 0.75. The mode exists, but at small mass
  it is bound by very little: δ_d < 10⁻⁶ for k = 0.75 at a₀ = 0.01 and for k = 1 up to
  a₀ = 0.1.
- **The threshold ε₀(k), which the theorem does not estimate.** Interpolating λ_edge = 1:
  a₀ ≈ 0.53 for k = 1.25 and a₀ ≈ 0.93 for k = 1.5. For k = 2 and 3 there is no discrete
  mode even with a shell as massive as the central mass (λ_edge = 0.82 and 0.74 at a₀ = 1).
  In the isochrone, with a narrow band, the threshold for the Wilson model (g = 2) was
  a₀ ≈ 0.091: the wide band of the point mass mixes much more strongly.
- **Their hypotheses.** The frequency Ω(J) is monotonic in all 36 equilibria. The
  single-gap condition, Ω_max/Ω_min ≥ 2, holds up to a₀ = 0.3. It is lost at a₀ = 0.6 for
  k = 0.75 (1.99) and at a₀ = 1 for k ≤ 1.5 (1.87–1.98).

## 2. Linear response in time at a₀ = 0.01

The linearized solver (1600 × 32 nodes, Δt = 0.5, up to t = 8000 ≈ 95 τ₁) gives an
algebraic decay for every k. The slope of the envelope against ln t over 1000 ≤ t ≤ 8000:

| k | \|h₁\| | ‖δΦ‖ |
|---|---|---|
| 0.75 | −1.67 | −1.60 |
| 1 | −2.00 | −1.93 |
| 1.5 | −2.55 | −2.50 |
| 2 | −3.00 | −2.59 |

The decay is close to t^{−(k+1)}, the tail of an edge at which the weight vanishes as
(J_t − J)ᵏ (battery η, Section 1.4). The sharper the edge, the slower the decay. For k ≤ 1
the discrete mode is too weakly bound to show up in this time. Their bound for pure
transport, ∂ₜU ≲ t^{−min(2,k)} [HRSS24], is an upper bound; the data used here vanish at the
circular orbit (s = (J/J_t)^{3/2}), which removes the limitation of the elliptic point.

## 3. PIC runs at a₀ = 0.01

Eight runs, 400 × 25 particles, Δr = 0.1, Δt = 0.1 (courant = 2, as in 11_landau),
BGtype = "sphere", up to t = 8000: for each k = 0.75, 1, 1.5 and 2, a run with ε = 0.1 and
its reference with ε = 0. The perturbation is (h₁, δΦ) of the D run minus that of the Z run,
divided by ε; h₁ uses B(J) = J² e^{−(J−0.35)²/0.2²} and the angle–action map of the equilibrium.

- **The PIC runs follow linear theory** in |h₁| to 1 % up to t ≈ 1000 (12 τ₁), except k = 2
  (5 %). At t ≈ 2000 the difference is 1 % for k = 0.75 and 8–50 % for the others. The energy
  is conserved to |ΔE/E| ≤ 8×10⁻⁸.
- **After that the runs reach a floor.** The D and Z runs decorrelate: at t ≈ 8000,
  |h₁,D − h₁,Z| is of the order of the noise of the unperturbed run itself, and ‖δΦ_ε‖ stays at
  2–4×10⁻⁷. The late signal is noise, not physics.
- At a₀ = 0.01 the PIC runs therefore show the damped edge response of every k, and cannot
  show the non-damping of k ≤ 1: the mode is bound by too little.

## 4. Next steps

1. **PIC runs where the mode of k ≤ 1 is bound visibly**: k = 0.75 at a₀ = 0.6 (δ_d = 0.027) or
   1.0 (δ_d = 0.075), and k = 1 at a₀ = 1.0 (δ_d = 0.044).
2. **A direct test of ε₀(k)**: at a₀ = 1.0, k = 1.25 has a discrete mode (δ_d = 0.019) and k = 2
   does not (λ_edge = 0.82).
3. **A lower floor**: larger ε (0.3–0.5, as in the demo) or more particles, to follow the decay
   beyond 12–24 τ₁.
