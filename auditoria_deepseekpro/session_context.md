# Conversation Summary and Context for Next Sessions

## 1. Conversation Overview

The conversation progressed through three distinct phases:

**Phase 1 — Codebase audit of `vlasov-poisson_PIC`:** The user requested an analysis of the PIC (particle-in-cell) Fortran codebase for the Vlasov–Poisson equation in spherical symmetry at fixed angular momentum. The assistant examined core source files (`main.f90`, `Makefile`, `parameters.f90`, `paramfile.f90`, `arrays.f90`, `utils.f90`, `functions.f90`, `distribution.f90`, `initial_data.f90`, `density.f90`, `grav_force.f90`, `poisson_rk.f90`, `energy.f90`, `analysish.f90`, `hdf5_io.f90`, `raw_io.f90`), internal documentation (`AUDITORIA_L0_2026-09-21.md`, `BUGS_TODO.md`), and provided a comprehensive summary of the architecture, physics, and state of the project.

**Phase 2 — Peer review of experimental documents:** The user asked for expert referee criticism of two LaTeX documents (`docs/demo_eta/demo_eta.tex` and `docs/experimento_eta/bateria_eta.tex`) describing numerical experiments on self-gravity versus phase mixing in the fixed-L reduction. The assistant verified physical/numerical consistency of key figures, evaluated scientific relevance, identified strengths, and produced detailed referee-level criticisms (major and minor).

**Phase 3 — File creation:** The user requested creation of a directory `auditoria_deepseekpro` with an `.md` file capturing all observations. This was done, with a path issue resolved (file initially created at `/home/erik/Documentos/Vlasov/auditoria_deepseekpro/` instead of the repository's subdirectory; corrected via `mv` and `rmdir`).

## 2. Active Development

**Most recent task completed:** Creation of `auditoria_deepseekpro/observaciones.md` inside the repository `/home/erik/Documentos/Vlasov/vlasov-poisson_PIC/`. The file (283 lines, 12,583 bytes) consolidates:
- Code architecture description
- Physics consistency verification results
- Full referee judgment of the eta study
- Additional code observations

No in-progress implementation or debugging remains open from the assistant's side. The work is a static audit and review, not active code modification.

## 3. Technical Stack

- **Fortran 90** with OpenMP (`!$OMP PARALLEL DO SCHEDULE(STATIC)`) — the PIC solver
- **HDF5** via `h5fc` wrapper (fallback to Debian/Ubuntu `libhdf5-dev` paths) and `raw` binary I/O
- **GCC/gfortran** with flags `-O3 -funroll-loops`, `-ffree-form`, `-fopenmp`, `-fallow-argument-mismatch`, `-J$(OBJDIR)`
- **Makefile** with automatic wildcard compilation of all `src/*.f90`
- **Python** post-processing scripts referenced but not directly inspected: `equilibrio.py`, `lineal.py`, `demo_eta.py`, `landau_analisis.py`, `colas_libres.py`, `aa_numerico.py`, `barrido_lineal.py`, `eta.py`, `rawgraph_io.py`
- **Numerical methods:**
  - Integrators: Euler, leapfrog (kick-drift-kick), Yoshida 4th/6th order (table-driven compositions), analytic action-angle advance
  - `rk4` deliberately NOT implemented (non-symplectic)
  - Poisson solved via enclosed mass integration (`poisson_rk.f90`) with linear ρ interpolation between nodes
  - B-spline deposition weights `W_n`, `Sn`, images at `(−r_j, −p_j)`
  - Matrix pencil method for pole extraction (frequency/damping)
- **Physical model:** `f(r,p_r,L,t) = F(r,p_r,t)δ(L−L₀)`; phase-space measure `8π²L dr dp dL`; Hamiltonian `H = p²/2 + L₀²/(2r²) + Φ_iso + Φ_self`; fixed `L₀=2` throughout the eta experiments
- **Data formats:** HDF5 (`output_format="hdf5"`), raw native binary (`output_format="raw"`), `.par` parameter files with `key = value` syntax

## 4. File Operations

### Files created/modified during this conversation

1. **`auditoria_deepseekpro/observaciones.md`** (created, 283 lines)
   - Location: `/home/erik/Documentos/Vlasov/vlasov-poisson_PIC/auditoria_deepseekpro/observaciones.md`
   - Purpose: Comprehensive audit/review record
   - Content: see sections 1–4 of the file (code description, consistency verification, referee judgment, additional observations)
   - Final content confirmed via `head -12`, `sed -n '140,175p'`, and `wc -l`

2. **Directory `auditoria_deepseekpro/`** (created)
   - Initially at wrong location: `/home/erik/Documentos/Vlasov/auditoria_deepseekpro/`
   - Corrected via: `mv /home/erik/Documentos/Vlasov/auditoria_deepseekpro/observaciones.md /home/erik/Documentos/Vlasov/vlasov-poisson_PIC/auditoria_deepseekpro/observaciones.md`
   - Then: `rmdir /home/erik/Documentos/Vlasov/auditoria_deepseekpro`

### Key existing files referenced/read during the audit

| File | Purpose |
|---|---|
| `Makefile` | Build configuration; automatic wildcard of `src/*.f90` |
| `src/main.f90` | Main loop, all integrators (euler, leapfrog, yoshida4/6, analytic, rk4 abort) |
| `src/parameters.f90` | Global parameters |
| `src/paramfile.f90` | Parameter reader with validation; aborts on `Lfix=0` |
| `src/arrays.f90` | Particle and grid arrays |
| `src/utils.f90` | Grid, memory, RNG, pericenter, `kepler_eta`, `reduce_arrays` |
| `src/functions.f90` | B-spline kernels `W_n`, `Sn` |
| `src/distribution.f90` | `F0(Q,J)` in action-angle variables |
| `src/initial_data.f90` | Initial states (`gaussian`, `aa`, `aa_halton`, `aa_quad`, `aa_random`, `checkpoint`) |
| `src/density.f90` | Deposition with images, `avg_rho`, `curr`, `rhomix` |
| `src/grav_force.f90` | Backgrounds + centrifugal + self-gravity |
| `src/poisson_rk.f90` | Poisson via enclosed mass, mirror exterior interpolation |
| `src/energy.f90` | Energies |
| `src/analysish.f90` | `h_k` projection; factorized Q quadrature; complex output |
| `src/hdf5_io.f90` | HDF5 output |
| `src/raw_io.f90` | Raw binary output |
| `AUDITORIA_L0_2026-09-21.md` | Exhaustive audit (E1–E22 findings, most corrected) |
| `BUGS_TODO.md` | Bug log with corrections, pending items |
| `docs/demo_eta/demo_eta.tex` | 850-line demo document (14 PIC runs) |
| `docs/experimento_eta/bateria_eta.tex` | 2432-line battery document (30 proposed runs) |

### Key quantitative checks documented in observaciones.md

- Isócrono: `c = 2.4142`, `Ω(0)=c^{-3}=0.07106`; epicíclica `κ² = 4πρ(r_c) + GM(<r)/r_c³` verified numerically
- Free decay: `σ_Ω` inferred values 2.73–3.07×10⁻³; `ΔΩ/σ_Ω ≈ 6.2–7.0` vs expected 6.1
- O'Neil parameter: `ν ≈ 11√ε` matches 3.0 (ε=0.075) and 10 (ε=0.84)
- Edge resonance positions: L5 (J_r≈0.140), L6 (J_r≈0.164) internally consistent
- Timescales: `τ₁ = 2π/ΔΩ = 277–556`; D5 rebound period 5700 ≈ 14τ₁

## 5. Solutions & Troubleshooting

### Issue 1: File created at incorrect path

**Problem:** The `create_new_file` tool created the audit file at `/home/erik/Documentos/Vlasov/auditoria_deepseekpro/observaciones.md` (wrong location: parent Vlasov directory instead of the PIC repo subdirectory).

**Detection:** `ls -la auditoria_deepseekpro/` inside the repo showed empty directory; `find` located the orphaned file at the wrong path.

**Solution:**
```
mv /home/erik/Documentos/Vlasov/auditoria_deepseekpro/observaciones.md /home/erik/Documentos/Vlasov/vlasov-poisson_PIC/auditoria_deepseekpro/observaciones.md
rmdir /home/erik/Documentos/Vlasov/auditoria_deepseekpro
```
Verified final file at correct location (283 lines, 12,583 bytes).

### Issue 2: No resolvable bug in conversation flow

A system message about reversion appeared in the transcript, but the conversation state was already valid and the assistant continued normally. No recovery action required.

### Background: prior work already completed (from codebase docs)

The repository contains a prior exhaustive audit (`AUDITORIA_L0_2026-09-21.md`, `BUGS_TODO.md`) with findings E1–E22, most corrected and verified bit-for-bit:
- E1: L₀=0 aborts (mitigated)
- E2: `state=gaussian` fixed
- E3: Reflection changes force sign (fixed)
- E4: Backgrounds `iso`, `isotrun`, `nfw`, `burkert` act on particles (fixed)
- E5: `null`/`sphere` reset/sum forces correctly (fixed)
- E6: `reduce_arrays` recomputes force (fixed)
- E7: Unbound nodes/NaN excluded from `h_k` (fixed)
- E8: Poisson via enclosed mass (fixed)
- E9: Mirror images in deposition (fixed)
- E10: Mirror exterior interpolation (fixed)
- E11: Courant pericenter constraint (fixed)
- E12/E13: Kepler Newton safeguarding + radicand clamping (fixed)
- E14: `vlasov_rhomix` shell density (fixed)
- E15: Midpoint nodes (`aa`, `aa_halton`) (fixed)
- E16: Grid size rounding (fixed)
- E17: Deterministic reductions (fixed)
- E19: Dead states removed (fixed)
- E20: Documentation updated (fixed)
- E21: Output density uses W_n volume (fixed)
- E22: `Sn(1)` half-weight at cell faces (fixed)
- E18: `eps` softening of centrifugal term (documented, kept 0 by default)

## 6. Outstanding Work

### Regarding the eta-study review (acceptance conditions for publication)

The referee verdict was **"major revisions"**. Four conditions must be met before the non-linear results (D3, D5, D6, D8) can be considered publishable:

1. **Firmly bound γ of the linear solver for D3** — currently γ_linear "<10⁻⁵", γ_PIC≈5×10⁻⁵; if true γ_linear is 4–5×10⁻⁵, the observed 19% loss is linear, not novel. Must demonstrate matrix pencil resolution and perform ε-sweep (D3 at ε→0.003) to show loss disappears when `ω_b < Ω_min − ω`.

2. **Convergence tests for D3/D5/D6** — split Δt, quadruplicate N; show fine component scales as N^(-1/2), smooth component invariant.

3. **Extend D8 with optimal ε (≈0.3–0.5)** and subtract static component to escape measurement floor (current window only ~13× above 2×10⁻⁴ floor, lost at ~15τ₁).

4. **Correct the conceptual error** in `bateria_eta.tex` about the Fortran `h_k` Q-factor `e^{-sin²(Q/2)/s_Q²}`: it is NOT "a constant per harmonic" — it convolves harmonics; only constant in the s_Q→∞ limit.

### Pending code observations (not part of current request but documented)

- Three dead files in `src/`: `poisson`, `poisson_ps`, `reduce_arrays` (no `.f90` extension, don't compile) — candidates for cleanup
- `aa_random` has no attempt limit when F≈0 in sampled region
- Documentation needed for output density using W_n volume vs. Poisson using `avg_rho`

### Pending scientific questions in repository (`PREGUNTAS_ABIERTAS.md`)

- Growing discreteness noise mechanism
- Extension of Landau measurements to `a0 ≥ 1e-2`
- Additional runs for the eta battery (block V to confirm η_c robustness, block B for convergence)

### What the user may request next

The user asked for the summary to enable seamless continuation. Logical next steps could include:
- Implementing/fixing the four referee conditions in the eta documents near the scripts `reproducir/scripts/demo_eta.py` or `lineal.py`
- Cleaning up dead `src/` files
- Running the convergence tests specified as outstanding work
- Addressing the `bateria_eta.tex` h_k factor error