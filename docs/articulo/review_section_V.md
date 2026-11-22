# Internal review of Section V (Initial data far from equilibrium)

22 November 2026. Copy of the Claude Doc https://claude.ai/artifact/KYMxCYMMaT6Xzfg28JreLV (revision 13).

## Scope and method

This review covers Section V, "Initial data far from equilibrium" (file sec6_far.tex, pages 13 to 20 of the PDF), in the version of commit a71b935 on the branch notas/mezcla-fuera-equilibrio (18 November 2026). No file of the manuscript was modified.

I read Section V in full, and the parts of Sections II to IV, VI to VIII, Appendix B, the abstract, the introduction and the conclusions that define or reuse its quantities. Each claim was checked against the equations, tables and figures of the manuscript. Where a point can only be settled with the run data, or with results that are not in the manuscript, the observation says so.

Importance is **high** when the issue affects a claim of the section or how a reader understands its main result; **medium** when it weakens an argument or makes it hard to follow; **low** for wording, notation and layout. Observations are numbered by block (1.1, 1.2, ...) so that they can be cited in comments.

## Summary

Section V supports its main observation: within 30 τ₁, the change between mixing and a dense core at the circular orbit is steep and continuous. Seven observations are of high importance, and only one of them (5.2) needs new simulations.

The weak points are what the window means of |h₁| measure, the definitions of μ and of the core, the convergence in N_Q near the change, and the place where Γ enters.

| # | Finding | Importance | Observations | New runs |
| --- | --- | --- | --- | --- |
| 1 | The means over 10 ≤ t/τ₁ ≤ 20 are called a "level", but Section VI shows that they sample a decaying tail. The steepness of the change and the absence of a threshold therefore refer to a fixed window. | High | 2.1, 2.3 | No |
| 2 | μ takes three values for the same run (0.52, 0.54 and 0.49 for M_gas/M = 0.03). For M_gas/M ≥ 0.05 it is partly an outcome of the run, and its square root is no longer motivated. | High | 3.1 to 3.5 | No |
| 3 | The two criteria for a core disagree near the change (μ = 0.37; σ_J = 0.2 with μ = 0.38). Without self-gravity, m_c = 0.051 exceeds the bound m_mix = 0.042. | High | 3.6, 3.7 | No (existing data) |
| 4 | Near the change, going from N_Q = 25 to 100 changes the mean modulus of h₁/h₀ by a factor 1.75 at μ = 0.43, as much as one step of 0.01 in μ. No finer N_Q was tested there. | High | 5.2, 5.3 | Yes |
| 5 | Γ is named in V.E and used in VI before its definition in VII.B, while V.A announces the main result of VII.E. | High | 6.1, 1.2 | No |
| 6 | V.B repeats the decay analysis of VI on coarser lattices. "Steep but continuous" is stated five times, and Tables I and II are cited before m_c and a_c are defined. | Medium | 1.1, 1.3 to 1.7 | No |
| 7 | The same values appear in Tables I to III, Table VIII and Figs. 5 and 9. Fig. 4(b) shows the same curves as Fig. 10. | Medium | 4.1 to 4.3 | No |

Suggested order of work:

1. Decide the parameter question (3.1 to 3.5, 6.1, 6.2), since it touches Sections V to IX.
2. Correct the wording of "level", "not damped" and "no critical value" (2.1, 2.3, 2.6).
3. Make the definition of the core consistent with the data already computed (3.6 to 3.9).
4. Reorganize the section and consolidate its figures and tables (block 1 and block 4).
5. Decide whether to run the tests of 5.2 to 5.6.

## 1. Structure and narrative

The order of the subsections is adequate: definitions, evolution, late state, core and comparison between families. The problems are a preview of Section VII in V.A, an overlap with Section VI, and repetitions.

| Subsection | Question it answers | Main evidence |
| --- | --- | --- |
| V.A Simulations and μ | Which runs are used, and which parameter compares self-gravity with the spread of the frequencies? | List of runs; τ₁, Eq. (44); μ, Eq. (45) |
| V.B Evolution of the cloud | Does the cloud still mix as its mass grows, and does its field relax? | Figs. 2 and 3 |
| V.C State at late times | Which state does the gas reach, and when is part of it not mixed? | Fig. 4(a), Table I, Eq. (46) |
| V.D Oscillation of the core | Is the radial oscillation of the core damped? | Fig. 4(b), Table II |
| V.E Shape and external potential | Does μ locate the change for other widths and for the Kepler potential? | Fig. 5, Table III |

### 1.1 Balance between mixing and the core

- **Location.** V.C (paragraphs 3 and 4), V.D and V.E.
- **Problem.** Three of the five subsections describe the core: its mass, its oscillation and where it appears. Mixing and the relaxation of the field take two paragraphs of V.B and one of V.C.
- **Why it matters.** The title and Section VI center the article on phase mixing. In Section V the core reads as the subject, not as the place where relaxation fails. The scheme of 16 November, cited in the note of VI, puts the core with the epicyclic model, in shorter form.
- **Importance.** Medium.
- **Recommendation.** Keep in V what answers three questions for each mass: does the gas mix, does the field relax, and which state is reached. Move the detailed dynamics of the oscillation (V.D paragraphs 2 to 4, Fig. 4(b), the a_c columns of Table II) next to VII.C and VII.D, where the model explains it.

### 1.2 V.A announces the main result of Section VII

- **Location.** V.A, paragraph 3, from "In section VII we derive a model" to "measured directly in the simulations"; V.E, paragraph 2, last four sentences.
- **Problem.** Before any result, V.A describes the mechanism of the core and states that the parameter of the model, not μ, decides whether a core forms. V.E repeats this with Γ, which is not yet defined.
- **Why it matters.** The reader meets the conclusion of VII.E before the evidence of V. Section V then uses μ while saying that μ is not the right parameter.
- **Importance.** Medium.
- **Recommendation.** In V.A keep one sentence: μ assumes that the width of the cloud sets both the spread of the frequencies and the distance from the circular orbit (see 3.5). Leave the comparison of μ and Γ to VII.E, or introduce Γ as proposed in 6.1.

### 1.3 V.B repeats the decay analysis of Section VI

- **Location.** V.B, paragraphs 2 and 3; Fig. 3.
- **Problem.** V.B gives the exponent without self-gravity (−3.03 over 10 ≤ t/τ₁ ≤ 30) and the relaxation of the field with self-gravity. VI.B gives −3.05 over 10 ≤ t/τ₁ ≤ 20 on 3200 × 100 points, and the slopes with self-gravity.
- **Why it matters.** The same quantity appears twice, with different windows and lattices. VI.A states that on the lattices of V the higher harmonics are dominated by discreteness from t ≈ 8τ₁ on.
- **Importance.** Medium.
- **Recommendation.** In V.B keep what depends on μ: the outcome (mixing or core) and the time at which ‖δΦ‖ falls below 10⁻³. Refer to VI.B for the decay laws.

### 1.4 "Steep but continuous" is stated five times

- **Location.** Opening paragraph; V.B, end of paragraphs 2 and 3; V.D, paragraphs 3 and 4; V.E, paragraph 2.
- **Problem.** The same conclusion is drawn five times, each time from a different observable.
- **Why it matters.** The evidence sits in one place, Table II, and the repetitions dilute it.
- **Importance.** Low.
- **Recommendation.** One statement after Table II, which collects |h₁|, ‖δΦ‖, m_c and a_c, and short pointers elsewhere.

### 1.5 Tables cited before their columns are defined

- **Location.** V.B, paragraphs 2 and 3 (Table II); V.C, paragraph 2 (Table I). m_c is defined by Eq. (46) and a_c in V.D, paragraph 1.
- **Problem.** Tables I and II are first cited before m_c and a_c are defined. V.B also names the core before V.C defines it.
- **Why it matters.** A reader who follows the citation meets columns that have not been explained.
- **Importance.** Low.
- **Recommendation.** Cite Table I after Eq. (46) and Table II from V.D on. In V.B describe the core as "a dense part of the cloud at the circular orbit (section V.C)".

### 1.6 The statement on resolution sits in V.E

- **Location.** V.E, paragraph 3 ("The results do not depend on the lattice, except close to the steep increase").
- **Problem.** A statement about all the simulations appears in the subsection on the families. It is also stronger than Appendix B.4, where a_c differs by up to a factor four between lattices at μ = 0.48.
- **Why it matters.** The limits of the resolution should be known before the results are read.
- **Importance.** Low.
- **Recommendation.** Move it to V.A, after the list of convergence runs, with the qualifications of B.4 (see 5.2 and 5.6).

### 1.7 Opening and closing paragraphs

- **Location.** Opening paragraph of V; last paragraph of V.E and its \EJ note.
- **Problem.** The opening gives the values 0.43 and 0.48 and the minimum and regrowth of the oscillation before any definition, and it does not mention Section VI. The section ends with the Kepler cloud and a note that asks for a test which VII.E already reports (PIC in the isochrone potential, model in the Kepler potential).
- **Why it matters.** The reader gets no statement of what V establishes, and the note is partly outdated.
- **Importance.** Low.
- **Recommendation.** Keep the opening qualitative and name Section VI. Replace the note by a pointer to VII.E. Close V with three or four sentences: the outcomes as functions of μ, the continuous change within 30τ₁, what μ does not capture, and the questions that VI and VII take up.

## 2. Physical consistency

Most interpretations are supported by the figures and tables of V. Two claims go beyond what a finite window can show: the "level" of |h₁| and the absence of a critical value. Others merge an observation with its interpretation.

| Claim (location) | Type | Support in the manuscript | Assessment |
| --- | --- | --- | --- |
| The cloud mixes for μ = 0.31, leaves a core for μ = 0.43, and the core oscillates for μ = 0.89 (V.B) | Observation | Fig. 2, Table I | Supported |
| Without self-gravity \|h₁\| ∝ t⁻³, the tail of Eq. (33) for J = 0 (V.B) | Observation and theory | Fig. 3(a), Section III; VI.B (exact) | Supported |
| The level reached at late times increases with μ (V.B) | Interpretation | Means over 10 ≤ t/τ₁ ≤ 20 | Contradicted by VI.A: a tail (2.1) |
| The field approaches a stationary profile while the DF mixes, the analogue of Landau damping (V.B) | Observation and analogy | Fig. 3(b) | Supported up to about 15τ₁; floors not separated (2.7) |
| The late state differs from F̂₀: 6.1 instead of 0.5 percent of the mass at J < 0.2σ_J (V.C) | Observation | Table I, three lattices | Supported |
| The redistribution happens because the actions change while the potential oscillates (V.C) | Interpretation | Timing: 1.0 to 5.8 percent within 3τ₁ | Plausible; the model of VII reproduces 0.060 |
| Values above F̂_max show transport to the circular orbit before mixing (V.C) | Interpretation | Bound F̂_max, Fig. 2 | The bound alone shows less (2.5) |
| The oscillation of the core is not damped within 30τ₁ (V.D) | Observation | Fig. 4(b), Table II | Wording (2.6) |
| The minimum and the regrowth at μ = 0.43 are not numerical artifacts (V.D) | Observation | N_Q, N_J, Δt, Δr | Existence supported; amplitude depends on Δr (5.4) |
| The change between μ = 0.43 and 0.48 is steep but continuous; no critical value (V.B, V.D, V.E) | Observation in a window; conjecture as t → ∞ | Table II | Valid within 30τ₁ (2.3) |
| For σ_Q ≤ 0.2 all families change at the same μ (V.E) | Observation | Table III | Only at the resolution of the brackets (2.4) |
| The core is self-bound for μ = 0.89 (V.D) | Inference | Shear time, estimate not shown | Needs the estimate (2.8) |
| Resemblance with trapped plasma waves and core-halo states (V.D) | Analogy | References in V.D | Hedged correctly |

### 2.1 A decaying tail called a level

- **Location.** V.B, paragraph 2 ("the level which is reached at late times increases with μ"); the means over 10 ≤ t/τ₁ ≤ 20 in Tables I to III and Fig. 5(a).
- **Problem.** For μ ≤ 0.43 the mean samples a tail that keeps decaying; VI.A says so and carries a note to change the wording of V.B. For μ ≥ 0.52 the same mean measures a persistent oscillation.
- **Why it matters.** The factor 1.4 to 1.7 per step of 0.01 in μ compares a tail at a chosen time with a plateau. A later window would make the change steeper, an earlier one smoother.
- **Importance.** High.
- **Recommendation.** Write "the mean value over 10 ≤ t/τ₁ ≤ 20" everywhere. Say once that below the change it keeps decaying (VI.B) and above it does not, so that the steepness depends on the window.

### 2.2 What |h₁| measures when a core is present

- **Location.** V.B, paragraph 2 (values 0.19 to 0.85 for μ ≥ 0.52); V.E, paragraph 1; Fig. 5(a).
- **Problem.** With B = 1, |h₁|/h₀ is the mass average of e^(−iQ). Above the change it is dominated by the coherent oscillation of the core, not by gas that is not mixed. The manuscript's numbers fit this reading: at μ = 0.43, a_c ≈ 0.05 against a core size of 0.15 goes with |h₁|/h₀ = 0.012 and m_c = 0.10; at μ = 0.89, an amplitude of 1.0 against 0.17 goes with 0.85 and m_c = 0.52.
- **Why it matters.** The growth of |h₁| by one to two orders of magnitude can be read as a failure of the mixing of the whole cloud. This reading would also explain the non-monotonic values above the change (0.32, 0.19 and 0.58 for μ = 0.60, 0.66 and 0.76), which the text does not discuss.
- **Importance.** Medium.
- **Recommendation.** State what |h₁| measures in the presence of a core. If the data allow it, give |h₁| for the gas outside the core. This reading is ours and should be checked with the data.

### 2.3 "No critical value" refers to a finite window

- **Location.** V.B, end of paragraph 3; V.D, paragraphs 3 and 4; the opening paragraph. The same claim is in the abstract and in Section IX.
- **Problem.** m_c is fixed during the first units of τ₁, and its continuity does not depend on the window. |h₁| and ‖δΦ‖ (10 ≤ t/τ₁ ≤ 20) and a_c (t ≥ 25τ₁) are taken in fixed windows, and V.D says that 30τ₁ do not decide whether the oscillation is eventually damped.
- **Why it matters.** If the oscillation decayed below some μ and persisted above it, the late-time limit would have a threshold, while every finite window is continuous. The text does not separate the two statements.
- **Importance.** High.
- **Recommendation.** Write "within 30τ₁ all the quantities change continuously with μ; the mass of the core does so independently of the window". State that a threshold for the late oscillation remains open, citing VII.D (the model up to 100τ₁).

### 2.4 Same location for all families, and m_c as a function of μ alone

- **Location.** V.E, paragraph 1 ("the location of this interval is the same for all the families") and paragraph 2 ("m_c is nearly a function of μ alone"); Table III.
- **Problem.** The brackets of Table III are 0.14 to 0.16 wide in μ for the four other families with σ_Q ≤ 0.2 (for instance 0.35 to 0.49 for σ_J = 0.05). Their common part is 0.43 to 0.48: they are consistent with one location but do not establish it. The values of m_c differ by −12 to +28 percent, or by up to 0.07 in μ, more than the width 0.05 of the change in the reference family.
- **Why it matters.** The first sentence claims more than the brackets show, and the second describes a collapse whose spread exceeds the feature it locates.
- **Importance.** Medium.
- **Recommendation.** "Consistent with the same location, within brackets of width 0.14 to 0.16." Give the residual differences and attribute them, as VII.D does, to the eccentricity of the orbits (ratio 0.21, 0.29 and 0.41 for σ_J = 0.05, 0.1 and 0.2).

### 2.5 From the bound F̂_max to transport before mixing

- **Location.** V.C, paragraph 3 ("Values above F̂_max cannot be obtained by mixing the cloud at fixed action. They show that part of the dense cloud has been transported to the circular orbit before it was mixed.").
- **Problem.** The bound shows that the gas at that action is not a mixture at fixed action. It does not by itself show transport of an unmixed piece, since V.C, paragraph 2, says that the actions change during the first units of τ₁.
- **Why it matters.** An observation (the bound is exceeded) and an interpretation (a compact part arrives before it is wound up) are merged. The interpretation rests on Fig. 2 and on Section VII.
- **Importance.** Medium.
- **Recommendation.** "Values above F̂_max show that the gas there is not a mixture at fixed action. Figure 2 shows that it is a compact part of the cloud, which reaches the circular orbit before it is wound up."

### 2.6 "Not damped" when the amplitude falls by two orders of magnitude

- **Location.** V.D, paragraph 1 ("For none of the masses which we have simulated is this oscillation damped within 30τ₁"); the abstract and Section IX.
- **Problem.** At 2τ₁ every core has a_c between 0.2 and 1.1. At μ = 0.43, a_c then falls to 0.007 at 15τ₁ and grows to 0.05 or 0.06, about one third of the size of the core (0.15).
- **Why it matters.** "Not damped" can be read as "does not decay".
- **Importance.** Medium.
- **Recommendation.** "Does not decay to zero within 30τ₁; for μ = 0.43 it decays by two orders of magnitude and then grows again to about one third of the size of the core."

### 2.7 Relaxation of the field: threshold, floors and recurrence

- **Location.** V.B, paragraph 3; Fig. 3(b).
- **Problem.** The relaxation times use the level 10⁻³, which is not justified. The means over 10 ≤ t/τ₁ ≤ 20 for μ ≤ 0.43 (7 × 10⁻⁵ to 2.4 × 10⁻⁴) overlap the interval in which the text attributes the growth of the curve without self-gravity to the recurrence of the harmonic n = 10 (about 15τ₁, Eq. (B5)). The reference profile is the average of the same run over its second half.
- **Why it matters.** The comparison with 4 × 10⁻⁵ without self-gravity is partly a comparison of numerical floors. The manuscript does not show whether the values for μ ≤ 0.31 lie above them.
- **Importance.** Medium.
- **Recommendation.** Show that the order of the times does not depend on the level 10⁻³. Compare the means only up to 15τ₁, or cite VI.B (3200 × 100) for small μ. Name the curve without self-gravity as the floor.

### 2.8 Self-binding of the core for μ = 0.89

- **Location.** V.D, paragraph 4, last two sentences.
- **Problem.** The times "about 10τ₁" and "exceeds 40τ₁" are given without the estimate behind them.
- **Why it matters.** The conclusion that the core is held together by its own gravity rests only on these numbers.
- **Importance.** Low.
- **Recommendation.** Give the estimate in one sentence: the spread of Ω across the core for its measured amplitude and size, and the time in which the spread of the phases reaches 2π.

### 2.9 Frequencies quoted for the core

- **Location.** V.B, paragraph 1 (0.074 for μ = 0.89, "close to" Ω(0) = 0.0711); V.D, paragraph 2 (0.0730 for μ = 0.43).
- **Problem.** 0.074 is 4 percent above Ω(0); the relevant comparison is the frequency shifted by the gas (VII.C). For the same core at μ = 0.43, VII.D gives 0.07284 and 0.07275, not 0.0730.
- **Why it matters.** Two sections quote different values for the same quantity.
- **Importance.** Low.
- **Recommendation.** Quote the frequency once, with its method (fit of r_c or slope of Q), and compare it with ω_c of VII.C.

### 2.10 Three values for the initial mass fraction at J < 0.2σ_J

- **Location.** V.C, paragraph 2 ("grows from 1.0 to 5.8 percent"); Table I (0.005 without self-gravity); VI.C ("0.6 percent in the initial datum").
- **Problem.** The fraction at t = 0 is 1.0, 0.5 or 0.6 percent depending on the place. Our check: the exact value for the datum is 0.59 percent; with 400 rows the interval holds 13 rows, which gives 0.55 percent; 1.0 percent is the value in the variables of the instantaneous potential.
- **Why it matters.** A reader may suspect that part of the growth comes from the change of variables.
- **Importance.** Low.
- **Recommendation.** One sentence: at t = 0 the actions of the instantaneous potential differ from those of the external potential at first order in M_gas/M (Appendix B.2), which nearly doubles this fraction.

### 2.11 Two short logical points

- **Location.** V.E, paragraph 2 ("the core consists of gas which has not been mixed, such that its DF cannot exceed F_max"); V.E, last paragraph ("which is the behavior described by Ref. [sCjL26] for small data").
- **Problem.** The DF never exceeds F_max, mixed or not, because it is constant along the characteristics. The Kepler clouds differ from the assumptions of Chaturvedi and Luk, as V.A says.
- **Why it matters.** The first reason is not the right one, and the second sentence suggests more agreement than there is.
- **Importance.** Low.
- **Recommendation.** "Since the DF is constant along the characteristics, it never exceeds F_max, which for a given mass decreases as σ_Q grows." And "which is consistent with the result of Ref. [sCjL26] for small data, although the assumptions differ (section V.A)".

## 3. Parameters and definitions

μ is a useful heuristic for clouds of one shape, but in its present form it is partly an outcome of the runs, it has three numerical conventions in the manuscript, and its absolute value depends on arbitrary choices. The core has two criteria which disagree near the change.

| Phenomenon | Captured by μ? | Evidence in the manuscript |
| --- | --- | --- |
| Location of the change for σ_J = 0.05, 0.1, 0.2 and for the Kepler potential (σ_Q = 0.1) | Yes, within brackets 0.14 to 0.16 wide | Table III, Fig. 5 |
| Mass of the core at a given μ for these families | Approximately: −12 to +28 percent | V.E, paragraph 2 |
| Width in angle | No: with σ_Q = 0.4 the change lies between μ = 0.52 and 0.69, and the core is smaller | Table III, V.E |
| Other dependences of the DF on J | No: μ at m_c = 0.1 ranges from 0.45 to 1.44 | VII.E, Table VIII |
| Distribution of the angular momenta | Not tested in V; σ_L/L₀ = 1.41 percent removes the core at μ = 0.52 | Section VIII |
| Persistence of the oscillation beyond 30τ₁ | Not decidable | V.D, VII.D |
| Rate of mixing below the change | Partly: at fixed Γ the exponent depends on the shape | VI.B |

### 3.1 μ is partly an outcome of the run

- **Location.** V.A, paragraph 4 (measurement of A); Table I; Fig. 5.
- **Problem.** A is measured in the self-gravitating run, during its first radial period. For M_gas/M ≤ 0.025 it equals the amplitude of the passive cloud within 5 percent. For larger masses it includes the contraction of the cloud: μ = 0.89 for M_gas/M = 0.1, where the law μ = 3.05(M_gas/M)^(1/2) gives 0.96, and 1.16 for 0.3 instead of 1.67.
- **Why it matters.** A control parameter should be fixed by the initial datum. With Γ = 11.43 M_gas/M, μ²/Γ falls from 0.81 for M_gas/M ≤ 0.025 to 0.69 at 0.1 and 0.39 at 0.3. "Γ proportional to μ²" (V.E, VII.B, IX) holds only for the smaller masses, and the axis of Fig. 5 is compressed at large μ.
- **Importance.** High.
- **Recommendation.** Compute μ for all the runs with the amplitude of the passive cloud, as VII.E already does, or with Eq. (58). State the range in which μ² is proportional to M_gas.

### 3.2 Three values of μ for the same run

- **Location.** V.A and Table I (0.43, 0.52 and 0.66 for M_gas/M = 0.02, 0.03 and 0.05); Table VIII (0.44, 0.54 and 0.69); Eq. (59) with Γ = 11.43 M_gas/M (0.40, 0.49 and 0.63). Sections VI and VIII use the values of V.
- **Problem.** Three conventions for A and ΔΩ coexist: measured in the run, passive cloud, and leading order in J.
- **Why it matters.** A reader who compares Table I with Table VIII sees the same simulation with different μ. The differences, 2 to 8 percent, are comparable to the width of the change (0.43 to 0.48).
- **Importance.** High.
- **Recommendation.** One definition for the whole article. If two are needed, give both in one table and say which one each figure uses.

### 3.3 The square root in μ is not motivated

- **Location.** V.A, paragraph 3, Eq. (45).
- **Problem.** The text says that the square of μ compares the change of the frequencies with their spread, but gives no reason for the root. A draft before commit ef178cd (15 November) gave it: μ = ω_b/ΔΩ, with the trapping frequency ω_b = (|Ω′|A)^(1/2), as in the criterion of O'Neil [tO65].
- **Why it matters.** Without it the natural parameter is μ², which is proportional to the mass. The root also changes how a reader takes "μ ≳ 0.4" and the factor per step of 0.01.
- **Importance.** Medium.
- **Recommendation.** Restore one sentence on ω_b, which also links V.D with Eq. (64) of VII.C, or use μ² and say that it is proportional to M_gas.

### 3.4 The value of μ depends on the one-percent convention

- **Location.** V.A, paragraph 2 (definition of the interval of actions whose ends set ΔJ and ΔΩ).
- **Problem.** ΔJ, ΔΩ, τ₁ and μ depend on the level at which the band is cut. Our computation for the reference cloud: with ten percent instead of one, ΔJ falls from 0.270 to 0.202 and ΔΩ from 0.0192 to 0.0147, so τ₁ grows from 327 to 426 and μ by 32 percent.
- **Why it matters.** Statements such as "μ ≳ 0.4" or "the change lies between 0.43 and 0.48" hold only for this convention. Comparisons between families of the same shape are not affected.
- **Importance.** Medium.
- **Recommendation.** Say so in one sentence where the band is defined.

### 3.5 The heuristic uses ΔJ both as a width and as a distance

- **Location.** V.A, paragraph 3 ("the frequencies change by an amount of the order of A/ΔJ across the cloud").
- **Problem.** ΔJ plays two roles: the spread of the frequencies, and the range of actions over which the potential energy of the cloud varies. This assumes that the cloud reaches the circular orbit, so that its width is also its distance from it.
- **Why it matters.** This is the assumption that VII.B removes. Stating it in V.A prepares Section VII without announcing its result.
- **Importance.** Medium.
- **Recommendation.** One sentence: "For the cloud (42), which extends to the circular orbit, its width in J is also its distance from that orbit; for other shapes the two differ (section VII)."

### 3.6 The two criteria for the core disagree near the change

- **Location.** V.C, paragraphs 3 and 4 (a peak above F̂_max from μ = 0.37 on; m_c > m_mix); Table I; V.E, paragraph 1.
- **Problem.** At μ = 0.37 the angle average of the DF exceeds F̂_max at the circular orbit, while m_c = 0.040 < m_mix. For σ_J = 0.2 and μ = 0.38 the peak is narrower than the ellipse (m_c = 0.045), and yet a_c behaves as for a core (0.04 to 0.07). VI.C notes the case μ = 0.37 as well.
- **Why it matters.** The central object of the section has two definitions which disagree exactly where the change occurs, so its location depends on the criterion.
- **Importance.** High.
- **Recommendation.** Say which criterion defines the core and use the other as a check. Note that m_c > m_mix is sufficient but not necessary: a peak narrower than the ellipse is missed. Near the change, report both the height of the peak relative to F̂_max and m_c.

### 3.7 The value of m_c without self-gravity exceeds m_mix

- **Location.** Table I (rows M_gas/M = 0 and 10⁻³: m_c = 0.051 and 0.048); V.C, paragraph 4 ("compatible with m_mix within the fluctuations").
- **Problem.** Without self-gravity, where the bound holds exactly, m_c = 0.051 exceeds m_mix = 0.042 by 20 percent. Centering the ellipse at the maximum of a noisy density is expected to bias m_c upwards.
- **Why it matters.** The reference for a core should be the spread of m_c for a mixed gas, which Table I gives (0.031 to 0.051), not m_mix. The value 0.059 for σ_Q = 0.2 and μ = 0.42, cited in V.E as intermediate, is only 16 percent above the value without self-gravity.
- **Importance.** Medium.
- **Recommendation.** Draw the value without self-gravity next to m_mix in Fig. 5(b). Quantify the fluctuation with the existing data, for instance by the spread of m_c over t ≥ 20τ₁.

### 3.8 J_∗ and the density estimator

- **Location.** Eq. (46); V.C, paragraph 4.
- **Problem.** J_∗ = 0.05σ_J is not justified, and the method used to estimate the density and its maximum (bins or kernel, and their size) is not given. For σ_J = 0.2 and μ = 0.38 the core is narrower than the ellipse, so the core does not scale with σ_J as the ellipse does.
- **Why it matters.** m_c, m_mix and the location of the change depend on these choices, and a reader cannot reproduce m_c.
- **Importance.** Medium.
- **Recommendation.** Describe the density estimator. For the runs near the change, give m_c for two or three values of J_∗, which needs only the existing data.

### 3.9 a_c follows particles selected at the final time

- **Location.** V.D, paragraph 1 (definition of r_c and a_c).
- **Problem.** r_c(t) uses the particles which lie in the ellipse at the end of the run. At early times a_c therefore describes the part of the cloud which ends in the core, not the core.
- **Why it matters.** If the core exchanges particles with the halo, the time of the minimum and the regrowth could depend on the time of the selection. The text does not say whether this was checked.
- **Importance.** Medium.
- **Recommendation.** State the choice and its meaning. Check it with the existing data, with a selection at 20τ₁ or with the motion of the maximum (r_∗, p_∗).

### 3.10 Other conventions

- **Location.** V.C (mass with J < 0.2σ_J); captions of Table I and Fig. 5 (averages over t ≥ 20τ₁).
- **Problem.** The level 0.2σ_J is not justified. All the runs extend to t = 10⁴ (2 × 10⁴ for σ_J = 0.05), so "t ≥ 20τ₁" means 20 to 30.6τ₁ for the reference family, up to 50τ₁ for σ_J = 0.2 and up to 62τ₁ for the Kepler potential.
- **Why it matters.** Fig. 5(b) compares averages over different windows in units of τ₁.
- **Importance.** Low.
- **Recommendation.** Justify 0.2σ_J in one clause. Use 20 ≤ t/τ₁ ≤ 30 for all the families, or state the windows.

## 4. Figures and tables

Figs. 2, 4(a) and 5 and Table II carry the argument and belong in the main text. Fig. 4(b) duplicates Fig. 10, Table III repeats Fig. 5, and the values at μ = 0.43 and 0.48 appear in six places.

| Element | Question it answers | Assessment | Place |
| --- | --- | --- | --- |
| Fig. 2: particles in the (r, p_r)-plane for μ = 0.31, 0.43, 0.89 | What does the core look like, and when does it appear? | Essential: the only direct picture of the core | Main text |
| Fig. 3(a): envelope of \|h₁\|/h₀ | How does the decay change with μ? | Overlaps with Fig. 6; crowded | Main text, reduced |
| Fig. 3(b): deviation of the potential | Does the field relax? | Supports the analogue of Landau damping | Main text |
| Fig. 4(a): angle average of the DF at late times | Which state is reached, and where is F̂_max exceeded? | Essential for the bound | Main text |
| Fig. 4(b): a_c(t) for 0.43 ≤ μ ≤ 0.48 | Is the oscillation damped? | Same PIC curves as Fig. 10 | Once, with the model |
| Fig. 5: all the families as functions of μ | Does μ locate the change? | Main summary of V | Main text |
| Table I: reference family at late times | \|h₁\|, mass at J < 0.2σ_J and m_c | Only the J < 0.2σ_J column is in no figure | Main text, reduced |
| Table II: the six masses of the change | Is the change continuous? | Core of the claim | Main text |
| Table III: brackets for each family | Where does each family change? | Repeats points of Fig. 5 | Supplement, or brackets drawn in Fig. 5 |

### 4.1 Fig. 3 overlaps with Fig. 6, and panel (a) is crowded

- **Location.** Fig. 3.
- **Problem.** Panel (a) shows eight values of μ in close shades of blue and green, and μ = 0.43, 0.48 and 0.52 are hard to tell apart. Dotted lines for N_Q = 100 double some of the curves. Section VI analyzes the same decay on finer lattices (Fig. 6).
- **Why it matters.** The figure tries to show the outcome, the decay law and the dependence on the lattice at once.
- **Importance.** Medium.
- **Recommendation.** Keep panel (b) and reduce panel (a) to four or five values of μ with distinct colors. Move the comparison of lattices to Appendix B.4.

### 4.2 Fig. 4(b) duplicates Fig. 10

- **Location.** Fig. 4(b); Fig. 10 (VII.D).
- **Problem.** The six PIC curves of a_c for 0.43 ≤ μ ≤ 0.48 appear in both figures, and Fig. 10 adds the model. Fig. 4 also joins a panel on the late state with one on the oscillation, which belong to different subsections.
- **Why it matters.** The material is duplicated, and one figure serves two subsections.
- **Importance.** Medium.
- **Recommendation.** Show the curves once, in Fig. 10 with the model. If observation 1.1 is adopted, V.D keeps one sentence and Table II.

### 4.3 The same numbers in several places

- **Location.** Tables I, II and III; Table VIII; Figs. 5 and 9.
- **Problem.** The mean of |h₁|/h₀ and m_c for μ = 0.43 and 0.48 appear in Tables I, II and III, again in Table VIII (there with μ = 0.44) and in Figs. 5 and 9. Table III repeats points of Fig. 5.
- **Why it matters.** Every copy must be kept consistent, and a reader cannot tell which table is the reference.
- **Importance.** Medium.
- **Recommendation.** Keep Table II as the table of the change. Reduce Table I to what no figure shows (the mass with J < 0.2σ_J), or merge it into Fig. 9. Move Table III to the supplement and draw its brackets in Fig. 5.

### 4.4 Readability of Figs. 2, 4(a) and 5(b)

- **Location.** Figs. 2, 4(a) and 5(b).
- **Problem.** Fig. 2 marks neither the circular orbit nor the ellipse (46), and its zoom panels use different axes in each row. Fig. 4(a) uses J instead of J/σ_J, marks neither 0.2σ_J nor J_∗, and has no error bar for the innermost interval, where the text says that μ = 0.66 reaches F_max "within the sampling error". Fig. 5(b) shows m_mix but not the value without self-gravity.
- **Why it matters.** The reader has to reconstruct the quantities the text refers to.
- **Importance.** Low.
- **Recommendation.** Mark the circular orbit (R₀ = 5.742, given in VII.A) and the ellipse in the zoom panels. Use J/σ_J in Fig. 4(a), with marks at 0.2σ_J and J_∗ and an error bar for the innermost interval. Add the value of m_c without self-gravity to Fig. 5(b) (see 3.7).

Candidates for the supplementary material: Table III; the comparison of lattices in Fig. 3(a); the details of a_c(t) for μ ≥ 0.47 (second minimum at 14 to 15τ₁, maximum at 22 to 25τ₁), which Appendix B.4 shows to depend on the lattice; and the full list of convergence runs of V.A.

## 5. Numerical validation

The code is well validated without self-gravity, and the reference family was repeated with N_Q, N_J, Δt and Δr changed. Missing are a test beyond N_Q = 100 near the change, finer lattices for the other families, and the window in which each observable is reliable.

| Test | Runs (M_gas/M) | Result reported | Where |
| --- | --- | --- | --- |
| Free mixing against the exact solution | No self-gravity | 6 × 10⁻⁹ relative error over 10⁴ time units | B.1 |
| Action-angle variables of the total potential | 10⁻⁴, 10⁻³, 10⁻² | Residue 2.2 × 10⁻³ with the external variables, 1.5 × 10⁻⁶ with the instantaneous ones (10⁻³) | B.2 |
| Cold collapse | No external potential | Second order in 1/N | B.3 |
| Energy | All the runs of V | At most 1.7 × 10⁻⁴, median 1.0 × 10⁻⁵ | B.3 |
| N_Q from 25 to 100 | 0.01, 0.02, 0.1 | Mean \|h₁\|/h₀: 2.62 to 2.44 × 10⁻³ (μ = 0.31), 2.1 to 1.2 × 10⁻² (0.43), 0.80 to 0.85 (0.89); m_c: 0.54 to 0.52 (0.89) | V.A, B.4 |
| N_J from 400 to 800 | 0.01, 0.02, 0.025 | Same \|h₁\| as 400 × 25 at μ = 0.43; 2.44 × 10⁻³ at 0.31 | B.4 |
| 3200 × 100 points | 0.02, 0.025 | Same as 400 × 100 at μ = 0.43; at 0.48, m_c within 6 percent and \|h₁\| within 13 percent, a_c up to a factor 4 apart | B.4 |
| Δt halved | 0.02, 0.025 | Unchanged | V.D, B.4 |
| Δr halved | 0.02, 0.025 | Final a_c at μ = 0.43 from 0.045–0.067 to 0.036–0.050; m_c from 0.103 to 0.106 on all lattices | V.D, B.4 |

### 5.1 Energy conservation is not a convergence test

- **Location.** IV.C, last sentence; Appendix B.3, paragraph 2.
- **Problem.** With a symplectic integrator and a fixed Δt, the discrete energy is conserved even when Δr or the number of particles is too small. The test checks the time integration and the implementation, not the convergence of the observables.
- **Why it matters.** The manuscript keeps the two apart (B.3 and B.4), and Section V does not use the energy as evidence. IV.C, however, lists the energy among the validations without saying what it does not test.
- **Importance.** Low.
- **Recommendation.** Keep the separation, and say so in one clause in IV.C.

### 5.2 No convergence in N_Q near the change

- **Location.** V.A (list of convergence runs); Appendix B.4; Table II.
- **Problem.** Only N_Q = 25 and 100 were used. At μ = 0.43 the mean of |h₁|/h₀ changes by a factor 1.75 between them, as large as the largest factor between consecutive masses in Table II (1.37 to 1.75 per 0.01 in μ); ‖δΦ‖ changes from 2.9 to 2.4 × 10⁻⁴. The test of VI.A with 6400 × 200 points is at Γ = 0.1 (μ ≈ 0.29), outside the change.
- **Why it matters.** The values of |h₁| in Table II could shift by about one step in μ. The continuity of the change would not be affected, its location could be.
- **Importance.** High.
- **Recommendation.** One run with N_Q = 200 at μ = 0.43 or 0.45 (a new simulation). Until then, state the uncertainty of the values in Table II.

### 5.3 The other families only on coarse lattices

- **Location.** V.A (lattices of the families); Table III; Fig. 5.
- **Problem.** σ_J = 0.05, σ_Q = 0.4 and the Kepler potential use 400 × 25 points (N_eff = 243); σ_Q = 0.2 uses N_Q = 50 and σ_J = 0.2 uses 800 × 25 (N_eff = 485). Near the change |h₁| depends on the number of particles (B.4), and the brackets of Table III rest on these runs.
- **Why it matters.** The comparison between families, which supports the use of μ, is made at a resolution at which the reference family is not converged.
- **Importance.** Medium.
- **Recommendation.** Repeat the two bracketing runs of each family with N_Q = 100 (new simulations), or state the caveat in V.E.

### 5.4 Δr and Δt

- **Location.** V.D, paragraph 2; Appendix B.3 (note) and B.4.
- **Problem.** Δr was tested with two values only, and the final amplitude at μ = 0.43 changes by 20 to 25 percent. The note in B.3 asks for a test of Δt with M_gas/M = 0.01 and 0.1, while V.A and B.4 already report Δt halved for 0.02 and 0.025.
- **Why it matters.** With two values of Δr the convergence of a_c cannot be assessed, and the note and the text disagree on what has been done.
- **Importance.** Medium.
- **Recommendation.** Give the dependence on Δr as an uncertainty of a_c. Update the note in B.3 to say what was done and what remains.

### 5.5 The window in which each observable is reliable

- **Location.** V.A and V.B; Appendices B.1 and B.2; VI.A.
- **Problem.** The text does not say until when each observable is reliable on each lattice. VI.A says that on the lattices of V the higher harmonics are dominated by discreteness from about 8τ₁ on. V.B attributes the late growth of ‖δΦ‖ without self-gravity to the recurrence of the harmonic n = 10 at about 15τ₁, and B.2 shows fluctuations from discreteness which grow with time and with the mass.
- **Why it matters.** The windows of the means (10 to 20τ₁, t ≥ 20τ₁, t ≥ 25τ₁) lie partly beyond these times.
- **Importance.** Medium.
- **Recommendation.** One sentence, or a column in a table of runs: the reliable window of |h₁|, ‖δΦ‖ and a_c for each lattice.

### 5.6 Sensitivity to small perturbations near the change

- **Location.** Appendix B.4 ("the outcome is sensitive to small perturbations"); V.D, paragraph 3.
- **Problem.** At μ = 0.48, a_c differs by up to a factor four between lattices, and there is no estimate of the spread for a fixed lattice. V.D still describes details for μ ≥ 0.47: a second minimum at 14 to 15τ₁ and a maximum at 22 to 25τ₁.
- **Why it matters.** These details lie below the reproducibility shown in B.4.
- **Importance.** Medium.
- **Recommendation.** Restrict the statements on a_c(t) for μ ≥ 0.46 to coarse ones. For an error bar, a small ensemble would do, for instance with the lattice shifted by half a cell in Q (new simulations).

### 5.7 Action variables near the circular orbit

- **Location.** Appendix A; Fig. 4(a), innermost interval (J below 10⁻³).
- **Problem.** Appendix A reports the accuracy of the action at the points of the lattice, in the external potential. The accuracy in the instantaneous potential near the circular orbit, where the core sits and the turning points merge, is not reported.
- **Why it matters.** The peak above F̂_max, and the value F_max in the innermost interval, depend on it.
- **Importance.** Low.
- **Recommendation.** One sentence with the error of J below 10⁻³ in the instantaneous potential, if it has been measured.

## 6. Connection with Sections VI to VIII

Section V works with μ, Section VI labels its runs with Γ, and Section VII defines Γ and relates it to μ. Γ can be introduced in V.A with quantities that V already uses, which removes the forward references and shows at once how Γ differs from μ.

### 6.1 Γ is used before its definition

- **Location.** V.A, paragraph 3 (Γ described, not named); V.E, paragraph 2 (named); VI.A (labels of the runs, with a note that Γ will move to V.A); VII.B, Eq. (57).
- **Problem.** Γ appears in V.E and throughout VI before Eq. (57) defines it. Its definition uses g, which comes from the averaged interaction of the shells in VII.A, and the median action J_med.
- **Why it matters.** A reader of VI cannot interpret the labels Γ = 0.025 to 0.15, and V.E cites a parameter that it has not given.
- **Importance.** High.
- **Recommendation.** Define Γ in V.A next to μ. Eqs. (57) and (58) give an identity that needs only the amplitude of the potential of a compact cloud with action J, A(J) = π g M_gas (2J)^(1/2), which is Eq. (58) with σ_J replaced by J:

```latex
\Gamma = \frac{A(J_{med})}{\sqrt{2}\,\pi\,|\Omega_0'|\,J_{med}^2},
\qquad
\mu^2 = \frac{A}{\Delta\Omega\,\Delta J} \approx \frac{A(\sigma_J)}{|\Omega_0'|\,\Delta J^2}.
```

Both compare the amplitude of the potential with the differential rotation over a squared range of actions. μ uses the width ΔJ of the band; Γ uses the median action J_med, that is, the distance of the cloud from the circular orbit. The approximation for μ² is ΔΩ ≈ |Ω₀′|ΔJ, as in VII.B.

g can be given with a pointer to VII.A. Section V can then state its results in Γ and keep μ for one remark (its relation with τ₁ and with trapping, see 3.3), and VI needs no forward reference. The identity follows from Eqs. (57) and (58) and changes no definition; it was already given in the chat of 17 November.

### 6.2 One value of μ for the whole article

- **Location.** VI.A (μ from 0.14 to 0.35), VI.C (μ = 0.37 with Γ = 0.17), VIII (μ = 0.31, 0.43, 0.52 and 0.66), VII.B and VII.E, the abstract and Section IX ("μ ≳ 0.4", "enter only through μ").
- **Problem.** VI and VIII use the values of V, VII.E those of the passive cloud, and Eq. (59) a third set (see 3.2). Section IX says that the mass, the width and the external potential "enter only through μ", which holds to leading order and, with the present definition, for M_gas/M ≤ 0.025 (see 3.1).
- **Why it matters.** Any change of the definition in V propagates to five places.
- **Importance.** Medium.
- **Recommendation.** After deciding 3.1 and 3.2, update VI, VII.B, VII.E, VIII, the abstract and IX together. Add "to leading order" in IX.

### 6.3 The transition from V to VI

- **Location.** End of V.E; opening paragraph of VI and VI.A.
- **Problem.** V ends with the Kepler clouds. VI opens with runs labeled by Γ and by four shapes that VII.E defines. The question of VI, how fast the clouds that relax mix, is not posed at the end of V.
- **Why it matters.** The reader moves from the location of the change to tail exponents without a bridge, and with two undefined ingredients.
- **Importance.** Medium.
- **Recommendation.** Close V with the summary of 1.7, ending with the two questions that follow: the rate of mixing below the change (VI) and the mechanism of the core (VII). Introduce the shapes where VI first needs them, moving the list of VII.E to VI.A or to V.A. Keep the comparison at equal μ and at equal Γ as the result of VII.E.

### 6.4 The range of Γ for a core is quoted in two ways

- **Location.** VI.A ("a core forms when Γ is between 0.23 and 0.29"); VII.E, the introduction and Section IX ("m_c reaches 0.1 at values of Γ between 0.24 and 0.29").
- **Problem.** The value 0.23 matches the PIC run with M_gas/M = 0.02 (Γ = 0.229, m_c = 0.103); 0.24 to 0.29 comes from the model.
- **Why it matters.** One statement has two ranges.
- **Importance.** Low.
- **Recommendation.** Use one range, and say whether it comes from the model or from the PIC runs.

## Questions the manuscript cannot settle

Nine questions raised above cannot be answered from the manuscript. Six of them need only the data already computed.

| Question | Why the manuscript cannot settle it | What would settle it | New runs |
| --- | --- | --- | --- |
| Does the oscillation of the core decay beyond 30τ₁? | PIC runs end at 30τ₁; the model is limited by the number of particles beyond 60τ₁ (VII.D) | Longer runs with more particles | Yes |
| Is N_Q = 100 converged near the change? | Only N_Q = 25 and 100 were used (5.2) | N_Q = 200 at μ = 0.43 or 0.45 | Yes |
| Do the brackets of the other families hold at N_Q = 100? | They rest on 400 × 25 or 800 × 25 points (5.3) | The two bracketing runs of each family | Yes |
| Which values of m_c does a mixed gas give with this estimator? | Only one run without self-gravity is reported (3.7) | Spread of m_c over time, or over shifted lattices | No |
| How does m_c depend on J_∗? | A single value, 0.05σ_J, is used (3.8) | m_c for two or three values of J_∗ near the change | No |
| Is the peak above F̂_max at μ = 0.37 a compact core or a narrow accumulation? | Fig. 2 shows μ = 0.31 and 0.43, not 0.37 (3.6) | A zoom of the (r, p_r)-plane at μ = 0.37 | No |
| Does a_c depend on when the core particles are selected? | Not reported (3.9) | a_c with a selection at 20τ₁, or the motion of the density maximum | No |
| Why is the mean of \|h₁\|/h₀ non-monotonic above the change? | Not discussed (2.2) | \|h₁\| of the gas outside the core, compared with a_c | No |
| How accurate are the action variables of the instantaneous potential near J = 0? | Appendix A reports the error only at the lattice points of the external potential (5.7) | A direct test below J = 10⁻³ | No |
