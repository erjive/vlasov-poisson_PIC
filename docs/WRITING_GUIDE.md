# Writing guide for reports

### Language

* Write all documents in English.

### General approach

Write scientific reports in a clear, direct, and pedagogical style, similar to a good research paper in computational physics.

The goal is to explain the mathematics, numerical method, and physical interpretation carefully, but without unnecessary rhetorical language.

### Writing style

* Prefer simple, direct sentences.
* Explain ideas step by step, but do not announce every step with rhetorical phrases.
* Use standard scientific prose rather than essay-like or literary language.
* Be detailed when the mathematics or numerical method requires it, but do not expand simple ideas unnecessarily.
* Avoid excessive sectioning and avoid turning every small idea into a subsection.
* Unless they are genuinely necessary, do not use phrases such as:

  * "This section builds the necessary machinery..."
  * "The ideas, step by step..."
  * "We now embark on..."
  * "Let us unpack..."
  * "The key conceptual point..."
  * "This provides the foundation for..."

* Avoid dramatic, promotional, or overly polished language.
* Do not restate the same idea in several different ways.
* Do not introduce a section with a paragraph explaining what the section is going to do if the title already makes that clear.

### Preferred style

Instead of:

"1. The ideas, step by step
This section builds, in six steps, what is necessary to understand the demo."

Write simply:

"## 1. Model and numerical method"

Then start directly with the content:

"We consider the spherically symmetric Vlasov–Poisson system with fixed angular momentum \(L_0\). The distribution function is written as ..."

The report should read like a research paper: mathematically precise, easy to follow, and sufficiently detailed for a reader to reproduce the argument or understand the simulation.

### Level of detail

Be comprehensive about:

* assumptions,
* definitions,
* equations,
* derivations that are important for correctness,
* numerical discretization,
* boundary conditions,
* diagnostics,
* validation tests,
* interpretation of the results.

Be concise about:

* obvious transitions,
* summaries of material already explained,
* basic observations that do not require justification.

When explaining an equation, prioritize:

1. what the equation means,
2. where it comes from,
3. how it is used in the numerical method.

Do not add prose merely to make the text sound sophisticated.

### Presenting numerical results

Added 22 November 2026. The reader is a physicist with a general knowledge of classical mechanics, statistical mechanics and numerical methods. This reader must be able to understand each experiment, identify its relevant results and grasp their physical meaning without opening the data files or the code.

Physical interpretation has priority over the enumeration of numbers. Do not present long lists of values, matrices, vectors or statistics, and do not give a number before explaining what it represents and why it matters.

Each result must answer:

1. What was computed?
2. What behavior is observed?
3. What does it mean physically?
4. What can be concluded?
5. What are the limits of this conclusion?

**Organization of an experiment.** Present each experiment in this order:

* **Objective.** The phenomenon or property under study.
* **Numerical setup.** Initial data, physical parameters, spatial and temporal resolution, number of particles and method. Use a table only when it helps to compare parameters.
* **Main results.** Prefer figures to tables. For each figure, say what the axes show, describe the behavior, identify the relevant trends, compare times or configurations where appropriate, and state what physical information it provides.
* **Physical interpretation.** Use the theory of the Vlasov-Poisson system. Separate numerical observations, physical interpretations and hypotheses.
* **Numerical validation.** Give the relevant conservation, convergence and stability tests. Say what the errors indicate and whether they are compatible with the expected accuracy of the method.
* **Conclusions.** The main findings, the limitations and the next experiments.

In a stand-alone report these six items are headings. In a manuscript section they fix the order of the argument and do not become headings, so that the text is not oversectioned.

**Use of numbers.**

* Use figures to show trends, small tables to compare configurations, and representative values to support a claim.
* Put complete data and detailed results in appendices.
* Give physical units where they apply. In the dimensionless variables of the article, state the units once.
* Do not invent data, interpretations or explanations that the available results do not support.

**Writing.** Avoid very long paragraphs, unnecessary enumerations of results, repetition, convoluted terminology, explanations without physical content and grandiose claims. Explain technical concepts when needed, without simplifying the physics incorrectly.

**Scientific judgment.** A visually attractive figure is not evidence of a physical result. Always distinguish:

1. what the data show directly,
2. what can be interpreted physically,
3. what still requires a numerical or theoretical check.

If the results do not allow a conclusion, say so explicitly.

### Overall rule

Write as an expert researcher explaining the work to another researcher who is intelligent but has not seen the code before.

Clarity is more important than elegance.
Precision is more important than rhetorical sophistication.
Detail is welcome when it improves understanding; verbosity is not.

Before finalizing a report, remove unnecessary introductory phrases, redundant explanations, and decorative language.
