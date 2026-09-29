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

### Overall rule

Write as an expert researcher explaining the work to another researcher who is intelligent but has not seen the code before.

Clarity is more important than elegance.
Precision is more important than rhetorical sophistication.
Detail is welcome when it improves understanding; verbosity is not.

Before finalizing a report, remove unnecessary introductory phrases, redundant explanations, and decorative language.
