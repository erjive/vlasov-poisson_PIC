# Introducción al problema

`vlasov_intro.tex`: notas didácticas sobre mecánica hamiltoniana, variables
ángulo-acción, phase mixing y amortiguamiento de Landau en Vlasov–Poisson con
simetría esférica, para alguien con formación en física o matemáticas que se
inicie en el tema. 42 páginas, 20 figuras.

Se compila con `latexmk -pdf vlasov_intro.tex` (solo paquetes de TeX Live
estándar).

## Figuras

Las figuras en `figuras/*.pdf` se versionan para que el documento compile sin
nada más. Se regeneran con

    cd figuras && python3 generar_figuras.py          # todas
    cd figuras && python3 generar_figuras.py fase     # solo una

Las seis primeras (Liouville, potencial, frecuencias, órbita, enrollamiento, h_k exacto)
son analíticas. Las demás leen las corridas de `exe/dfstudy`, `exe/dfsnap` y
`exe/sg`, que no se versionan; se regeneran con `paper_runs/scripts/run_*.sh`
(las capturas de la sección 8.2, con `run_dfsnap.sh`).

Cada figura imprime en la terminal las cifras que cita el texto, para poder
comprobar que documento y datos coinciden.
