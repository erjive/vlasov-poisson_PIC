"""Componente de h_1 que gira, a tiempos largos, frente a la resolución de la cuadratura.

La parte estática S de h_1 es un efecto del mapa del isócrono (aa_meseta.py) y no
depende de la resolución; la parte que gira es la misma con cualquier mapa, así que
se mide directamente sobre el h_1 que escribe el código: R = h_1 - <h_1> por ventana.

Si R fuera un piso de discreción de la rejilla, bajaría al duplicar las filas en J
(Nrc) o los nodos en Q (Npc). Si es una propiedad resuelta de la distribución, no.

Uso:  python3 giro_resolucion.py                      (a0=1e-3, las corridas de abajo)
      python3 giro_resolucion.py base otra [otra ...]  (la primera es la referencia)
"""
import os, sys, numpy as np

SG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'sg')
CORRIDAS = [('Nrc=400, Npc=25 (N=1e4)', 'long20k_snap'),
            ('Nrc=800, Npc=25 (N=2e4)', 'long20k_nrc800'),
            ('Nrc=400, Npc=50 (N=2e4)', 'long20k_npc50'),
            ('sin autogravedad, N=1e4', 'long20k_nosg')]
VENTANAS = [(1600, 2000), (2000, 6000), (6000, 12000), (12000, 20000), (15000, 20000)]
if len(sys.argv) > 2:
    CORRIDAS = [(d, d) for d in sys.argv[1:]]


def leer(d):
    a = np.loadtxt(os.path.join(SG, d, 'hk1_complex.tl'))
    return a[:, 0], (a[:, 3] + 1j*a[:, 4])/a[0, 1]


base = {}
print(f"{'corrida':>26} {'ventana':>15} {'|S|':>10} {'rms giro':>10} {'giro/giro(N=1e4)':>17} {'corr. con N=1e4':>15}")
for nombre, d in CORRIDAS:
    if not os.path.exists(os.path.join(SG, d, 'hk1_complex.tl')):
        print(f'{nombre:>26}   (no existe {d})')
        continue
    t, z = leer(d)
    for lo, hi in VENTANAS:
        v = (t >= lo) & (t <= hi)
        R = z[v] - z[v].mean()
        r = np.sqrt(np.mean(np.abs(R)**2))
        if d == CORRIDAS[0][1]:
            base[(lo, hi)] = (r, R)
            cociente, corr = 1.0, 1.0
        else:
            rb, Rb = base[(lo, hi)]
            cociente = r/rb
            corr = abs(np.vdot(R, Rb))/np.linalg.norm(R)/np.linalg.norm(Rb)
        print(f'{nombre:>26} {f"[{lo},{hi}]":>15} {abs(z[v].mean()):>10.3e} {r:>10.3e}'
              f' {cociente:>17.3f} {corr:>15.3f}')
