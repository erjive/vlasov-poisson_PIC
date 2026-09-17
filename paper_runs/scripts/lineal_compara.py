"""Compara la teoría lineal (lineal.py) con la simulación y con el phase mixing libre.

Uso:  python3 lineal_compara.py <prefijo PIC> <nrc PIC> <lineal.npz> [<lineal.npz> ...]
      (el PIC se lee de exe/landau/<prefijo>_n<nrc>_e0 y _e0.1)
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from landau_cola import cargar, ajustar

prefijo, nrc = sys.argv[1], int(sys.argv[2])
t, h1_pic, _, libre, _, _, _ = cargar(prefijo, nrc)
lin = [(f, np.load(f)) for f in sys.argv[3:]]

print(f"{'t':>6} {'PIC':>11} {'libre':>11}" + ''.join(f' {os.path.basename(f)[:22]:>24}' for f, _ in lin))
for tt in [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]:
    n = np.argmin(np.abs(t - tt))
    fila = f'{t[n]:>6.0f} {abs(h1_pic[n]):>11.4e} {abs(libre[n]):>11.4e}'
    for _, d in lin:
        m = np.argmin(np.abs(d['t'] - tt))
        fila += f' {abs(d["h1"][m]):>24.4e}'
    print(fila)

print('\nDiferencia relativa |h1_lineal - h1_PIC| / |h1_PIC| (complejos, fase incluida):')
for f, d in lin:
    hl = np.interp(t, d['t'], d['h1'].real) + 1j*np.interp(t, d['t'], d['h1'].imag)
    for lo, hi in [(0, 800), (800, 1400), (1400, 2000)]:
        v = (t >= lo) & (t <= hi)
        print(f'  {os.path.basename(f)}: t in [{lo},{hi}]: mediana {np.median(np.abs(hl[v]-h1_pic[v])/np.abs(h1_pic[v])):.2e}')

print('\nAjuste de la cola (ventana [1100,1800]):')
r = ajustar(t, h1_pic, 1100, 1800)
print(f"  PIC Nrc={nrc}: " + '   '.join(f'{k}: w={w:.5f} g={g:.3e}' for k, (w, g) in r.items()))
for f, d in lin:
    r = ajustar(d['t'], d['h1'], 1100, 1800)
    print(f"  {os.path.basename(f)}: " + '   '.join(f'{k}: w={w:.5f} g={g:.3e}' for k, (w, g) in r.items()))
