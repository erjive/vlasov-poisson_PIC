"""lambda(0) y los autovalores siguientes del lazo autoconsistente en los equilibrios de L fijo.

Dos comprobaciones del criterio lambda_borde > 1 (lambda_borde.py):

  - lambda(0) < 1. Con frecuencia imaginaria, omega = i sigma, |R_k| decrece con sigma, así que
    hay un modo creciente si y solo si lambda(0) > 1: lambda(0) < 1 es la estabilidad lineal.
  - El segundo autovalor de |R|^{1/2} (-M) |R|^{1/2} queda por debajo de 1 hasta el borde de la
    banda, de modo que hay a lo sumo un modo discreto.

Recorre los equilibrios de masa puntual del informe de Hadžić (exe/hadzic/lineal) y los de
isócrono de la batería y la demo de eta (exe/eta_lineal) con F_eq monótona.

    python3 lambda_cero.py

Escribe exe/eta_lineal/lambda_cero.txt. Tarda unos 10 minutos con cuatro hilos.
"""
import os, sys, glob, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from lambda_borde import Lazo

BASE = os.path.join(AQUI, '..', '..', 'exe')


def tres(z, w):
    """Los tres mayores autovalores del operador simétrico a la frecuencia w."""
    h = np.sqrt(np.abs(z.R(w)))
    ev = np.linalg.eigvalsh(h[:, None]*(-z.M)*h[None, :])
    return ev[-1], ev[-2], ev[-3]


if __name__ == '__main__':
    casos = [os.path.join(BASE, 'hadzic', 'lineal', f'P_k{k:g}_a{a0:g}_equilibrio.npz')
             for k in (0.75, 1, 1.25, 1.5, 2, 3) for a0 in (0.01, 0.1, 0.6, 1)]
    casos += sorted(glob.glob(os.path.join(BASE, 'eta_lineal', '[ADLMG]*_equilibrio.npz')))
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('delta en anchos de banda bajo Omega_min; l1, l2, l3: los tres mayores autovalores.')
    w(f'{"caso":>22} {"a0":>6} {"lam(0)":>8} {"l1(1e-3)":>9} {"l2(1e-3)":>9} {"l3(1e-3)":>9} '
      f'{"l1(borde)":>9} {"l2(borde)":>9}')
    for npz in casos:
        if not os.path.exists(npz):
            w(f'falta {npz}'); continue
        z = Lazo(npz)
        if not z.monotona:
            continue
        an = z.Om_max - z.Om_min
        a, b = tres(z, z.Om_min - 1e-3*an), tres(z, z.Om_min)
        w(f'{os.path.basename(npz)[:-15]:>22} {z.a0:6g} {z.lam(0.0):8.4f} {a[0]:9.4f} {a[1]:9.4f} '
          f'{a[2]:9.4f} {b[0]:9.4f} {b[1]:9.4f}')
    open(os.path.join(BASE, 'eta_lineal', 'lambda_cero.txt'), 'w').write('\n'.join(lineas) + '\n')
