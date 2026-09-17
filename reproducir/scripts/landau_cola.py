"""Frecuencia y tasa de amortiguamiento de la cola colectiva de h_1.

La respuesta a la perturbación se aísla restando la corrida sin perturbación
(misma rejilla, mismo desajuste entre el Poisson discreto y el continuo):

    h1_pert(t)   = [h1(eps) - h1(0)]/eps
    dPhi_pert    = [dPhi(eps) - dPhi(0)]/eps

y se compara con el phase mixing puro en el potencial de equilibrio
(landau_libre.py). En la ventana donde la cola domina, se ajusta un modo
amortiguado h ~ A exp(-i w t - g t) de dos maneras independientes:

  1. pendientes: ln|h| frente a t (da -g) y fase desenrollada frente a t (da -w);
  2. matrix pencil (Hua y Sarkar 1990) con pocos términos: polos complejos
     z = exp(s dt), s = -g - i w, del modo de mayor amplitud.

Uso:  python3 landau_cola.py <prefijo> <nrc> [<nrc> ...]
      (lee exe/landau/<prefijo>_n<nrc>_e0 y _e0.1, con landau.npz y libre.npz)
"""
import os, sys, numpy as np

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'landau')
EPS = 0.1


def matrix_pencil(y, dt, M):
    """Devuelve polos s (y amplitudes) de y(t_n) ~ sum_m a_m exp(s_m t_n)."""
    N = len(y); L = N//2
    Y = np.array([y[i:i+L+1] for i in range(N-L)])
    U, S, Vh = np.linalg.svd(Y, full_matrices=False)
    # Las filas de Y son combinaciones de [1, z, z^2, ...]: el subespacio de la
    # señal lo generan las filas de Vh tal cual, sin conjugar.
    V = Vh[:M].T
    V1, V2 = V[:-1], V[1:]
    z = np.linalg.eigvals(np.linalg.pinv(V1) @ V2)
    s = np.log(z)/dt
    E = np.exp(np.outer(np.arange(N)*dt, s))
    a, *_ = np.linalg.lstsq(E, y, rcond=None)
    return s, a


def cargar(prefijo, nrc):
    d0 = np.load(os.path.join(BASE, f'{prefijo}_n{nrc}_e0', 'landau.npz'))
    d1 = np.load(os.path.join(BASE, f'{prefijo}_n{nrc}_e{EPS}', 'landau.npz'))
    lb = np.load(os.path.join(BASE, f'{prefijo}_n{nrc}_e{EPS}', 'libre.npz'))
    t = d1['t']
    h1 = (d1['hk'][:, 1] - d0['hk'][:, 1])/EPS
    h1_sin_restar = d1['hk'][:, 1]/EPS
    libre = lb['hk'][:, 1]/EPS
    dphi = (d1['dphi'] - d0['dphi'])/EPS
    return t, h1, h1_sin_restar, libre, dphi, d1['r'], d0['hk'][:, 1]


def ajustar(t, h, lo, hi):
    v = (t >= lo) & (t <= hi)
    tv, hv = t[v], h[v]
    g = -np.polyfit(tv, np.log(np.abs(hv)), 1)[0]
    w = -np.polyfit(tv, np.unwrap(np.angle(hv)), 1)[0]
    res = {'pendientes': (w, g)}
    for M in (1, 2, 3):
        s, a = matrix_pencil(hv, tv[1]-tv[0], M)
        k = np.argmax(np.abs(a))
        res[f'pencil M={M}'] = (-s[k].imag, -s[k].real)
    return res


if __name__ == '__main__':
    prefijo = sys.argv[1]
    for nrc in [int(x) for x in sys.argv[2:]]:
        t, h1, h1s, libre, dphi, r, h1_0 = cargar(prefijo, nrc)
        print(f'\n=== {prefijo}, Nrc={nrc}')
        print(f"{'t':>6} {'|h1_pert|':>11} {'|h1 libre|':>11} {'cociente':>9} {'|h1(e=0)|/eps':>14}")
        for tt in [0, 400, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000]:
            n = np.argmin(np.abs(t - tt))
            print(f'{t[n]:>6.0f} {abs(h1[n]):>11.3e} {abs(libre[n]):>11.3e}'
                  f' {abs(h1[n])/max(abs(libre[n]),1e-300):>9.2f} {abs(h1_0[n])/EPS:>14.2e}')
        for lo, hi in [(1000, 1600), (1100, 1800), (1200, 2000)]:
            res = ajustar(t, h1, lo, hi)
            print(f'  ventana [{lo},{hi}]: ' + '   '.join(f'{k}: w={w:.5f} g={g:.2e}' for k, (w, g) in res.items()))
