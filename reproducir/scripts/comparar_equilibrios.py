"""Compara los equilibrios de la demo eta construidos de dos maneras: la iteración
de punto fijo de equilibrio.py (los _equilibrio.npz de exe/demo_eta/ic/) y la
ecuación diferencial de equilibrio_edo.py.

    python3 comparar_equilibrios.py            tabla de diferencias por equilibrio
    python3 comparar_equilibrios.py --refinar  además, A4 con la iteración en una
                                               malla más fina, para ver cuál converge
"""
import os, sys, time, argparse, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from equilibrio import Equilibrio
from equilibrio_edo import EquilibrioEDO

IC = os.path.join(AQUI, '..', '..', 'exe', 'demo_eta', 'ic')
CASOS = [('A1', 'D1'), ('A3', 'D2'), ('A4', 'D5'), ('L5', 'D3'), ('L6', 'D4'),
         ('G1a', 'D8'), ('M5', 'D9')]


def picard(npz):
    d = np.load(npz)
    jt, w0 = float(d['J_borde']), float(d['w0'])
    E_t = float(np.interp(jt, d['J_t'], d['E_t']))
    return dict(r=d['r'], phi=d['phi_self'], rho=d['rho'], A=float(d['A']), a0=float(d['a0']),
                jt=jt, g=float(d['k_borde']), w0=w0, E_t=E_t, T=(E_t - float(d['E_t'][0]))/w0)


def diferencias(p, e):
    phi_e, rho_e = e.phi_self(p['r']), e.rho(p['r'])
    escala = np.max(np.abs(p['phi']))
    return dict(dphi=np.max(np.abs(phi_e - p['phi']))/escala,
                drho=np.max(np.abs(rho_e - p['rho']))/np.max(p['rho']),
                dEt=abs(e.E_t - p['E_t'])/abs(p['E_t']), dT=abs(e.T - p['T'])/p['T'],
                dA=abs(e.A - p['A'])/p['A'])


class EquilibrioFino(Equilibrio):
    """La misma iteración con más puntos en E y en p_r; la malla radial se pasa aparte."""
    def tabla_J_de_E(self, m, n=8000):
        return super().tabla_J_de_E(m, n)

    def densidad(self, m, E_t, J_t, npm=2401):
        return super().densidad(m, E_t, J_t, npm)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--refinar', action='store_true')
    a = ap.parse_args()
    print('Diferencias relativas, EDO frente a la iteración de punto fijo:')
    print(f'{"caso":5} {"a0":>7}  {"Phi_self":>8} {"rho":>8} {"E_t":>8} {"T":>8} {"A":>8}   tiempo EDO')
    for caso, corrida in CASOS:
        p = picard(os.path.join(IC, f'{corrida}_equilibrio.npz'))
        t0 = time.time()
        e = EquilibrioEDO(p['a0'], p['jt'], p['g'], p['w0']).resolver()
        dt = time.time() - t0
        d = diferencias(p, e)
        print(f'{caso:5} {p["a0"]:7.4f}  {d["dphi"]:8.1e} {d["drho"]:8.1e} {d["dEt"]:8.1e} '
              f'{d["dT"]:8.1e} {d["dA"]:8.1e}   {dt:.1f} s  (residuo {e.residuo:.0e})')
    # La EDO contra sí misma, con tolerancias y nodos más laxos.
    p = picard(os.path.join(IC, 'D5_equilibrio.npz'))
    e1 = EquilibrioEDO(p['a0'], p['jt'], p['g'], p['w0']).resolver()
    e2 = EquilibrioEDO(p['a0'], p['jt'], p['g'], p['w0'], rtol=1e-10, nodos=100).resolver()
    print(f'EDO en A4 con rtol 1e-10 y 100 nodos frente a 1e-12 y 200: Phi_self '
          f'{np.max(np.abs(e1.phi_self(p["r"]) - e2.phi_self(p["r"])))/np.max(np.abs(p["phi"])):.1e}')
    if a.refinar:
        print('A4 con la iteración en mallas más finas (r cada h, 8000 energías, 2401 puntos en p_r):')
        for h in (0.01, 0.005, 0.0025):
            r = np.arange(h, 25.0 + 1e-9, h)
            t0 = time.time()
            q = EquilibrioFino(p['a0'], r_malla=r, jmax=p['jt'], forma='maxwell', jt=p['jt'],
                               k=p['g'], w0=p['w0']).iterar(verboso=False)
            dif = np.max(np.abs(e1.phi_self(r) - q.phi_self))/np.max(np.abs(q.phi_self))
            print(f'  h = {h:<7} Phi_self: {dif:.1e}   ({time.time()-t0:.0f} s, q = {q.q:.3f})')
