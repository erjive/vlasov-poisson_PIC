"""Colas de la mezcla de fases libre para cada forma de F_eq, sin simular.

Sin autogravedad la solución lineal es exacta: con F = F_eq(J) [1 + eps s(J) cos Q],

    dF(Q, J, t) = eps F_eq(J) s(J) cos(Q - Omega(J) t),

así que los dos observables son cuadraturas en J, sin evolución en el tiempo:

    h_1(t)    ∝ ∫ dJ B(J) F_eq s e^{-i Omega t}
    dPhi(r,t) = 8 pi^2 L0 eps Re ∫ dJ F_eq s C(r,J) e^{-i Omega t},
    C(r, J)   = ∫ dQ e^{iQ} G(r, r(Q,J)),   G(r,r') = -1/max(r,r'),

con G la función de Green de capas esféricas. Solo hace falta el mapa inverso
(Q,J) -> r, una vez. La forma de la cola la dicta la regularidad de cada
integrando en J (un extremo o punto donde se comporta como x^a da t^-(a+1)):

  * en J = 0, C ~ sqrt(J) (amplitud epicíclica), así que dPhi va como
    F_eq s J^{1/2}: si F_eq(0) > 0, t^{-3/2} con s = 1, t^{-2} con
    s = sqrt(J/J_max) y t^{-3} con s = (J/J_max)^{3/2}; si F_eq ~ J^2, t^{-7/2};
  * donde un punto de retorno de la órbita cruza el radio r, C tiene una
    singularidad (J - J*)^{3/2}: t^{-5/2}, para cualquier F_eq. h_1 no la tiene
    porque B(J) es suave.

Se usa el isócrono desnudo (a0 -> 0): la cola es propiedad de la forma de F_eq.
Reporta |h_1(t)|/|h_1(0)| y el rms en r in [3, 15] de |dPhi| relativo a t = 0,
en t = n tau_1 con tau_1 = 2 pi/(ancho de la banda); la pendiente logarítmica
entre 20 y 40 tau_1 (la asintótica); y dPhi(0) por unidad de masa y de eps
relativo al primer caso, que recalibra el parámetro de atrapamiento mu.
La convergencia se comprueba con media malla en J.

Uso:  python3 colas_libres.py [--nu 2400] [--nq 128] [--ntau 1 2 5 10]
"""
import os, sys, argparse, time, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from equilibrio import Equilibrio, invertir, F_perfil, F_maxwell
from landau_libre import omega_de_J

warnings.filterwarnings('ignore', category=RuntimeWarning)
L0 = 2.0
C_ISO = 0.5*(L0 + np.sqrt(L0**2 + 4))
RMED = np.linspace(3, 15, 121)                  # los radios de lineal.py


def prueba(j1, sj1):
    return lambda J: np.exp(-(J - j1)**2/sj1**2)*J**2


def trapecio(u):
    w = np.gradient(u)*0 + (u[1] - u[0])
    w[0] = w[-1] = 0.5*(u[1] - u[0])
    return w


def validar():
    """Contra landau_libre (suma sobre las partículas de 11_landau, sigma = 0.10)."""
    raiz = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'landau')
    eqf = os.path.join(raiz, 'L_a1e-2_n400_e0.1_equilibrio.npz')
    if not os.path.exists(eqf):
        print('validación omitida: falta', eqf)
        return
    eq = np.load(eqf)
    om = omega_de_J(eq['E_t'], eq['J_t'])
    u = np.linspace(0, 1, 40001); J = 0.6*u**2; wJ = trapecio(u)*1.2*u
    g = wJ*prueba(0.10, 0.10)(J)*F_perfil(J, 'gauss', 0.10)
    t = np.array([400, 800, 1000, 1200, 1400.0])
    h = np.abs(np.exp(-1j*np.outer(t, om(J))) @ g)/np.sum(g)
    print('validación con 11_landau (sigma_J = 0.10): |h_1|/|h_1(0)| libre')
    print('   t      cuadratura   landau_libre')
    for tt, hh, ref in zip(t, h, (5.5e-1, 8.0e-2, 1.6e-2, 1.5e-3, 9.9e-5)):
        print(f'{tt:6.0f}   {hh:10.2e}   {ref:10.1e}')
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--nu', type=int, default=2400, help='nodos en u, con J = J_end u^2 (par)')
    ap.add_argument('--nq', type=int, default=128, help='nodos en Q para C(r,J)')
    ap.add_argument('--ntau', type=float, nargs='+', default=[1, 2, 5, 10],
                    help='tiempos que se muestran, en unidades de tau_1')
    arg = ap.parse_args()
    NTAU = list(arg.ntau)
    NT = NTAU + [20, 40]                             # la pendiente, siempre en 20-40 tau_1
    validar()

    # Una sola malla (Q, J) para todos los casos: J = J_end u^2 resuelve sqrt(J) en 0.
    J_end = 0.30
    u = np.linspace(0, 1, arg.nu + 1)
    J = J_end*u**2
    Qn = (np.arange(arg.nq) + 0.5)*2*np.pi/arg.nq
    eq = Equilibrio(1e-12, jmax=J_end)
    m = eq.mapa()
    E_t, J_t = eq.tabla_J_de_E(m)
    JJ, QQ = np.meshgrid(J[1:], Qn, indexing='ij')
    t0 = time.time()
    r = np.empty(JJ.size)
    for s in range(0, JJ.size, 8000):
        r[s:s+8000], _ = invertir(m, E_t, J_t, QQ.ravel()[s:s+8000], JJ.ravel()[s:s+8000])
    r = r.reshape(JJ.shape)
    print(f'mapa inverso: {JJ.size} nodos en {time.time()-t0:.0f} s')
    eQ = np.exp(1j*Qn)*2*np.pi/arg.nq
    C = np.zeros((len(RMED), len(J)), complex)       # C(r, J = 0) = 0
    for i, rr in enumerate(RMED):
        C[i, 1:] = (-1.0/np.maximum(rr, r)) @ eQ
    Om = 1.0/(J + C_ISO)**3

    Ecl = -0.5/(J + C_ISO)**2                        # E(J) del isócrono desnudo

    def maxwell(g, w0, jt=0.138):
        Eb, Ec = -0.5/(jt + C_ISO)**2, -0.5/C_ISO**2
        return F_maxwell(Ecl, Eb, (Eb - Ec)/w0, g)

    casos = [   # nombre, F_eq(J), perturbación, J del soporte (None: gaussiana)
        ('gauss sigma=0.05',     F_perfil(J, 'gauss', 0.05), 'plana', None),
        ('J^2 (J_t-J)^2',        F_perfil(J, 'politropo', jt=0.138, k=2.0, m=2.0), 'plana', 0.138),
        ('Wilson W0=3',          maxwell(2.0, 3.0), 'plana', 0.138),
        ('Wilson W0=3',          maxwell(2.0, 3.0), 'suave', 0.138),
        ('Wilson W0=3',          maxwell(2.0, 3.0), 'suave3', 0.138),
        ('Wilson W0=6',          maxwell(2.0, 6.0), 'suave3', 0.138),
    ]
    EXP = {'plana': 0.0, 'suave': 0.5, 'suave3': 1.5}
    B = prueba(0.05, 0.05)                           # j1 = sj1 = 0.36 J_t, como 11_landau
    print(f'isócrono L0 = {L0:g}; B(J) con j1 = sj1 = 0.05; rms de dPhi en r in [3, 15]\n')
    cab = (f'{"F_eq":>20} {"pert":>6} {"tau_1":>6} | ' +
           ' '.join(f'{"h1 %gtau" % n:>9}' for n in NTAU) + f' {"pend.":>6} | ' +
           ' '.join(f'{"dPhi %gtau" % n:>9}' for n in NTAU) + f' {"pend.":>6} {"conv.":>6}'
           f' {"dPhi0/M":>8}')
    print(cab); print('-'*len(cab))
    ref0 = None
    w = trapecio(u)*2*J_end*u
    for nombre, F, pert, jsop in casos:
        s = (J/(jsop if jsop else 1.0))**EXP[pert]
        if jsop is None:                             # gaussiana: F > 1 % del máximo
            sop = F >= 0.01*F.max()
            ancho = Om[sop].max() - Om[sop].min()
        else:
            ancho = 1.0/C_ISO**3 - 1.0/(jsop + C_ISO)**3
        tau1 = 2*np.pi/ancho
        t = np.array([0.0] + [n*tau1 for n in NT])
        fase = np.exp(-1j*np.outer(t, Om))
        res = []
        for paso in (1, 2):                          # malla completa y media malla
            uu = u[::paso]; ww = trapecio(uu)*2*J_end*uu
            g = (F*s)[::paso]*ww
            h = np.abs(fase[:, ::paso] @ (B(J[::paso])*g))
            Z = (C[:, ::paso]*g[None, :]) @ fase[:, ::paso].T      # (radio, tiempo)
            ph = np.sqrt(np.mean(np.abs(Z)**2, axis=0))
            res.append((h/h[0], ph/ph[0], ph[0]))
        (h, ph, ph0), (_, ph2, _) = res
        conv = np.max(np.abs(ph2[1:]/ph[1:] - 1))
        porM = ph0/np.sum(F*w)                       # dPhi(0) por unidad de masa y de eps
        ref0 = porM if ref0 is None else ref0
        pend = lambda x: np.log(x[-1]/x[-2])/np.log(2.0)
        n = len(NTAU)
        print(f'{nombre:>20} {pert:>6} {tau1:6.0f} | ' +
              ' '.join(f'{x:9.1e}' for x in h[1:n+1]) + f' {pend(h):6.1f} | ' +
              ' '.join(f'{x:9.1e}' for x in ph[1:n+1]) + f' {pend(ph):6.1f} {conv:6.0e}'
              f' {porM/ref0:8.2f}')
    print('\npend. = d log(amplitud)/d log(t) entre 20 y 40 tau_1; conv. = cambio relativo '
          'máximo de dPhi con media malla en J (hasta 40 tau_1);\ndPhi0/M = dPhi(0) rms por '
          'unidad de masa y de eps, relativo al primer caso.')


if __name__ == '__main__':
    main()
