"""Ganancia del lazo autoconsistente, lambda(omega), de los modos radiales con L fijo, y la
hipótesis de que hay modo discreto bajo la banda si y solo si lambda(Omega_min) > 1.

Para un modo dF = a(Q,J) e^{-i omega t}, con a = sum_k a_k(J) e^{ikQ} y phi = dPhi(r(Q,J)),
la ecuación linealizada da (k Omega - omega) a_k = k F_eq' phi_k. La densidad solo ve la
parte par en Q, porque r(-Q) = r(Q); sumando k y -k,

    b_k = a_k + a_{-k} = R_k phi_k,     R_k = 2 k^2 Omega F_eq' / (k^2 Omega^2 - omega^2),

y Poisson con la función de Green de las cáscaras, G = -1/max(r, r'), cierra el lazo:

    phi_k(J) = sum_k' int M_kk'(J, J') b_k'(J') dJ',
    M_kk'(J, J') = 4 pi L0 int int cos(kQ) cos(k'Q') G(r(Q,J), r(Q',J')) dQ dQ'.

Hay modo a la frecuencia omega si K(omega) = M R(omega) tiene autovalor 1. M es negativa
(la gravedad atrae); con F_eq' < 0 y omega < Omega_min también lo es R, y K es semejante a
la matriz simétrica y positiva |R|^{1/2} (-M) |R|^{1/2}. lambda(omega) es su mayor autovalor;
crece con omega, y lambda_borde = lambda(Omega_min). Con F_eq' de signo variable (familia
hueca) se usa el mayor autovalor real de K.

Discretización: J = J_t (1 - s^p) con Gauss-Legendre en s, Q en nq puntos medios, r(Q,J)
con el mapa inverso del equilibrio (equilibrio.invertir), y las integrales en Q depositando
cada punto de la órbita en una malla radial de nr nodos con pesos lineales; entonces
M = 4 pi L0 P G P^T, con G = -1/max(r_m, r_n) en la malla. Con un borde (E_t - E)^g, en
omega = Omega_min el integrando va como s^(p (g-1) - 1) ds: p = 2 lo deja regular si
g >= 1.5, y con 1 < g < 1.5 se toma p = 1/(g-1). Con g <= 1 la integral diverge en el borde
y p = 4 acerca los nodos a él, para resolver lambda hasta 1e-12 anchos de banda del borde.
E_t - E y Omega - Omega_min se evalúan en los nodos en función de la distancia al borde,
J_t - J = J_t s^p, sin restar números casi iguales (borde='spline'). borde='tabla' es la
aritmética anterior al 2026-10-07: E_t con la interpolación lineal de equilibrio.borde_E,
que queda ~1e-9 por debajo del spline que da E en los nodos; F se anula entonces ~1e-8
anchos de banda antes de Omega_min, y ese corte redondea la singularidad del borde. El
efecto en lambda_borde es 0.025 con g = 1.25 (masa puntual, a0 = 1), 3e-4 con g = 1.5 y
menor que 1e-6 con g >= 2; omega_d no cambia (menos de 1e-8).

    python3 lambda_borde.py [--casos A4 L5 ...] [--nj 200] [--kmax 6] [--nq 128] [--nr 2000]
                            [--borde] [--debil]

Lee exe/eta_lineal/<caso>_equilibrio.npz, compara con exe/eta_lineal/tabla.txt (el polo
del solucionador lineal, barrido_lineal.py) y escribe exe/eta_lineal/lambda_borde.txt.
Con --borde, cómo se acerca lambda a su valor en el borde según g; con --debil, lambda
cerca del borde en los bordes abruptos con acoplamiento débil.
"""
import os, sys, argparse, numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from aa_numerico import MapaAA
from equilibrio import invertir, FORMAS_E, borde_E, F_E, dFdE_E, F_perfil, dF_perfil

DIR = os.path.join(AQUI, '..', '..', 'exe', 'eta_lineal')


class Lazo:
    """El operador del lazo de un equilibrio, listo para evaluar lambda(omega)."""

    def __init__(self, npz, nj=200, kmax=6, nq=128, nr=2000, pot=None, borde='spline'):
        d = np.load(npz)
        self.a0, L0, jt = float(d['a0']), float(d['L0']), float(d['J_borde'])
        forma, g = str(d['forma']), float(d['k_borde'])
        if pot is None:                 # J_t - J = J_t s^pot: |R_1| ds ~ s^(pot (g-1) - 1) ds en el borde
            pot = 2 if forma not in FORMAS_E or borde == 'tabla' or g >= 1.5 else 4 if g <= 1 else 1/(g - 1)
        self.pot = pot
        sp = CubicSpline(d['J_t'], d['E_t'])
        Omf = sp.derivative()
        x, w = np.polynomial.legendre.leggauss(nj)
        s, ws = 0.5*(x + 1), 0.5*w
        self.Om_min, self.Om_max = float(Omf(jt)), float(Omf(0.0))
        A = float(d['A'])
        self.dOm = None
        if borde == 'tabla':            # la aritmética anterior al 2026-10-07, bit a bit
            assert pot == 2
            J = jt*(1 - s**2)
            self.wJ = ws*2*jt*s
            E, self.Om = sp(J), Omf(J)
        else:
            xb = jt*s**pot              # distancia al borde, J_t - J, sin resta
            J = jt - xb
            self.wJ = ws*pot*jt*s**(pot - 1)
            E, self.Om = sp(J), Omf(J)
            # E_t - E y Omega - Omega_min. Dentro del último tramo del spline se desarrolla su
            # polinomio en J_t, y así no hay cancelación en los nodos pegados al borde.
            i = int(np.clip(np.searchsorted(d['J_t'], jt, side='left') - 1, 0, len(d['J_t']) - 2))
            a, (c0, c1, c2) = jt - float(d['J_t'][i]), (float(c) for c in sp.c[:3, i])
            dentro = xb <= a
            self.uE = np.where(dentro, (c0*(3*a**2 - 3*a*xb + xb**2) + c1*(2*a - xb) + c2)*xb, float(sp(jt)) - E)
            self.dOm = np.where(dentro, (3*c0*xb - (6*c0*a + 2*c1))*xb, self.Om - self.Om_min)
        if forma in FORMAS_E and borde == 'tabla':
            Et, T = borde_E(forma, d['E_t'], d['J_t'], jt, float(d['w0']))
            F = A*F_E(forma, E, Et, T, g)
            self.dF = A*dFdE_E(forma, E, Et, T, g)*self.Om
        elif forma in FORMAS_E:
            # La energía del borde es la del spline en J_t, el mismo que da E en los nodos.
            # equilibrio.borde_E la interpola linealmente en la tabla y queda ~1e-9 por debajo:
            # con ella F se anula ~1e-8 anchos de banda antes de Omega_min, y ese corte redondea
            # la singularidad del borde (con k = 1.25 y a0 = 1, lambda_borde baja de 1.583 a 1.558).
            T = (float(sp(jt)) - float(d['E_t'][0]))/float(d['w0']) if forma == 'maxwell' else None
            F = A*F_E(forma, -self.uE, 0.0, T, g)
            self.dF = A*dFdE_E(forma, -self.uE, 0.0, T, g)*self.Om
        else:
            m = float(d['m_borde'])
            F = A*F_perfil(J, forma, None if forma != 'gauss' else float(d['sigma_J']), jt, g, m)
            self.dF = A*dF_perfil(J, forma, None if forma != 'gauss' else float(d['sigma_J']), jt, g, m)
        self.masa = 16*np.pi**3*L0*np.sum(self.wJ*F)          # debe ser a0
        self.monotona = bool(np.all(self.dF <= 0))
        # r(Q, J) en los nodos, con el mapa inverso del equilibrio.
        self.fondo = str(d['fondo']) if 'fondo' in d.files else 'isocrono'
        mapa = MapaAA(d['r'], d['phi_self'], L=L0, fondo=self.fondo)
        Q = (np.arange(nq) + 0.5)*2*np.pi/nq
        JJ, QQ = np.meshgrid(J, Q, indexing='ij')
        r = np.empty(JJ.size)
        for i in range(0, JJ.size, 4000):
            r[i:i+4000], _ = invertir(mapa, d['E_t'], d['J_t'], QQ.ravel()[i:i+4000], JJ.ravel()[i:i+4000])
        r = r.reshape(nj, nq)
        # P[(k-1) nj + i, m] = int dQ cos(kQ) w_m(r(Q, J_i)), con pesos lineales en la malla.
        rg = np.linspace(r.min() - 1e-3, r.max() + 1e-3, nr)
        dr = rg[1] - rg[0]
        u = (r - rg[0])/dr
        i0 = np.minimum(np.floor(u).astype(int), nr - 2)
        f = u - i0
        self.kmax, self.nj = kmax, nj
        P = np.zeros((kmax*nj, nr))
        filas = np.repeat(np.arange(nj)[:, None], nq, axis=1)
        for k in range(1, kmax + 1):
            c = np.cos(k*Q)[None, :]*2*np.pi/nq
            fk = (k - 1)*nj + filas
            np.add.at(P, (fk, i0), c*(1 - f))
            np.add.at(P, (fk, i0 + 1), c*f)
        G = -1.0/np.maximum(rg[:, None], rg[None, :])
        self.M = 4*np.pi*L0*(P @ G @ P.T)
        self.k2 = (np.arange(1, kmax + 1)**2)[:, None]*np.ones((1, nj))

    def R(self, omega):
        """R_k(J; omega) en el orden (k-1) nj + i, por los pesos de la cuadratura en J."""
        Om = self.Om[None, :]
        if self.dOm is None:
            return (2*self.k2*Om*self.dF[None, :]/(self.k2*Om**2 - omega**2)*self.wJ[None, :]).ravel()
        n = np.sqrt(self.k2)
        menos = n*Om - omega                        # n Omega - omega;
        menos[0] = self.dOm + (self.Om_min - omega)  # con n = 1, sin cancelación junto al borde
        return (2*self.k2*Om*self.dF[None, :]/(menos*(n*Om + omega))*self.wJ[None, :]).ravel()

    def lam(self, omega):
        D = self.R(omega)
        if self.monotona:
            h = np.sqrt(np.abs(D))
            return float(np.linalg.eigvalsh(h[:, None]*(-self.M)*h[None, :])[-1])
        ev = np.linalg.eigvals(self.M*D[None, :])
        return float(np.max(ev[np.abs(ev.imag) < 1e-9*np.abs(ev).max()].real))

    def delta_modo(self, dmin=1e-12):
        """Distancia del modo al borde, delta_d = (Omega_min - omega_d)/(Omega_max - Omega_min),
        de lambda(omega_d) = 1 con F_eq monótona. La raíz se busca en log(delta), para resolver
        los modos pegados al borde (g <= 1 con masa pequeña). None si lambda < 1 hasta dmin
        anchos de banda del borde, o si lambda(0) >= 1."""
        an = self.Om_max - self.Om_min
        f = lambda x: self.lam(self.Om_min - 10**x*an) - 1
        a, b = np.log10(dmin), np.log10(self.Om_min/an)        # hasta omega = 0
        if f(a) <= 0 or self.lam(0.0) >= 1:
            return None
        return 10**brentq(f, a, b, xtol=1e-9)

    def omega_modo(self):
        """Frecuencia del modo discreto, lambda(omega) = 1 bajo la banda, si lo hay. Se busca
        el primer cruce desde abajo en 2 anchos de banda; si F_eq' cambia de signo lambda no
        tiene por qué ser monótona, y por eso se recorre una malla antes de afinar."""
        ancho = self.Om_max - self.Om_min
        if self.lam(self.Om_min) <= 1:
            return None
        # En [0, Omega_min): con una banda ancha (masa puntual) 2 anchos llegarían a omega < 0.
        ws = np.maximum(self.Om_min - ancho*np.concatenate([np.linspace(2, 0.1, 20), [0.03, 0.01, 0.0]]),
                        0.0)
        ls = [self.lam(w) for w in ws]
        for i in range(len(ws) - 1):
            if ls[i] < 1 <= ls[i+1]:
                return brentq(lambda w: self.lam(w) - 1, ws[i], ws[i+1], xtol=1e-9)
        return None


def tabla_lineal():
    """omega y x del solucionador lineal (exe/eta_lineal/tabla.txt)."""
    res = {}
    for l in open(os.path.join(DIR, 'tabla.txt')):
        p = l.split()
        if len(p) > 9 and p[0] not in ('caso',) and p[1][0].isdigit():
            res[p[0]] = (float(p[6]), float(p[8]))
    return res


def borde(w, casos=(('G3', 3), ('A4', 2), ('G15', 1.5), ('G1', 1)), nj=400):
    """Cómo se acerca lambda a su valor en el borde: lambda_borde - lambda(Omega_min - delta)
    ~ delta^min(g-1, 1) (con logaritmo si g = 2); para g = 1, lambda ~ ln(1/delta)."""
    ds = np.array([3e-2, 1e-2, 3e-3, 1e-3, 3e-4])
    w('Aproximación al borde, delta en anchos de banda: ' + ' '.join(f'{d:.0e}' for d in ds))
    for n, g in casos:
        z = Lazo(os.path.join(DIR, f'{n}_equilibrio.npz'), nj=nj)
        ancho = z.Om_max - z.Om_min
        lb = z.lam(z.Om_min)
        l = np.array([z.lam(z.Om_min - d*ancho) for d in ds])
        if g > 1:
            p = np.diff(np.log(lb - l))/np.diff(np.log(ds))
            w(f'  {n:4} g = {g:g}: lambda_borde = {lb:.5f}; pendiente de ln(lambda_borde - lambda) '
              f'en ln delta: ' + ' '.join(f'{x:.2f}' for x in p))
        else:
            p = np.diff(l)/np.diff(np.log(1/ds))
            w(f'  {n:4} g = {g:g}: d lambda/d ln(1/delta): ' + ' '.join(f'{x:.3f}' for x in p))


def debil(w, casos=('G075a', 'G075b', 'G1a', 'G1b'), nj=400):
    """Bordes abruptos con acoplamiento débil: lambda muy cerca del borde, y la distancia
    delta_d del modo bajo el borde, en anchos de banda, si cae por encima de 1e-8. Separar
    el modo de la respuesta del borde lleva un tiempo del orden de tau_1/delta_d."""
    ds = 10.0**-np.arange(1, 9)
    w('Bordes abruptos, acoplamiento débil: lambda en delta/ancho = 1e-1 ... 1e-8, y delta_d')
    for n in casos:
        z = Lazo(os.path.join(DIR, f'{n}_equilibrio.npz'), nj=nj)
        ancho = z.Om_max - z.Om_min
        f = lambda ld: z.lam(z.Om_min - np.exp(ld)*ancho) - 1
        l = [z.lam(z.Om_min - d*ancho) for d in ds]
        txt = (f'delta_d = {np.exp(brentq(f, np.log(1e-8), np.log(1e-1), xtol=1e-6)):.1e}'
               if l[-1] > 1 > l[0] else 'lambda < 1 hasta 1e-8')
        w(f'  {n:5} a0 = {z.a0:.4f}: ' + ' '.join(f'{x:.3f}' for x in l) + f'   {txt}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--casos', nargs='+', default=['L1', 'A1', 'L2', 'A2', 'L3', 'A3', 'L4', 'A4',
                                                    'L5', 'A5', 'L6', 'M5', 'A6', 'A7', 'A8', 'A9',
                                                    'D1', 'D2', 'E3', 'G3', 'G15', 'G1', 'G075',
                                                    'E1', 'H9', 'H3'])
    ap.add_argument('--nj', type=int, default=200)
    ap.add_argument('--kmax', type=int, default=6)
    ap.add_argument('--nq', type=int, default=128)
    ap.add_argument('--nr', type=int, default=2000)
    ap.add_argument('--borde', action='store_true', help='además, la aproximación al borde')
    ap.add_argument('--debil', action='store_true', help='además, bordes abruptos con acoplamiento débil')
    a = ap.parse_args()
    lin = tabla_lineal()
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w(f'nj = {a.nj}, kmax = {a.kmax}, nq = {a.nq}, nr = {a.nr}')
    w(f'{"caso":6} {"a0":>7} {"masa/a0-1":>10} {"monót.":>6} {"lambda_borde":>12} '
      f'{"omega_modo":>11} {"x_modo":>7} | lineal: {"omega":>8} {"x":>6}')
    for n in a.casos:
        z = Lazo(os.path.join(DIR, f'{n}_equilibrio.npz'), a.nj, a.kmax, a.nq, a.nr)
        lb = z.lam(z.Om_min)
        wm = z.omega_modo()
        wl, xl = lin.get(n, (float('nan'), float('nan')))
        txt = (f'{wm:11.5f} {(wm - z.Om_min)/(z.Om_max - z.Om_min):+7.3f}' if wm is not None
               else f'{"--":>11} {"--":>7}')
        w(f'{n:6} {z.a0:7.4f} {z.masa/z.a0 - 1:10.1e} {"sí" if z.monotona else "no":>6} {lb:12.4f} '
          f'{txt} | lineal: {wl:8.5f} {xl:+6.2f}')
    if a.borde:
        w(); borde(w)
    if a.debil:
        w(); debil(w)
    open(os.path.join(DIR, 'lambda_borde.txt'), 'w').write('\n'.join(lineas) + '\n')
