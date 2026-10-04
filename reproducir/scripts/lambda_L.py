"""Ganancia del lazo, lambda(omega), de los modos radiales de un equilibrio esférico f0(E, L)
con distribución en el momento angular: la norma del operador de Mathur.

Es la generalización de lambda_borde.py (L fijo) a dos acciones. Con dF = a(Q, J, L) e^{-i omega t}
y perturbaciones radiales, L se conserva y la ecuación linealizada es la misma en cada L,

    (k Omega - omega) a_k = k (df0/dJ) phi_k,     df0/dJ = Omega df0/dE,
    b_k = a_k + a_{-k} = R_k phi_k,               R_k = 2 k^2 Omega df0/dJ / (k^2 Omega^2 - omega^2),

con Omega = Omega(E, L). La masa es 8 pi^2 int L dL int dJ int dQ f, de modo que

    phi_k(J, L) = sum_k' int L' dL' dJ' M_kk'(J, L; J', L') b_k'(J', L'),
    M_kk' = 4 pi int int cos(kQ) cos(k'Q') G(r, r') dQ dQ',      G = -1/max(r, r'):

respecto de L fijo solo cambia la medida, L0 dJ -> L dL dJ. En las variables (E, L),
dJ = dE/Omega y el peso de cada nodo es D = L dL dE 2 k^2 Omega df0/dE / (k^2 Omega^2 - omega^2).

Con df0/dE < 0 y omega bajo la banda, el operador es semejante a |D|^{1/2} (-M) |D|^{1/2},
positivo, con -M = 4 pi P (-G) P^T y P[n k, a] = int dQ cos(kQ) w_a(r(Q)) (pesos lineales en
una malla radial). Sus autovalores no nulos son los de 4 pi C^T (P^T |D| P) C, con -G = C C^T:
una matriz del tamaño de la malla radial, sea cual sea el número de órbitas. Es el operador
de Mathur escrito en el espacio de configuración (Hadžić, Rein y Straub 2022).

Órbitas: con s = r^2, dt = ds / (2 sqrt(g(s))), g(s) = 2 (E - Phi) s - L^2, que tiene ceros
simples en s_- y s_+ también cuando L -> 0; con s = s_m + s_a sin(theta), dt = dtheta/(2 sqrt(h)),
h = g/((s - s_-)(s_+ - s)), suave. La regla del punto medio en theta converge muy rápido.

    python3 lambda_L.py validar        el límite de L fijo frente a lambda_borde.py
    python3 lambda_L.py politropos     lambda_borde(k) de los politropos isótropos (E0 - E)^k y k*
    python3 lambda_L.py king           lambda_borde(kappa) de los modelos de King y kappa*
    python3 lambda_L.py convergencia   k* frente a cada parámetro numérico
    python3 lambda_L.py dispersion     L fijo frente a una dispersión en L, con k <= 1
Las tablas quedan en exe/lambda_L/.
"""
import os, sys, argparse, numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import beta
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)


def orbitas(phi, dphi, E, L, r_lo, r_hi, nth=96, tabla=4000):
    """Órbitas (E, L) del potencial phi(r) entre r_lo y r_hi. Devuelve Omega y, en los nth
    puntos medios en theta de media órbita (del pericentro al apocentro), r, Q y dQ.

    El radio de la órbita circular de cada L sale de invertir r^3 phi'(r) = L^2 (creciente);
    los puntos de retorno, por bisección de w^2 = 2 (E - phi) - L^2/r^2 a cada lado."""
    E, L = np.asarray(E, float), np.asarray(L, float)
    rt = np.linspace(r_lo, r_hi, tabla)
    c = rt**3*dphi(rt)
    if np.any(np.diff(c) <= 0):
        raise SystemExit('r^3 phi\' no es creciente: hay órbitas circulares inestables')
    rc = np.interp(L**2, c, rt)
    w2 = lambda r: 2*(E - phi(r)) - L**2/r**2
    if np.any(w2(rc) <= 0):
        raise SystemExit('nodo sin órbita: E por debajo de la órbita circular')
    if np.any(w2(np.full_like(rc, r_lo)) >= 0) or np.any(w2(np.full_like(rc, r_hi)) >= 0):
        raise SystemExit('órbita que se sale de [r_lo, r_hi]')

    def raiz(a, b):                       # w2(a) < 0 < w2(b), o al revés
        fa = w2(a)
        for _ in range(60):
            m = 0.5*(a + b); fm = w2(m)
            iz = np.sign(fm) == np.sign(fa)
            a, b = np.where(iz, m, a), np.where(iz, b, m)
        return 0.5*(a + b)
    rm, rp = raiz(np.full_like(rc, r_lo), rc.copy()), raiz(np.full_like(rc, r_hi), rc.copy())
    sm, sa = 0.5*(rp**2 + rm**2), 0.5*(rp**2 - rm**2)
    th = (np.arange(nth) + 0.5)*np.pi/nth - 0.5*np.pi
    s = sm[:, None] + sa[:, None]*np.sin(th)[None, :]
    r = np.sqrt(s)
    g = 2*(E[:, None] - phi(r))*s - (L**2)[:, None]
    h = g/((s - (rm**2)[:, None])*((rp**2)[:, None] - s))
    dt = (np.pi/nth)/(2*np.sqrt(h))
    T = 2*dt.sum(axis=1)
    Om = 2*np.pi/T
    t = np.cumsum(dt, axis=1) - 0.5*dt
    return Om, r, Om[:, None]*t, Om[:, None]*dt, rm, rp


class LazoL:
    """Operador del lazo de un conjunto de órbitas (E, L) con pesos w = L dL dE y df0/dE."""

    def __init__(self, phi, dphi, E, L, w, dfdE, r_lo, r_hi, kmax=6, nth=96, nr=400):
        self.Om, r, Q, dQ, self.rm, self.rp = orbitas(phi, dphi, E, L, r_lo, r_hi, nth)
        n = len(E)
        self.rg = (np.arange(nr) + 0.5)*(r.max()*(1 + 1e-9))/nr          # centros de celda en [0, r_max]
        u = r/(self.rg[1] - self.rg[0]) - 0.5
        i0 = np.clip(np.floor(u).astype(int), 0, nr - 2)
        f = u - i0
        filas, cols, vals = [], [], []
        for k in range(1, kmax + 1):
            c = 2*np.cos(k*Q)*dQ                                  # las dos mitades de la órbita
            fk = (k - 1)*n + np.repeat(np.arange(n)[:, None], nth, axis=1)
            filas += [fk.ravel(), fk.ravel()]; cols += [i0.ravel(), (i0 + 1).ravel()]
            vals += [(c*(1 - f)).ravel(), (c*f).ravel()]
        self.P = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(filas), np.concatenate(cols))),
                                   shape=(kmax*n, nr))
        self.C = np.linalg.cholesky(1.0/np.maximum(self.rg[:, None], self.rg[None, :]))
        self.k2 = np.repeat(np.arange(1, kmax + 1)**2, n)
        self.Omk, self.wk, self.dfk = np.tile(self.Om, kmax), np.tile(w, kmax), np.tile(dfdE, kmax)
        self.Om_nodos = float(self.Om.min())
        if np.any(dfdE > 0):
            raise SystemExit('df0/dE > 0: el operador simétrico no vale')

    def lam(self, omega):
        D = np.abs(self.wk*2*self.k2*self.Omk*self.dfk/(self.k2*self.Omk**2 - omega**2))
        K = (self.P.T @ sparse.diags(D) @ self.P).toarray()
        return float(np.linalg.eigvalsh(4*np.pi*self.C.T @ K @ self.C)[-1])


# ------------------------------------------------------------------ equilibrios isótropos
class Isotropo:
    """f0 = Phi(E0 - E), autogravitante. Con y = E0 - U0, y'' + 2 y'/r = -4 pi rho(y), y(0) = kappa,
    rho(y) = 4 pi sqrt(2) int_0^y Phi(eta) (y - eta)^{1/2} d eta. Se toma E0 = 0, así que
    phi = -y y la energía de ligadura es e = -E, entre 0 (el borde) y kappa (el centro)."""

    def __init__(self, kappa, Phi, dPhi, rho, nombre=''):
        self.kappa, self.Phi, self.dPhi, self.nombre = kappa, Phi, dPhi, nombre
        f = lambda r, z: [z[1], -4*np.pi*rho(max(z[0], 0.0)) - 2*z[1]/r]
        cero = lambda r, z: z[0]
        cero.terminal, cero.direction = True, -1
        r0, rc = 1e-6, rho(kappa)
        sol = solve_ivp(f, [r0, 1e5], [kappa - 4*np.pi*rc*r0**2/6, -4*np.pi*rc*r0/3], events=cero,
                        rtol=1e-12, atol=1e-14, dense_output=True)
        self.R = float(sol.t_events[0][0])
        r = np.concatenate([[0.0], np.linspace(r0, self.R, 6000)])
        y = np.concatenate([[kappa], sol.sol(r[1:])[0]]); y[-1] = 0.0
        self.y = CubicSpline(r, y, bc_type=((1, 0.0), 'not-a-knot'))
        self.dy = self.y.derivative()
        self.M = -self.R**2*float(self.dy(self.R))
        self.phi = lambda r: -self.y(r)
        self.dphi = lambda r: -self.dy(r)

    def Omega_borde(self, nth=4000):
        """Frecuencia de la órbita radial de energía E0 (L = 0, de 0 a R): el periodo más largo."""
        th = (np.arange(nth) + 0.5)*np.pi/nth - 0.5*np.pi
        s = 0.5*self.R**2*(1 + np.sin(th))
        h = 2*self.y(np.sqrt(s))/(self.R**2 - s)
        return 2*np.pi/(2*np.sum((np.pi/nth)/(2*np.sqrt(h))))

    def nodos(self, ne=48, nx=24, q=2.0):
        """Nodos (E, L): e = kappa s^q con s de Gauss-Legendre en (0, 1) (se acumulan en el
        borde e -> 0), y L = L_c(E) x, x de Gauss-Legendre en (0, 1)."""
        xs, ws = np.polynomial.legendre.leggauss(ne); s, ws = 0.5*(xs + 1), 0.5*ws
        xx, wx = np.polynomial.legendre.leggauss(nx); x, wx = 0.5*(xx + 1), 0.5*wx
        e, we = self.kappa*s**q, self.kappa*q*s**(q - 1)*ws
        rt = np.linspace(1e-6*self.R, self.R, 20000)
        Ec = self.phi(rt) + 0.5*rt*self.dphi(rt)                 # energía de la órbita circular
        rcirc = np.interp(-e, Ec, rt)
        Lc = np.sqrt(rcirc**3*self.dphi(rcirc))
        E = np.repeat(-e, nx)
        L = (Lc[:, None]*x[None, :]).ravel()
        w = (we[:, None]*Lc[:, None]**2*x[None, :]*wx[None, :]).ravel()      # L dL dE
        return E, L, w, -self.dPhi(np.repeat(e, nx))

    def lazo(self, ne=48, nx=24, kmax=16, nth=128, nr=400, q=2.0):
        E, L, w, dfdE = self.nodos(ne, nx, q)
        return LazoL(self.phi, self.dphi, E, L, w, dfdE, 1e-7*self.R, self.R, kmax, nth, nr)

    def masa(self, ne=48, nx=24):
        """8 pi^2 int L dL dE T f0, que debe dar M."""
        E, L, w, _ = self.nodos(ne, nx)
        Om = orbitas(self.phi, self.dphi, E, L, 1e-7*self.R, self.R)[0]
        return 8*np.pi**2*np.sum(w*(2*np.pi/Om)*self.Phi(-E))


def politropo(k):
    """f0 = (E0 - E)^k: Lane-Emden de índice k + 3/2. La familia es invariante de escala; kappa = 1."""
    ck = 4*np.pi*np.sqrt(2)*beta(k + 1, 1.5)
    return Isotropo(1.0, lambda e: e**k, lambda e: k*e**(k - 1), lambda y: ck*y**(k + 1.5), f'politropo k = {k:g}')


def king(kappa):
    """f0 = e^{E0 - E} - 1, con y(0) = kappa (la convención de Straub 2024)."""
    xu, wu = np.polynomial.legendre.leggauss(64); u, wu = 0.5*(xu + 1), 0.5*wu
    rho = lambda y: 8*np.pi*np.sqrt(2)*y**1.5*np.sum(wu*u**2*np.expm1(y*(1 - u**2)))
    return Isotropo(kappa, np.expm1, np.exp, rho, f'King kappa = {kappa:g}')


def resultado(p, **kw):
    """lambda(0), lambda_borde, y la frecuencia del modo si lo hay."""
    z = p.lazo(**kw)
    ob = p.Omega_borde()
    lb = z.lam(ob)
    wd = brentq(lambda w: z.lam(w) - 1, 0.0, ob, xtol=1e-12*ob) if lb > 1 else None
    return dict(R=p.R, M=p.M, Om_b=ob, Om_max=float(z.Om.max()), Om_nodos=z.Om_nodos, lam0=z.lam(0.0),
                lam_b=lb, omega_d=wd, z=z)


# ------------------------------------------------------------------ pasos
def validar():
    """Límite de L fijo: f0 = F(E) g(L), con g estrecha alrededor de L0 e int g L dL = L0, en el
    potencial total de un equilibrio de L fijo. Debe dar lambda de lambda_borde.Lazo."""
    from lambda_borde import Lazo
    from aa_numerico import FONDOS
    from equilibrio import borde_E, dFdE_E
    base = os.path.join(AQUI, '..', '..', 'exe', 'hadzic', 'lineal')
    print('Límite de L fijo (masa puntual, polE): lambda(omega) de lambda_borde.Lazo y de LazoL')
    for nombre in ('P_k2_a1', 'P_k3_a1', 'P_k2_a0.1'):
        npz = os.path.join(base, nombre + '_equilibrio.npz')
        d = np.load(npz)
        L0, jt, k, A = float(d['L0']), float(d['J_borde']), float(d['k_borde']), float(d['A'])
        ps = CubicSpline(d['r'], d['phi_self'])
        fondo = FONDOS[str(d['fondo'])]
        phi = lambda r: fondo(r) + ps(r)
        dphi = lambda r: 1.0/r**2 + ps(r, 1)                       # fondo puntual
        Et, T = borde_E('polE', d['E_t'], d['J_t'], jt, None)
        z = Lazo(npz)
        for sig, nL in ((1e-3, 3), (1e-2, 5)):
            xl, wl = np.polynomial.hermite_e.hermegauss(nL)
            Ls, gL = L0*(1 + sig*xl), wl/wl.sum()*L0               # sum gL = L0, pesos de L dL
            xs, ws = np.polynomial.legendre.leggauss(160); s, ws = 0.5*(xs + 1), 0.5*ws
            E, L, w, dfdE = [], [], [], []
            rt = np.linspace(float(d['r'][1]), float(d['r'][-1]), 20000)
            for Lm, gm in zip(Ls, gL):
                rc = np.interp(Lm**2, rt**3*dphi(rt), rt)
                emax = Et - (phi(rc) + Lm**2/(2*rc**2))
                e = emax*s**2
                E.append(Et - e); L.append(np.full_like(e, Lm)); w.append(gm*2*emax*s*ws)
                dfdE.append(A*dFdE_E('polE', Et - e, Et, T, k))
            E, L, w, dfdE = map(np.concatenate, (E, L, w, dfdE))
            zl = LazoL(phi, dphi, E, L, w, dfdE, 1.0001, 19.9, kmax=6, nth=128, nr=1200)
            fila = []
            for x in (0.5, 0.9, 1.0):
                om = x*z.Om_min
                fila.append(f'{z.lam(om):.6f} / {zl.lam(om):.6f}')
            print(f'  {nombre:10} sigma_L/L0 = {sig:g}, {nL} nodos en L:  omega/Omega_min = 0.5, 0.9, 1:  ' + '   '.join(fila)
                  + f'   Omega_min {z.Om_min:.6f} / {zl.Om_nodos:.6f}')


SAL = os.path.join(AQUI, '..', '..', 'exe', 'lambda_L')


def _tabla(nombre, modelos, etiqueta):
    os.makedirs(SAL, exist_ok=True)
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Omega_b: frecuencia de la órbita radial de energía E0 (el borde inferior de la banda). x_d = '
      '1 - omega_d/Omega_b; d_d = (Omega_b - omega_d)/(Omega_max - Omega_b).')
    w(f'{etiqueta:>6} {"R":>9} {"M":>9} {"masa/M-1":>9} {"Omega_b":>9} {"Omax/Omb":>8} {"lam(0)":>7} {"lam_borde":>9} '
      f'{"omega_d":>9} {"x_d":>9} {"d_d":>9}')
    for par, p in modelos:
        r = resultado(p)
        if r['omega_d'] is None:
            fin = f'{"--":>9} {"--":>9} {"--":>9}'
        else:
            fin = (f'{r["omega_d"]:9.5f} {1 - r["omega_d"]/r["Om_b"]:9.2e} '
                   f'{(r["Om_b"] - r["omega_d"])/(r["Om_max"] - r["Om_b"]):9.2e}')
        w(f'{par:6g} {r["R"]:9.4f} {r["M"]:9.5f} {p.masa()/p.M - 1:+9.1e} {r["Om_b"]:9.5f} {r["Om_max"]/r["Om_b"]:8.3f} '
          f'{r["lam0"]:7.4f} {r["lam_b"]:9.5f} {fin}')
    return lineas, w


def politropos(ks=(0.25, 0.5, 0.75, 1.0, 1.1, 1.15, 1.2, 1.22, 1.24, 1.25, 1.3, 1.5, 2.0, 2.5, 3.0)):
    """lambda_borde(k) de los politropos isótropos, el modo si lo hay, y el umbral k*."""
    lineas, w = _tabla('politropos', [(k, politropo(k)) for k in ks], 'k')
    ks_ = brentq(lambda k: (lambda p: p.lazo().lam(p.Omega_borde()) - 1)(politropo(k)), 1.2, 1.3, xtol=1e-8)
    w(f'umbral: lambda_borde(k*) = 1 en k* = {ks_:.5f};   12/pi^2 = {12/np.pi**2:.5f}')
    open(os.path.join(SAL, 'politropos.txt'), 'w').write('\n'.join(lineas) + '\n')


def modelos_king(kappas=(0.5, 1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0)):
    """lambda_borde(kappa) de los modelos de King y el umbral kappa*."""
    lineas, w = _tabla('king', [(c, king(c)) for c in kappas], 'kappa')
    g = lambda c: (lambda p: p.lazo().lam(p.Omega_borde()) - 1)(king(c))
    if g(kappas[0])*g(kappas[-1]) < 0:
        w(f'umbral: lambda_borde(kappa*) = 1 en kappa* = {brentq(g, kappas[0], kappas[-1], xtol=1e-7):.5f}')
    else:
        w('lambda_borde - 1 no cambia de signo en el intervalo')
    open(os.path.join(SAL, 'king.txt'), 'w').write('\n'.join(lineas) + '\n')


def dispersion(casos=('P_k0.75_a1', 'P_k1_a1', 'P_k0.75_a0.1'), anchos=(0.0025, 0.01, 0.03, 0.1), nL=24, ne=96, q=4.0):
    """¿Sobrevive a una dispersión en L la divergencia de lambda en el borde que da L fijo con
    k <= 1? f0 = F(E) g(L), con g parabólica de semiancho w L0 alrededor de L0 e int g L dL = L0,
    en el potencial total del equilibrio de L fijo (no se recalcula el potencial con la
    dispersión). Se evalúa lambda en omega = Omega_b (1 - d), con Omega_b la menor frecuencia
    de la órbita de energía E_t en el intervalo de L."""
    from lambda_borde import Lazo
    from equilibrio import borde_E, dFdE_E
    os.makedirs(SAL, exist_ok=True)
    base = os.path.join(AQUI, '..', '..', 'exe', 'hadzic', 'lineal')
    ds = (1e-2, 1e-3, 1e-4, 1e-5)
    lineas = []
    def w_(x=''):
        lineas.append(x); print(x, flush=True)
    w_('lambda(Omega_b (1 - d)) con d = 1e-2, 1e-3, 1e-4, 1e-5. Masa puntual, polE, L0 = 2.')
    xs, ws = np.polynomial.legendre.leggauss(ne); s, ws = 0.5*(xs + 1), 0.5*ws
    xl, wl = np.polynomial.legendre.leggauss(nL)
    for nombre in casos:
        npz = os.path.join(base, nombre + '_equilibrio.npz'); d = np.load(npz)
        L0, jt, k, A = float(d['L0']), float(d['J_borde']), float(d['k_borde']), float(d['A'])
        ps = CubicSpline(d['r'], d['phi_self'])
        phi = lambda r: -1.0/r + ps(r)
        dphi = lambda r: 1.0/r**2 + ps(r, 1)
        Et, T = borde_E('polE', d['E_t'], d['J_t'], jt, None)
        z = Lazo(npz, nj=400)
        w_(f'{nombre}: k = {k:g}, a0 = {float(d["a0"]):g}')
        w_(f'   L fijo                                              ' + ' '.join(f'{z.lam(z.Om_min*(1 - x)):7.3f}' for x in ds))
        rt = np.linspace(1.0001, 19.9, 20000); c = rt**3*dphi(rt)
        for wrel in anchos:
            Ls = L0*(1 + wrel*xl); gL = wl*(1 - xl**2)*Ls; gL = gL/gL.sum()*L0
            E, L, w, dfdE = [], [], [], []
            for Lm, gm in zip(Ls, gL):
                rc = np.interp(Lm**2, c, rt); emax = Et - (phi(rc) + Lm**2/(2*rc**2)); e = emax*s**q
                E.append(Et - e); L.append(np.full_like(e, Lm)); w.append(gm*q*emax*s**(q - 1)*ws)
                dfdE.append(A*dFdE_E('polE', Et - e, Et, T, k))
            E, L, w, dfdE = map(np.concatenate, (E, L, w, dfdE))
            zl = LazoL(phi, dphi, E, L, w, dfdE, 1.0001, 19.9, kmax=6, nth=128, nr=1200)
            Lb = L0*(1 + wrel*np.linspace(-1, 1, 41))
            Ob = orbitas(phi, dphi, Et - np.full_like(Lb, 1e-10), Lb, 1.0001, 19.9, 256)[0]
            w_(f'   semiancho {wrel:6.4f} L0, Omega_b(L) en [{Ob.min():.6f}, {Ob.max():.6f}]  '
               + ' '.join(f'{zl.lam(Ob.min()*(1 - x)):7.3f}' for x in ds))
    open(os.path.join(SAL, 'dispersion.txt'), 'w').write('\n'.join(lineas) + '\n')


def convergencia():
    """k* frente a cada parámetro numérico."""
    os.makedirs(SAL, exist_ok=True)
    lineas = [f'{"ne":>4} {"nx":>4} {"nth":>4} {"nr":>5} {"kmax":>4} {"q":>3} {"lam(1.2)":>9} {"lam(1.25)":>9} {"k*":>9}']
    print(lineas[0])
    pol = {k: politropo(k) for k in (1.2, 1.25)}
    for ne, nx, nth, nr, kmax, q in [(48, 24, 128, 400, 16, 2), (96, 48, 128, 400, 16, 2), (48, 24, 256, 400, 16, 2),
                                     (48, 24, 128, 800, 16, 2), (48, 24, 128, 1600, 16, 2), (48, 24, 128, 400, 6, 2),
                                     (48, 24, 128, 400, 10, 2), (48, 24, 128, 400, 24, 2), (48, 24, 128, 400, 16, 3),
                                     (96, 48, 256, 1600, 24, 2)]:
        kw = dict(ne=ne, nx=nx, nth=nth, nr=nr, kmax=kmax, q=q)
        l = [pol[k].lazo(**kw).lam(pol[k].Omega_borde()) for k in (1.2, 1.25)]
        ks_ = brentq(lambda k: (lambda p: p.lazo(**kw).lam(p.Omega_borde()) - 1)(politropo(k)), 1.2, 1.3, xtol=1e-8)
        lineas.append(f'{ne:4d} {nx:4d} {nth:4d} {nr:5d} {kmax:4d} {q:3g} {l[0]:9.6f} {l[1]:9.6f} {ks_:9.6f}')
        print(lineas[-1], flush=True)
    open(os.path.join(SAL, 'convergencia.txt'), 'w').write('\n'.join(lineas) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['validar', 'politropos', 'king', 'convergencia', 'dispersion'])
    a = ap.parse_args()
    {'king': modelos_king}.get(a.paso, globals().get(a.paso))()
