"""Equilibrio autoconsistente de la Maxwelliana rebajada como ecuación diferencial
ordinaria: una construcción independiente de la iteración de equilibrio.py.

Con F = A E_gamma(g, (E_t - E)/T) para E < E_t, la densidad en un radio depende
solo del potencial en ese radio. Integrando en p_r, con E = p_r^2/2 + Phi_ef(r),

    4 pi r^2 rho(r) = C E_gamma(g + 1/2, psi(r)),   C = 8 pi^2 L0 A sqrt(2 pi T),
    psi(r) = (E_t - Phi_ef(r))/T,                  rho = 0 donde psi <= 0,

porque la integral de E_gamma(g, psi - p^2/2T) sobre |p| <= sqrt(2 T psi) vale
sqrt(2 pi T) E_gamma(g + 1/2, psi): término a término de la serie, con la función
beta B(g + n + 1, 1/2). La ecuación de Poisson queda como el sistema

    dM/dr = C E_gamma(g + 1/2, psi(r)),   du/dr = M/r^2,

con u = Phi_self - Phi_self(0) y psi = (Et - Phi_iso - L0^2/2r^2 - u)/T, donde
Et = E_t - Phi_self(0). La barrera centrífuga deja vacío el centro: en el vacío
interior M = u = 0, y fuera del soporte Phi_self = -M/r, lo que fija Phi_self(0).
Solo se integra sobre el soporte [r_-, r_+], donde psi > 0 y el lado derecho es
liso, así que el integrador nunca cruza el quiebre del borde.

Las incógnitas son tres números (Et, T, C), y las condiciones también:
    masa total M(r_+) = a0;
    T = (E_t - E_c)/W0, es decir psi = W0 en la órbita circular (el mínimo de
    Phi_ef), porque E_t - Phi_ef = T psi;
    J(E_t) = J_t: la órbita del borde, cuyos puntos de retorno son r_- y r_+, tiene
    acción (1/pi) int sqrt(2 T psi) dr = J_t.
Se resuelven con scipy.optimize.root (Powell híbrido). No hay iteración sobre el
campo: cada evaluación es una integración, con el error fijado por sus tolerancias.

Uso como módulo:
    eq = EquilibrioEDO(a0=0.075, jt=0.138, g=2, w0=3).resolver()
    eq.phi_self(r), eq.rho(r), eq.E_t, eq.T, eq.A
Solo vale para la Maxwelliana rebajada (F función de E). La familia hueca y la
gaussiana dependen de J, que no es local en el potencial.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, brentq


def E_gamma_escalar(s, x):
    """sum_n x^(s+n)/Gamma(s+n+1) para x > 0: la serie de equilibrio.E_gamma."""
    if x <= 0.0:
        return 0.0
    t = x**s/math.gamma(s + 1)
    suma = t
    n = 1
    while True:
        t *= x/(s + n)
        suma += t
        if t <= 1e-17*suma:
            return suma
        n += 1


def dphi_iso(r):
    """d Phi_iso/dr para Phi_iso = -1/(1 + sqrt(1 + r^2))."""
    s = np.sqrt(1.0 + r*r)
    return r/(s*(1.0 + s)**2)


class EquilibrioEDO:
    def __init__(self, a0, jt, g, w0, l0=2.0, rtol=1e-12, atol=1e-16, nodos=200):
        self.a0, self.jt, self.g, self.w0, self.l0 = a0, jt, g, w0, l0
        self.rtol, self.atol = rtol, atol
        self.xg, self.wg = np.polynomial.legendre.leggauss(nodos)

    # -------------------------------------------------------------- integración
    def _psi(self, r, u, Et, T):
        return (Et + 1.0/(1.0 + np.sqrt(1.0 + r*r)) - 0.5*self.l0**2/(r*r) - u)/T

    def integrar(self, Et, T, C):
        """Integra el soporte; devuelve un dict con la solución y sus derivados."""
        g5 = self.g + 0.5
        psi0 = lambda r: self._psi(r, 0.0, Et, T)
        rmin = 0.05
        # Vacío interior (M = u = 0): el soporte empieza donde psi(r, 0) cruza 0.
        rr = np.geomspace(rmin, 50.0, 4000)
        pos = np.flatnonzero(psi0(rr) > 0)
        if len(pos) == 0:
            return None
        r_menos = brentq(psi0, rr[pos[0] - 1], rr[pos[0]], xtol=1e-15, rtol=1e-15)

        def rhs(r, y):
            p = self._psi(r, y[1], Et, T)
            return [C*E_gamma_escalar(g5, p) if p > 0 else 0.0, y[0]/(r*r)]

        def salida(r, y):
            return self._psi(r, y[1], Et, T)
        salida.terminal, salida.direction = True, -1
        sol = solve_ivp(rhs, (r_menos, 60.0), [0.0, 0.0], method='DOP853', rtol=self.rtol,
                        atol=self.atol, dense_output=True, events=salida)
        if sol.status != 1:
            return None
        r_mas = float(sol.t_events[0][0])
        M, u_mas = (float(v) for v in sol.y_events[0][0])
        # Órbita circular: dPhi_ef/dr = dPhi_iso/dr + M(r)/r^2 - L0^2/r^3 = 0.
        dphief = lambda r: dphi_iso(r) + sol.sol(r)[0]/r**2 - self.l0**2/r**3
        r_c = brentq(dphief, r_menos, r_mas, xtol=1e-15, rtol=1e-15)
        psi_c = float(self._psi(r_c, sol.sol(r_c)[1], Et, T))
        # Acción de la órbita del borde: r = rm + ra sin(th) quita la singularidad.
        rm, ra = 0.5*(r_mas + r_menos), 0.5*(r_mas - r_menos)
        th = 0.5*np.pi*self.xg
        r = rm + ra*np.sin(th)
        psi = np.maximum(self._psi(r, sol.sol(r)[1], Et, T), 0.0)
        # (1/pi) int dr = (1/pi)(pi/2) sum_k w_k ra cos(th_k) (...)
        J = 0.5*np.sum(self.wg*np.sqrt(2*T*psi)*ra*np.cos(th))
        return dict(sol=sol, r_menos=r_menos, r_mas=r_mas, M=M, u_mas=u_mas, r_c=r_c,
                    psi_c=psi_c, J=J, Et=Et, T=T, C=C)

    def residuos(self, z):
        s = self.integrar(z[0], math.exp(z[1]), math.exp(z[2]))
        if s is None:
            return [1e3, 1e3, 1e3]
        return [s['M']/self.a0 - 1.0, s['psi_c']/self.w0 - 1.0, s['J']/self.jt - 1.0]

    # -------------------------------------------------------------- solución
    def inicial(self):
        """Punto de partida: el isócrono solo (u = 0) con la masa de a0."""
        c = 0.5*(self.l0 + math.sqrt(self.l0**2 + 4))
        E_t = -0.5/(self.jt + c)**2
        T = (E_t + 0.5/c**2)/self.w0
        s = self.integrar(E_t, T, 1e-30)
        return np.array([E_t, math.log(T), math.log(self.a0/(s['M']/1e-30))])

    def resolver(self, z0=None, verboso=False):
        z0 = self.inicial() if z0 is None else z0
        # El criterio es el residuo: con xtol tan chico MINPACK suele terminar con
        # 'no further improvement', que solo dice que llegó al redondeo.
        res = root(self.residuos, z0, method='hybr', options=dict(xtol=1e-14, maxfev=400))
        if np.max(np.abs(res.fun)) > 1e-10:
            # Continuación en la masa: partir de la solución con la mitad.
            if self.a0 < 1e-4:
                raise RuntimeError(f'EDO sin converger: {res.message}, residuos {res.fun}')
            medio = EquilibrioEDO(0.5*self.a0, self.jt, self.g, self.w0, self.l0,
                                  self.rtol, self.atol, len(self.xg)).resolver()
            res = root(self.residuos, medio.z, method='hybr', options=dict(xtol=1e-14, maxfev=400))
            if np.max(np.abs(res.fun)) > 1e-10:
                raise RuntimeError(f'EDO sin converger: {res.message}, residuos {res.fun}')
        self.z, self.residuo = res.x, np.max(np.abs(res.fun))
        s = self.integrar(res.x[0], math.exp(res.x[1]), math.exp(res.x[2]))
        self.__dict__.update(s)
        self.phi0 = -self.M/self.r_mas - self.u_mas      # Phi_self en el vacío interior
        self.E_t = self.Et + self.phi0
        self.A = self.C/(8*np.pi**2*self.l0*math.sqrt(2*np.pi*self.T))
        self.E_c = self.E_t - self.T*self.psi_c
        if verboso:
            print(f'  EDO: E_t = {self.E_t:.12e}, T = {self.T:.12e}, A = {self.A:.12e}, '
                  f'masa = {self.M:.12e}, soporte r = [{self.r_menos:.6f}, {self.r_mas:.6f}], '
                  f'r_c = {self.r_c:.6f}, residuo = {self.residuo:.1e}')
        return self

    def phi_self(self, r):
        r = np.asarray(r, float)
        out = np.full_like(r, self.phi0)
        sop = (r >= self.r_menos) & (r <= self.r_mas)
        out[sop] = self.phi0 + self.sol.sol(r[sop])[1]
        fuera = r > self.r_mas
        out[fuera] = -self.M/r[fuera]
        return out

    def rho(self, r):
        r = np.asarray(r, float)
        out = np.zeros_like(r)
        sop = (r > self.r_menos) & (r < self.r_mas)
        psi = self._psi(r[sop], self.sol.sol(r[sop])[1], self.Et, self.T)
        g5 = self.g + 0.5
        out[sop] = self.C*np.array([E_gamma_escalar(g5, p) for p in psi])/(4*np.pi*r[sop]**2)
        return out


def adoptar(eq, edo):
    """Pasa la solución de la EDO a un Equilibrio de equilibrio.py, para generar con
    él la condición inicial: Phi_self en su malla, la tabla J(E) de ese potencial
    (con la que se invierte el mapa) y los parámetros (E_t, T, A) de la EDO."""
    eq.phi_self = edo.phi_self(eq.r)
    E_tab, J_tab = eq.tabla_J_de_E(eq.mapa())
    eq.E_tab, eq.J_tab = E_tab, J_tab
    eq.E_borde, eq.T, eq.A = edo.E_t, edo.T, edo.A
    eq.E_t, eq.J_t = E_tab, J_tab
    eq.rho, eq.masa = edo.rho(eq.r), edo.M
    return eq
