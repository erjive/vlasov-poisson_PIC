"""Equilibrio autoconsistente con L fijo y condición inicial para medir Landau.

Un equilibrio es F_eq(J), función solo de la acción VERDADERA: la del potencial
total Phi = Phi_iso + Phi_self[F_eq]. Como la acción depende del potencial y el
potencial de la distribución, se itera:

    Phi_self -> mapa (r,p_r) -> J -> F_eq(J) -> rho(r) -> Poisson -> Phi_self

La masa no depende del potencial: dr dp_r = dQ dJ (transformación canónica), así
que M = 8 pi^2 L0 * 2 pi * Int F(J) dJ y la amplitud de F_eq queda fija de entrada.
Solo cambia la forma del potencial.

Con el equilibrio convergido se colocan los nodos de cuadratura en una rejilla
regular de las (Q,J) verdaderas y se invierte el mapa numérico para obtener (r,p_r).
La perturbación es F = F_eq(J) [1 + eps cos Q], que no cambia la masa.

Uso:
    python3 equilibrio.py --a0 1e-2 --eps 0.1 --nrc 400 --npc 25 --salida exe/landau/ic.dat
"""
import os, sys, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico import MapaAA, phi_iso, L0

SIGMA_J = 0.10
J_MAX = 0.60                       # F_eq ~ J^2 exp(-J^2/0.01): e^-36 fuera


def F_forma(J):
    return J**2*np.exp(-J**2/SIGMA_J**2)


def amplitud(a0):
    Jq = np.linspace(0, J_MAX, 200001)
    return a0/(16*np.pi**3*L0*np.trapezoid(F_forma(Jq), Jq))


class Equilibrio:
    def __init__(self, a0, r_malla=np.arange(0.01, 25.0 + 1e-9, 0.01)):
        self.a0 = a0
        self.A = amplitud(a0)
        self.r = r_malla
        self.phi_self = np.zeros_like(r_malla)

    def mapa(self):
        return MapaAA(self.r, self.phi_self)

    def tabla_J_de_E(self, m, n=4000):
        """J(E) a L fijo: J solo depende de E. Se evalúa en el pericentro de
        órbitas que pasan por r_c con distintos p_r."""
        E0 = m.phi_ef(m.rc)
        Emax = self.E_de_J_aprox(m, J_MAX*1.05)
        E = np.linspace(E0, Emax, n)
        p = np.sqrt(2*np.maximum(E - E0, 0))
        _, J, _ = m(np.full_like(E, m.rc), p)
        J[0] = 0.0
        return E, J

    def E_de_J_aprox(self, m, J):
        # Isócrono como cota de partida; se amplía hasta que la tabla cubra J.
        c = 0.5*(L0 + np.sqrt(L0**2 + 4))
        E = -0.5/(J + c)**2
        for _ in range(20):
            _, Jt, _ = m(np.array([m.rc]), np.array([np.sqrt(2*max(E - m.phi_ef(m.rc), 0))]))
            if Jt[0] >= J:
                return E
            E = 0.5*E                  # E < 0: acercarse a cero sube J
        return E

    def densidad(self, m, E_t, J_t, npm=1201):
        """rho(r) 4 pi r^2 = 8 pi^2 L0 Int F(J(E(r,p))) dp."""
        Emax = E_t[-1]
        phief = m.phi_ef(self.r)
        pmax = np.sqrt(2*np.maximum(Emax - phief, 0.0))
        u = np.linspace(-1, 1, npm)
        P = pmax[:, None]*u[None, :]
        E = 0.5*P**2 + phief[:, None]
        J = np.interp(E, E_t, J_t, right=J_MAX*10)
        F = self.A*F_forma(J)
        integral = np.trapezoid(F, u, axis=1)*pmax
        return 8*np.pi**2*L0*integral/(4*np.pi*self.r**2)

    def poisson(self, rho):
        r = self.r
        dM = 4*np.pi*r**2*rho
        M = np.concatenate([[0], np.cumsum(0.5*(dM[1:] + dM[:-1])*np.diff(r))])
        g = 4*np.pi*r*rho
        ext = np.concatenate([np.cumsum((0.5*(g[1:] + g[:-1])*np.diff(r))[::-1])[::-1], [0]])
        return -M/r - ext, M[-1]

    def iterar(self, tol=1e-13, alfa=1.0, maxit=60, verboso=True):
        for it in range(maxit):
            m = self.mapa()
            E_t, J_t = self.tabla_J_de_E(m)
            rho = self.densidad(m, E_t, J_t)
            nuevo, M = self.poisson(rho)
            cambio = np.max(np.abs(nuevo - self.phi_self))
            self.phi_self = (1 - alfa)*self.phi_self + alfa*nuevo
            if verboso:
                print(f'  iteración {it:2d}: max|dPhi| = {cambio:.2e}   masa = {M:.10e}'
                      f'   Phi_self(r_c) = {np.interp(m.rc, self.r, self.phi_self):.6e}')
            if cambio < tol:
                break
        self.E_t, self.J_t, self.rho, self.masa = E_t, J_t, rho, M
        return self


def invertir(m, E_t, J_t, Q, J, nb=60, ng=32):
    """(Q,J) verdaderas -> (r,p_r) para el potencial del mapa m."""
    E = np.interp(J, J_t, E_t)
    r1 = m._raiz(E, 1e-2, m.rc)
    r2 = m._raiz(E, m.rc, 60.0)
    rm, ra = 0.5*(r1 + r2), 0.5*(r2 - r1)
    x, w = np.polynomial.legendre.leggauss(ng)

    def t_parcial(th_hi):
        lo = -0.5*np.pi
        mid, half = 0.5*(th_hi + lo), 0.5*(th_hi - lo)
        th = mid[:, None] + half[:, None]*x[None, :]
        rr = rm[:, None] + ra[:, None]*np.sin(th)
        v = np.sqrt(np.maximum(2*(E[:, None] - m.phi_ef(rr)), 0))
        jac = ra[:, None]*np.cos(th)
        coc = np.where(v > 0, jac/np.where(v > 0, v, 1), 0)
        return half*np.sum(w*coc, axis=1)

    T = 2*t_parcial(np.full_like(E, 0.5*np.pi))
    ida = Q <= np.pi
    objetivo = np.where(ida, Q, 2*np.pi - Q)*T/(2*np.pi)
    a = np.full_like(E, -0.5*np.pi); b = np.full_like(E, 0.5*np.pi)
    for _ in range(nb):
        c = 0.5*(a + b)
        bajo = t_parcial(c) < objetivo
        a = np.where(bajo, c, a); b = np.where(bajo, b, c)
    th = 0.5*(a + b)
    r = rm + ra*np.sin(th)
    p = np.sqrt(np.maximum(2*(E - m.phi_ef(r)), 0))
    return r, np.where(ida, p, -p)


def condicion_inicial(eq, eps, nrc, npc):
    m = eq.mapa()
    dJ = J_MAX/nrc; dQ = 2*np.pi/npc
    Jn = (np.arange(nrc) + 0.5)*dJ
    Qn = (np.arange(npc) + 0.5)*dQ
    JJ, QQ = np.meshgrid(Jn, Qn, indexing='ij')        # orden (i-1)*Npc+j del código
    JJ, QQ = JJ.ravel(), QQ.ravel()
    r = np.empty_like(JJ); p = np.empty_like(JJ)
    for s in range(0, len(JJ), 4000):
        r[s:s+4000], p[s:s+4000] = invertir(m, eq.E_t, eq.J_t, QQ[s:s+4000], JJ[s:s+4000])
    F = eq.A*F_forma(JJ)*(1 + eps*np.cos(QQ))
    return r, p, F, QQ, JJ


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--a0', type=float, default=1e-2)
    ap.add_argument('--eps', type=float, default=0.1)
    ap.add_argument('--nrc', type=int, default=400)
    ap.add_argument('--npc', type=int, default=25)
    ap.add_argument('--salida', required=True)
    arg = ap.parse_args()
    print(f'equilibrio: a0={arg.a0:g}, F_eq = A J^2 exp(-J^2/{SIGMA_J}^2), L0={L0}')
    eq = Equilibrio(arg.a0).iterar()
    r, p, F, Q, J = condicion_inicial(eq, arg.eps, arg.nrc, arg.npc)
    m = eq.mapa()
    Qc, Jc, _ = m(r, p)
    dQ = np.abs(np.angle(np.exp(1j*(Qc - Q))))
    peso = F > 1e-6*F.max()
    print(f'inversión: max|J - J_nodo| = {np.max(np.abs(Jc - J)[peso]):.1e},'
          f' max|Q - Q_nodo| = {dQ[peso].max():.1e}  (nodos con peso)')
    os.makedirs(os.path.dirname(os.path.abspath(arg.salida)), exist_ok=True)
    np.savetxt(arg.salida, np.column_stack([r, p, F]), fmt='%.17e')
    base = os.path.splitext(arg.salida)[0]
    np.savez(base + '_equilibrio.npz', r=eq.r, phi_self=eq.phi_self, rho=eq.rho,
             E_t=eq.E_t, J_t=eq.J_t, A=eq.A, a0=arg.a0, eps=arg.eps,
             nrc=arg.nrc, npc=arg.npc, J_max=J_MAX, sigma_J=SIGMA_J)
    print(f'escrito {arg.salida} ({len(r)} partículas) y {base}_equilibrio.npz')
