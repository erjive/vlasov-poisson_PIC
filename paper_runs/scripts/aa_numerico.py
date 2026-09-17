"""Variables ángulo-acción numéricas para un potencial esférico cualquiera, con L fijo.

El potencial es el isócrono analítico más una corrección tabulada en una malla
radial (el potencial propio que guarda el código), interpolada con un spline
cúbico natural. Para cada partícula (r, p_r):

    E  = p_r^2/2 + Phi_ef(r),          Phi_ef = Phi + L^2/(2 r^2)
    r_-, r_+ : raíces de Phi_ef(r) = E, por bisección a cada lado del mínimo
    J  = (1/pi) Int_{r_-}^{r_+} sqrt(2(E - Phi_ef)) dr
    T  = 2 Int_{r_-}^{r_+} dr / sqrt(2(E - Phi_ef))
    Q  = (2 pi/T) Int_{r_-}^{r} dr'/sqrt(2(E - Phi_ef))      (p_r >= 0)
    Q  = 2 pi - (lo anterior)                                  (p_r < 0)

Las integrales se hacen con r = rm + ra sin(theta), que cancela la singularidad
inversa de raíz cuadrada en los puntos de retorno, y Gauss-Legendre.

Uso como módulo:  mapa = MapaAA(r_malla, phi_self_malla); Q, J, E = mapa(r, p)
"""
import numpy as np

L0 = 2.0


def phi_iso(r):
    return -1.0/(1.0 + np.sqrt(1.0 + r**2))


class SplineCubico:
    """Spline cúbico natural en una malla uniforme (sin scipy)."""

    def __init__(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        n = len(x); h = x[1] - x[0]
        assert np.allclose(np.diff(x), h)
        # Segundas derivadas M_i: sistema tridiagonal con M_0 = M_{n-1} = 0.
        a = np.full(n-2, 1.0); b = np.full(n-2, 4.0); c = np.full(n-2, 1.0)
        d = 6.0*(y[2:] - 2*y[1:-1] + y[:-2])/h**2
        for i in range(1, n-2):                      # Thomas
            m = a[i]/b[i-1]
            b[i] -= m*c[i-1]; d[i] -= m*d[i-1]
        M = np.zeros(n)
        M[n-2] = d[-1]/b[-1]
        for i in range(n-4, -1, -1):
            M[i+1] = (d[i] - c[i]*M[i+2])/b[i]
        self.x0, self.h, self.y, self.M, self.n = x[0], h, y, M, n

    def __call__(self, xq):
        xq = np.asarray(xq, float)
        i = np.clip(((xq - self.x0)//self.h).astype(int), 0, self.n-2)
        xi = self.x0 + i*self.h
        A = (xi + self.h - xq)/self.h; B = 1.0 - A
        return (A*self.y[i] + B*self.y[i+1]
                + ((A**3 - A)*self.M[i] + (B**3 - B)*self.M[i+1])*self.h**2/6.0)


class MapaAA:
    def __init__(self, r_malla=None, phi_self=None, L=L0, nodos=64):
        self.L = L
        self.self_ = None if phi_self is None else SplineCubico(r_malla, phi_self)
        if phi_self is not None:
            # Fuera de la tabla el potencial propio es kepleriano: -M/r.
            self.r_ult = float(r_malla[-1])
            self.M_self = -float(phi_self[-1])*self.r_ult
        x, w = np.polynomial.legendre.leggauss(nodos)
        self.x, self.w = x, w
        # Radio de la órbita circular: mínimo de Phi_ef, en malla fina y refinado.
        rr = np.linspace(0.5, 19.0, 200001)
        self.rc = rr[np.argmin(self.phi_ef(rr))]

    def phi(self, r):
        p = phi_iso(r)
        if self.self_ is None:
            return p
        r = np.asarray(r, float)
        return p + np.where(r <= self.r_ult, self.self_(np.minimum(r, self.r_ult)),
                            -self.M_self/np.maximum(r, 1e-300))

    def phi_ef(self, r):
        return self.phi(r) + 0.5*self.L**2/r**2

    def _raiz(self, E, a, b):
        """Bisección vectorizada de Phi_ef(r) = E en [a, b] (un cambio de signo)."""
        a = np.full_like(E, a); b = np.full_like(E, b)
        ga = self.phi_ef(a) - E
        for _ in range(80):
            m = 0.5*(a + b)
            gm = self.phi_ef(m) - E
            izq = np.sign(gm) == np.sign(ga)
            a = np.where(izq, m, a); ga = np.where(izq, gm, ga)
            b = np.where(izq, b, m)
        return 0.5*(a + b)

    def __call__(self, r, p, lote=4000):
        r = np.asarray(r, float); p = np.asarray(p, float)
        Q = np.empty_like(r); J = np.empty_like(r); E = np.empty_like(r)
        for s in range(0, len(r), lote):
            Q[s:s+lote], J[s:s+lote], E[s:s+lote] = self._lote(r[s:s+lote], p[s:s+lote])
        return Q, J, E

    def _lote(self, r, p):
        E = 0.5*p**2 + self.phi_ef(r)
        r1 = self._raiz(E, 1e-2, self.rc)
        r2 = self._raiz(E, self.rc, 60.0)
        rm, ra = 0.5*(r1 + r2), 0.5*(r2 - r1)
        x, w = self.x, self.w

        def integrales(th_hi):
            """Int_{-pi/2}^{th_hi} de sqrt(2(E-Phi)) ra cos y de ra cos/sqrt(2(E-Phi))."""
            lo = -0.5*np.pi
            mid, half = 0.5*(th_hi + lo), 0.5*(th_hi - lo)
            th = mid[:, None] + half[:, None]*x[None, :]
            rr = rm[:, None] + ra[:, None]*np.sin(th)
            v2 = 2.0*(E[:, None] - self.phi_ef(rr))
            v = np.sqrt(np.maximum(v2, 0.0))
            jac = ra[:, None]*np.cos(th)
            # El cociente jac/v es finito en los extremos; se evita 0/0 numérico.
            cociente = np.where(v > 0, jac/np.where(v > 0, v, 1.0), 0.0)
            I_p = half*np.sum(w*v*jac, axis=1)
            I_t = half*np.sum(w*cociente, axis=1)
            return I_p, I_t

        Ip, It = integrales(np.full_like(r, 0.5*np.pi))
        J = Ip/np.pi
        T = 2.0*It
        s = np.clip((r - rm)/np.where(ra > 0, ra, 1.0), -1.0, 1.0)
        _, t_r = integrales(np.arcsin(s))
        Q = 2.0*np.pi*t_r/T
        Q = np.where(p >= 0, Q, 2.0*np.pi - Q)
        return np.mod(Q, 2.0*np.pi), J, E
