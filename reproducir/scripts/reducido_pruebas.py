"""Pruebas del modelo reducido con momento angular distribuido (docs/modelo_reducido, Sección IV).

    python3 reducido_pruebas.py            todas (un minuto)
    python3 reducido_pruebas.py P7         solo las que se nombren

P1  El núcleo G(rho, d) es (pi/2) por el promedio de |d - rho cos(theta)| sobre theta.
P2  dG/drho = sqrt(1 - d^2/rho^2) si |d| < rho y 0 si no (el factor de la ecuación de movimiento).
P3  La suma campofl de reducidof.c (precisión simple) frente a la suma directa en doble precisión.
P4  Con todos los L iguales, campofl y energial dan lo mismo que las rutinas de L fija.
P5  La suma de la ecuación de movimiento es el gradiente de la energía de interacción.
P6  Una corrida corta con L distribuido conserva H y P, y el error baja al reducir el paso.
P7  Dos cascarones: las ecuaciones exactas (EOM de las notas, con el cruce resuelto como evento) frente
    al modelo reducido de los mismos dos cascarones, con L igual, con L distinto y cruce, y sin cruce.

P7 no depende del código PIC ni de la malla radial: pone a prueba el promedio sobre la oscilación
radial y el núcleo G con la dinámica exacta de un par. Resultado (2026-11-12): con L igual el modelo da
la acción intercambiada con un 2 % de error. Con L distinto la fase relativa de los dos cascarones gira
un 2 a 6 % más rápido en el modelo que en las ecuaciones exactas, y el error baja al 0.4 a 1.5 % si se
añade a la frecuencia de cada cascarón el corrimiento por la masa que encierra en promedio,
0.1454 m_otro [1/2 + asin(d/rho)/pi]: con L igual ese corrimiento es el mismo para todos (Eq.
UniformShift de las notas) y no afecta al modelo, con L distinto no lo es. Es el primer término que
el modelo de la Sección IV deja fuera. Cada prueba escribe PASA o FALLA; el código de salida es el
número de pruebas que fallan.
"""
import sys, ctypes
import numpy as np
from scipy.integrate import quad, solve_ivp
import reducido as R

_dp, _fp = R._dp, R._fp
FALLOS = []


def anota(nombre, pasa, texto):
    print(f'{"PASA" if pasa else "FALLA"} {nombre} {texto}', flush=True)
    if not pasa:
        FALLOS.append(nombre)


def G(rho, d):
    """Núcleo promediado, Eq. (Kernel) de las notas."""
    rho, d = np.asarray(rho, float), np.asarray(d, float)
    c = np.abs(d) < rho
    q = np.sqrt(np.where(c, rho*rho - d*d, 0.0))
    return np.where(c, q + d*np.arcsin(np.clip(d/np.where(rho > 0, rho, 1.0), -1, 1)), 0.5*np.pi*np.abs(d))


def suma_directa(x, y, u, m):
    """s_i = sum_j m_j sqrt(max(rho_ij^2 - d_ij^2, 0)) (Z_i - Z_j)/rho_ij^2, en doble precisión."""
    dx, dy, du = x[:, None] - x[None, :], y[:, None] - y[None, :], u[:, None] - u[None, :]
    r2 = dx*dx + dy*dy
    w = np.sqrt(np.maximum(r2 - du*du, 0.0))/np.where(r2 > 0, r2, 1.0)*m[None, :]
    return (dx*w).sum(1), (dy*w).sum(1)


def nube(n, semilla=1, ancho_u=0.05):
    g = np.random.default_rng(semilla)
    x, y = 0.3*g.standard_normal(n), 0.3*g.standard_normal(n)
    return x, y, ancho_u*g.standard_normal(n), g.uniform(0.5, 1.5, n)/n


def p1():
    g = np.random.default_rng(0); peor = 0.0
    for _ in range(200):
        rho, d = g.uniform(0.01, 2.0), g.uniform(-2.5, 2.5)
        cortes = [np.arccos(d/rho)] if abs(d) < rho else None
        prom = quad(lambda th: abs(d - rho*np.cos(th)), 0, np.pi, points=cortes, epsabs=1e-13, epsrel=1e-13, limit=200)[0]/np.pi
        peor = max(peor, abs(0.5*np.pi*prom - G(rho, d))/max(rho, abs(d)))
    anota('P1', peor < 1e-10, f'núcleo G frente al promedio directo: error relativo máximo {peor:.1e} en 200 pares (rho, d)')


def p2():
    g = np.random.default_rng(1); peor = 0.0; h = 1e-5
    for _ in range(200):
        rho, d = g.uniform(0.05, 2.0), g.uniform(-2.5, 2.5)
        if abs(abs(d) - rho) < 1e-3:
            continue
        num = (G(rho + h, d) - G(rho - h, d))/(2*h)
        peor = max(peor, abs(num - np.sqrt(max(1 - d*d/rho**2, 0.0))))
    anota('P2', peor < 1e-7, f'dG/drho frente a sqrt(1 - d^2/rho^2): error máximo {peor:.1e}')


def p3():
    lib, libf = R.biblioteca(); n = 1500
    x, y, u, m = nube(n)
    sx, sy = suma_directa(x, y, u, m)
    f = lambda a: a.ctypes.data_as(_fp)
    xf, yf, uf, mf = (a.astype(np.float32) for a in (x, y, u, m))
    ax, ay = np.empty(n, np.float32), np.empty(n, np.float32)
    libf.campofl(n, f(xf), f(yf), f(uf), f(mf), f(ax), f(ay))
    e = np.hypot(ax - sx, ay - sy).max()/np.hypot(sx, sy).max()
    anota('P3', e < 1e-5, f'campofl (simple) frente a la suma directa (doble), {n} partículas: error relativo máximo {e:.1e}')


def p4():
    lib, libf = R.biblioteca(); n = 1500
    x, y, _, m = nube(n); u = np.full(n, 1.2345)
    f = lambda a: a.ctypes.data_as(_fp); p = lambda a: a.ctypes.data_as(_dp)
    xf, yf, uf, mf = (a.astype(np.float32) for a in (x, y, u, m))
    ax, ay, bx, by = (np.empty(n, np.float32) for _ in range(4))
    libf.campofl(n, f(xf), f(yf), f(uf), f(mf), f(ax), f(ay))
    libf.campof(n, f(xf), f(yf), f(mf), f(bx), f(by))
    e1 = np.hypot(ax - bx, ay - by).max()/np.hypot(bx, by).max()
    ea, eb = np.empty(n), np.empty(n)
    lib.energial(n, p(x), p(y), p(u), p(m), p(ea)); lib.energia(n, p(x), p(y), p(m), p(eb))
    e2 = np.abs(ea - eb).max()/np.abs(eb).max()
    anota('P4', e1 < 2e-6 and e2 < 1e-12, f'L iguales: campofl frente a campof {e1:.1e}; energial frente a energia {e2:.1e}')


def p5():
    lib, _ = R.biblioteca(); n = 300
    x, y, u, m = nube(n, 2); p = lambda a: a.ctypes.data_as(_dp)

    def U(x, y):
        e = np.empty(n)
        lib.energial(n, p(np.ascontiguousarray(x)), p(np.ascontiguousarray(y)), p(u), p(m), p(e))
        return 0.5*(m*e).sum()

    sx, sy = suma_directa(x, y, u, m); h = 1e-6; peor = 0.0
    for i in range(0, n, 7):
        dx = np.zeros(n); dx[i] = h
        gx = (U(x + dx, y) - U(x - dx, y))/(2*h); gy = (U(x, y + dx) - U(x, y - dx))/(2*h)
        peor = max(peor, np.hypot(gx - m[i]*sx[i], gy - m[i]*sy[i])/(m[i]*np.hypot(sx, sy).max()))
    anota('P5', peor < 1e-6, f'gradiente de la energía frente a la suma de la ecuación de movimiento: error relativo máximo {peor:.1e}')


def p6():
    res = []
    for dt in (5.0, 2.5):
        t, X, Y, m, inv = R.correr('prueba', 0.03, nj=60, nq=16, dt=dt, tfin=1000.0, guardar=False, verboso=False,
                                   esquema='strang', nl=4, sl=0.01)
        ok = np.isfinite(inv[:, 0])
        H = inv[ok, 0] + inv[ok, 1]
        res.append((np.abs(H - H[0]).max()/abs(inv[0, 1]), np.abs(inv[ok, 2]/inv[0, 2] - 1).max(), X.shape[1]))
    (h1, q1, n), (h2, q2, _) = res
    anota('P6', h1 < 2e-3 and h2 < h1 and q1 < 1e-4,
          f'{n} partículas, N_L = 4, hasta t = 1000: |dH|/U = {h1:.1e} (dt = 5) y {h2:.1e} (dt = 2.5); |dP|/P = {q1:.1e} y {q2:.1e}')


# --- P7: dos cascarones, ecuaciones exactas frente al modelo reducido -------------------------------

def _dphi(r):
    s = np.sqrt(1 + r*r)
    return r/(s*(1 + s)**2)


def aa_isocrono(r, p, L):
    """Variables de ángulo-acción del isócrono (G = M = b = 1), con Q = 0 en el pericentro."""
    E = -1/(1 + np.sqrt(1 + r*r)) + 0.5*L*L/(r*r) + 0.5*p*p
    c = -0.5/E - 1.0
    e = np.sqrt(np.clip(1 - (L*L/c)*(1 + 1/c), 0, None))
    s = 1 + np.sqrt(1 + r*r)                      # s = 2 + c (1 - e cos eta)
    eta = np.arccos(np.clip((1 - (s - 2)/c)/np.where(e > 0, e, 1.0), -1, 1))
    eta = np.where(p >= 0, eta, 2*np.pi - eta)
    return np.mod(eta - e*c/(c + 1)*np.sin(eta), 2*np.pi), 1/np.sqrt(-2*E) - R.c_de_L(L)


def isocrono_rp(Q, J, L):
    E = -0.5/(J + R.c_de_L(L))**2
    c = -0.5/E - 1.0
    e = np.sqrt(max(1 - (L*L/c)*(1 + 1/c), 0.0)); k = e*c/(c + 1)
    eta = Q
    for _ in range(80):
        eta = eta - (eta - k*np.sin(eta) - Q)/(1 - k*np.cos(eta))
    s = 2 + c*(1 - e*np.cos(eta)); r = np.sqrt(s*(s - 2))
    p = np.sqrt(max(2*(E + 1/(1 + np.sqrt(1 + r*r)) - 0.5*L*L/(r*r)), 0.0))*(1.0 if np.sin(eta) >= 0 else -1.0)
    return r, p


def exacto(J, Q, L, m, tfin, dts):
    """Dos cascarones con las ecuaciones exactas: cada uno siente la masa del otro cuando lo tiene
    dentro. La fuerza salta en el cruce, así que cada tramo entre cruces se integra por separado
    (DOP853) y el cruce r_1 = r_2 se localiza como evento. Devuelve t, r, p (2 x nt) y el error
    relativo de la energía."""
    y = np.array([*isocrono_rp(Q[0], J[0], L[0]), *isocrono_rp(Q[1], J[1], L[1])])   # r1, p1, r2, p2

    def energia(y):
        r1, p1, r2, p2 = y
        return (m[0]*(0.5*p1*p1 + 0.5*L[0]**2/r1**2 - 1/(1 + np.sqrt(1 + r1*r1)))
                + m[1]*(0.5*p2*p2 + 0.5*L[1]**2/r2**2 - 1/(1 + np.sqrt(1 + r2*r2))) - m[0]*m[1]/max(r1, r2))

    def campo(fuera):                             # fuera: índice del cascarón exterior en este tramo
        def f(t, y):
            r1, p1, r2, p2 = y
            return [p1, L[0]**2/r1**3 - _dphi(r1) - (m[1]/r1**2 if fuera == 0 else 0.0),
                    p2, L[1]**2/r2**3 - _dphi(r2) - (m[0]/r2**2 if fuera == 1 else 0.0)]
        return f

    cruce = lambda t, y: y[0] - y[2]
    cruce.terminal = True
    ts = np.arange(0.0, tfin + 0.5*dts, dts); sal = np.empty((4, len(ts))); k = 0
    t = 0.0; fuera = 0 if y[0] > y[2] else 1; E0 = energia(y); ncruces = 0
    while t < tfin:
        cruce.direction = -1.0 if fuera == 0 else 1.0
        s = solve_ivp(campo(fuera), (t, tfin), y, method='DOP853', rtol=1e-11, atol=1e-13, events=cruce, dense_output=True)
        k1 = np.searchsorted(ts, s.t[-1], side='right')
        sal[:, k:k1] = s.sol(ts[k:k1]); k = k1
        t, y = s.t[-1], s.y[:, -1].copy()
        if s.status == 1:                         # cruce: cambia cuál es el exterior
            fuera = 1 - fuera; ncruces += 1
            y[0 if fuera == 0 else 2] += 1e-13    # un paso mínimo al otro lado para no repetir el evento
    return ts, sal[[0, 2]], sal[[1, 3]], abs(energia(y)/E0 - 1), ncruces


def reducido_par(J, Q, L, m, tfin, dts, corr=0.0):
    """El modelo reducido para los mismos dos cascarones (RK de orden 8, sin partículas vecinas).
    corr > 0 añade a la frecuencia de cada cascarón el corrimiento por la masa que encierra en
    promedio, corr m_otro [1/2 + asin(d/rho)/pi] (el término que el modelo no contiene)."""
    F = R.fondo(); kap, g = F['kappa'], F['g']
    c = R.c_de_L(np.asarray(L)); u = np.sqrt(kap)*R_de_L(np.asarray(L)); d = u[0] - u[1]

    def f(t, y):
        z = y[:2] + 1j*y[2:]
        w = (0.5*np.abs(z)**2 + c)**-3 - kap
        dz = z[0] - z[1]; rho2 = abs(dz)**2
        K = np.sqrt(max(rho2 - d*d, 0.0))*dz/rho2
        if corr:
            a = np.arcsin(np.clip(d/np.sqrt(rho2), -1, 1))/np.pi
            w = w + corr*np.array([m[1]*(0.5 + a), m[0]*(0.5 - a)])
        v = 1j*w*z + 1j*g*np.array([m[1]*K, -m[0]*K])
        return [*v.real, *v.imag]

    z0 = np.sqrt(2*np.asarray(J))*np.exp(1j*np.asarray(Q))
    ts = np.arange(0.0, tfin + 0.5*dts, dts)
    s = solve_ivp(f, (0, tfin), [*z0.real, *z0.imag], method='DOP853', rtol=1e-10, atol=1e-12, t_eval=ts)
    return ts, s.y[:2] + 1j*s.y[2:]


R_de_L = R.R_de_L


def p7():
    """Pares con masas iguales m = 0.002 y acciones J = 0.004 y 0.001 (4 y 1 % de sigma_J), desfasados
    un cuarto de vuelta, hasta t = 9000 (unas cien vueltas radiales). d = sqrt(kappa_0) [R(L_1) - R(L_2)]
    se da en unidades de rho_0 = |zeta_1 - zeta_2| inicial:
      a   L igual: la acción que intercambian y su curva;
      b   L distinto con cruce (d = 0.25 y 0.5 rho_0): acción intercambiada y giro de la fase relativa,
          con el modelo y con el modelo más el corrimiento por la masa encerrada;
      c   sin cruce (d mayor que |zeta_1| + |zeta_2|): ningún intercambio, ni exacto ni en el modelo.
    La acción y la fase exactas se promedian sobre una vuelta radial."""
    import reducido_frecuencias as RF
    F = R.fondo(); kap = F['kappa']; dts = 0.5; nv = int(round(2*np.pi/kap/dts)); tfin = 9000.0
    h = 1e-4; dR = (R_de_L(2 + h) - R_de_L(2 - h))/(2*h)
    m = [0.002, 0.002]; J = [0.004, 0.001]; Q = [0.3, 0.3 + 0.5*np.pi]
    z0 = np.sqrt(2*np.array(J))*np.exp(1j*np.array(Q)); rho0 = abs(z0[0] - z0[1])
    caja = lambda a: np.convolve(a, np.ones(nv)/nv, mode='valid')

    def par(x):
        dL = x*rho0/(np.sqrt(kap)*dR); L = [2.0 + 0.5*dL, 2.0 - 0.5*dL]
        t, r, p, dE, nc = exacto(J, Q, L, m, tfin, dts)
        Q1, J1 = aa_isocrono(r[0], p[0], L[0]); Q2, J2 = aa_isocrono(r[1], p[1], L[1])
        tm = caja(t); Jm = caja(J1); giro = np.ptp(caja(np.unwrap(Q1 - Q2))[[0, -1]])
        mod = []
        for corr in (0.0, RF.COEF_U):
            tr, z = reducido_par(J, Q, L, m, tfin, dts, corr)
            Jr = np.interp(tm, tr, 0.5*np.abs(z[0])**2)
            pr = np.interp(tm, tr, np.unwrap(np.angle(z[0]*np.conj(z[1]))))
            mod.append(dict(amp=np.ptp(Jr)/np.ptp(Jm) - 1, rms=np.sqrt(np.mean((Jm - Jr)**2))/np.ptp(Jm),
                            giro=np.ptp(pr[[0, -1]])/giro - 1, cte=np.ptp(0.5*np.abs(z)**2, axis=1).max()))
        return dict(dL=dL, dE=dE, nc=nc, amp=np.ptp(Jm), giro=giro, J2=np.abs(J2/J2[0] - 1).max(), mod=mod)

    e = par(0.0); a = e['mod'][0]
    anota('P7a', abs(a['amp']) < 0.05 and a['rms'] < 0.10,
          f'L igual: {e["nc"]} cruces, |dE/E| = {e["dE"]:.0e}; J_1 varía {e["amp"]:.2e}; modelo: amplitud {a["amp"]:+.1%}, '
          f'diferencia rms/amplitud {a["rms"]:.3f}')
    for x in (0.25, 0.5):
        e = par(x); a, b = e['mod']
        anota('P7b', abs(a['giro']) < 0.08 and abs(b['giro']) < 0.025 and abs(b['amp']) < 0.05,
              f'd = {x:g} rho_0 (dL = {e["dL"]:.4f}): {e["nc"]} cruces; fase relativa gira {e["giro"]:.1f} rad; modelo: giro {a["giro"]:+.1%}, '
              f'amplitud {a["amp"]:+.1%}; con la masa encerrada: giro {b["giro"]:+.1%}, amplitud {b["amp"]:+.1%}')
    x = 1.25*np.abs(z0).sum()/rho0
    e = par(x); a = e['mod'][0]
    anota('P7c', e['nc'] == 0 and e['J2'] < 1e-8 and a['cte'] < 1e-10,
          f'd = {x:.2f} rho_0 (dL = {e["dL"]:.4f}): {e["nc"]} cruces; la acción del cascarón interior cambia {e["J2"]:.0e} (exacto) '
          f'y las del modelo {a["cte"]:.0e}')


if __name__ == '__main__':
    pruebas = dict(P1=p1, P2=p2, P3=p3, P4=p4, P5=p5, P6=p6, P7=p7)
    for k in (sys.argv[1:] or pruebas):
        pruebas[k]()
    print(f'\n{len(FALLOS)} pruebas fallan' + (': ' + ', '.join(FALLOS) if FALLOS else ''))
    sys.exit(len(FALLOS))
