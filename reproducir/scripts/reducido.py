"""Modelo reducido de la nube cerca de la órbita circular (notas: docs/modelo_reducido).

Cerca de la órbita circular (radio R_0, frecuencia radial kappa_0) dos cascarones se atraen con la
energía m m' |xi - xi'|/(2 R_0^2), con xi el radio medido desde R_0. Al promediar sobre la oscilación
radial, con zeta = sqrt(2J) exp(i(Q - kappa_0 t)),

    d zeta_i/dt = i (Omega(J_i, L_i) - kappa_0) zeta_i + i g sum_j m_j K(zeta_i - zeta_j, d_ij),
    g = 1/(pi R_0^2 sqrt(kappa_0)),    K(z, d) = sqrt(max(|z|^2 - d^2, 0)) z/|z|^2,
    d_ij = sqrt(kappa_0) (R(L_i) - R(L_j)),    Omega(J, L) = (J + c(L))^-3.

Con L fija d = 0 y K = z/|z|. Se conservan H = sum m (E(J, L) - kappa_0 J) + g sum_{i<j} m m' G(|z|, d),
con G = sqrt(|z|^2 - d^2) + d asin(d/|z|) si |d| < |z| y (pi/2)|d| si no, y P = sum m J.
El dato inicial es el del PIC: la retícula de puntos medios en (Q, J) recortada, con el peso
exp(-sin^2(Q/2)/s_Q^2) J^2 exp(-J^2/s_J^2), por exp(-(L - L_0)^2/s_L^2) L si L está distribuido.

    python3 reducido.py serie referencia|familias|largas|L|Lb|control   corre lo que falte de una serie
    python3 reducido.py correr <nombre> <M_gas> [nj nq dt tfin esquema s_J s_Q isocrono|kepler N_L s_L]
    python3 reducido.py comparar [<reducido>:<corrida del PIC> ...]   tablas frente al PIC (bulto.py)
    python3 reducido.py tablas | figuras                        tablas y figuras de las notas
    python3 reducido.py largas                                  corridas hasta 100 tau_1, por ventanas

Las salidas quedan en exe/demo_eta/reducido/ (o en la carpeta de la variable REDUCIDO_DATOS). Todo
corre en serie con cuatro hilos fijados a los núcleos físicos.
"""
import os, sys, time, ctypes, subprocess, warnings
for _k, _v in (('OMP_NUM_THREADS', '4'), ('OMP_PLACES', 'cores'), ('OMP_PROC_BIND', 'close')):
    os.environ.setdefault(_k, _v)               # antes de cargar la biblioteca: los cuatro núcleos físicos
import numpy as np
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
warnings.filterwarnings('ignore', category=RuntimeWarning)
from demo_eta import envolvente, RAIZ, EXE
DATOS = os.environ.get('REDUCIDO_DATOS', os.path.join(EXE, 'demo_eta', 'reducido'))
DOC = os.path.join(RAIZ, 'docs', 'modelo_reducido')
C_L = 0.5*(2.0 + np.sqrt(8.0)); KAPPA = C_L**-3; LFIX = 2.0
_dphi = lambda r: r/(np.sqrt(1 + r*r)*(1 + np.sqrt(1 + r*r))**2)
RC = brentq(lambda r: r**3*_dphi(r) - LFIX**2, 1.0, 20.0)
G = 1.0/(np.pi*RC**2*np.sqrt(KAPPA))


def c_de_L(L, nombre='isocrono'):
    """c(L) de Omega(J, L) = (J + c(L))^-3: (L + sqrt(L^2 + 4))/2 en el isócrono, L en Kepler."""
    return 0.5*(L + np.sqrt(L*L + 4.0)) if nombre == 'isocrono' else 1.0*L


def R_de_L(L, nombre='isocrono'):
    """Radio de la órbita circular: sqrt(c^4 - 1) en el isócrono, L^2 en Kepler."""
    c = c_de_L(L, nombre)
    return np.sqrt(c**4 - 1.0) if nombre == 'isocrono' else c*c


def fondo(nombre='isocrono', l0=LFIX):
    """Constantes del potencial externo en L = l0: c de Omega(J) = (J + c)^-3 y h_0 = -1/(2 (J + c)^2),
    frecuencia radial kappa y radio r_c de la órbita circular, y acoplamiento g = 1/(pi r_c^2 sqrt(kappa))."""
    c, rc = c_de_L(l0, nombre), R_de_L(l0, nombre)
    kappa = c**-3
    return dict(c=c, kappa=kappa, rc=rc, g=1.0/(np.pi*rc**2*np.sqrt(kappa)))


def tau1(sj=0.1, nombre='isocrono'):
    """tau_1 = 2 pi/Delta Omega con la banda al 1 % del factor en J, como bulto.banda()."""
    import bulto as B
    if sj == 0.1 and nombre == 'isocrono':
        return B.TAU1
    f = lambda x: x*np.exp(1 - x) - 0.01
    jlo, jhi = sj*np.sqrt(brentq(f, 1e-9, 1.0)), sj*np.sqrt(brentq(f, 1.0, 50.0))
    c = fondo(nombre)['c']
    return 2*np.pi/((jlo + c)**-3 - (jhi + c)**-3)


SALIDA = 10.0                                   # como el PIC: una instantánea cada 10 unidades de tiempo
_dp = ctypes.POINTER(ctypes.c_double)
_fp = ctypes.POINTER(ctypes.c_float)


def biblioteca():
    so, src = os.path.join(AQUI, 'libreducido.so'), os.path.join(AQUI, 'reducido.c')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        subprocess.check_call(['gcc', '-O3', '-march=native', '-fopenmp', '-shared', '-fPIC', '-o', so, src, '-lm'])
    lib = ctypes.CDLL(so)
    lib.campo.argtypes = [ctypes.c_int] + [_dp]*5
    lib.energia.argtypes = [ctypes.c_int] + [_dp]*4
    lib.energial.argtypes = [ctypes.c_int] + [_dp]*5
    so, src = os.path.join(AQUI, 'libreducidof.so'), os.path.join(AQUI, 'reducidof.c')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        subprocess.check_call(['gcc', '-Ofast', '-march=native', '-fopenmp', '-shared', '-fPIC', '-o', so, src, '-lm'])
    libf = ctypes.CDLL(so)
    libf.campof.argtypes = [ctypes.c_int] + [_fp]*5
    libf.campofl.argtypes = [ctypes.c_int] + [_fp]*6
    return lib, libf


def reticula(nj, nq, umbral=1e-7, sj=0.1, sq=0.1):
    """Partículas (Q, J, peso) de la retícula recortada, en el orden del PIC (J fuera, Q dentro)."""
    import bulto as B
    J, Q = B.recorte_reticula(nj, nq, umbral, sj, sq)
    JJ, QQ = np.meshgrid(J, Q, indexing='ij')
    F = np.exp(-np.sin(0.5*QQ)**2/sq**2)*np.exp(-JJ**2/sj**2)*JJ**2
    return QQ.ravel(), JJ.ravel(), F.ravel()


def reticula_L(nj, nq, nl, sl, l0=LFIX, ancho=2.0, umbral=1e-7, sj=0.1, sq=0.1):
    """Retícula del bulto con momento angular distribuido, como el estado aa_quad de
    VlasovPoisson_PIC_sp con dftype = gauss: el peso de (Q, J) por exp(-(L - l0)^2/sl^2), en nl puntos
    medios de [l0 - ancho sl, l0 + ancho sl]. La masa de cada partícula lleva el factor L del
    elemento de volumen 8 pi^2 dQ dJ L dL. Devuelve Q, J, L y el peso (masa relativa)."""
    Q, J, F = reticula(nj, nq, umbral, sj, sq)
    Ln = l0 + ancho*sl*(2*(np.arange(nl) + 0.5)/nl - 1.0)
    gl = np.exp(-(Ln - l0)**2/sl**2)*Ln
    return (np.tile(Q, nl), np.tile(J, nl), np.repeat(Ln, len(Q)), (gl[:, None]*F[None, :]).ravel())


class Modelo:
    """Integra el modelo reducido. esquema = 'strang': rotación exacta con Omega(J) - kappa durante
    medio paso, flujo de la interacción con el punto medio, y otro medio paso de rotación (segundo
    orden, dos evaluaciones de la suma por paso). esquema = 'rk4': Runge-Kutta 4 sobre todo el campo
    (referencia). simple = True usa la suma sobre pares en precisión simple (error ~1e-6).
    Con L distribuido, efectos elige qué entra: 'ambos' (el modelo), 'frecuencia' (Omega(J, L_i) con el
    núcleo de interacción de L igual) o 'nucleo' (el núcleo con d_ij y Omega(J, L_0)). Los dos últimos
    son controles para ver cuál de los dos efectos de la dispersión de L deshace el núcleo de la nube."""

    def __init__(self, Q, J, w, masa, g=None, esquema='strang', simple=True, pot='isocrono', L=None, l0=LFIX,
                 efectos='ambos'):
        assert efectos in ('ambos', 'frecuencia', 'nucleo')
        self.lib, self.libf = biblioteca()
        F = fondo(pot, l0); self.c, self.kappa = F['c'], F['kappa']
        g = F['g'] if g is None else g
        self.L = L
        self.con_d = L is not None and efectos != 'frecuencia'      # ¿el núcleo de interacción lleva d_ij?
        if L is not None and efectos != 'nucleo':                   # frecuencia con la L de cada partícula
            self.c = c_de_L(np.asarray(L, float), pot)
        if L is not None:                       # u_i = sqrt(kappa_0) R(L_i)
            self.u = np.ascontiguousarray(np.sqrt(self.kappa)*R_de_L(np.asarray(L, float), pot))
            self.uf = self.u.astype(np.float32)
        self.n = len(Q)
        self.m = np.ascontiguousarray(masa*w/w.sum())
        self.mf = self.m.astype(np.float32)
        self.g, self.esquema, self.simple = g, esquema, simple
        self.x = np.ascontiguousarray(np.sqrt(2*J)*np.cos(Q))
        self.y = np.ascontiguousarray(np.sqrt(2*J)*np.sin(Q))
        self.sx, self.sy = np.empty(self.n), np.empty(self.n)
        self.sxf, self.syf = np.empty(self.n, np.float32), np.empty(self.n, np.float32)

    def _p(self, a):
        return a.ctypes.data_as(_dp)

    def suma(self, x, y):
        """s_i = sum_j m_j (Z_i - Z_j)/|Z_i - Z_j|."""
        if self.simple:
            xf, yf = x.astype(np.float32), y.astype(np.float32)
            f = lambda a: a.ctypes.data_as(_fp)
            if not self.con_d:
                self.libf.campof(self.n, f(xf), f(yf), f(self.mf), f(self.sxf), f(self.syf))
            else:
                self.libf.campofl(self.n, f(xf), f(yf), f(self.uf), f(self.mf), f(self.sxf), f(self.syf))
            return self.sxf.astype(float), self.syf.astype(float)
        assert not self.con_d, 'con L distribuido solo está la suma en precisión simple'
        x, y = np.ascontiguousarray(x), np.ascontiguousarray(y)
        self.lib.campo(self.n, self._p(x), self._p(y), self._p(self.m), self._p(self.sx), self._p(self.sy))
        return self.sx.copy(), self.sy.copy()

    def derivada(self, x, y):
        sx, sy = self.suma(x, y)
        w = (0.5*(x*x + y*y) + self.c)**-3 - self.kappa
        return -w*y - self.g*sy, w*x + self.g*sx

    def rota(self, tau):
        a = ((0.5*(self.x**2 + self.y**2) + self.c)**-3 - self.kappa)*tau
        c, s = np.cos(a), np.sin(a)
        self.x, self.y = c*self.x - s*self.y, s*self.x + c*self.y

    def inter(self, dt):
        sx, sy = self.suma(self.x, self.y)
        sx, sy = self.suma(self.x - 0.5*dt*self.g*sy, self.y + 0.5*dt*self.g*sx)
        self.x, self.y = self.x - dt*self.g*sy, self.y + dt*self.g*sx

    def paso(self, dt):
        if self.esquema == 'strang':
            self.rota(0.5*dt); self.inter(dt); self.rota(0.5*dt)
            return
        x, y = self.x, self.y
        k1x, k1y = self.derivada(x, y)
        k2x, k2y = self.derivada(x + 0.5*dt*k1x, y + 0.5*dt*k1y)
        k3x, k3y = self.derivada(x + 0.5*dt*k2x, y + 0.5*dt*k2y)
        k4x, k4y = self.derivada(x + dt*k3x, y + dt*k3y)
        self.x = x + dt*(k1x + 2*k2x + 2*k3x + k4x)/6
        self.y = y + dt*(k1y + 2*k2y + 2*k3y + k4y)/6

    def invariantes(self):
        u = np.empty(self.n)
        x, y = np.ascontiguousarray(self.x), np.ascontiguousarray(self.y)
        if not self.con_d:
            self.lib.energia(self.n, self._p(x), self._p(y), self._p(self.m), self._p(u))
        else:
            self.lib.energial(self.n, self._p(x), self._p(y), self._p(self.u), self._p(self.m), self._p(u))
        J = 0.5*(x**2 + y**2)
        h0 = -0.5*(J + self.c)**-2 - self.kappa*J
        return (self.m*h0).sum(), 0.5*self.g*(self.m*u).sum(), (self.m*J).sum()


def correr(nombre, masa, nj=400, nq=100, dt=5.0, tfin=10000.0, sj=0.1, sq=0.1, g=None, guardar=True, verboso=True,
           esquema='strang', simple=True, pot='isocrono', nl=1, sl=0.0, l0=LFIX, ancho=2.0, efectos='ambos'):
    if nl > 1:
        Q, J, L, w = reticula_L(nj, nq, nl, sl, l0, ancho, sj=sj, sq=sq)
        M = Modelo(Q, J, w, masa, g, esquema, simple, pot, L, l0, efectos)
    else:
        Q, J, w = reticula(nj, nq, sj=sj, sq=sq)
        M = Modelo(Q, J, w, masa, g, esquema, simple, pot, None, l0)
    npas = int(round(SALIDA/dt)); nsal = int(round(tfin/SALIDA))
    X = np.empty((nsal + 1, M.n), np.float32); Y = np.empty((nsal + 1, M.n), np.float32)
    inv = np.empty((nsal + 1, 3))
    X[0], Y[0] = M.x, M.y; inv[0] = M.invariantes()
    t0 = time.time()
    for k in range(1, nsal + 1):
        for _ in range(npas):
            M.paso(dt)
        X[k], Y[k] = M.x, M.y
        if k % 50 == 0 or k == nsal:
            inv[k] = M.invariantes()
            if verboso and k % 250 == 0:
                print(f'  {nombre}: t = {k*SALIDA:7.0f}  dH/U = {(inv[k,0]+inv[k,1]-inv[0,0]-inv[0,1])/inv[0,1]:+.2e}  dP/P = {inv[k,2]/inv[0,2]-1:+.2e}  ({time.time()-t0:.0f} s)', flush=True)
        else:
            inv[k] = np.nan
    t = SALIDA*np.arange(nsal + 1)
    if guardar:
        os.makedirs(DATOS, exist_ok=True)
        np.savez(os.path.join(DATOS, f'{nombre}.npz'), t=t, X=X, Y=Y, m=M.m, inv=inv,
                 par=np.array([masa, nj, nq, dt, tfin, sj, sq, M.g, fondo(pot, l0)['c'], M.kappa, nl, sl, l0, ancho,
                               ('ambos', 'frecuencia', 'nucleo').index(efectos)]),
                 L=(M.L if M.L is not None else np.array([l0])))
    return t, X, Y, M.m, inv


def medidas_plano(t, X, Y, m, sj, tau, kappa):
    """Medidas del estado de la nube en el plano de zeta = sqrt(2J) exp(i Q): X e Y son las partes
    real e imaginaria de zeta de cada partícula en cada instante (instantes x partículas) y m sus masas.
    Son las medidas que bulto.nucleo() hace sobre el PIC de L fija: h1 (peso 1 y peso
    J^2 exp(-J^2/s_J^2)), m02 (masa con J < 0.2 s_J), rho_c (densidad en J < 0.025 s_J en unidades de
    Fbar_max), mnuc (masa en el disco de acción 0.05 s_J centrado en el máximo de densidad, media en
    t >= 20 tau_1), ac (amplitud radial del centro de las partículas que están en ese disco en la última
    instantánea, |Z_c|/sqrt(kappa)) y dr (su ancho radial rms). Ninguna depende de una rotación común del
    plano, así que valen igual con zeta en el marco que rota con kappa_0 (el modelo reducido) y sin él
    (el PIC, también con L distribuido: VlasovPoisson_PIC_sp/reproducir/scripts/nube_L.py)."""
    W = m.sum()
    n0max = np.exp(-1.0)/(sj*np.sqrt(np.pi)/4)
    rad = np.sqrt(2*0.05*sj); b = np.arange(-1.8, 1.8 + rad, 0.5*rad)
    nt = len(t)
    h1 = np.empty(nt, complex); hB = np.empty(nt, complex); m02 = np.empty(nt); rho_c = np.empty(nt)

    def disco(k):
        x, y = X[k].astype(float), Y[k].astype(float)
        h = np.histogram2d(x, y, bins=(b, b), weights=m)[0]
        q = h.copy(); q[1:] += h[:-1]; q[:-1] += h[1:]
        q2 = q.copy(); q2[:, 1:] += q[:, :-1]; q2[:, :-1] += q[:, 1:]
        i, j = np.unravel_index(np.argmax(q2), q2.shape)
        x0, y0 = 0.5*(b[i] + b[i + 1]), 0.5*(b[j] + b[j + 1])
        for _ in range(8):
            sel = (x - x0)**2 + (y - y0)**2 < rad**2
            if m[sel].sum() == 0:
                break
            x0, y0 = np.average(x[sel], weights=m[sel]), np.average(y[sel], weights=m[sel])
        return (x - x0)**2 + (y - y0)**2 < rad**2, (x0, y0)

    for k in range(nt):
        x, y = X[k].astype(float), Y[k].astype(float)
        J = 0.5*(x*x + y*y); rr = np.sqrt(2*J) + 1e-300
        e1 = (x - 1j*y)/rr
        h1[k] = (m*e1).sum()/W
        Bw = m*J**2*np.exp(-J**2/sj**2)
        hB[k] = (Bw*e1).sum()/Bw.sum()
        m02[k] = m[J < 0.2*sj].sum()/W
        rho_c[k] = m[J < 0.025*sj].sum()/W/(0.025*sj)/n0max
    tarde = np.where(t >= 20*tau)[0]
    mnuc = np.mean([m[disco(k)[0]].sum() for k in tarde[::5]])/W
    sel, centro = disco(nt - 1)
    ms = m[sel]
    Zc = ((X[:, sel].astype(float) + 1j*Y[:, sel].astype(float)) @ ms)/ms.sum()
    d2 = ((np.abs(X[:, sel].astype(float) + 1j*Y[:, sel].astype(float) - Zc[:, None])**2) @ ms)/ms.sum()
    return dict(t=t, tau=tau, h1=h1, hB=hB, m02=m02, rho_c=rho_c, mnuc=mnuc, rho_p=mnuc/(0.05*sj*n0max),
                msel=ms.sum()/W, Zc=Zc, ac=np.abs(Zc)/np.sqrt(kappa), dr=np.sqrt(d2/(2*kappa)), centro=np.array(centro))


def diagnosticos(nombre, rehacer=False):
    """Las medidas de medidas_plano() sobre la corrida <nombre> del modelo reducido; se guardan en
    <nombre>_diag.npz."""
    sal = os.path.join(DATOS, f'{nombre}_diag.npz')
    if os.path.exists(sal) and not rehacer:
        return np.load(sal)
    d = np.load(os.path.join(DATOS, f'{nombre}.npz'))
    par = d['par']; sj = float(par[5])
    c, kappa = (float(par[8]), float(par[9])) if len(par) > 8 else (C_L, KAPPA)
    l0 = float(par[12]) if len(par) > 12 else LFIX
    tau = tau1(sj, 'kepler' if abs(c - l0) < 1e-12 else 'isocrono')
    np.savez(sal, **medidas_plano(d['t'], d['X'], d['Y'], d['m'], sj, tau, kappa))
    return np.load(sal)


TS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 20, 25, 30]


def comparar(pares):
    """pares: lista de (corrida del modelo reducido, corrida del PIC)."""
    import bulto as B
    print('amplitud radial del núcleo a_c en t/tau_1 =', TS)
    filas = []
    for red, pic in pares:
        D = diagnosticos(red); t, tau = D['t'], float(D['tau'])
        print(f'{red:10} ' + ' '.join(f'{D["ac"][np.argmin(np.abs(t - x*tau))]:.3f}' for x in TS))
        u = t >= 20*tau; v = (t >= 10*tau) & (t <= 20*tau)
        f = [red, np.abs(D['h1'][v]).mean(), np.abs(D['h1'][u]).mean(), np.abs(D['hB'][v]).mean(), D['rho_c'][u].mean(), D['m02'][u].mean(),
             float(D['mnuc']), D['dr'][u].mean(), D['ac'][u].mean()]
        filas.append(f)
        if pic and os.path.exists(B.ruta(pic, 'nucleo.npz')):
            nu = B.nucleo(pic); tp = nu['t']
            k = int(round(45/(tp[1] - tp[0])))
            media = np.array([nu['rc'][max(0, i - k):i + k + 1].mean() for i in range(len(tp))])
            env = envolvente(np.abs(nu['rc'] - media), tp)
            print(f'{"PIC " + pic:10} ' + ' '.join(f'{env[np.argmin(np.abs(tp - x*tau))]:.3f}' for x in TS))
            hk = np.load(B.ruta(pic, 'bulto.npz'))['hk']
            up = tp >= 20*tau; vp = (tp >= 10*tau) & (tp <= 20*tau)
            filas.append(['PIC ' + pic, np.abs(hk[vp, 1]).mean(), np.abs(hk[up, 1]).mean(), np.nan, nu['rho_c'][up].mean(), nu['m02'][up].mean(),
                          float(nu['mnuc']), nu['dr'][up].mean(), env[up].mean()])
    print('\nestado tardío: |h1| en 10-20 y >= 20 tau_1, |h1| con peso J^2 en 10-20, rho_c, m02, m_c, ancho dr y amplitud a_c en >= 20 tau_1')
    for f in filas:
        print(f'{f[0]:10} {f[1]:.2e} {f[2]:.2e} {f[3]:.2e}  rho_c {f[4]:6.2f}  m02 {f[5]:.3f}  m_c {f[6]:.3f}  dr {f[7]:.3f}  a_c {f[8]:.3f}')


def _ref(m):
    return f'R{int(round(m*1000)):03d}'


# Series: (nombre, M_gas, N_J, N_Q, dt, t_fin, opciones de correr(), corrida del PIC con la que se compara).
# La retícula es la de la corrida del PIC (recortada donde el peso baja de 1e-7 de su máximo).
SERIES = {
    'referencia': [(_ref(m), m, 400, 100, 2.5, 1e4, {}, p) for m, p in (
        (0.001, 'B1'), (0.01, 'B2q'), (0.015, 'B015'), (0.02, 'B02q'), (0.021, 'B021'), (0.022, 'B022'), (0.023, 'B023'),
        (0.024, 'B024'), (0.025, 'B025h'), (0.03, 'B3'), (0.04, 'B04'), (0.05, 'B05'), (0.07, 'B07'), (0.1, 'B4q'))],
    'familias': ([(f'rK_{m:g}', m, 400, 25, 1.25, 1e4, dict(pot='kepler'), f'K_{m:g}') for m in (0.01, 0.02, 0.04, 0.08)]
                 + [(f'rJ05_{m:g}', m, 400, 25, 2.5, 2e4, dict(sj=0.05), f'J05_{m:g}') for m in (0.005, 0.01, 0.02)]
                 + [(f'rJ20_{m:g}', m, 800, 25, 1.25, 1e4, dict(sj=0.2), f'J20_{m:g}') for m in (0.02, 0.04, 0.08)]
                 + [(f'rQ4_{m:g}', m, 400, 25, 2.5, 1e4, dict(sq=0.4), f'Q4_{m:g}') for m in (0.01, 0.02, 0.04, 0.08)]
                 + [(f'rQ2_{m:g}', m, 400, 50, 2.5, 1e4, dict(sq=0.2), f'Q2_{m:g}') for m in (0.01, 0.02, 0.04, 0.08)]),
    # hasta 100 tau_1: ¿muere la oscilación del núcleo?
    'largas': ([(_ref(m) + '_100', m, 400, 100, 2.5, 32730.0, {}, None) for m in (0.02, 0.025, 0.022, 0.03)]
               # las mismas con retículas más gruesas (200 x 50 y 300 x 75): dependencia con el número de partículas
               + [(_ref(m) + '_100' + x, m, nj, nq, 2.5, 32730.0, {}, None)
                  for x, nj, nq in (('c', 200, 50), ('m', 300, 75)) for m in (0.02, 0.025)]),
    # momento angular distribuido, exp(-(L - 2)^2/s_L^2), con N_L >= 480 s_L puntos en L (recurrencia en L
    # fuera de la corrida); retícula 200 x 50 y esquema de segundo orden para que el costo sea de minutos.
    'L': ([('rL0_0.03', 0.03, 200, 50, 5.0, 1e4, dict(esquema='strang'), None)]
          + [(f'rL{int(round(s*1000)):03d}_0.03', 0.03, 200, 50, 5.0, 1e4, dict(esquema='strang', nl=n, sl=s), None)
             for s, n in ((0.005, 4), (0.01, 6), (0.02, 10))]
          + [('rL020_0.01', 0.01, 200, 50, 5.0, 1e4, dict(esquema='strang', nl=10, sl=0.02), None),
             ('rL0_0.01', 0.01, 200, 50, 5.0, 1e4, dict(esquema='strang'), None),
             ('rL040_0.03', 0.03, 200, 50, 5.0, 1e4, dict(esquema='strang', nl=20, sl=0.04), None)]),
    # La serie 'L' usa N_L puntos en |L - L_0| <= 2 s_L con paso 0.005 a 0.008, y sus capas vuelven a estar en
    # fase en t_rec = 2 pi/(|dOmega/dL| dL) = 10417 a 16667. La recurrencia se nota desde t_rec - 70/s_L: los
    # valores tardíos (t >= 20 tau_1) de esa serie no son fiables con s_L >= 0.01 (2026-11-12; el paso
    # `reticulas` de VlasovPoisson_PIC_sp/reproducir/scripts/nube_L.py lo calcula). 'Lb' es la misma serie con
    # |L - L_0| <= 2.5 s_L y N_L = 8, 12, 20 y 30, sin recurrencia dentro de la corrida. Costo estimado de la
    # suma directa: 10 min, 21 min, 1 h y 2.2 h. No se ha corrido.
    'Lb': [(f'rL{int(round(s*1000)):03d}b_{m:g}', m, 200, 50, 5.0, 1e4, dict(esquema='strang', nl=n, sl=s, ancho=2.5), None)
           for m, s, n in ((0.03, 0.02, 20), (0.03, 0.01, 12), (0.03, 0.005, 8), (0.01, 0.02, 20), (0.03, 0.04, 30))],
    # Controles con la retícula de 'L' (basta para t <= 20 tau_1): solo el efecto de L en la frecuencia (f) y
    # solo el del núcleo de interacción (k). 15 min cada uno. No se han corrido.
    'control': [(f'rL020{x}_0.03', 0.03, 200, 50, 5.0, 1e4, dict(esquema='strang', nl=10, sl=0.02, efectos=e), None)
                for x, e in (('f', 'frecuencia'), ('k', 'nucleo'))],
}


def serie(nombre):
    for n, m, nj, nq, dt, tfin, op, _ in SERIES[nombre]:
        if os.path.exists(os.path.join(DATOS, n + '.npz')):
            continue
        t0 = time.time()
        correr(n, m, nj, nq, dt, tfin, **dict(dict(esquema='rk4'), **op))
        print(f'OK {n} ({time.time()-t0:.0f} s)', flush=True)


def medidas_pic(n):
    """Las medidas de la corrida n del PIC con las que se compara el modelo (bulto.nucleo)."""
    import bulto as B
    nu = B.nucleo(n); t = nu['t']; tau = B.banda(n)[3]
    k = int(round(45/(t[1] - t[0])))
    media = np.array([nu['rc'][max(0, i - k):i + k + 1].mean() for i in range(len(t))])
    env = envolvente(np.abs(nu['rc'] - media), t)
    hk = np.load(B.ruta(n, 'bulto.npz'))['hk']
    u = t >= 20*tau; v = (t >= 10*tau) & (t <= 20*tau)
    return dict(t=t/tau, ac=env, mu=B.escalas(n)[2], mc=float(nu['mnuc']), h1=np.abs(hk[v, 1]).mean(), m02=nu['m02'][u].mean(),
                dr=nu['dr'][u].mean())


def medidas(r):
    D = diagnosticos(r); t, tau = D['t'], float(D['tau'])
    u = t >= 20*tau; v = (t >= 10*tau) & (t <= 20*tau)
    return dict(t=t/tau, ac=envolvente(D['ac'], t), mc=float(D['mnuc']), h1=np.abs(D['h1'][v]).mean(), m02=D['m02'][u].mean(),
                dr=D['dr'][u].mean())


def centro_instantaneo(n, cada=2):
    """Amplitud radial del centro del disco de máxima densidad en cada instante, |centro|/sqrt(kappa_0), y
    masa en ese disco. No fija un conjunto de partículas, así que sirve para corridas largas, en las que
    el núcleo intercambia partículas con el halo. Devuelve t, amplitud, masa y N_ef = (sum m)^2/sum m^2."""
    sal = os.path.join(DATOS, f'{n}_centro.npz')
    if os.path.exists(sal):
        d = np.load(sal)
        return d['t'], d['a'], d['mc'], float(d['neff'])
    d = np.load(os.path.join(DATOS, f'{n}.npz')); t, m = d['t'], d['m']; X, Y = d['X'], d['Y']; par = d['par']
    sj = float(par[5]); kappa = float(par[9]); rad = np.sqrt(2*0.05*sj); b = np.arange(-1.8, 1.8 + rad, 0.5*rad)
    ks = np.arange(0, len(t), cada); a = np.empty(len(ks)); mc = np.empty(len(ks))
    for i, k in enumerate(ks):
        x, y = X[k].astype(float), Y[k].astype(float)
        h = np.histogram2d(x, y, bins=(b, b), weights=m)[0]
        q = h.copy(); q[1:] += h[:-1]; q[:-1] += h[1:]
        q2 = q.copy(); q2[:, 1:] += q[:, :-1]; q2[:, :-1] += q[:, 1:]
        ii, jj = np.unravel_index(np.argmax(q2), q2.shape)
        x0, y0 = 0.5*(b[ii] + b[ii + 1]), 0.5*(b[jj] + b[jj + 1])
        for _ in range(8):
            sel = (x - x0)**2 + (y - y0)**2 < rad**2
            if m[sel].sum() == 0:
                break
            x0, y0 = np.average(x[sel], weights=m[sel]), np.average(y[sel], weights=m[sel])
        a[i] = np.hypot(x0, y0)/np.sqrt(kappa); mc[i] = m[sel].sum()/m.sum()
    neff = m.sum()**2/(m**2).sum()
    np.savez(sal, t=t[ks], a=a, mc=mc, neff=neff)
    return t[ks], a, mc, neff


def largas():
    """Tabla de las corridas hasta 100 tau_1: amplitud del centro instantáneo del núcleo y |h_1| por ventanas."""
    tau = tau1()
    V = [(5, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100)]
    print('amplitud radial del centro instantáneo del núcleo y |h_1|, medias en ventanas de t/tau_1:', V)
    for n, *_ in sorted(SERIES['largas'], key=lambda x: (x[1], -x[2])):
        if not os.path.exists(os.path.join(DATOS, n + '.npz')):
            continue
        t, a, mc, neff = centro_instantaneo(n)
        D = diagnosticos(n); td = D['t']; h = np.abs(D['h1'])
        print(f'{n:10} N_ef = {neff:5.0f}  a_c: ' + ' '.join(f'{a[(t >= x*tau) & (t < y*tau)].mean():.3f}' for x, y in V))
        print(f'{"":10}             |h_1|: ' + ' '.join(f'{h[(td >= x*tau) & (td < y*tau)].mean():.3f}' for x, y in V)
              + f'   m_c en 20-30 y 90-100: {mc[(t >= 20*tau) & (t < 30*tau)].mean():.3f} {mc[t >= 90*tau].mean():.3f}')


def tabla_largas():
    """Filas de la tabla de las corridas hasta 100 tau_1 (tab_largas.tex): |h_1| medio por ventanas."""
    tau = tau1(); V = [(10, 20), (20, 30), (30, 40), (40, 60), (60, 80), (80, 100)]; filas = []
    for n, m, nj, nq, *_ in sorted(SERIES['largas'], key=lambda x: (x[1], x[2])):
        if not os.path.exists(os.path.join(DATOS, n + '.npz')):
            continue
        D = diagnosticos(n); t = D['t']; h = np.abs(D['h1'])
        mm = np.load(os.path.join(DATOS, n + '.npz'))['m']; neff = mm.sum()**2/(mm**2).sum()
        filas.append(f'${m:g}$ & ${nj}\\times {nq}$ & ${neff:.0f}$ & ' + ' & '.join(f'${h[(t >= a*tau) & (t < b*tau)].mean():.3f}$' for a, b in V) + ' \\\\')
    _macro('TabLargas', 'tab_largas.tex', filas)
    print('\n'.join(filas))


def tabla_L():
    """Filas de la tabla de las notas con momento angular distribuido (tab_L.tex) y la misma tabla en pantalla."""
    tau = tau1(); filas = []
    for n, m, nj, nq, dt, tfin, op, _ in sorted(SERIES['L'], key=lambda x: (-x[1], x[6].get('sl', 0.0))):
        if not os.path.exists(os.path.join(DATOS, n + '.npz')):
            continue
        D = diagnosticos(n); t = D['t']; h = np.abs(D['h1'])
        u = t >= 20*tau; v = (t >= 10*tau) & (t <= 20*tau)
        sl, nl = op.get('sl', 0.0), op.get('nl', 1)
        mc = float(D['mnuc'])
        print(f'{n:11} s_L = {sl:5.3f} N_L = {nl:2d} | |h1| {h[v].mean():.2e} {h[u].mean():.2e} | rho_c {D["rho_c"][u].mean():5.2f} '
              f'm02 {D["m02"][u].mean():.3f} m_c {mc:.3f} dr {D["dr"][u].mean():.3f}')
        filas.append(f'${m:g}$ & ${sl:g}$ & ${100*sl/np.sqrt(2)/LFIX:.2f}$ & ${nl}$ & {_sci(h[v].mean())} & ${D["rho_c"][u].mean():.2f}$ & '
                     f'${D["m02"][u].mean():.3f}$ & ${mc:.3f}$ \\\\')
    os.makedirs(DOC, exist_ok=True)
    _macro('TabL', 'tab_L.tex', filas)


def _sci(x):
    if x >= 0.0995:
        return f'${x:.2f}$'
    e = int(np.floor(np.log10(x))); a = x/10.0**e
    if round(a, 1) >= 10:
        a, e = a/10, e + 1
    return rf'${a:.1f}\times 10^{{{e}}}$'


def _macro(nombre, archivo, filas):
    """Escribe las filas de una tabla como la macro \\<nombre> (un \\input dentro de tabular deja una fila vacía)."""
    os.makedirs(DOC, exist_ok=True)
    cuerpo = '%\n'.join(filas).replace('$--$', '--')
    open(os.path.join(DOC, archivo), 'w').write('\\newcommand{\\' + nombre + '}{%\n' + cuerpo + '%\n\\hline}\n')


def tablas():
    """Filas de las tablas de las notas (docs/modelo_reducido/tab_*.tex): PIC y modelo reducido."""
    os.makedirs(DOC, exist_ok=True)
    f = lambda x: f'{x:.3f}'
    filas = []
    for r, m, *_, p in SERIES['referencia']:
        P, Q = medidas_pic(p), medidas(r)
        nuc = (f(P['mc']), f(Q['mc'])) if m >= 0.02 else ('--', '--')
        filas.append(f'${m:g}$ & ${P["mu"]:.2f}$ & {_sci(P["h1"])} & {_sci(Q["h1"])} & ${f(P["m02"])}$ & ${f(Q["m02"])}$ & ${nuc[0]}$ & ${nuc[1]}$ \\\\')
    _macro('TabReferencia', 'tab_referencia.tex', filas)
    filas = []; grupos = {'rK': 'Kepler', 'rJ05': r'$\sigma_J = 0.05$', 'rJ20': r'$\sigma_J = 0.2$', 'rQ2': r'$\sigma_Q = 0.2$',
                          'rQ4': r'$\sigma_Q = 0.4$'}
    orden = ['rK', 'rJ05', 'rJ20', 'rQ2', 'rQ4']
    for g in orden:
        primero = True
        for r, m, *_, p in SERIES['familias']:
            if r.split('_')[0] != g:
                continue
            P, Q = medidas_pic(p), medidas(r)
            hay = P['dr'] < 0.3 and Q['dr'] < 0.3                # hay un núcleo compacto en los dos
            nuc = (f(P['mc']), f(Q['mc'])) if hay else ('--', '--')
            filas.append(f'{grupos[g] if primero else ""} & ${m:g}$ & ${P["mu"]:.2f}$ & {_sci(P["h1"])} & {_sci(Q["h1"])} & ${f(P["m02"])}$ & ${f(Q["m02"])}$ & ${nuc[0]}$ & ${nuc[1]}$ \\\\')
            primero = False
        if g != orden[-1]:
            filas.append('\\hline')
    _macro('TabFamilias', 'tab_familias.tex', filas)
    print('tablas escritas en', DOC)


def figuras():
    """Figuras de las notas: oscilación del núcleo en el cruce y estado tardío de la familia de referencia.
    Dos series con los dos primeros colores de la paleta de referencia (validados con validate_palette.js
    sobre el fondo claro) y marcadores distintos."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sup, tinta, tinta2, reja, c_pic, c_red = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df', '#2a78d6', '#eb6834'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    os.makedirs(os.path.join(DOC, 'figuras'), exist_ok=True)

    def ejes(ax):
        ax.grid(True, axis='y', color=reja, lw=0.6); ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=3, width=0.6); ax.tick_params(which='minor', length=2, width=0.4)

    ref = SERIES['referencia']
    fig, ax = plt.subplots(2, 3, figsize=(10, 5.6), sharex=True, sharey=True)
    for a, (r, m, *_, p) in zip(ax.ravel(), ref[3:9]):
        P, Q = medidas_pic(p), medidas(r)
        a.plot(P['t'], P['ac'], color=c_pic, lw=1.6, solid_capstyle='round', label='PIC')
        a.plot(Q['t'], Q['ac'], color=c_red, lw=1.6, solid_capstyle='round', label='reduced model')
        a.set_yscale('log'); a.set_ylim(3e-3, 1.0); a.set_xlim(0, 30.6)
        a.set_xticks([0, 10, 20, 30]); a.set_yticks([0.01, 0.1, 1.0]); a.set_yticklabels(['0.01', '0.1', '1'])
        ejes(a)
        a.set_title(rf'$M_{{gas}}/M = {m:g}$   ($\mu = {P["mu"]:.2f}$)', loc='left', fontsize=9.5, color=tinta, pad=6)
    for a in ax[-1]:
        a.set_xlabel(r'$t/\tau_1$')
    for a in ax[:, 0]:
        a.set_ylabel('amplitude of the core oscillation')
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=1.6, h_pad=1.4)
    fig.legend(h, l, loc='upper right', ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.99, 1.0), handlelength=2.2)
    fig.savefig(os.path.join(DOC, 'figuras', 'reducido_series.pdf')); plt.close(fig)

    P = [medidas_pic(p) for *_, p in ref]; Q = [medidas(r) for r, *_ in ref]
    mu = np.array([p['mu'] for p in P])
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.3))
    for a, (k, titulo, log) in zip(ax, (('m02', r'mass fraction with $J < 0.2\,\sigma_J$', False),
                                        ('mc', 'mass fraction of the core', False),
                                        ('h1', r'$|h_1|/h_0$, mean over $10\leq t/\tau_1\leq 20$', True))):
        a.plot(mu, [p[k] for p in P], color=c_pic, lw=1.2, marker='o', ms=5.5, mec=sup, mew=0.8, label='PIC')
        a.plot(mu, [q[k] for q in Q], color=c_red, lw=1.2, marker='D', ms=4.8, mec=sup, mew=0.8, label='reduced model')
        if log:
            a.set_yscale('log')
        else:
            a.set_ylim(0, None)
        a.set_xlim(0, 0.95); a.set_xticks([0, 0.2, 0.4, 0.6, 0.8]); ejes(a)
        a.set_title(titulo, loc='left', fontsize=9.5, color=tinta, pad=6); a.set_xlabel(r'$\mu$')
    h, l = ax[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=1.8)
    fig.legend(h, l, loc='upper right', ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.99, 1.0))
    fig.savefig(os.path.join(DOC, 'figuras', 'reducido_estado.pdf')); plt.close(fig)
    print('figuras escritas en', os.path.join(DOC, 'figuras'))


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == 'correr':
        a = sys.argv[2:]
        nombre, masa = a[0], float(a[1])
        nj = int(a[2]) if len(a) > 2 else 400
        nq = int(a[3]) if len(a) > 3 else 100
        dt = float(a[4]) if len(a) > 4 else 2.5
        tfin = float(a[5]) if len(a) > 5 else 10000.0
        esquema = a[6] if len(a) > 6 else 'rk4'
        sj = float(a[7]) if len(a) > 7 else 0.1
        sq = float(a[8]) if len(a) > 8 else 0.1
        pot = a[9] if len(a) > 9 else 'isocrono'
        nl = int(a[10]) if len(a) > 10 else 1
        sl = float(a[11]) if len(a) > 11 else 0.0
        t0 = time.time()
        correr(nombre, masa, nj, nq, dt, tfin, sj, sq, esquema=esquema, pot=pot, nl=nl, sl=sl)
        print(f'OK {nombre} ({time.time()-t0:.0f} s)', flush=True)
    elif len(sys.argv) > 2 and sys.argv[1] == 'serie':
        serie(sys.argv[2])
    elif len(sys.argv) > 1 and sys.argv[1] == 'comparar':
        pares = [tuple(a.split(':')) if ':' in a else (a, None) for a in sys.argv[2:]]
        comparar(pares or [(r, p) for r, *_, p in SERIES['referencia'] + SERIES['familias']])
    elif len(sys.argv) > 1 and sys.argv[1] == 'tablas':
        tablas(); tabla_L(); tabla_largas()
    elif len(sys.argv) > 1 and sys.argv[1] == 'largas':
        largas()
    elif len(sys.argv) > 1 and sys.argv[1] == 'figuras':
        figuras()
    else:
        print(__doc__)
