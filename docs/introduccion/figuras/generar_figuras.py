"""Genera todas las figuras de vlasov_intro.tex.

Uso, desde este directorio:

    python3 generar_figuras.py            # todas
    python3 generar_figuras.py orbita     # solo las que contengan 'orbita'

Las figuras analíticas (potencial, frecuencias, órbita, enrollamiento, h_k
exacto) no necesitan nada más. Las demás leen las corridas de exe/dfstudy y
exe/sg, que no se versionan; se regeneran con paper_runs/scripts/run_*.sh.

Cada figura imprime las cifras que el texto cita, para que documento y datos no
puedan divergir sin que se note.
"""
import os, sys, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..', '..'))
sys.path.insert(0, os.path.join(RAIZ, 'paper_runs', 'scripts'))
from df0 import df0, omega, rp_to_QJ, L0, c as C_ISO     # noqa: E402
from exact import make_hk                               # noqa: E402

EXE = os.path.join(RAIZ, 'exe')
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Estilo: tipografía que case con Computer Modern, trazos finos, rejilla tenue.
# Tres colores categóricos en orden fijo; negro para exacto/control.
AZUL, NARANJA, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
TINTA, GRIS, GRIS_CLARO = '#0b0b0b', '#52514e', '#b9b8b3'
plt.rcParams.update({
    'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 8.5,
    'axes.labelsize': 9, 'axes.titlesize': 9, 'legend.fontsize': 7.5,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.edgecolor': GRIS, 'axes.labelcolor': TINTA, 'xtick.color': GRIS,
    'ytick.color': GRIS, 'text.color': TINTA,
    'axes.grid': True, 'grid.color': '#e4e3df', 'grid.linewidth': 0.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'lines.linewidth': 1.2, 'legend.frameon': False,
    'figure.dpi': 150, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
})
ANCHO = 6.1   # pulgadas: el ancho de texto del documento


def guardar(fig, nombre):
    fig.savefig(os.path.join(AQUI, nombre + '.pdf'))
    plt.close(fig)
    print(f'  -> {nombre}.pdf')


def QTICKS(ax, eje='x'):
    t, l = [0, np.pi, 2*np.pi], ['0', r'$\pi$', r'$2\pi$']
    (ax.set_xticks if eje == 'x' else ax.set_yticks)(t)
    (ax.set_xticklabels if eje == 'x' else ax.set_yticklabels)(l)


# ---------------------------------------------------------------------------
# Isócrono con L fijo: todo lo analítico que usan las primeras figuras.
def phi_iso(r):
    return -1.0/(1.0 + np.sqrt(1.0 + r**2))


def phi_ef(r, L=L0):
    return phi_iso(r) + 0.5*L**2/r**2


def E_de_J(J):
    return -0.5/(J + C_ISO)**2


def domega(J):
    return -3.0/(J + C_ISO)**4


def radio_circular(L=L0):
    """Mínimo de Phi_ef, por bisección sobre la derivada."""
    d = lambda r: r/(np.sqrt(1+r**2)*(1+np.sqrt(1+r**2))**2) - L**2/r**3
    a, b = 0.1, 50.0
    for _ in range(200):
        m = 0.5*(a+b)
        a, b = (m, b) if d(m) < 0 else (a, m)
    return 0.5*(a+b)


def retorno(E, L=L0):
    """Puntos de retorno r-<r+ por bisección a cada lado del mínimo."""
    rc = radio_circular(L)
    g = lambda r: phi_ef(r, L) - E

    def bis(a, b):
        for _ in range(200):
            m = 0.5*(a+b)
            if np.sign(g(m)) == np.sign(g(a)):
                a = m
            else:
                b = m
        return 0.5*(a+b)
    return bis(1e-3, rc), bis(rc, 1e4)


def accion_y_periodo(E, L=L0, n=64):
    """J y T_r por cuadratura, con r = rm + ra sin(theta): la sustitución de
    la sección 9, que deja integrandos suaves en los puntos de retorno."""
    r1, r2 = retorno(E, L)
    rm, ra = 0.5*(r1+r2), 0.5*(r2-r1)
    x, w = np.polynomial.legendre.leggauss(n)
    th = 0.5*np.pi*x
    w = 0.5*np.pi*w
    r = rm + ra*np.sin(th)
    v = np.sqrt(np.maximum(2*(E - phi_ef(r, L)), 0.0))
    jac = ra*np.cos(th)
    J = np.sum(w*v*jac)/np.pi
    T = 2*np.sum(w*jac/v)
    return J, T


# ===========================================================================
def fig_potencial():
    print('[potencial] potencial efectivo y retrato de fases')
    rc = radio_circular()
    Js = [0.05, 0.15, 0.30]
    cols = [AZUL, NARANJA, AQUA]
    fig, (a, b) = plt.subplots(1, 2, figsize=(ANCHO, 2.5))

    r = np.linspace(0.9, 12, 800)
    a.plot(r, phi_ef(r), color=TINTA, label=r'$\Phi_{\rm ef}(r)$')
    a.plot(r, phi_iso(r), color=GRIS_CLARO, ls='--', lw=1,
           label=r'$\Phi_{\rm iso}(r)$')
    for J, col in zip(Js, cols):
        E = E_de_J(J)
        r1, r2 = retorno(E)
        a.hlines(E, r1, r2, color=col, lw=1.4)
        a.plot([r1, r2], [E, E], 'o', ms=3, color=col)
    a.plot(rc, phi_ef(rc), 'o', ms=3.5, color=TINTA)
    a.annotate('órbita circular ($J=0$)', xy=(rc, phi_ef(rc)),
               xytext=(rc+1.2, phi_ef(rc)-0.0065), fontsize=7.5, color=GRIS,
               arrowprops=dict(arrowstyle='-', lw=0.6, color=GRIS))
    a.annotate(r'$r_-$', xy=retorno(E_de_J(0.30))[0:1] + (E_de_J(0.30),),
               xytext=(-13, -3), textcoords='offset points', fontsize=8)
    a.annotate(r'$r_+$', xy=(retorno(E_de_J(0.30))[1], E_de_J(0.30)),
               xytext=(4, -9), textcoords='offset points', fontsize=8)
    a.set_xlim(0.9, 12)
    a.set_ylim(-0.098, -0.035)
    a.set_xlabel('$r$')
    a.set_ylabel('energía')
    a.legend(loc='upper right')
    a.set_title('(a) potencial efectivo, $L_0=2$', loc='left')

    for J, col in zip(Js, cols):
        E = E_de_J(J)
        r1, r2 = retorno(E)
        rr = r1 + (r2-r1)*0.5*(1 - np.cos(np.linspace(0, np.pi, 400)))
        pp = np.sqrt(np.maximum(2*(E - phi_ef(rr)), 0))
        b.fill_between(rr, -pp, pp, color=col, alpha=0.10, lw=0)
        b.plot(np.r_[rr, rr[::-1]], np.r_[pp, -pp[::-1]], color=col,
               label=f'$J={J:.2f}$')
        Jn, _ = accion_y_periodo(E)
        print(f'   J={J:.2f}: area/2pi por cuadratura = {Jn:.12f}')
    b.plot(rc, 0, 'o', ms=3.5, color=TINTA)
    b.annotate('', xy=(rc+1.2, 0.165), xytext=(rc-0.4, 0.165),
               arrowprops=dict(arrowstyle='->', lw=0.8, color=GRIS))
    b.text(rc+2.0, 0.158, 'sentido del flujo', fontsize=7, color=GRIS,
           ha='left')
    b.set_xlabel('$r$')
    b.set_ylabel('$p_r$')
    b.set_ylim(-0.19, 0.19)
    b.set_xlim(3.2, 13.5)
    b.legend(loc='lower right', ncol=1)
    b.set_title(r'(b) órbitas: área encerrada $=2\pi J$', loc='left')
    fig.tight_layout(w_pad=1.5)
    print(f'   r_circular = {rc:.6f}   Phi_ef(rc) = {phi_ef(rc):.8f}'
          f'   E(J=0) = {E_de_J(0):.8f}')
    guardar(fig, 'potencial')


# ===========================================================================
def area_poligono(x, y):
    """Fórmula del cordón de zapato (shoelace) para un polígono cerrado."""
    return 0.5*abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def fig_liouville():
    print('[liouville] una mancha de condiciones iniciales: el área se conserva')
    n = 20000
    s = np.linspace(0, 4, n, endpoint=False)
    def borde(x0, x1, y0, y1):          # contorno de un rectángulo, antihorario
        u = np.mod(s, 1)
        lado = s.astype(int)
        x = np.select([lado == 0, lado == 1, lado == 2, lado == 3],
                      [x0 + u*(x1-x0), np.full(n, x1), x1 - u*(x1-x0), np.full(n, x0)])
        y = np.select([lado == 0, lado == 1, lado == 2, lado == 3],
                      [np.full(n, y0), y0 + u*(y1-y0), np.full(n, y1), y1 - u*(y1-y0)])
        return x, y

    def yoshida(q, p, fuerza, dt, pasos):
        w1 = 1/(2 - 2**(1/3))
        w0 = -2**(1/3)*w1
        for _ in range(pasos):
            for w in (w1, w0, w1):
                p = p + 0.5*w*dt*fuerza(q)
                q = q + w*dt*p
                p = p + 0.5*w*dt*fuerza(q)
        return q, p

    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.6))
    cols = [AZUL, NARANJA, AQUA]
    # (a) oscilador armónico, omega=1: rotación rígida.
    q, p = borde(0.6, 1.2, -0.3, 0.3)
    A0 = area_poligono(q, p)
    ts = [0.0, np.pi/2, 3*np.pi/4]
    tprev = 0.0
    for t, col, lab in zip(ts, cols, ['$t=0$', r'$t=\pi/2$', r'$t=3\pi/4$']):
        q, p = yoshida(q, p, lambda x: -x, 0.001, int(round((t - tprev)/0.001)))
        tprev = t
        ax[0].fill(q, p, color=col, alpha=0.25, lw=0)
        ax[0].plot(np.r_[q, q[0]], np.r_[p, p[0]], color=col, lw=0.9,
                   label=lab)
        print(f'   armonico t={t:.3f}: area/area0 - 1 = {area_poligono(q, p)/A0-1:+.1e}')
    th = np.linspace(0, 2*np.pi, 400)
    for R in (0.6, 1.2*np.sqrt(1.0625)):
        ax[0].plot(R*np.cos(th), R*np.sin(th), color=GRIS_CLARO, lw=0.6, ls=':')
    ax[0].set_aspect('equal')
    ax[0].set_xlim(-1.45, 1.45)
    ax[0].set_ylim(-1.45, 1.45)
    ax[0].set_xlabel('$q$')
    ax[0].set_ylabel('$p$')
    ax[0].legend(loc='upper left', fontsize=6.5)
    ax[0].set_title('(a) oscilador armónico', loc='left')

    # (b) potencial efectivo del isócrono: se cizalla.
    F = lambda r: -r/(np.sqrt(1+r**2)*(1+np.sqrt(1+r**2))**2) + L0**2/r**3
    r, p = borde(6.3, 7.3, -0.025, 0.025)
    A0 = area_poligono(r, p)
    ts = [0, 150, 600]
    tprev = 0
    for t, col in zip(ts, cols):
        r, p = yoshida(r, p, F, 0.05, int(round((t - tprev)/0.05)))
        tprev = t
        ax[1].fill(r, p, color=col, alpha=0.25, lw=0)
        ax[1].plot(np.r_[r, r[0]], np.r_[p, p[0]], color=col, lw=0.8,
                   label=f'$t={t}$')
        print(f'   isocrono t={t}: area/area0 - 1 = {area_poligono(r, p)/A0-1:+.1e}')
    ax[1].plot(radio_circular(), 0, 'o', ms=3, color=TINTA)
    ax[1].set_xlabel('$r$')
    ax[1].set_ylabel('$p_r$')
    ax[1].legend(loc='lower right', fontsize=6.5)
    ax[1].set_title(r'(b) potencial efectivo del isócrono', loc='left')
    fig.tight_layout(w_pad=1.5)
    guardar(fig, 'liouville')


# ===========================================================================
def fig_frecuencias():
    print('[frecuencias] E(J), omega(J), omega\'(J) y validación por cuadratura')
    J = np.linspace(0, 1.0, 400)
    Jn = np.linspace(0.02, 0.98, 13)
    En = E_de_J(Jn)
    JT = np.array([accion_y_periodo(E) for E in En])
    errJ = np.max(np.abs(JT[:, 0] - Jn))
    errw = np.max(np.abs(2*np.pi/JT[:, 1] - omega(Jn))/omega(Jn))
    print(f'   max |J_cuadratura - J|        = {errJ:.2e}')
    print(f'   max |omega_cuad/omega - 1|    = {errw:.2e}')
    # Frecuencia epicíclica en la órbita circular: omega(0) debe ser kappa.
    rc, h = radio_circular(), 1e-3
    kap = np.sqrt((phi_ef(rc+h) - 2*phi_ef(rc) + phi_ef(rc-h))/h**2)
    print(f'   kappa(r_c) = {kap:.6f}   omega(0) = {omega(0.0):.6f}')
    # Isócrono: omega = (-2E)^{3/2} para todo L.
    print(f'   omega - (-2E)^1.5 = {np.max(np.abs(omega(J)-(-2*E_de_J(J))**1.5)):.1e}')
    print(f'   banda J in [0,0.35]: omega in [{omega(0.35):.4f}, {omega(0):.4f}]')

    fig, ax = plt.subplots(1, 3, figsize=(ANCHO, 2.15))
    for a in ax:
        a.axvspan(0, 0.35, color='#eef3fa', lw=0, zorder=0)
        a.set_xlabel('$J$')
        a.set_xlim(0, 1)
    ax[0].plot(J, E_de_J(J), color=AZUL, label='cerrada')
    ax[0].plot(Jn, En, 'o', ms=3.2, mfc='white', mec=TINTA, mew=0.8,
               label='cuadratura')
    ax[0].set_title('(a) $E(J)$', loc='left')
    ax[0].legend(loc='lower right', handlelength=1.5)
    ax[1].plot(J, omega(J), color=AZUL)
    ax[1].plot(JT[:, 0], 2*np.pi/JT[:, 1], 'o', ms=3.2, mfc='white', mec=TINTA,
               mew=0.8)
    ax[1].axhline(kap, color=GRIS_CLARO, ls='--', lw=0.9)
    ax[1].text(0.97, kap-0.0015, r'$\kappa(r_c)$', ha='right', va='top',
               fontsize=7.5, color=GRIS)
    ax[1].set_title(r'(b) $\omega(J)=\partial E/\partial J$', loc='left')
    ax[2].plot(J, domega(J), color=AZUL)
    ax[2].set_title(r"(c) $\omega'(J)$", loc='left')
    ax[2].text(0.175, -0.03, 'soporte\nde $F_0$', ha='center', fontsize=7,
               color=GRIS)
    fig.tight_layout(w_pad=0.8)
    guardar(fig, 'frecuencias')


# ===========================================================================
def fig_orbita():
    print('[orbita] r(t), Q(t), J(t) a lo largo de una órbita')
    J0 = 0.30
    E0 = E_de_J(J0)
    r1, _ = retorno(E0)
    fuerza = lambda r: (-r/(np.sqrt(1+r**2)*(1+np.sqrt(1+r**2))**2)
                        + L0**2/r**3)
    dt, T = 0.01, 2*2*np.pi/omega(J0)
    n = int(T/dt)
    r, p = r1, 0.0
    R, P = np.empty(n), np.empty(n)
    # Yoshida de cuarto orden, el mismo integrador de las corridas.
    w1 = 1/(2 - 2**(1/3))
    w0 = -2**(1/3)*w1
    for i in range(n):
        R[i], P[i] = r, p
        for w in (w1, w0, w1):
            p += 0.5*w*dt*fuerza(r)
            r += w*dt*p
            p += 0.5*w*dt*fuerza(r)
    t = dt*np.arange(n)
    Q, Jt = rp_to_QJ(R, P)
    print(f'   max |J(t)-J0| = {np.max(np.abs(Jt-J0)):.1e}')
    Qlin = np.mod(omega(J0)*t, 2*np.pi)
    dQ = np.angle(np.exp(1j*(Q - Qlin)))
    print(f'   max |Q(t) - omega t| = {np.max(np.abs(dQ)):.1e}')

    fig, ax = plt.subplots(3, 1, figsize=(ANCHO, 3.6), sharex=True,
                           gridspec_kw=dict(height_ratios=[1.2, 1.2, 0.8]))
    ax[0].plot(t, R, color=AZUL)
    ax[0].set_ylabel('$r$')
    Tr = 2*np.pi/omega(J0)
    for k in range(3):
        ax[0].axvline(k*Tr, color=GRIS_CLARO, lw=0.7, ls=':')
    ax[0].set_title(f'(a) radio: movimiento no uniforme ($J={J0}$, '
                    f'$T_r=2\\pi/\\omega={Tr:.1f}$)', loc='left')
    Qp = Q.copy()
    Qp[np.abs(np.diff(Q, append=Q[-1])) > np.pi] = np.nan
    ax[1].plot(t, Qp, color=AZUL, label='$Q$ desde $(r,p_r)$')
    tl = np.linspace(0, T, 2000)
    ql = np.mod(omega(J0)*tl, 2*np.pi)
    ql[np.abs(np.diff(ql, append=ql[-1])) > np.pi] = np.nan
    ax[1].plot(tl, ql, color=NARANJA, ls='--', lw=1, label=r'$\omega(J)\,t$')
    QTICKS(ax[1], 'y')
    ax[1].set_ylabel('$Q$')
    ax[1].legend(loc='upper left', ncol=2)
    ax[1].set_title('(b) ángulo: avanza uniformemente', loc='left')
    ax[2].plot(t, 1e14*(Jt - J0), color=AZUL)
    ax[2].set_ylabel(r'$10^{14}(J-J_0)$')
    ax[2].set_title('(c) acción: constante', loc='left')
    ax[2].set_xlabel('$t$')
    ax[2].set_xlim(0, T)
    fig.tight_layout(h_pad=0.4)
    guardar(fig, 'orbita')


# ===========================================================================
def fig_enrollamiento():
    print('[enrollamiento] F(Q,J,t) exacto en (Q,J) y en (r,p_r)')
    ts = [0, 400, 2000]
    Q = np.linspace(0, 2*np.pi, 700)
    J = np.linspace(0, 0.32, 500)
    QQ, JJ = np.meshgrid(Q, J)
    E = E_de_J(0.32)
    r1, r2 = retorno(E)
    r = np.linspace(r1*0.97, r2*1.02, 700)
    pm = np.sqrt(2*(E - phi_ef(radio_circular())))*1.05
    p = np.linspace(-pm, pm, 600)
    RR, PP = np.meshgrid(r, p)
    Qr, Jr = rp_to_QJ(RR, PP)
    fuera = ~np.isfinite(Jr) | (Jr > 0.32)

    fig, ax = plt.subplots(2, 3, figsize=(ANCHO, 3.9))
    for j, t in enumerate(ts):
        F = df0(np.mod(QQ - omega(JJ)*t, 2*np.pi), JJ, 'gauss')
        ax[0, j].pcolormesh(Q, J, F/F.max(), cmap='Blues', shading='auto',
                            rasterized=True, vmin=0, vmax=1)
        ax[0, j].set_title(f'$t={t}$', loc='left')
        QTICKS(ax[0, j])
        ax[0, j].set_xlabel('$Q$')
        ax[0, j].grid(False)
        G = df0(np.mod(Qr - omega(Jr)*t, 2*np.pi), Jr, 'gauss')
        G = np.where(fuera, np.nan, G)
        ax[1, j].pcolormesh(r, p, G/np.nanmax(G), cmap='Blues',
                            shading='auto', rasterized=True, vmin=0, vmax=1)
        ax[1, j].set_xlabel('$r$')
        ax[1, j].grid(False)
        if t > 0:
            dJ = 2*np.pi/(abs(domega(0.1))*t)
            print(f'   t={t}: ancho de franja 2pi/(|omega\'|t) = {dJ:.3f}')
    ax[0, 0].set_ylabel('$J$')
    ax[1, 0].set_ylabel('$p_r$')
    for j in (1, 2):
        ax[0, j].set_yticklabels([])
        ax[1, j].set_yticklabels([])
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    guardar(fig, 'enrollamiento')


# ===========================================================================
PRUEBA = dict(j0=0.10, sj=0.10, sq=0.40)    # la función de prueba Phi_1


def fig_hk_exacto():
    print('[hk_exacto] envolvente gaussiana del phase mixing, solución exacta')
    hk = make_hk('gauss', 1e-3, PRUEBA['j0'], PRUEBA['sj'], PRUEBA['sq'],
                 nJ=20001)
    t = np.linspace(0, 1000, 501)
    h0 = hk(0, [0.0])[0].real
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.4))
    for k, col in zip([1, 2, 3], [AZUL, NARANJA, AQUA]):
        h = np.abs(hk(k, t))/h0
        ax[0].semilogy(t, h, color=col, label=f'$k={k}$')
        m = t <= 250
        ax[1].plot(t[m]**2, np.log(h[m]/h[0]), color=col, label=f'$k={k}$')
        mm = (t > 0) & (t <= 150)
        b = np.polyfit(t[mm]**2, np.log(h[mm]/h[0]), 1)[0]
        print(f'   k={k}: pendiente frente a t^2 = {b:.3e}   '
              f'(razon a k=1: {b/b1 if k > 1 else 1:.2f})')
        if k == 1:
            b1 = b
    ax[0].set_ylim(1e-8, 2)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel('$|h_k(t)|/h_0$')
    ax[0].legend(loc='lower left')
    ax[0].set_title('(a) escala logarítmica en $t$', loc='left')
    ax[1].set_xlabel('$t^2$')
    ax[1].set_ylabel(r'$\ln\,|h_k(t)/h_k(0)|$')
    ax[1].ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax[1].set_title('(b) frente a $t^2$: rectas', loc='left')
    fig.tight_layout(w_pad=1.5)
    guardar(fig, 'hk_exacto')


# ===========================================================================
# Estudio agnóstico, fondo fijo
DFS = ['bimodal', 'spiral', 'king']
NS = [500, 5000, 50000]
COL_DF = dict(bimodal=AZUL, spiral=NARANJA, king=AQUA)
NOMBRE_DF = dict(bimodal='bimodal', spiral='espiral', king='King')


def leer_hk(ruta):
    a = np.loadtxt(os.path.join(ruta, 'hk1_complex.tl'))
    return a[:, 0], np.array([a[:, 1+2*k] + 1j*a[:, 2+2*k] for k in range(5)])


_EXACTOS = {}


def exacto_df(dft, t):
    if dft not in _EXACTOS:
        hk = make_hk(dft, 1e-4, PRUEBA['j0'], PRUEBA['sj'], PRUEBA['sq'],
                     nJ=40001)
        _EXACTOS[dft] = np.array([hk(k, t) for k in range(5)])
    return _EXACTOS[dft]


def fig_convergencia():
    print('[convergencia] Monte Carlo frente a cuadratura')
    D = os.path.join(EXE, 'dfstudy')
    t, _ = leer_hk(os.path.join(D, 'bimodal_quad_50000'))
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.5), sharey=True)
    for dft in DFS:
        EX = exacto_df(dft, t)
        h0 = EX[0].real[0]
        err = lambda d: np.mean(np.abs(leer_hk(os.path.join(D, d))[1][1]
                                       - EX[1]))/h0
        m = np.array([[err(f'{dft}_mcs{s}_{N}') for s in range(1, 6)]
                      for N in NS])
        mu, sd = m.mean(1), m.std(1, ddof=1)
        ax[0].errorbar(NS, mu, yerr=sd, fmt='-o', ms=3.5, capsize=2.5, lw=1,
                       color=COL_DF[dft], label=NOMBRE_DF[dft])
        pend = np.polyfit(np.log(NS), np.log(mu), 1)[0]
        q = [err(f'{dft}_quad_{N}') for N in NS]
        ax[1].loglog(NS, q, '-s', ms=3.5, lw=1, color=COL_DF[dft],
                     label=NOMBRE_DF[dft])
        print(f'   {dft:>8}: MC media {mu}  pendiente {pend:+.2f}')
        print(f'   {dft:>8}: cuadratura {np.array(q)}')
    ref = np.array(NS, float)
    ax[0].plot(ref, 6e-2*np.sqrt(NS[0]/ref), ':', color=TINTA, lw=1,
               label=r'$\propto N^{-1/2}$')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_title('(a) Monte Carlo, 5 semillas por punto', loc='left')
    ax[0].set_ylabel(r'error medio de $h_1$ $/\,h_0$')
    ax[1].axhline(1.2e-9, color=GRIS_CLARO, ls='--', lw=0.9)
    ax[1].text(9000, 4.5e-10, 'error del integrador', fontsize=7, color=GRIS,
               ha='center')
    ax[1].set_ylim(2e-10, 0.3)
    ax[1].set_title('(b) cuadratura en $(Q,J)$', loc='left')
    ax[0].legend(loc='lower left')
    ax[1].legend(loc='upper right')
    for a in ax:
        a.set_xlabel('$N$')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'convergencia')


def fig_agnosticas():
    print('[agnosticas] las tres predicciones falsables')
    D = os.path.join(EXE, 'dfstudy')
    fig, ax = plt.subplots(1, 3, figsize=(ANCHO, 2.3))

    t, h = leer_hk(os.path.join(D, 'bimodal_quad_50000'))
    EX = exacto_df('bimodal', t)
    h0 = EX[0].real[0]
    ax[0].semilogy(t, np.abs(h[2])/h0, color=AZUL, label='$k=2$ (señal)')
    ax[0].semilogy(t, np.abs(h[3])/h0, color=NARANJA, lw=0.7, label='$k=3$')
    ax[0].semilogy(t, np.abs(h[4])/h0, color=AQUA, lw=0.7, label='$k=4$')
    print(f'   bimodal: mediana |h3|/h0 = {np.median(np.abs(h[3]))/h0:.1e}'
          f'  |h4|/h0 = {np.median(np.abs(h[4]))/h0:.1e}')
    ax[0].set_ylim(1e-14, 1)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel('$|h_k|/h_0$')
    ax[0].legend(loc='center right', fontsize=6.5)
    ax[0].set_title('(a) bimodal: modos prohibidos', loc='left')

    t, h = leer_hk(os.path.join(D, 'spiral_quad_50000'))
    EX = exacto_df('spiral', t)
    h0 = EX[0].real[0]
    tst = 50.0/abs(domega(0.15))
    for k, col in zip([1, 2], [AZUL, NARANJA]):
        ax[1].plot(t, np.abs(EX[k])/h0, color=col, lw=1, label=f'exacto $k={k}$')
        ax[1].plot(t[::25], np.abs(h[k][::25])/h0, 'o', ms=2.2, color=col,
                   mfc='white', mew=0.7)
    ax[1].axvline(tst, color=GRIS, ls=':', lw=0.9)
    ax[1].text(tst+110, 0.29, f'$t^*={tst:.0f}$', fontsize=7.5, color=GRIS)
    tmax = t[np.argmax(np.abs(h[1]))]
    print(f'   espiral: t* predicho = {tst:.1f}   maximo medido de |h1| = {tmax:.0f}')
    ax[1].set_xlabel('$t$')
    ax[1].set_ylim(-0.02, 0.72)
    ax[1].legend(loc='upper right', fontsize=6.2, handlelength=1.4)
    ax[1].set_title('(b) espiral: se desenrolla', loc='left')

    t, _ = leer_hk(os.path.join(D, 'king_quad_500'))
    EX = exacto_df('king', t)
    h0 = EX[0].real[0]
    e = [np.mean(np.abs(leer_hk(os.path.join(D, f'king_quad_{N}'))[1][1]
                        - EX[1]))/h0 for N in NS]
    ax[2].loglog(NS, e, '-s', ms=3.5, color=AQUA, label='King, cuadratura')
    ax[2].loglog(NS, e[0]*(NS[0]/np.array(NS, float))**2, ':', color=TINTA,
                 lw=1, label=r'$\propto N^{-2}$')
    print(f'   King: factores {e[0]/e[1]:.0f} y {e[1]/e[2]:.0f}')
    ax[2].set_xlabel('$N$')
    ax[2].set_ylabel(r'error de $h_1$ $/\,h_0$')
    ax[2].legend(loc='lower left', fontsize=6.5)
    ax[2].set_title('(c) King: orden $h^2$', loc='left')
    fig.tight_layout(w_pad=1.4)
    guardar(fig, 'agnosticas')


# ===========================================================================
# Autogravedad
SG = os.path.join(EXE, 'sg')
CORRIDAS = [('sin autogravedad', 'quad_nosg', TINTA),
            ('con autogravedad, cuadratura', 'quad', AZUL),
            ('con autogravedad, Monte Carlo', 'mc', NARANJA)]


def fig_salud():
    import h5py
    print('[salud] energía (tal como la reporta el código y con 1/2 W_self) y h_0')
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.35))
    f = h5py.File(os.path.join(SG, 'long_fino', 'vlasov_output.h5'))
    rg = f['grid']['r'][:]
    st = sorted([s for s in f if s.startswith('step_')],
                key=lambda s: int(s.split('_')[1]))
    tt, Ed, Eb, Ec = [], [], [], []
    for s in st:
        x = f[s]
        r, p, w = x['r_part'][:], x['p_part'][:], x['f'][:]
        # El grid 'potential' es el total; la parte propia es lo que sobra del
        # isócrono, que se anula en el infinito como debe.
        ps = np.interp(r, rg, x['potential'][:] - phi_iso(rg))
        base = 0.5*p**2 + phi_ef(r)
        tt.append(x.attrs['time'])
        Ed.append(x.attrs['total_energy'])
        Eb.append(np.sum((base + ps)*w))       # lo que suma energy.f90
        Ec.append(np.sum((base + 0.5*ps)*w))   # la energía que se conserva
    tt, Ed, Eb, Ec = map(np.array, (tt, Ed, Eb, Ec))
    print(f'   recalculo del diagnostico del codigo: {Eb[-1]/Eb[0]-1:+.4e}'
          f'  (atributo {Ed[-1]/Ed[0]-1:+.4e})')
    print(f'   con 1/2 W_self: final {Ec[-1]/Ec[0]-1:+.2e}  '
          f'max {np.max(np.abs(Ec/Ec[0]-1)):.2e}')
    m = tt > 0
    ax[0].semilogy(tt[m], np.abs(Ed[m]/Ed[0] - 1), color=NARANJA, lw=1,
                   label=r'$\sum f\,(p_r^2/2+\Phi)$, como en el código')
    ax[0].semilogy(tt[m], np.abs(Ec[m]/Ec[0] - 1), color=AZUL, lw=1,
                   label=r'con $\frac{1}{2}$ en la autoenergía')
    ax[0].set_ylim(1e-9, 3e-3)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel(r'$|E(t)/E(0)-1|$')
    ax[0].legend(loc='center right', fontsize=6.3)
    ax[0].set_title('(a) energía, con autogravedad', loc='left')
    for nom, d, col in CORRIDAS:
        a = np.loadtxt(os.path.join(SG, d, 'hk1.tl'))
        ls = '-' if d != 'mc' else (0, (4, 2))
        ax[1].plot(a[:, 0], a[:, 1]/a[0, 1] - 1, color=col, ls=ls, lw=1.1,
                   label=nom)
        print(f'   {d:>10}: cambio de h0 = {a[-1,1]/a[0,1]-1:+.4e}')
    ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel(r'$h_0(t)/h_0(0)-1$')
    ax[1].legend(loc='center right', fontsize=6.3)
    ax[1].set_title('(b) cambio de $h_0$', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'salud')


def fig_envolvente():
    print('[envolvente] test gaussiana frente a exponencial')
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.5))
    for nom, d, col in CORRIDAS:
        t, h = leer_hk(os.path.join(SG, d))
        y = np.abs(h[1])/h[0].real[0]
        ls = '-' if d != 'mc' else (0, (4, 2))
        ax[0].semilogy(t, y, color=col, ls=ls, lw=1, label=nom)
        m = (t > 0) & (t < 250)
        r2 = np.corrcoef(t[m]**2, np.log(y[m]))[0, 1]**2
        r1 = np.corrcoef(t[m], np.log(y[m]))[0, 1]**2
        print(f'   {d:>10}: R2(t^2) = {r2:.5f}   R2(t) = {r1:.3f}')
        mm = t < 400
        ax[1].plot(t[mm]**2, np.log(y[mm]), color=col, ls=ls, lw=1.1, label=nom)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel('$|h_1|/h_0$')
    ax[0].set_title('(a) decaimiento y meseta', loc='left')
    ax[0].legend(loc='upper right', fontsize=6.5)
    ax[0].axhline(1.77e-3, color=GRIS_CLARO, ls='--', lw=0.8)
    ax[0].text(30, 1.0e-3, r'$1.77\times10^{-3}$', ha='left', va='top', fontsize=7,
               color=GRIS)
    ax[1].ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax[1].set_xlabel('$t^2$')
    ax[1].set_ylabel(r'$\ln\,(|h_1|/h_0)$')
    ax[1].set_title(r'(b) $t<400$ frente a $t^2$', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'envolvente')


def meseta(d):
    a = np.loadtxt(os.path.join(SG, d, 'hk1.tl'))
    return np.median(a[a[:, 0] > 1600, 2])/a[0, 1]


def fig_barridos():
    print('[barridos] la meseta frente a la discretización')
    fam = [('(a) malla $\\Delta r$', [('0.1', 'quad'), ('0.2', 'scan_dr_0.200'),
                                      ('0.05', 'scan_dr_0.050'),
                                      ('0.025', 'scan_dr_0.025')]),
           ('(b) orden del B-spline', [('1', 'quad'), ('2', 'scan_bspl_2'),
                                       ('3', 'scan_bspl_3')]),
           ('(c) nodos en $Q$', [('25', 'quad'), ('100', 'scan_npc_100'),
                                 ('200', 'scan_npc_200')]),
           ('(d) partículas $N$', [('$10^4$', 'quad'), ('$10^3$', 'scan_N_1000'),
                                   ('$10^5$', 'scan_N_100000')])]
    fig, ax = plt.subplots(1, 4, figsize=(ANCHO, 2.2), sharey=True)
    ctrl = np.loadtxt(os.path.join(SG, 'quad_nosg', 'hk1.tl'))
    for a, (tit, runs) in zip(ax, fam):
        a.semilogy(ctrl[:, 0], ctrl[:, 2]/ctrl[0, 1], color=GRIS_CLARO, lw=1,
                   label='sin autograv.')
        for (lab, d), col, ls in zip(runs, [TINTA, AZUL, NARANJA, AQUA],
                                     ['-', (0, (5, 2)), (0, (2, 1.5)),
                                      (0, (1, 1))]):
            x = np.loadtxt(os.path.join(SG, d, 'hk1.tl'))
            a.semilogy(x[:, 0], x[:, 2]/x[0, 1], color=col, ls=ls, lw=1,
                       label=lab)
            print(f'   {tit:<28} {lab:>7}: meseta {meseta(d):.4e}')
        a.set_title(tit, loc='left', fontsize=8)
        a.set_xlabel('$t$')
        a.set_ylim(1e-6, 1.5)
        a.set_xticks([0, 1000, 2000])
        a.legend(loc='lower left', fontsize=6, handlelength=2.2)
    ax[0].set_ylabel('$|h_1|/h_0$')
    fig.tight_layout(w_pad=0.3)
    guardar(fig, 'barridos')


LARGAS = [(1e-4, 'long_a0_1e-4', AQUA), (1e-3, 'long_a0_1e-3', AZUL),
          (1e-2, 'long_a0_1e-2', NARANJA)]


def mediana_en(t, y, a, b):
    m = (t >= a) & (t <= b)
    return np.median(y[m])


def fig_masa():
    print('[masa] la meseta frente a a_0, y el ruido del Monte Carlo')
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.5))
    A = np.array([x[0] for x in LARGAS])
    for (a, b), col, mk, nom in [((1600, 2000), AZUL, 'o', r'$t\in[1600,2000]$'),
                                 ((15000, 20000), NARANJA, 's',
                                  r'$t\in[15000,20000]$')]:
        ms = []
        for a0, d, _ in LARGAS:
            t, h = leer_hk(os.path.join(SG, d))
            m_ = (t >= a) & (t <= b)
            # Parte estática: el promedio complejo. La mediana de |h_1| la infla
            # la componente que gira (ruido de discreción) con a0=1e-2.
            ms.append(abs((h[1][m_]/h[0].real[0]).mean()))
        ms = np.array(ms)
        p = np.polyfit(np.log(A), np.log(ms), 1)[0]
        ax[0].loglog(A, ms, '-', marker=mk, ms=4, lw=1, color=col,
                     label=f'{nom}, pendiente ${p:+.2f}$')
        print(f'   ventana [{a},{b}]: mesetas {ms}  pendiente {p:+.3f}')
    ax[0].loglog(A, 1.77e-3*A/1e-3, ':', color=TINTA, lw=1,
                 label=r'$\propto a_0$')
    ax[0].set_xlabel('$a_0$')
    ax[0].set_ylabel(r'parte estática $|\langle h_1\rangle|/h_0$')
    ax[0].legend(loc='upper left', fontsize=6.5)
    ax[0].set_title('(a) escala con la masa', loc='left')

    # Monte Carlo frente a cuadratura: la señal queda bajo el ruido del MC.
    for d, N, ls in [('mc_1000', '10^3', (0, (1, 1))),
                     ('mc', '10^4', (0, (4, 2))),
                     ('mc_100000', '10^5', '-')]:
        x = np.loadtxt(os.path.join(SG, d, 'hk1.tl'))
        ax[1].semilogy(x[:, 0], x[:, 2], color=NARANJA, ls=ls, lw=0.9,
                       label=f'Monte Carlo, $N={N}$')
        print(f'   {d:>10}: meseta absoluta {np.median(x[x[:,0]>1600,2]):.2e}')
    x = np.loadtxt(os.path.join(SG, 'quad', 'hk1.tl'))
    ax[1].semilogy(x[:, 0], x[:, 2], color=AZUL, lw=1.3,
                   label='cuadratura, $N=10^4$')
    print(f'   cuadratura: meseta absoluta {np.median(x[x[:,0]>1600,2]):.2e}')
    ax[1].set_ylim(1e-9, 5e-6)
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel(r'$|h_1|$ sin normalizar')
    ax[1].legend(loc='upper right', fontsize=6.3)
    ax[1].set_title('(b) el Monte Carlo mide su ruido', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'masa')


def fig_fase():
    print('[fase] 200 periodos: amplitud, fase y plano complejo de h_1')
    fig, ax = plt.subplots(1, 3, figsize=(ANCHO, 2.35),
                           gridspec_kw=dict(width_ratios=[1.15, 1.15, 0.9]))
    for a0, d, col in LARGAS:
        t, h = leer_hk(os.path.join(SG, d))
        z = h[1]/h[0].real[0]
        y = np.abs(z)
        fase = np.unwrap(np.angle(z))
        lab = f'$a_0=10^{{{int(np.log10(a0))}}}$'
        ax[0].semilogy(t, y, color=col, lw=0.6, label=lab)
        ax[1].plot(t, fase, color=col, lw=1.1, label=lab)
        print(f'   a0={a0:.0e}:')
        for a, b in [(100, 400), (1300, 2000), (2000, 15000), (15000, 20000)]:
            m = (t >= a) & (t <= b)
            w = -np.polyfit(t[m], fase[m], 1)[0]
            S = z[m].mean()
            giro = np.sqrt(2)*np.std(z[m].real)
            sp = np.abs(np.fft.rfft(z[m].real - z[m].real.mean()))
            wp = (np.fft.rfftfreq(m.sum(), t[1]-t[0])*2*np.pi)[np.argmax(sp)]
            print(f'      [{a},{b}]: omega efectiva {w:+.2e}  fase recorre '
                  f'{np.ptp(fase[m]):6.2f} rad  mediana |h1|/h0 '
                  f'{np.median(y[m]):.4e}  parte estatica {abs(S):.4e} '
                  f'(arg {np.angle(S):+.3f})  parte que gira {giro:.2e} '
                  f'= {giro/abs(S):.1%} a omega {wp:.4f}')
    # Plano complejo, normalizado por la parte estática de cada masa.
    for a0, d, col, (a, b), lw in [
            (1e-2, 'long_a0_1e-2', NARANJA, (15000, 20000), 0.35),
            (1e-3, 'long_a0_1e-3', AZUL, (2000, 20000), 0.35)]:
        t, h = leer_hk(os.path.join(SG, d))
        z = h[1]/h[0].real[0]
        m = (t >= a) & (t <= b)
        S = np.abs(z[(t >= 2000) & (t <= 15000)].mean())
        ax[2].plot(z[m].real/S, z[m].imag/S, color=col, lw=lw,
                   label=f'$a_0=10^{{{int(np.log10(a0))}}}$, $t\\geq{a}$')
    ax[2].plot(0, 0, '+', color=TINTA, ms=6, mew=1)
    ax[2].plot(1, 0, 'o', color=TINTA, ms=2.5)
    ax[2].set_aspect('equal')
    ax[2].set_xlim(-1.5, 2.5)
    ax[2].set_ylim(-1.7, 2.7)
    ax[2].set_xlabel(r'${\rm Re}\,h_1/S$')
    ax[2].set_ylabel(r'${\rm Im}\,h_1/S$')
    ax[2].legend(loc='upper center', fontsize=5.8, handlelength=1.2)
    ax[2].set_title('(c) plano complejo', loc='left')
    w0 = omega(0.1)
    tl = np.linspace(0, 20000, 50)
    ax[1].plot(tl, -w0*tl, color=GRIS, ls='--', lw=0.9,
               label=rf'modo, $\omega={w0:.3f}$')
    ax[0].set_ylim(5e-5, 1.5)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel('$|h_1|/h_0$')
    ax[0].legend(loc='upper right', fontsize=6)
    ax[0].set_title('(a) amplitud', loc='left')
    ax[1].set_ylim(-260, 15)
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel(r'$\arg h_1$ (rad)')
    ax[1].legend(loc='lower left', fontsize=6)
    ax[1].set_title('(b) fase desenrollada', loc='left')
    for x in ax[:2]:
        x.set_xticks([0, 10000, 20000])
        x.set_xticklabels(['0', '10\u2009000', '20\u2009000'])
    fig.tight_layout(w_pad=0.8)
    guardar(fig, 'fase')


def fig_accion():
    import h5py
    print('[accion] J isócrona de partículas individuales, con autogravedad')
    f = h5py.File(os.path.join(SG, 'long_fino', 'vlasov_output.h5'))
    st = sorted([s for s in f if s.startswith('step_')],
                key=lambda s: int(s.split('_')[1]))
    idx = [4978, 7500, 9900]
    t = np.array([f[s].attrs['time'] for s in st])
    R = np.array([f[s]['r_part'][idx] for s in st])
    P = np.array([f[s]['p_part'][idx] for s in st])
    _, J = rp_to_QJ(R, P)
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.4),
                           gridspec_kw=dict(width_ratios=[1.6, 1]))
    for i, col in zip(range(3), [AZUL, NARANJA, AQUA]):
        x = J[:, i]
        Jm = x.mean()
        ax[0].plot(t, 100*(x/Jm - 1), color=col, lw=0.9,
                   label=rf'$\langle J\rangle={Jm:.3f}$')
        sp = np.abs(np.fft.rfft(x - Jm))
        fr = np.fft.rfftfreq(len(x), t[1]-t[0])*2*np.pi
        ax[1].plot(fr/omega(Jm), sp/sp.max(), color=col, lw=1)
        tr = t > 500
        d1 = x[(t > 500) & (t <= 1250)].mean()
        d2 = x[t > 1250].mean()
        print(f'   <J>={Jm:.4f}: pico a pico (t>500) {100*np.ptp(x[tr])/Jm:.2f}%  '
              f'deriva de la media {100*(d2-d1)/Jm:+.3f}%  '
              f'omega osc = {fr[np.argmax(sp)]:.4f}  orbital = {omega(Jm):.4f}')
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel(r'$J_{\rm iso}/\langle J\rangle-1$ (%)')
    ax[0].legend(loc='upper right', fontsize=6.5, ncol=3)
    ax[0].set_ylim(-0.75, 0.95)
    ax[0].set_title('(a) acción isócrona por partícula', loc='left')
    ax[1].axvline(1, color=GRIS_CLARO, ls='--', lw=0.9)
    ax[1].set_xlim(0, 3)
    ax[1].set_xlabel(r'frecuencia $/\,\omega(\langle J\rangle)$')
    ax[1].set_ylabel('espectro (normalizado)')
    ax[1].set_title('(b) espectro de la oscilación', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'accion')


# ===========================================================================
# Capturas de la evolución de las tres distribuciones de prueba (PIC real).
TIEMPOS_DF = dict(bimodal=[0, 200, 800, 2000], spiral=[0, 360, 720, 2000],
                  king=[0, 200, 800, 2000])


def fig_capturas():
    import h5py
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('peso', ['#dce9f7', AZUL, '#0b2e5c'])
    for dft in DFS:
        print(f'[capturas] {dft}')
        f = h5py.File(os.path.join(EXE, 'dfsnap', dft, 'vlasov_output.h5'))
        pasos = sorted([k for k in f if k.startswith('step_')],
                       key=lambda k: int(k.split('_')[1]))
        tiempos = np.array([f[k].attrs['time'] for k in pasos])
        w = f[pasos[0]]['f'][:]
        vis = w > 2e-3*w.max()             # los nodos de peso despreciable no se ven
        orden = np.argsort(w[vis])
        c = (w[vis]/w.max())[orden]
        datos = []
        for t in TIEMPOS_DF[dft]:
            k = pasos[int(np.argmin(np.abs(tiempos - t)))]
            r = f[k]['r_part'][:][vis][orden]
            p = f[k]['p_part'][:][vis][orden]
            Q, J = rp_to_QJ(r, p)
            datos.append((float(f[k].attrs['time']), r, p, Q, J))
            print(f'   t={datos[-1][0]:7.1f}: {vis.sum()} nodos visibles de {len(w)}')
        rr = np.concatenate([d[1] for d in datos])
        pp = np.concatenate([d[2] for d in datos])
        JJ = np.concatenate([d[4] for d in datos])
        fig, ax = plt.subplots(2, 4, figsize=(ANCHO, 3.35))
        kw = dict(c=c, cmap=cmap, vmin=0, vmax=1, s=0.25, lw=0, rasterized=True)
        for j, (t, r, p, Q, J) in enumerate(datos):
            ax[0, j].scatter(Q, J, **kw)
            ax[0, j].set_xlim(0, 2*np.pi)
            ax[0, j].set_ylim(0, JJ.max()*1.04)
            QTICKS(ax[0, j])
            ax[0, j].set_title(f'$t={t:.0f}$', loc='left')
            ax[0, j].set_xlabel('$Q$', labelpad=1)
            ax[1, j].scatter(r, p, **kw)
            ax[1, j].set_xlim(rr.min() - 0.2, rr.max() + 0.2)
            ax[1, j].set_ylim(pp.min()*1.08, pp.max()*1.08)
            ax[1, j].set_xlabel('$r$', labelpad=1)
            for x in ax[:, j]:
                x.grid(False)
                if j > 0:
                    x.set_yticklabels([])
        ax[0, 0].set_ylabel('$J$')
        ax[1, 0].set_ylabel('$p_r$')
        fig.tight_layout(h_pad=0.5, w_pad=0.3)
        guardar(fig, f'capturas_{dft}')


# ===========================================================================
def fig_mapa_numerico():
    print('[mapa_numerico] h_1 con el mapa del isócrono y con el del potencial real')
    d = np.load(os.path.join(SG, 'long20k_snap', 'aa_meseta.npz'))
    t = d['t']
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.5),
                           gridspec_kw=dict(width_ratios=[1.3, 1]))
    for k, col, lab in [('iso', AZUL, 'mapa del isócrono'),
                        ('promedio', NARANJA, 'mapa numérico')]:
        ax[0].semilogy(t, np.abs(d[f'h1_{k}']), color=col, lw=0.8, label=lab)
        v = t >= 2000
        z = d[f'h1_{k}'][v]
        print(f'   {k:>9}: |<h1>|/h0 = {abs(z.mean()):.2e}   rms parte que gira = '
              f'{np.sqrt(np.mean(np.abs(z - z.mean())**2)):.3e}')
    ax[0].set_ylim(1e-6, 1.5)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel('$|h_1|/h_0$')
    ax[0].legend(loc='upper right', fontsize=6.5)
    ax[0].set_xticks([0, 10000, 20000])
    ax[0].set_xticklabels(['0', '10\u2009000', '20\u2009000'])
    ax[0].set_title('(a) la meseta desaparece', loc='left')
    # (b) con long_fino, instantáneas cada 4 unidades: la oscilación orbital
    # (periodo ~125) queda resuelta. Mismo potencial promediado que en (a).
    import h5py
    from aa_numerico import MapaAA
    g = h5py.File(os.path.join(SG, 'long20k_snap', 'vlasov_output.h5'))
    pas = sorted([k for k in g if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    tg = np.array([g[k].attrs['time'] for k in pas])
    rg = g['grid']['r'][:]
    ps = np.mean([g[k]['potential'][:] - phi_iso(rg)
                  for k, tt in zip(pas, tg) if tt >= 2000], axis=0)
    mapa = MapaAA(rg, ps)
    h = h5py.File(os.path.join(SG, 'long_fino', 'vlasov_output.h5'))
    pas = sorted([k for k in h if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    tf = np.array([h[k].attrs['time'] for k in pas])
    sel = tf >= 1500
    r = np.array([h[k]['r_part'][4978] for k, s_ in zip(pas, sel) if s_])
    p = np.array([h[k]['p_part'][4978] for k, s_ in zip(pas, sel) if s_])
    _, Ji = rp_to_QJ(r, p)
    _, Jn, _ = mapa(r, p)
    for x, col, lab in [(Ji, AZUL, 'isócrona'), (Jn, NARANJA, 'numérica')]:
        ax[1].plot(tf[sel], 100*(x/x.mean() - 1), color=col, lw=0.9, label=lab)
        print(f'   partícula 4978, t in [1500,2000], J {lab}: pico a pico '
              f'{100*np.ptp(x)/x.mean():.3f}%')
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel(r'$J/\langle J\rangle-1$ (%)')
    ax[1].legend(loc='lower right', fontsize=6.5, ncol=2)
    ax[1].set_ylim(-0.32, 0.32)
    ax[1].set_title(r'(b) acción de una partícula, $\langle J\rangle=0.30$', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'mapa_numerico')


# ===========================================================================
# Sección 10: el mapa ángulo-acción numérico.
def fig_mapa_validacion():
    import h5py
    from aa_numerico import MapaAA
    print('[mapa_validacion] contra el isócrono y convergencia de la cuadratura')
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.6))

    # (a) numérico sin potencial propio frente al analítico, partícula por partícula.
    f = h5py.File(os.path.join(SG, 'long_fino_nosg', 'vlasov_output.h5'))
    pas = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    r = f[pas[250]]['r_part'][:]; p = f[pas[250]]['p_part'][:]
    Qa, Ja = rp_to_QJ(r, p)
    Qn, Jn, _ = MapaAA()(r, p)
    dQ = np.abs(np.angle(np.exp(1j*(Qn - Qa))))
    dJ = np.abs(Jn - Ja)
    dist = np.minimum.reduce([Qa, np.abs(Qa - np.pi), 2*np.pi - Qa])
    piso = 1e-17
    ax[0].loglog(np.maximum(dist, 1e-6), np.maximum(dQ, piso), '.', ms=1.2, color=AZUL,
                 alpha=0.5, rasterized=True, label=r'$|Q_{\rm num}-Q_{\rm an}|$')
    ax[0].loglog(np.maximum(dist, 1e-6), np.maximum(dJ, piso), '.', ms=1.2, color=NARANJA,
                 alpha=0.5, rasterized=True, label=r'$|J_{\rm num}-J_{\rm an}|$')
    ax[0].set_ylim(1e-17, 1e-4)
    ax[0].set_xlabel('distancia angular al punto de retorno (rad)')
    ax[0].set_ylabel('diferencia')
    ax[0].legend(loc='upper right', fontsize=6.5, markerscale=5)
    ax[0].set_title('(a) numérico frente a analítico', loc='left')
    print(f'   max|dJ| = {dJ.max():.1e}; |dQ|: mediana {np.median(dQ):.1e}, máx {dQ.max():.1e};'
          f' máx |dQ| con distancia > 0.01: {dQ[dist > 0.01].max():.1e}')

    # (b) convergencia del periodo radial con y sin la sustitución r = rm + ra sin(th).
    g = h5py.File(os.path.join(SG, 'long20k_snap', 'vlasov_output.h5'))
    pas = sorted([k for k in g if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    tg = np.array([g[k].attrs['time'] for k in pas])
    rg = g['grid']['r'][:]
    ps = np.mean([g[k]['potential'][:] - phi_iso(rg) for k, tt in zip(pas, tg) if tt >= 2000], axis=0)
    m = MapaAA(rg, ps)
    nodos = np.array([4, 8, 16, 32, 64, 128])
    for Jobj, col in [(0.05, AZUL), (0.30, NARANJA)]:
        E = E_de_J(Jobj)
        r1 = m._raiz(np.array([E]), 1e-2, m.rc)[0]
        r2 = m._raiz(np.array([E]), m.rc, 60.0)[0]
        rm, ra = 0.5*(r1 + r2), 0.5*(r2 - r1)

        def T_sust(n):
            x, w = np.polynomial.legendre.leggauss(n)
            th = 0.5*np.pi*x
            rr = rm + ra*np.sin(th)
            v = np.sqrt(np.maximum(2*(E - m.phi_ef(rr)), 0))
            return 2*np.sum(0.5*np.pi*w*ra*np.cos(th)/v)

        def T_ingenuo(n):
            x, w = np.polynomial.legendre.leggauss(n)
            rr = rm + ra*x
            v = np.sqrt(np.maximum(2*(E - m.phi_ef(rr)), 0))
            return 2*np.sum(ra*w/v)

        ref = T_sust(1024)
        es = [abs(T_sust(n)/ref - 1) for n in nodos]
        ei = [abs(T_ingenuo(n)/ref - 1) for n in nodos]
        ax[1].loglog(nodos, np.maximum(es, 1e-16), '-o', ms=3, color=col,
                     label=f'con sustitución, $J\\approx{Jobj}$')
        ax[1].loglog(nodos, ei, '--s', ms=3, color=col, mfc='white',
                     label=f'sin sustitución, $J\\approx{Jobj}$')
        print(f'   J~{Jobj}: error relativo de T con sustitución {np.array(es)}')
        print(f'   J~{Jobj}: sin sustitución {np.array(ei)}')
    ax[1].loglog(nodos, 0.5*nodos**-1.0, ':', color=TINTA, lw=0.9, label=r'$\propto n^{-1}$')
    ax[1].set_ylim(1e-16, 2)
    ax[1].set_xlabel('nodos de Gauss--Legendre')
    ax[1].set_ylabel(r'error relativo de $T_r$')
    ax[1].legend(loc='lower left', fontsize=5.8, ncol=1)
    ax[1].set_title('(b) cuadratura del periodo', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'mapa_validacion')


def fig_mapa_marco():
    import h5py
    print('[mapa_marco] el potencial propio en el tiempo y h_0 en tres marcos')
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.5))
    g = h5py.File(os.path.join(SG, 'long20k_snap', 'vlasov_output.h5'))
    pas = sorted([k for k in g if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    tg = np.array([g[k].attrs['time'] for k in pas])
    rg = g['grid']['r'][:]
    cols = [GRIS_CLARO, AQUA, NARANJA, AZUL, TINTA]
    for tt, col in zip([0, 400, 1000, 2000, 20000], cols):
        n = np.argmin(np.abs(tg - tt))
        ps = g[pas[n]]['potential'][:] - phi_iso(rg)
        ax[0].plot(rg, 1e4*ps, color=col, lw=1.1 if tt < 20000 else 0.9,
                   ls='-' if tt < 20000 else (0, (4, 2)), label=f'$t={tt}$')
    i5 = np.argmin(np.abs(rg - 5.0))
    p0 = g[pas[0]]['potential'][i5] - phi_iso(rg[i5])
    p2 = g[pas[np.argmin(np.abs(tg - 2000))]]['potential'][i5] - phi_iso(rg[i5])
    tarde = np.array([g[pas[n]]['potential'][i5] - phi_iso(rg[i5]) for n in np.where(tg >= 2000)[0]])
    print(f'   Phi_self(r=5): t=0 {p0:.4e}, t=2000 {p2:.4e} (cambio {p2/p0-1:+.1%});'
          f' fluctuación para t>=2000: {tarde.std()/abs(tarde.mean()):.1e}')
    ax[0].set_xlim(0, 20)
    ax[0].set_xlabel('$r$')
    ax[0].set_ylabel(r'$10^4\,\Phi_{\rm self}(r)$')
    ax[0].legend(loc='lower right', fontsize=6.5)
    ax[0].set_title('(a) potencial propio, $a_0=10^{-3}$', loc='left')

    d = np.load(os.path.join(SG, 'long20k_snap', 'aa_meseta.npz'))
    t = d['t']
    for k, col, lab in [('iso', AZUL, 'mapa del isócrono'),
                        ('promedio', NARANJA, 'potencial promediado'),
                        ('instante', AQUA, 'potencial instantáneo')]:
        h0 = d[f'h0_{k}']
        ax[1].plot(t, 100*(h0/h0[0] - 1), color=col, lw=1, label=lab)
        print(f'   {k:>9}: h0(20000)/h0(0)-1 = {h0[-1]/h0[0]-1:+.4e}')
    ax[1].set_xlim(0, 4000)
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel(r'$h_0(t)/h_0(0)-1$ (%)')
    ax[1].legend(loc='center right', fontsize=6.5)
    ax[1].set_title(r'(b) $h_0$ según el marco, $t\leq4000$', loc='left')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'mapa_marco')


# ===========================================================================
# Amortiguamiento de Landau en un equilibrio autoconsistente.
LAN = os.path.join(EXE, 'landau')
VENT_COLA = [(1000, 1600), (1100, 1800), (1200, 2000)]


def _pert(sufijo, nrc=400):
    d0 = np.load(os.path.join(LAN, f'L_a1e-2_n{nrc}_e0{sufijo}', 'landau.npz'))
    d1 = np.load(os.path.join(LAN, f'L_a1e-2_n{nrc}_e0.1{sufijo}', 'landau.npz'))
    return d1['t'], (d1['hk'][:, 1] - d0['hk'][:, 1])/0.1, d0


def fig_landau_respuesta():
    from landau_cola import ajustar
    print('[landau_respuesta] simulación, phase mixing libre y teoría lineal')
    t, h, d0 = _pert('', 800)
    lb = np.load(os.path.join(LAN, 'L_a1e-2_n800_e0.1', 'libre.npz'))
    lin = np.load(os.path.join(LAN, 'lineal', 'lin_1600x64.npz'))
    fig, ax = plt.subplots(1, 2, figsize=(ANCHO, 2.6),
                           gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax[0].semilogy(lb['t'], np.abs(lb['hk'][:, 1])/0.1, color=GRIS, lw=1,
                   ls=(0, (4, 2)), label='phase mixing libre')
    ax[0].semilogy(t, np.abs(d0['hk'][:, 1])/0.1, color=GRIS_CLARO, lw=0.8,
                   label=r'piso: corrida sin perturbación$/\varepsilon$')
    ax[0].semilogy(t, np.abs(h), color=AZUL, lw=1.2, label='simulación')
    ax[0].semilogy(lin['t'], np.abs(lin['h1']), color=NARANJA, lw=0.9,
                   ls=(0, (1.5, 1.5)), label='teoría lineal')
    ax[0].set_ylim(1e-6, 1)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel(r'$|h_1|/\varepsilon$')
    ax[0].legend(loc='upper right', fontsize=6.2)
    ax[0].set_title(r'(a) respuesta, $a_0=10^{-2}$', loc='left')
    r = ajustar(t, h, 1100, 1800)
    w, g = r['pencil M=2']
    v = (t >= 900) & (t <= 2000)
    ax[1].semilogy(t[v], np.abs(h[v]), color=AZUL, lw=1.2, label='simulación')
    m = (lin['t'] >= 900) & (lin['t'] <= 2000)
    ax[1].semilogy(lin['t'][m], np.abs(lin['h1'][m]), color=NARANJA, lw=0.9,
                   ls=(0, (1.5, 1.5)), label='teoría lineal')
    vv = (t >= 1100) & (t <= 1800)
    A = np.exp(np.mean(np.log(np.abs(h[vv])) + g*t[vv]))
    tt = np.linspace(1000, 2000, 50)
    ax[1].semilogy(tt, A*np.exp(-g*tt), color=TINTA, lw=0.8, ls=':',
                   label=rf'$\propto e^{{-\gamma t}}$, $\gamma={g*1e3:.2f}\times10^{{-3}}$')
    ax[1].set_xlabel('$t$')
    ax[1].legend(loc='upper right', fontsize=6.2)
    ax[1].set_title('(b) la cola colectiva', loc='left')
    print(f'   ajuste [1100,1800], pencil M=2: w={w:.5f} g={g:.4e}')
    fig.tight_layout(w_pad=1.2)
    guardar(fig, 'landau_respuesta')


def fig_landau_controles():
    from landau_cola import ajustar
    print('[landau_controles] equilibrio, linealidad y robustez de gamma')
    fig, ax = plt.subplots(1, 3, figsize=(ANCHO, 2.35),
                           gridspec_kw=dict(width_ratios=[1, 1, 1.15]))
    # (a) prueba de equilibrio: corrida sin perturbación
    d0 = np.load(os.path.join(LAN, 'L_a1e-2_n400_e0', 'landau.npz'))
    t0 = d0['t']; r = d0['r']; z = (r >= 3) & (r <= 15)
    rms = np.sqrt(np.mean(d0['dphi'][:, z]**2, axis=1))/d0['escala']
    ax[0].semilogy(t0[1:], np.abs(d0['hk'][1:, 1]), color=AZUL, lw=0.9, label=r'$|h_1|/h_0$')
    ax[0].semilogy(t0, rms, color=NARANJA, lw=0.9, label=r'rms $\delta\Phi/|\Phi_{\rm self}|$')
    ax[0].semilogy(t0[1:], np.abs(d0['hk'][1:, 0].real - 1), color=AQUA, lw=0.9, label=r'$|h_0-1|$')
    ax[0].set_ylim(1e-7, 1e-3)
    ax[0].set_xlabel('$t$')
    ax[0].legend(loc='lower right', fontsize=5.8)
    ax[0].set_title(r'(a) equilibrio, $\varepsilon=0$', loc='left')
    # (b) linealidad
    for e, col, ls in [('0.1', AZUL, '-'), ('0.05', NARANJA, (0, (3, 2)))]:
        d = np.load(os.path.join(EXE, 'landau', f'eq_a1e-2_eps{e}', 'landau.npz'))
        ax[1].semilogy(d['t'], np.abs(d['hk'][:, 1])/float(e), color=col, ls=ls, lw=1,
                       label=rf'$\varepsilon={e}$')
    ax[1].set_ylim(1e-5, 1)
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel(r'$|h_1|/\varepsilon$')
    ax[1].legend(loc='lower left', fontsize=6)
    ax[1].set_title('(b) linealidad', loc='left')
    # (c) gamma con cada variante
    casos = [('PIC $N_J$=400', lambda: _pert('', 400)[:2]),
             ('PIC $N_J$=800', lambda: _pert('', 800)[:2]),
             (r'PIC $\Delta t$=0.05', lambda: _pert('_c1', 400)[:2]),
             (r'PIC $\Delta t$=0.2', lambda: _pert('_c4', 400)[:2])]
    for nom in ['lin_400x25', 'lin_800x25', 'lin_1600x32', 'lin_1600x32_dt025', 'lin_1600x64']:
        casos.append((nom.replace('lin_', 'lineal ').replace('_dt025', r' $\Delta t$=0.25'),
                      (lambda n=nom: (lambda d: (d['t'], d['h1']))(np.load(os.path.join(LAN, 'lineal', n + '.npz'))))))
    for i, (nom, carga) in enumerate(casos):
        t, h = carga()
        gs = [ajustar(t, h, lo, hi)['pencil M=2'][1]*1e3 for lo, hi in VENT_COLA]
        col = AZUL if nom.startswith('PIC') else NARANJA
        ax[2].errorbar(np.mean(gs), i, xerr=[[np.mean(gs) - min(gs)], [max(gs) - np.mean(gs)]],
                       fmt='o', ms=3, color=col, capsize=2, lw=0.9)
        print(f'   {nom:>28}: gamma*1e3 en las tres ventanas = ' + ', '.join(f'{x:.3f}' for x in gs))
    ax[2].set_yticks(range(len(casos)))
    ax[2].set_yticklabels([c[0] for c in casos], fontsize=6)
    ax[2].invert_yaxis()
    ax[2].set_xlabel(r'$\gamma\times10^{3}$')
    ax[2].set_xlim(4.95, 5.2)
    ax[2].set_title(r'(c) $\gamma$', loc='left')
    fig.tight_layout(w_pad=0.6)
    guardar(fig, 'landau_controles')


TODAS = [fig_liouville, fig_potencial, fig_frecuencias, fig_orbita, fig_enrollamiento,
         fig_hk_exacto, fig_convergencia, fig_agnosticas, fig_salud,
         fig_envolvente, fig_barridos, fig_masa, fig_fase, fig_accion,
         fig_capturas, fig_mapa_numerico, fig_mapa_validacion, fig_mapa_marco,
         fig_landau_respuesta, fig_landau_controles]

if __name__ == '__main__':
    filtro = sys.argv[1:]
    for fn in TODAS:
        if not filtro or any(s in fn.__name__ for s in filtro):
            fn()
