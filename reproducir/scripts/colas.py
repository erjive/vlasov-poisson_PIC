"""Leyes de decaimiento de la mezcla y del campo del gas, por forma de la FD y por masa.

Dos observables de las corridas ya hechas: |h_n(t)|/h_0, que mide la mezcla (variables de ángulo-acción del
potencial de cada instante, peso B = 1), y D(t) = ||Phi_gas(t) - <Phi_gas>||/max |<Phi_gas>| en 3 <= r <= 15, que mide
cuánto falta para que el campo del gas sea estacionario (bulto.campo_particulas). De cada uno se toma la envolvente
(máximo en una ventana de +-50) y su pendiente logarítmica en ventanas de t/tau_1.

Lo que se espera. Sin autogravedad la cola de h_n la fija la FD en los extremos de su soporte en J (Sección 3 del
artículo): F_n ~ (J - J_e)^k da t^-(k+1). Para J^p exp(-J^2/s_J^2) es t^-(p+1); para (J - a)_+^2 exp(...), t^-3; para
la meseta, t^-1; para el anillo gaussiano, más rápida que cualquier potencia. `libre` lo comprueba con la integral
exacta. Una FD suave y genérica en la órbita circular (suave en las coordenadas cartesianas del plano (r, p_r), sin
anularse allí) tiene F_n ~ J^(n/2): sus colas son t^-(1 + n/2), es decir t^-1.5, t^-2 y t^-2.5 para n = 1, 2, 3, y el
potencial, que pesa el primer armónico con la amplitud radial sqrt(J), decae como t^-2. La nube del artículo sin
autogravedad daría t^-3.5 para el potencial.

    python3 colas.py libre      colas exactas sin autogravedad de las seis formas
    python3 colas.py tabla      pendientes y niveles de las corridas en las que el gas se mezcla
    python3 colas.py figura     figura de trabajo (exe/demo_eta/colas/colas.png)
    python3 colas.py alta       las corridas de alta resolución (tanda `colas` de VlasovPoisson_PIC_sp): pendientes,
                                perfiles junto a la órbita circular, borde interior y mezcla sola frente a la corrida
    python3 colas.py figura_alta   figura de trabajo de esas corridas (exe/demo_eta/colas/colas_alta.png)
    python3 colas.py articulo   tabla (docs/articulo/tab_colas.tex), figura (docs/articulo/figuras/colas_alta.pdf) y
                                cifras de la sección del artículo sobre la rapidez de la mezcla (Sec:Rates)

No corre simulaciones: lee exe/demo_eta/bulto (este repositorio) y exe/rep/12_nube de VlasovPoisson_PIC_sp.
"""
import os, sys
import numpy as np
from scipy.special import ive
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import bulto as B
import reducido as R
import formas as F
from demo_eta import envolvente

SAL = os.path.join(R.EXE, 'demo_eta', 'colas')
SP = os.path.join(os.environ.get('VP_PIC_SP', os.path.join(os.path.dirname(R.RAIZ), 'VlasovPoisson_PIC_sp')), 'exe', 'rep', '12_nube')
SJ, SQ, L0 = 0.1, 0.1, 2.0
# Exponente k de la FD en el extremo que domina la cola sin autogravedad (None: más rápida que cualquier potencia).
K_LIBRE = {'p2': 2, 'p4': 4, 'h05': 2, 'anillo': None, 'meseta': 0, 'p0': 0}
# Corridas del PIC en las que el gas se mezcla (m_c compatible con m_mix), por forma: (corrida, origen).
MEZCLA_FORMAS = {'p2': ('B2q', 'fijo'), 'p4': ('F_p4_020', 'sp'), 'h05': ('F_h05_020', 'sp'),
                 'anillo': ('G_anillo_150', 'sp'), 'meseta': ('G_meseta_150', 'sp')}
# Familias de la Sección 6 (factor J^2 exp(-J^2/s_J^2)): corridas de L fija con mu <= 0.42, sin núcleo o casi.
FAMILIAS = [('reference', ['B1', 'B2q', 'B2j', 'B015']), (r'$\sigma_J=0.05$', ['J05_0.005']),
            (r'$\sigma_J=0.2$', ['J20_0.02', 'J20_0.04']), ('Kepler', ['K_0.01', 'K_0.02']),
            (r'$\sigma_Q=0.2$', ['Q2_0.01', 'Q2_0.02']), (r'$\sigma_Q=0.4$', ['Q4_0.01', 'Q4_0.02'])]
OTRAS_FORMAS = ['F_p4_020', 'F_p4_030', 'F_h05_020', 'F_h05_030', 'G_anillo_150', 'G_meseta_150']


def pendiente(x, y, a, b):
    """Pendiente de log y frente a log x en a <= x <= b."""
    z = (x >= a) & (x <= b) & (y > 0)
    return np.polyfit(np.log(x[z]), np.log(y[z]), 1)[0]


def h_libre(forma, n, ts, nj=200001):
    """|h_n(t)|/h_0 exacto sin autogravedad para la forma dada (isócrono, s_Q = s_J = 0.1, B = 1):
    (I_n(z)/I_0(z)) |int f(J) exp(-i n Omega(J) t) dJ| / int f(J) dJ, z = 1/(2 s_Q^2), por la regla del trapecio."""
    J = np.linspace(0.0, 8*SJ, nj)
    f = R.perfil_J(J, SJ, F.TODAS[forma][0]); f[0] *= 0.5; f[-1] *= 0.5
    om = (J + R.c_de_L(L0))**-3
    z = 0.5/SQ**2
    fac = ive(n, z)/ive(0, z)
    return fac*np.array([abs((f*np.exp(-1j*n*om*t)).sum()) for t in ts])/f.sum()


def curvas(n, origen):
    """t/tau_1, envolventes de |h_k| (k = 1, 2, 3; solo k = 1 para las corridas de VlasovPoisson_PIC_sp, cuyos
    armónicos superiores están en el cluster) y de D, y tau_1 de la corrida."""
    if origen == 'fijo':
        d = B.cargar(n); t = d['t']; tau = B.banda(n)[3]
        hk = np.abs(np.asarray(d['hk']))
        h = [envolvente(hk[:, k], t) for k in (1, 2, 3)]
        tD, D = B.campo_particulas(n)
    else:
        M = np.load(os.path.join(SP, n, 'medidas.npz')); t = M['t']; tau = float(M['tau'])
        h = [envolvente(np.abs(M['h1']), t)]
        c = np.load(os.path.join(SP, n, 'campo.npz')); tD, D = c['t'], c['D']
    return t/tau, h, tD/tau, envolvente(D, tD), tau


def gamma_fijo(n):
    """Gamma de una corrida de L fija de bulto.py (factor J^2 exp(-J^2/s_J^2))."""
    sj = B.parametro(n, 'sr', 0.1); pot = 'isocrono' if B.FONDO.get(n, 'isocrono') == 'isocrono' else 'kepler'
    Fo = R.fondo(pot, L0)
    return Fo['g']*B.info[n]['a0']/(3.0*Fo['c']**-4*R.J_promedio(sj)**1.5)


def libre():
    print('Colas exactas sin autogravedad (isócrono, s_Q = s_J = 0.1): pendiente de la envolvente de |h_n| en 8-25 tau_1')
    print(f'{"forma":8} {"esperada":>9}   n = 1   n = 2   n = 3    |h_1| en 10 y 20 tau_1')
    for forma in ('p2', 'p4', 'h05', 'anillo', 'meseta', 'p0'):
        tau = R.banda_perfil(F.TODAS[forma][0], SJ)[3]
        ts = np.arange(0.0, 26*tau, 10.0)
        e = [envolvente(h_libre(forma, n, ts), ts) for n in (1, 2, 3)]
        k = K_LIBRE[forma]
        print(f'{forma:8} {"-" if k is None else f"{-(k + 1):d}":>9}  ' + ' '.join(f'{pendiente(ts/tau, y, 8, 25):6.2f}' for y in e)
              + f'    {e[0][np.argmin(abs(ts - 10*tau))]:.1e}  {e[0][np.argmin(abs(ts - 20*tau))]:.1e}')


def tabla():
    print('Corridas en las que el gas se mezcla. Pendientes de las envolventes en 3-10 tau_1 (y en 10-20 tau_1 para |h_1|);')
    print('niveles de |h_1| en 10 y 20 tau_1 y de D en 10 tau_1. Con FD genérica en la órbita circular: -1.5, -2, -2.5 y D: -2.')
    cab = f'{"familia":14} {"corrida":13} {"M_gas":>7} {"Gamma":>6} {"mu":>5} | {"n=1":>6} {"n=2":>6} {"n=3":>6} {"D":>6} | {"n=1 tarde":>9} | {"|h_1|(10)":>9} {"|h_1|(20)":>9} {"D(10)":>8}'
    print(cab)
    filas = [('free', 'B0')] + [(f, n) for f, L in FAMILIAS for n in L]
    for fam, n in filas:
        x, h, xD, D, tau = curvas(n, 'fijo')
        mu = B.escalas(n)[2] if B.info[n]['sg'] else 0.0
        gam = gamma_fijo(n) if B.info[n]['sg'] else 0.0
        p = [pendiente(x, y, 3, 10) for y in h]
        print(f'{fam:14} {n:13} {B.info[n]["a0"] if B.info[n]["sg"] else 0:7g} {gam:6.3f} {mu:5.2f} | {p[0]:6.2f} {p[1]:6.2f} {p[2]:6.2f} {pendiente(xD, D, 3, 10):6.2f} | '
              f'{pendiente(x, h[0], 10, 20):9.2f} | {h[0][np.argmin(abs(x - 10))]:9.1e} {h[0][np.argmin(abs(x - 20))]:9.1e} {D[np.argmin(abs(xD - 10))]:8.1e}')
    import importlib
    for n in OTRAS_FORMAS:
        if not os.path.exists(os.path.join(SP, n, 'medidas.npz')):
            continue
        x, h, xD, D, tau = curvas(n, 'sp')
        r = np.load(os.path.join(SP, n, 'resumen.npz'))
        forma = n.split('_')[1]
        m = float(open(os.path.join(SP, n, 'params_usados.par')).read().split('a0')[1].split('=')[1].split()[0]) if os.path.exists(os.path.join(SP, n, 'params_usados.par')) else np.nan
        gam = R.G*m/(abs(F.OMEGA_P)*R.J_promedio(SJ, F.TODAS[forma][0])**1.5)
        print(f'{forma:14} {n:13} {m:7g} {gam:6.3f} {float(r["mu"]):5.2f} | {pendiente(x, h[0], 3, 10):6.2f} {"-":>6} {"-":>6} {pendiente(xD, D, 3, 10):6.2f} | '
              f'{pendiente(x, h[0], 10, 20):9.2f} | {h[0][np.argmin(abs(x - 10))]:9.1e} {h[0][np.argmin(abs(x - 20))]:9.1e} {D[np.argmin(abs(xD - 10))]:8.1e}')


def figura():
    """Figura de trabajo, tres paneles: (a) |h_1| de las seis formas sin autogravedad (exacta, líneas finas) y con
    ella en la corrida que se mezcla (gruesas); (b) pendiente de |h_1| en 3-10 tau_1 frente a mu para las familias de
    la Sección 6; (c) D(t) de las corridas con 0.3 <= mu <= 0.42 de todas las familias y formas."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 4.0))
    # (a) formas
    a = ax[0]
    for forma, (n, origen) in MEZCLA_FORMAS.items():
        pf, pot, nj, col, mk, rot = F.TODAS[forma]
        tau = R.banda_perfil(pf, SJ)[3]
        ts = np.arange(0.0, 30*tau, 10.0)
        a.plot(ts[1:]/tau, envolvente(h_libre(forma, 1, ts), ts)[1:], color=col, lw=0.9, ls=(0, (4, 2.5)))
        x, h, xD, D, _ = curvas(n, origen)
        a.plot(x[1:], h[0][1:], color=col, lw=1.6, label=rot.replace(' (sections 6-8)', ''))
    xx = np.array([4.0, 28.0])
    a.plot(xx, 0.12*xx**-1.5, color=tinta2, lw=0.8, ls=':'); a.text(13, 0.12*13**-1.5*1.5, r'$t^{-3/2}$', fontsize=8, color=tinta2)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1, 30); a.set_ylim(1e-7, 1.5)
    a.set_title(r'$|h_1|/h_0$: without (dashed) and with self-gravity', loc='left', fontsize=9.5, pad=6)
    a.set_xlabel(r'$t/\tau_1$'); a.legend(frameon=False, fontsize=7.5, loc='lower left')
    # (b) pendiente frente a mu
    a = ax[1]
    marcas = ['o', 'v', '^', '*', 's', 'D']; cols = ['#0b0b0b', '#2a78d6', '#1baf7a', '#6250d6', '#eb6834', '#e87ba4']
    x0, h0, _, _, _ = curvas('B0', 'fijo')
    a.plot([0.0], [pendiente(x0, h0[0], 3, 10)], marker='o', ms=6, color=cols[0], mfc='none', ls='none')
    for (fam, L), mk, col in zip(FAMILIAS, marcas, cols):
        mus, ps = [], []
        for n in L:
            x, h, _, _, _ = curvas(n, 'fijo'); mus.append(B.escalas(n)[2]); ps.append(pendiente(x, h[0], 3, 10))
        a.plot(mus, ps, marker=mk, ms=6, color=col, ls='none', mec=sup, mew=0.6, label=fam)
    for y, txt in ((-3.0, r'$-3$ (free)'), (-1.5, r'$-3/2$')):
        a.axhline(y, color=tinta2, lw=0.8, ls=':'); a.text(0.13, y + 0.06, txt, fontsize=8, color=tinta2)
    a.set_xlim(-0.02, 0.46); a.set_ylim(-3.7, -0.8)
    a.set_title(r'slope of $|h_1|$ over $3\leq t/\tau_1\leq 10$', loc='left', fontsize=9.5, pad=6)
    a.set_xlabel(r'$\mu$'); a.legend(frameon=False, fontsize=7.5, loc='lower right')
    # (c) campo
    a = ax[2]
    x, h, xD, D, _ = curvas('B0', 'fijo'); a.plot(xD[1:], D[1:], color=tinta2, lw=0.9, ls=(0, (4, 2.5)), label='reference cloud, free')
    for (fam, L), col in zip(FAMILIAS, cols):
        for n in L:
            mu = B.escalas(n)[2]
            if 0.3 <= mu <= 0.42 and n != 'B2j':
                x, h, xD, D, _ = curvas(n, 'fijo'); a.plot(xD[1:], D[1:], color=col, lw=1.2)
    for forma, (n, origen) in MEZCLA_FORMAS.items():
        if origen == 'sp':
            x, h, xD, D, _ = curvas(n, origen); a.plot(xD[1:], D[1:], color=F.TODAS[forma][3], lw=1.2, ls=(0, (1, 1)))
    a.plot(xx, 0.02*xx**-2.0, color=tinta2, lw=0.8, ls=':'); a.text(13, 0.02*13**-2.0*1.6, r'$t^{-2}$', fontsize=8, color=tinta2)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1, 30); a.set_ylim(2e-5, 0.2)
    a.set_title(r'deviation $D$ of the potential of the gas', loc='left', fontsize=9.5, pad=6)
    a.set_xlabel(r'$t/\tau_1$')
    a.legend(handles=[Line2D([], [], color=tinta2, lw=0.9, ls=(0, (4, 2.5)), label='reference cloud, free'),
                      Line2D([], [], color=tinta2, lw=1.2, label=r'families of section 6, $0.3\leq\mu\leq 0.42$'),
                      Line2D([], [], color=tinta2, lw=1.2, ls=(0, (1, 1)), label='other shapes (colors of the left panel)')],
             frameon=False, fontsize=7.5, loc='lower left')
    for a in ax:
        a.grid(True, color=reja, lw=0.5); a.set_axisbelow(True)
        for lado in ('top', 'right'):
            a.spines[lado].set_visible(False)
    fig.tight_layout(w_pad=1.5)
    os.makedirs(SAL, exist_ok=True)
    fig.savefig(os.path.join(SAL, 'colas.png'), dpi=120); plt.close(fig)
    print('escrita', os.path.join(SAL, 'colas.png'))


# Corridas de alta resolución (3200 x 100; nube_L.py, tanda `colas`): forma -> [(Gamma, corrida), ...].
ALTA = {'p2': [(0.025, 'H_p2_025'), (0.05, 'H_p2_050'), (0.10, 'H_p2_100'), (0.15, 'H_p2_150')],
        'p4': [(0.05, 'H_p4_050'), (0.10, 'H_p4_100')], 'h05': [(0.05, 'H_h05_050'), (0.10, 'H_h05_100')],
        'anillo': [(0.05, 'H_anillo_050'), (0.10, 'H_anillo_100')], 'meseta': [(0.05, 'H_meseta_050'), (0.10, 'H_meseta_100')]}
# El número del nombre es 1000 Gamma con la mediana de J, con la que se fijaron las masas. El artículo usa el promedio <J>.
ALTA = {f: [(g*(F.J_mediana(F.TODAS[f][0])/R.J_promedio(SJ, F.TODAS[f][0]))**1.5, n) for g, n in L] for f, L in ALTA.items()}


def carga_alta(n):
    """t/tau_1, envolventes de |h_k| (k = 1, 2, 3) y de D de una corrida de alta resolución."""
    A = np.load(os.path.join(SP, n, 'armonicos.npz')); tau = float(np.load(os.path.join(SP, n, 'medidas.npz'))['tau'])
    t = A['t']; hk = np.abs(A['hk'])
    c = np.load(os.path.join(SP, n, 'campo.npz'))
    return t/tau, [envolvente(hk[:, k], t) for k in (1, 2, 3)], c['t']/tau, envolvente(c['D'], c['t']), tau


def exponente_perfil(n, k=1, i=2, a=0.01, b=0.10):
    """Exponente q de a_k(J) ~ J^q junto a la órbita circular (perfiles.npz, tiempo de índice i: 2 es 10 tau_1), por
    ajuste en a <= J/s_J <= b; nan si no hay gas ahí."""
    P = np.load(os.path.join(SP, n, 'perfiles.npz')); x = 0.5*(P['x'][1:] + P['x'][:-1])
    y = P['f0'][i] if k == 0 else P['a'][i, k - 1]
    z = (x >= a) & (x <= b) & (y > 0)
    return np.polyfit(np.log(x[z]), np.log(y[z]), 1)[0] if z.sum() > 3 else np.nan


def borde_interior(n, i=2, nivel=0.01):
    """J/s_J de la primera franja en que la masa por unidad de J pasa de `nivel` de su máximo (perfiles_ancho.npz)."""
    P = np.load(os.path.join(SP, n, 'perfiles_ancho.npz')); x = 0.5*(P['x'][1:] + P['x'][:-1]); f = P['f0'][i]
    return x[np.nonzero(f >= nivel*f.max())[0][0]]


def cociente_cinematico(n, i0=1, k=0, ventana=(15.0, 25.0)):
    """Media geométrica, en la ventana de t/tau_1, del cociente entre la envolvente de |h_{k+1}| de la corrida y la de
    su continuación cinemática desde t0 (índice i0: 0, 1, 2 son 3, 5 y 10 tau_1). 1: después de t0 solo hay mezcla."""
    c = np.load(os.path.join(SP, n, 'cinematica.npz')); t = c['t']; x = t/float(c['tau'])
    real = envolvente(np.abs(c['hk'][:, k]), t)
    cin = np.abs(c['hk_cin'][i0, k]); cin = envolvente(np.where(np.isfinite(cin), cin, 0.0), t)
    z = (x >= ventana[0]) & (x <= ventana[1])
    return np.exp(np.mean(np.log(real[z]/cin[z])))


def alta():
    print('Corridas de alta resolución (3200 x 100, L igual). Pendientes de las envolventes en ventanas de t/tau_1;')
    print('q_1: exponente de a_1(J) en 0.01-0.10 s_J en t = 10 tau_1 (la cola cinemática sería t^-(q_1+1)); borde: J/s_J donde la FD')
    print('pasa del 1 % de su máximo en t = 10 tau_1; real/cin: |h_1| de la corrida entre el de la mezcla sola desde 5 y 10 tau_1 (15-25 tau_1).')
    print(f'{"forma":7} {"Gamma":>5} | {"|h_1| 3-10":>10} {"10-20":>6} {"20-30":>6} | {"|h_2| 10-20":>11} {"|h_3|":>6} | {"D 3-10":>7} {"10-20":>6} | '
          f'{"q_1":>5} {"-(q_1+1)":>8} | {"borde":>5} | {"real/cin 5":>10} {"10":>5} | {"|h_1|(20)":>9}')
    for forma, L in ALTA.items():
        tau = R.banda_perfil(F.TODAS[forma][0], SJ)[3]
        ts = np.arange(0.0, 31*tau, 10.0); e = envolvente(h_libre(forma, 1, ts), ts); xl = ts/tau
        k = K_LIBRE[forma]
        print(f'{forma:7} {"libre":>5} | {pendiente(xl, e, 3, 10):10.2f} {pendiente(xl, e, 10, 20):6.2f} {pendiente(xl, e, 20, 30):6.2f} | '
              f'{"":11} {"":6} | {"":7} {"":6} | {"-" if k is None else k:>5} {"-" if k is None else -(k + 1):>8} | {"":5} | {"":10} {"":5} | '
              f'{e[np.argmin(abs(xl - 20))]:9.1e}')
        for gam, n in L:
            x, h, xD, D, _ = carga_alta(n)
            q = exponente_perfil(n)
            print(f'{"":7} {gam:5.3f} | {pendiente(x, h[0], 3, 10):10.2f} {pendiente(x, h[0], 10, 20):6.2f} {pendiente(x, h[0], 20, 30):6.2f} | '
                  f'{pendiente(x, h[1], 10, 20):11.2f} {pendiente(x, h[2], 10, 20):6.2f} | {pendiente(xD, D, 3, 10):7.2f} {pendiente(xD, D, 10, 20):6.2f} | '
                  f'{q:5.2f} {-(q + 1):8.2f} | {borde_interior(n):5.2f} | {cociente_cinematico(n, 1):10.2f} {cociente_cinematico(n, 2):5.2f} | '
                  f'{h[0][np.argmin(abs(x - 20))]:9.1e}')
    x, h, xD, D, _ = carga_alta('H_p2_100'); x2, h2, xD2, D2, _ = carga_alta('H2_p2_100')
    z = (x >= 3) & (x <= 30)
    print('convergencia (3200 x 100 frente a 6400 x 200, Gamma = 0.10): mayor diferencia relativa de las envolventes en 3-30 tau_1: '
          + ', '.join(f'|h_{k+1}| {np.abs(h[k][z]/h2[k][z] - 1).max():.2f}' for k in range(3)) + f', D {np.abs(D[(xD >= 3) & (xD <= 30)]/D2[(xD2 >= 3) & (xD2 <= 30)] - 1).max():.2f}'
          + '; hasta 20 tau_1: ' + ', '.join(f'{np.abs(h[k][(x >= 3) & (x <= 20)]/h2[k][(x >= 3) & (x <= 20)] - 1).max():.2f}' for k in range(3)))


def figura_alta():
    """Figura de trabajo de las corridas de alta resolución, seis paneles: (a) |h_1| de la nube del artículo al
    crecer Gamma; (b) su potencial; (c) su FD junto a la órbita circular en 10 tau_1; (d) |h_1| de las otras formas
    con Gamma = 0.05 y sin autogravedad; (e) borde interior de la FD de dos formas huecas; (f) la corrida frente a la
    mezcla sola desde 5 tau_1."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    fig, ax = plt.subplots(2, 3, figsize=(12.5, 7.6))
    azules = plt.cm.Blues(np.linspace(0.4, 0.95, 4))
    tau = B.TAU1; ts = np.arange(0.0, 30.6*tau, 10.0)
    # (a) y (b): la nube del artículo
    a = ax[0, 0]
    a.plot(ts[1:]/tau, envolvente(h_libre('p2', 1, ts), ts)[1:], color=tinta2, lw=1.0, ls=(0, (4, 2.5)), label=r'$\Gamma=0$ (exact)')
    for (gam, n), col in zip(ALTA['p2'], azules):
        x, h, xD, D, _ = carga_alta(n)
        a.plot(x[1:], h[0][1:], color=col, lw=1.5, label=rf'$\Gamma={gam:.3f}$')
        ax[0, 1].plot(xD[1:], D[1:], color=col, lw=1.5)
    xx = np.array([6.0, 28.0])
    a.plot(xx, 3.0*xx**-3.0, color=tinta2, lw=0.8, ls=':'); a.text(20, 3.0*20**-3.0*0.35, r'$t^{-3}$', fontsize=8, color=tinta2)
    a.plot(xx, 0.25*xx**-1.2, color=tinta2, lw=0.8, ls=':'); a.text(20, 0.25*20**-1.2*1.5, r'$t^{-1.2}$', fontsize=8, color=tinta2)
    a.set_title(r'(a) $|h_1|/h_0$, cloud of the article', loc='left', fontsize=9.5, pad=6); a.legend(frameon=False, fontsize=7.5, loc='lower left')
    a.set_ylim(5e-6, 1.5)
    xB, hB, xDB, DB, _ = curvas('B0', 'fijo')
    a = ax[0, 1]
    a.plot(xx, 0.6*xx**-3.5, color=tinta2, lw=0.8, ls=':'); a.text(7, 0.6*7**-3.5*0.3, r'$t^{-3.5}$', fontsize=8, color=tinta2)
    a.plot(xx, 0.02*xx**-1.7, color=tinta2, lw=0.8, ls=':'); a.text(20, 0.02*20**-1.7*1.6, r'$t^{-1.7}$', fontsize=8, color=tinta2)
    a.set_title(r'(b) deviation $D$ of the potential of the gas', loc='left', fontsize=9.5, pad=6); a.set_ylim(1e-6, 0.2)
    # (c) perfiles junto a la órbita circular
    a = ax[0, 2]
    xp = np.linspace(0.005, 0.395, 40)
    a.plot(xp, xp**2*np.exp(-xp**2)/(SJ*np.sqrt(np.pi)/4), color=tinta2, lw=1.0, ls=(0, (4, 2.5)), label='initial datum')
    for (gam, n), col in zip(ALTA['p2'], azules):
        P = np.load(os.path.join(SP, n, 'perfiles.npz')); x = 0.5*(P['x'][1:] + P['x'][:-1])
        a.plot(x, P['f0'][2], color=col, lw=1.5, label=rf'$\Gamma={gam:.3f}$: $q_1={exponente_perfil(n):.2f}$')
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel(r'$J/\sigma_J$')
    a.set_title(r'(c) mass per unit $J$ at the circular orbit, $t=10\tau_1$', loc='left', fontsize=9.5, pad=6); a.legend(frameon=False, fontsize=7.5, loc='lower right')
    # (d) otras formas
    a = ax[1, 0]
    for forma in ('p4', 'h05', 'anillo', 'meseta'):
        pf, pot, nj, col, mk, rot = F.TODAS[forma]
        tauf = R.banda_perfil(pf, SJ)[3]; tsf = np.arange(0.0, 30.6*tauf, 10.0)
        a.plot(tsf[1:]/tauf, envolvente(h_libre(forma, 1, tsf), tsf)[1:], color=col, lw=0.9, ls=(0, (4, 2.5)))
        x, h, xD, D, _ = carga_alta(ALTA[forma][0][1])
        a.plot(x[1:], h[0][1:], color=col, lw=1.5, label=rot)
    a.set_ylim(5e-7, 1.5)
    a.set_title(r'(d) $|h_1|/h_0$, other shapes: free (dashed) and $\Gamma=0.05$', loc='left', fontsize=9.5, pad=6); a.legend(frameon=False, fontsize=7.5, loc='lower left')
    # (e) borde interior
    a = ax[1, 1]
    for forma, ls in (('h05', '-'), ('meseta', (0, (1, 1)))):
        pf, pot, nj, col, mk, rot = F.TODAS[forma]
        Jg = np.linspace(0.0, 3*SJ, 3001); f = R.perfil_J(Jg, SJ, pf); f = f/(f.sum()*(Jg[1] - Jg[0]))
        a.plot(Jg/SJ, f, color=tinta2, lw=0.9, ls=ls)
        for (gam, n), al in zip(ALTA[forma], (0.55, 1.0)):
            P = np.load(os.path.join(SP, n, 'perfiles_ancho.npz')); x = 0.5*(P['x'][1:] + P['x'][:-1])
            a.plot(x, P['f0'][2], color=col, lw=1.4, ls=ls, alpha=al)
    a.set_xlim(0, 1.2); a.set_xlabel(r'$J/\sigma_J$')
    a.set_title(r'(e) inner edge at $t=10\tau_1$: initial (gray), $\Gamma=0.05$, $0.10$', loc='left', fontsize=9.5, pad=6)
    a.legend(handles=[Line2D([], [], color=F.TODAS['h05'][3], lw=1.4, label=F.TODAS['h05'][5]),
                      Line2D([], [], color=F.TODAS['meseta'][3], lw=1.4, ls=(0, (1, 1)), label=F.TODAS['meseta'][5])], frameon=False, fontsize=7.5, loc='lower right')
    # (f) corrida frente a mezcla sola
    a = ax[1, 2]
    for n, col, rot in (('H_p2_050', azules[1], r'article, $\Gamma=0.05$'), ('H_p2_150', azules[3], r'article, $\Gamma=0.15$'),
                        ('H_h05_050', F.TODAS['h05'][3], r'hollow, $\Gamma=0.05$')):
        c = np.load(os.path.join(SP, n, 'cinematica.npz')); t = c['t']; x = t/float(c['tau'])
        a.plot(x[1:], envolvente(np.abs(c['hk'][:, 0]), t)[1:], color=col, lw=1.5, label=rot)
        cin = np.abs(c['hk_cin'][1, 0]); z = np.isfinite(cin)
        a.plot(x[z], envolvente(cin[z], t[z]), color=col, lw=1.0, ls=(0, (1, 1)))
    a.set_ylim(1e-5, 0.3)
    a.set_title(r'(f) run (solid) and phase mixing alone from $5\tau_1$ (dotted)', loc='left', fontsize=9.5, pad=6); a.legend(frameon=False, fontsize=7.5, loc='lower left')
    for a in (ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 2]):
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1, 31); a.set_xlabel(r'$t/\tau_1$')
    for a in ax.ravel():
        a.grid(True, color=reja, lw=0.5); a.set_axisbelow(True)
        for lado in ('top', 'right'):
            a.spines[lado].set_visible(False)
    fig.tight_layout(w_pad=1.5, h_pad=1.5)
    os.makedirs(SAL, exist_ok=True)
    fig.savefig(os.path.join(SAL, 'colas_alta.png'), dpi=120); plt.close(fig)
    print('escrita', os.path.join(SAL, 'colas_alta.png'))


def masa_alta(n):
    """M_gas/M de una corrida de alta resolución (a0 de params_usados.par)."""
    for linea in open(os.path.join(SP, n, 'params_usados.par')):
        c = linea.split()
        if len(c) >= 3 and c[0] == 'a0' and c[1] == '=':
            return float(c[2])
    raise ValueError(n)


def borde_inicial(forma, nivel=0.01):
    """J/s_J del borde interior del dato inicial, con las franjas de perfiles_ancho (0.01 s_J) y el mismo nivel."""
    x = np.arange(300)*0.01 + 0.005
    f = R.perfil_J(x*SJ, SJ, F.TODAS[forma][0])
    return x[np.nonzero(f >= nivel*f.max())[0][0]]


def num(v, cifras=2):
    """Número para las tablas del artículo: decimales fijos, o mantisa y potencia de 10 si es pequeño."""
    if not np.isfinite(v):
        return '--'
    if abs(v) < 1e-2:
        e = int(np.floor(np.log10(abs(v)))); m = v/10**e
        if round(m, 1) >= 10:
            m, e = m/10, e + 1
        return f'${m:.1f}\\times 10^{{{e}}}$'
    return f'${v:.{cifras}f}$'


ROTULO = {'p2': r'$J^2e^{-J^2/\sigma_J^2}$', 'p4': r'$J^4e^{-J^2/\sigma_J^2}$', 'h05': r'$(J-0.05)_+^2e^{-J^2/\sigma_J^2}$',
          'anillo': 'ring', 'meseta': 'water bag'}


def articulo():
    """Tabla tab_colas.tex y figura figuras/colas_alta.pdf de la sección sobre la rapidez de la mezcla, y las cifras
    que cita el texto. Seis paneles: (a) |h_1| de la nube de referencia al crecer Gamma; (b) su potencial; (c) primer
    armónico de su FD junto a la órbita circular en 10 tau_1; (d) |h_1| de las otras formas, libre y con Gamma = 0.05 y
    0.10; (e) borde interior de dos formas huecas en 10 tau_1; (f) pendiente de |h_1| en 10-20 tau_1 frente a Gamma."""
    art = os.path.join(R.RAIZ, 'docs', 'articulo')
    filas, pend = [], {}
    print('forma   Gamma  M_gas/M |  |h_1| 3-10 10-20 20-30 |  D 3-10 10-20 | |h_1|(20) |  q_1   q_0 | J_in/s_J')
    for forma, L in ALTA.items():
        tau = R.banda_perfil(F.TODAS[forma][0], SJ)[3]
        ts = np.arange(0.0, 31*tau, 10.0); e = envolvente(h_libre(forma, 1, ts), ts); xl = ts/tau
        k = K_LIBRE[forma]
        p1 = [pendiente(xl, e, 3, 10), pendiente(xl, e, 10, 20)] if k is not None else [np.nan, np.nan]
        h20 = e[np.argmin(abs(xl - 20))]
        jin = borde_inicial(forma) if forma in ('h05', 'anillo', 'meseta') else np.nan
        pend[forma] = [(0.0, p1[1])] if k is not None else []
        print(f'{forma:7} {"0":>5}  {"":7} | {p1[0]:10.2f} {p1[1]:5.2f} {"":5} | {"":7} {"":5} | {h20:9.1e} | {k if k is not None else "-":>4} {"":5} | {jin:.2f}')
        celdas = [ROTULO[forma], '$0$', '--', num(p1[0]), num(p1[1]), '--', '--', num(h20),
                  f'${k}$' if forma in ('p2', 'p4') else '--', num(jin)]
        filas.append(' & '.join(celdas) + r' \\%')
        for gam, n in L:
            x, h, xD, D, _ = carga_alta(n)
            s = [pendiente(x, h[0], 3, 10), pendiente(x, h[0], 10, 20), pendiente(x, h[0], 20, 30)]
            sD = [pendiente(xD, D, 3, 10), pendiente(xD, D, 10, 20)]
            h20 = h[0][np.argmin(abs(x - 20))]
            q1, q0 = (exponente_perfil(n, 1), exponente_perfil(n, 0)) if forma in ('p2', 'p4') else (np.nan, np.nan)
            jin = borde_interior(n) if forma in ('h05', 'anillo', 'meseta') else np.nan
            m = masa_alta(n)
            pend[forma].append((gam, s[1]))
            print(f'{"":7} {gam:5.3f}  {m:7.5f} | {s[0]:10.2f} {s[1]:5.2f} {s[2]:5.2f} | {sD[0]:7.2f} {sD[1]:5.2f} | {h20:9.1e} | {q1:4.2f} {q0:5.2f} | {jin:.2f}')
            celdas = ['', f'${gam:.3f}$', f'${m:g}$', num(s[0]), num(s[1]),
                      num(sD[0]), num(sD[1]), num(h20), num(q1), num(jin)]
            filas.append(' & '.join(celdas) + r' \\%')
        filas.append(r'\hline%')
    filas = filas[:-1]
    with open(os.path.join(art, 'tab_colas.tex'), 'w') as fh:
        fh.write('% Escrito por reproducir/scripts/colas.py articulo; no editar a mano.\n\\newcommand{\\TabColas}{%\n'
                 + '\n'.join(filas) + '\n}\n')
    print('escrita', os.path.join(art, 'tab_colas.tex'))
    # Convergencia: retícula de 3200 x 100 frente a 6400 x 200 (Gamma = 0.10).
    x, h, xD, D, _ = carga_alta('H_p2_100'); x2, h2, xD2, D2, _ = carga_alta('H2_p2_100')
    for b in (20, 30):
        z = (x >= 3) & (x <= b); zD = (xD >= 3) & (xD <= b)
        print(f'convergencia en 3-{b} tau_1, mayor diferencia relativa de las envolventes: '
              + ', '.join(f'|h_{k+1}| {np.abs(h[k][z]/h2[k][z] - 1).max():.3f}' for k in range(3))
              + f', D {np.abs(D[zD]/D2[zD] - 1).max():.3f}')
    # Las partículas de las corridas (con peso).
    for n in ('H_p2_100', 'H2_p2_100'):
        A = np.load(os.path.join(SP, n, 'medidas.npz'))
        print(n, 'tau_1 =', f'{float(A["tau"]):.1f}', '; t final / tau_1 =', f'{A["t"][-1]/float(A["tau"]):.1f}')
    for forma, L in ALTA.items():
        print(f'tau_1 de {forma}: {R.banda_perfil(F.TODAS[forma][0], SJ)[3]:.1f}')
    # Figura.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    fig, ax = plt.subplots(2, 3, figsize=(10, 6.6))
    azules = ['#86b6ef', '#4f93e3', '#2a78d6', '#11468f']
    tau = R.banda_perfil(F.TODAS['p2'][0], SJ)[3]; ts = np.arange(0.0, 30.6*tau, 10.0)
    a = ax[0, 0]
    a.plot(ts[1:]/tau, envolvente(h_libre('p2', 1, ts), ts)[1:], color=tinta2, lw=1.0, ls=(0, (4, 2.5)), label=r'$\Gamma=0$')
    for (gam, n), col in zip(ALTA['p2'], azules):
        x, h, xD, D, _ = carga_alta(n)
        a.plot(x[1:], h[0][1:], color=col, lw=1.4, label=rf'$\Gamma={gam:.3f}$')
        ax[0, 1].plot(xD[1:], D[1:], color=col, lw=1.4)
    a.set_ylim(5e-6, 1.5); a.legend(frameon=False, fontsize=7.5, loc='lower left')
    a.set_title(r'(a) $|h_1|/h_0$, reference shape', loc='left', fontsize=9.5, pad=6)
    a = ax[0, 1]; a.set_ylim(1e-6, 0.2)
    a.set_title(r'(b) deviation $\|\delta\Phi\|$ of the potential', loc='left', fontsize=9.5, pad=6)
    a = ax[0, 2]
    xp = np.linspace(0.005, 0.395, 40)
    a.plot(xp, ive(1, 0.5/SQ**2)/ive(0, 0.5/SQ**2)*xp**2*np.exp(-xp**2)/(SJ*np.sqrt(np.pi)/4), color=tinta2, lw=1.0,
           ls=(0, (4, 2.5)), label=r'initial datum: $q_1=2$')
    for (gam, n), col in zip(ALTA['p2'], azules):
        P = np.load(os.path.join(SP, n, 'perfiles.npz')); xm = 0.5*(P['x'][1:] + P['x'][:-1])
        a.plot(xm, P['a'][2, 0], color=col, lw=1.4, label=rf'$\Gamma={gam:.3f}$: $q_1={exponente_perfil(n):.2f}$')
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel(r'$J/\sigma_J$'); a.set_xlim(0.005, 0.4)
    a.xaxis.set_major_locator(FixedLocator((0.01, 0.03, 0.1, 0.3))); a.xaxis.set_minor_formatter(NullFormatter())
    a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    a.set_title(r'(c) first harmonic near the circular orbit', loc='left', fontsize=9.5, pad=6); a.set_ylabel(r'$a_1$')
    a.legend(frameon=False, fontsize=7.5, loc='lower right')
    a = ax[1, 0]
    for forma in ('p4', 'h05', 'anillo', 'meseta'):
        col = F.TODAS[forma][3]
        tauf = R.banda_perfil(F.TODAS[forma][0], SJ)[3]; tsf = np.arange(0.0, 30.6*tauf, 10.0)
        a.plot(tsf[1:]/tauf, envolvente(h_libre(forma, 1, tsf), tsf)[1:], color=col, lw=0.9, ls=(0, (4, 2.5)))
        for (gam, n), lw in zip(ALTA[forma], (0.8, 1.5)):
            x, h, xD, D, _ = carga_alta(n)
            a.plot(x[1:], h[0][1:], color=col, lw=lw, alpha=0.75 if lw < 1 else 1.0)
    a.set_ylim(5e-8, 1.5)
    a.set_title(r'(d) $|h_1|/h_0$, other shapes', loc='left', fontsize=9.5, pad=6)
    a.legend(handles=[Line2D([], [], color=F.TODAS[f][3], lw=1.5, label=ROTULO[f]) for f in ('p4', 'h05', 'anillo', 'meseta')]
             + [Line2D([], [], color=tinta2, lw=0.9, ls=(0, (4, 2.5)), label=r'$\Gamma=0$'),
                Line2D([], [], color=tinta2, lw=0.8, label=r'$\Gamma\approx 0.05$'), Line2D([], [], color=tinta2, lw=1.5, label=r'$\Gamma\approx 0.1$')],
             frameon=False, fontsize=7, loc='lower left', ncol=2, handlelength=2.0, columnspacing=1.0)
    a = ax[1, 1]
    for forma, ls in (('h05', '-'), ('meseta', (0, (1, 1)))):
        col = F.TODAS[forma][3]
        Jg = np.linspace(0.0, 3*SJ, 3001); f = R.perfil_J(Jg, SJ, F.TODAS[forma][0]); f = f/(f.sum()*(Jg[1] - Jg[0]))
        a.plot(Jg/SJ, f, color=tinta2, lw=0.9, ls=ls)
        for (gam, n), al in zip(ALTA[forma], (0.5, 1.0)):
            P = np.load(os.path.join(SP, n, 'perfiles_ancho.npz')); xm = 0.5*(P['x'][1:] + P['x'][:-1])
            a.plot(xm, P['f0'][2], color=col, lw=1.4, ls=ls, alpha=al)
    a.set_xlim(0, 1.2); a.set_ylim(0, 10); a.set_xlabel(r'$J/\sigma_J$')
    a.set_title(r'(e) inner edge of two hollow shapes', loc='left', fontsize=9.5, pad=6); a.set_ylabel(r'mass per unit $J$')
    a.legend(handles=[Line2D([], [], color=F.TODAS['h05'][3], lw=1.4, label=ROTULO['h05']),
                      Line2D([], [], color=F.TODAS['meseta'][3], lw=1.4, ls=(0, (1, 1)), label=ROTULO['meseta']),
                      Line2D([], [], color=tinta2, lw=0.9, label=r'$t=0$'),
                      Line2D([], [], color=tinta2, lw=1.4, alpha=0.5, label=r'$\Gamma\approx 0.05$'),
                      Line2D([], [], color=tinta2, lw=1.4, label=r'$\Gamma\approx 0.1$')], frameon=False, fontsize=7, loc='upper left', ncol=2,
             handlelength=2.0, columnspacing=1.0)
    a = ax[1, 2]
    for forma in ALTA:
        col, mk = F.TODAS[forma][3], F.TODAS[forma][4]
        g, s = np.array(pend[forma]).T
        a.plot(g, s, color=col, marker=mk, ms=5.5, mec=sup, mew=0.6, lw=1.0, label=ROTULO[forma])
    for y in (-1.0, -3.0):
        a.axhline(y, color=tinta2, lw=0.7, ls=':')
    a.set_xlim(-0.01, 0.16); a.set_ylim(-5.4, -0.6); a.set_xlabel(r'$\Gamma$')
    a.set_title(r'(f) slope of $|h_1|$ over $10\leq t/\tau_1\leq 20$', loc='left', fontsize=9.5, pad=6)
    a.legend(frameon=False, fontsize=7, loc='lower right')
    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1, 31); a.set_xlabel(r'$t/\tau_1$')
        a.xaxis.set_major_locator(FixedLocator((1, 2, 5, 10, 20, 30))); a.xaxis.set_minor_formatter(NullFormatter())
        a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    for a in ax.ravel():
        a.grid(True, which='major', color=reja, lw=0.6); a.set_axisbelow(True)
        a.tick_params(length=3, width=0.6); a.tick_params(which='minor', length=2, width=0.4)
        for lado in ('top', 'right'):
            a.spines[lado].set_visible(False)
    fig.tight_layout(w_pad=1.2, h_pad=1.4)
    sal = os.path.join(art, 'figuras', 'colas_alta.pdf')
    fig.savefig(sal); os.makedirs(SAL, exist_ok=True); fig.savefig(os.path.join(SAL, 'colas_alta_articulo.png'), dpi=110)
    plt.close(fig)
    print('escrita', sal)


if __name__ == '__main__':
    pasos = dict(libre=libre, tabla=tabla, figura=figura, alta=alta, figura_alta=figura_alta, articulo=articulo)
    if len(sys.argv) < 2 or sys.argv[1] not in pasos:
        print(__doc__)
    else:
        pasos[sys.argv[1]]()
