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
    return Fo['g']*B.info[n]['a0']/(3.0*Fo['c']**-4*(F.J_mediana(None)*sj/SJ)**1.5)


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
        gam = R.G*m/(abs(F.OMEGA_P)*F.J_mediana(F.TODAS[forma][0])**1.5)
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


if __name__ == '__main__':
    pasos = dict(libre=libre, tabla=tabla, figura=figura)
    if len(sys.argv) < 2 or sys.argv[1] not in pasos:
        print(__doc__)
    else:
        pasos[sys.argv[1]]()
