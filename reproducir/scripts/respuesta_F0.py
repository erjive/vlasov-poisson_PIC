"""Análisis de la prueba con F_0 lento (respuesta_mezcla.py: perfiles, lineal --perfiles --cuasilineal).

Pasos:
    python3 respuesta_F0.py tabla <nombre> <perf.npz> <ext.npz> <lin.npz> [<lin2.npz>]
    python3 respuesta_F0.py figura <salida.png> <nombre> <perf.npz> <lin.npz> [<lin2.npz>] [-- <nombre> <perf> <lin> ...]

tabla, por corrida:
1. Cambio de F_0: el medido (F0_lento, suavizado con el núcleo de la evolución) frente al cuasilineal, dG/dJ con G
   el flujo integrado, promediado en +-tau_1/2 como el medido. Amplitud en unidades del máximo de F_0 suavizado (el
   del histograma crudo es unas 4 veces mayor por las filas de la retícula), correlación y proyección.
2. Deriva del fondo: Omega(J) de los potenciales medios en ventanas de tau_1 frente a Omega_0(J), con la acción como
   invariante adiabático; dispersión, pesada con F_0, de la fase acumulada (una fase igual en todo J no cambia |h_1|).
3. Dónde está la respuesta en J: la corrida menos la mezcla libre exacta de las mismas partículas (armónico 1,
   complejo, franjas de 4 filas) frente a la respuesta de la rejilla de cada variante (lineal - libre), con ventanas
   triangulares en J que suman 1 (nodos en los cuantiles de masa 0.02, 0.25, 0.5, 0.75, 0.98). Con bordes duros no
   sirve: las filas de la retícula oscilan en J con la respuesta y cruzan los bordes de golpe, y en las bandas sin
   respuesta lineal la corrida tendría una "respuesta" de 6 a 50 veces la predicha que se cancela en el total.
figura: fila 1, F_0(J) en t0 suavizado; filas 2-4, Delta F_0 medido y cuasilineal en 10, 15 y 20 tau_1, frente a
J/J_med.
"""
import os, sys, time
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from respuesta_mezcla import F0_lento, tabla_E_J
from aa_numerico import MapaAA
from landau_libre import omega_de_J
from scipy.ndimage import gaussian_filter1d

TS = (8, 12, 16, 20, 25)


def cerca(x, v):
    return int(np.argmin(np.abs(x - v)))


def cargar(lins):
    """Junta los archivos de `lineal`; la cuasilineal sale del que tenga G, y solo_cl del de los armónicos."""
    d = {}
    for f in lins:
        z = np.load(f); d.update({k: z[k] for k in z.files if k not in d and k != 'solo_cl'})
        if 'G' in z.files:
            d.update({k: z[k] for k in ('G', 'tG', 'h1_corc', 'Hg_linc') if k in z.files})
        if 'Hg_lin' in z.files:
            d['solo_cl'] = bool(z['solo_cl']) if 'solo_cl' in z.files else False
    return d


def cambio_cuasilineal(d, s, dJ, tau):
    """Delta F_0 cuasilineal, dG/dJ suavizado y promediado en +-tau_1/2 (como F0_lento), desde t0."""
    q0 = gaussian_filter1d(d['G'], s, axis=1, order=1, mode='constant')/dJ
    tG = d['tG']
    dq = np.array([q0[np.abs(tG - x) <= 0.5*tau + 1e-9].mean(0) for x in tG])
    return tG, dq - dq[0]


def tabla(nombre, perf, ext, *lins):
    p = np.load(perf); e = np.load(ext); d = cargar(lins)
    tau, t0, L = float(p['tau']), float(p['t0']), float(p['L'])
    nj, dJ = int(p['nj']), float(p['dJ'])
    c = 8*np.pi**2*L*2*np.pi*dJ
    s = float(d.get('suave', 8.0))
    F0, Jn = d['F0'], d['Jn']
    mF = F0/F0.sum()
    Fmax = gaussian_filter1d(F0, 8.0, mode='nearest').max()
    Jmed = Jn[np.searchsorted(np.cumsum(mF), 0.5)]
    print(f'\n=== {nombre}: t0 = {t0/tau:.2f} tau_1, J_med = {Jmed:.4g}, suave = {s:g} franjas '
          f'({s*dJ/Jmed:.3f} J_med); err_F0 = {float(d.get("err_F0", np.nan)):.1e}')
    fuera = p['fuera']/float(e['mgas'])
    print(f'  masa fuera de jmax: máx {fuera.max():.1e};  sin mapa: máx {int(p["malos"].max())} partículas')

    # 1. Cambio de F_0
    tF, dlt = F0_lento(p, tau, t0, c)
    dm = gaussian_filter1d(dlt, s, axis=1, order=0, mode='nearest')
    tiene_cl = 'G' in d
    if tiene_cl:
        tG, dq = cambio_cuasilineal(d, s, dJ, tau)
    print('  1. Cambio de F_0 (suavizado): máx|dF_0|/máx F_0 medido' + (' y cuasilineal; correlación; medido/cuasilineal'
                                                                      ' (proyección); J del máximo/J_med' if tiene_cl else ''))
    for T in TS:
        a = dm[cerca(tF, T*tau)]
        txt = f'     {T:2d} tau_1: medido {np.abs(a).max()/Fmax:6.3f}'
        if tiene_cl:
            b = dq[cerca(tG, T*tau)]
            r = np.dot(a, b)/np.sqrt(np.dot(a, a)*np.dot(b, b)); pr = np.dot(a, b)/np.dot(b, b)
            txt += (f'  cuasilineal {np.abs(b).max()/Fmax:6.3f}  corr {r:5.2f}  proy {pr:5.2f}'
                    f'  J_máx {Jn[np.argmax(np.abs(a))]/Jmed:5.2f} / {Jn[np.argmax(np.abs(b))]/Jmed:5.2f}')
        print(txt)

    # 2. Deriva del fondo
    t1 = time.time()
    jmax = float(p['jmax'])
    sel = [k for k, T in enumerate(p['ventanas']) if abs(T - round(T)) < 1e-9]
    om = []
    for k in sel:
        m = MapaAA(p['rg'], p['phiw'][k], L=L)
        E_t, J_t = tabla_E_J(m, jmax)
        om.append(omega_de_J(E_t, J_t)(Jn))
    om = np.array(om); tw = p['tw'][sel]
    dfase = np.concatenate([np.zeros((1, nj)), np.cumsum(0.5*(om[1:] - om[0] + om[:-1] - om[0])*np.diff(tw)[:, None], 0)])
    dfase -= (dfase*mF).sum(1)[:, None]
    disp = np.sqrt((dfase**2*mF).sum(1))
    print(f'  2. Deriva del fondo ({time.time() - t1:.0f} s): dispersión de la fase acumulada (rad), pesada con F_0; '
          f'Phi_0 frente a la ventana 1: {np.abs(p["phiw"][0] - e["phi0"]).max():.1e}')
    print('     ' + '  '.join(f'{T:.0f}: {disp[cerca(tw, T*tau)]:.3f}' for T in (10, 15, 20, 25, 30)))
    dOm = (om[-1] - om[0])/np.abs(om[0]).max()
    print(f'     cambio de Omega en 30 tau_1 respecto de máx|Omega_0|: pesado {np.sqrt((dOm**2*mF).sum()):.1e}, '
          f'máx en el soporte {np.abs(dOm[mF > 1e-3*mF.max()]).max():.1e}')

    # 3. Respuesta por ventanas triangulares en J
    if 'Hg_lib' not in d:
        return
    ok = np.isfinite(e['Q0']) & np.isfinite(e['J0'])
    Q0, J0, mi = e['Q0'][ok], e['J0'][ok], e['m'][ok]
    mg = float(e['mgas'])
    Om = np.interp(J0, Jn, d['om'])
    z = np.clip((J0/dJ).astype(int), 0, nj - 1)//4
    tH = d['tHg']; ta = p['t_arm']
    nb = nj//4
    Hx = np.array([np.bincount(z, weights=mi*np.cos(Q0 + Om*(x - t0)), minlength=nb)
                   - 1j*np.bincount(z, weights=mi*np.sin(Q0 + Om*(x - t0)), minlength=nb) for x in ta])
    Jb = (np.arange(nb)*4 + 2)*dJ
    cm = np.cumsum(F0)/F0.sum()
    nodos = Jn[np.searchsorted(cm, [0.02, 0.25, 0.5, 0.75, 0.98])]
    ventanas = []
    for q in range(len(nodos)):
        y = np.zeros(len(nodos)); y[q] = 1.0
        ventanas.append(np.interp(Jb, nodos, y))
    ventanas.append(np.ones(nb))
    variantes = [(k, n) for k, n in (('Hg_lin', 'F_0 fijo'), ('Hg_linv', 'F_0 lento'), ('Hg_linc', 'cuasilineal'))
                 if k in d and not (k == 'Hg_lin' and d.get('solo_cl', False))]
    ir = [cerca(ta, T*tau) for T in TS]; ig = [cerca(tH, T*tau) for T in TS]
    Rrun = (p['H'][:, 0, :] - Hx)/mg
    tot = np.abs(p['H'][:, 0, :].sum(1))/mg
    print('  3. Respuesta (corrida - mezcla libre exacta) del armónico 1 con ventanas triangulares en J (nodos en los '
          'cuantiles de masa 0.02, 0.25, 0.5, 0.75, 0.98), en ' + ', '.join(f'{T}' for T in TS) + ' tau_1:')
    print('     |R| de la corrida, y por variante |R_corrida|/|R_curva| y fase de R_corrida/R_curva (grados)')
    print('     |h_1| total de la corrida: ' + ' '.join(f'{tot[i]:.1e}' for i in ir))
    for q, W in enumerate(ventanas):
        rr = Rrun @ W
        nombre = 'total' if q == len(ventanas) - 1 else f'nodo J/J_med {nodos[q]/Jmed:4.2f}'
        print(f'     {nombre:20} |R| ' + ' '.join(f'{abs(rr[i]):.1e}' for i in ir))
        for k, n in variantes:
            rv = ((d[k][:, 0, :] - d['Hg_lib'][:, 0, :])/mg) @ W
            print(f'        {n:12} ' + ' '.join(f'{abs(rr[i])/abs(rv[j]):5.2f}' for i, j in zip(ir, ig))
                  + '   fase ' + ' '.join(f'{np.degrees(np.angle(rr[i]/rv[j])):5.0f}' for i, j in zip(ir, ig)))


def figura(sal, *args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    grupos, g = [], []
    for x in args:
        if x == '--':
            grupos.append(g); g = []
        else:
            g.append(x)
    grupos.append(g)
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    MED, CL = '#1baf7a', '#eda100'
    plt.rcParams.update({'font.size': 8.5, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    TT = (10, 15, 20)
    nc = len(grupos)
    fig, ax = plt.subplots(1 + len(TT), nc, figsize=(2.9*nc + 0.6, 7.2), squeeze=False, sharex='col')
    for j, g in enumerate(grupos):
        nombre, perf, lins = g[0], g[1], g[2:]
        p = np.load(perf); d = cargar(lins)
        tau, t0, L = float(p['tau']), float(p['t0']), float(p['L'])
        dJ = float(p['dJ']); c = 8*np.pi**2*L*2*np.pi*dJ; s = float(d.get('suave', 8.0))
        F0, Jn = d['F0'], d['Jn']
        Fmax = gaussian_filter1d(F0, s, mode='nearest').max()
        cm = np.cumsum(F0)/F0.sum(); Jmed = Jn[np.searchsorted(cm, 0.5)]; x = Jn/Jmed
        Jhi = Jn[np.searchsorted(cm, 0.995)]/Jmed
        tF, dlt = F0_lento(p, tau, t0, c)
        dm = gaussian_filter1d(dlt, s, axis=1, order=0, mode='nearest')/Fmax
        tiene = 'G' in d
        if tiene:
            tG, dq = cambio_cuasilineal(d, s, dJ, tau); dq = dq/Fmax
        a = ax[0, j]
        a.plot(x, gaussian_filter1d(F0, s, mode='nearest')/Fmax, color=tinta, lw=1.2)
        a.set_title(nombre, loc='left', fontsize=9.5, pad=4)
        if j == 0:
            a.set_ylabel(r'$F_0/\max F_0$ at $t_0$')
        for i, T in enumerate(TT):
            a = ax[1 + i, j]
            a.axhline(0, color=reja, lw=0.8)
            a.plot(x, dm[cerca(tF, T*tau)], color=MED, lw=1.4)
            if tiene:
                a.plot(x, dq[cerca(tG, T*tau)], color=CL, lw=1.4, ls=(0, (4, 2)))
            a.text(0.98, 0.94, f'{T} $\\tau_1$', transform=a.transAxes, ha='right', va='top', color=tinta2)
            if j == 0:
                a.set_ylabel(r'$\Delta F_0/\max F_0$')
        for a in ax[:, j]:
            a.set_xlim(0, Jhi); a.grid(True, color=reja, lw=0.5); a.set_axisbelow(True)
            for lado in ('top', 'right'):
                a.spines[lado].set_visible(False)
        ax[-1, j].set_xlabel(r'$J/J_{\rm med}$')
    fig.legend(handles=[Line2D([], [], color=MED, lw=1.4, label="slow part of the run's $F_0$ change"),
                        Line2D([], [], color=CL, lw=1.4, ls=(0, (4, 2)), label='quasilinear prediction')],
               loc='upper center', ncol=2, frameon=False, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(sal, dpi=120); plt.close(fig)
    print('escrita', sal)


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('tabla', 'figura'):
        sys.exit(__doc__)
    (tabla if sys.argv[1] == 'tabla' else figura)(*sys.argv[2:])
