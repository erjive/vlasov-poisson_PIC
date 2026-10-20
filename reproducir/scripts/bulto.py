"""El bulto: el dato inicial de 09_autogravedad a varias masas.

F0 = exp(-sin^2(Q/2)/sp^2) J^2 exp(-J^2/sr^2), sp = sr = 0.1, en la malla regular
(Q, J) del isócrono desnudo (state = aa_quad, dftype = gauss: lo genera el propio
código). No es un equilibrio más una perturbación: toda la masa está concentrada
en una fase orbital, con armónicos hasta k ~ 1/sp. La pregunta es si, al crecer
la masa, la autogravedad del bulto lo mantiene unido contra la cizalla de la
mezcla de fases.

Corridas: B0 sin autogravedad (la mezcla libre), B1..B5 con a0 = 1e-3 ... 0.3 y
Npc = 25 como el original, y B4q (a0 = 0.1, Npc = 100) para la resolución en Q.

Pasos: preparar, correr, analizar, figuras, videos (como demo_eta.py). Todo
queda en exe/demo_eta/bulto/; figuras en docs/demo_eta/figuras/bulto_*; videos en
reproducir/videos/demo_eta/bulto_*.mp4.

El mapa (Q, J) de cada instantánea se calcula en el potencial de ESE instante
(el dato no es estacionario y no hay un equilibrio fijo que sirva de marco).
"""
import os, sys, subprocess, time, argparse, warnings, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
warnings.filterwarnings('ignore', category=RuntimeWarning)
from demo_eta import envolvente, norma, metadatos, RAIZ, EXE, FIGDIR, VIDDIR

BASE = os.path.join(EXE, 'demo_eta', 'bulto')
PARDIR = os.path.join(RAIZ, 'reproducir', 'corridas', '12_demo_eta')
PLANTILLA = os.path.join(RAIZ, 'reproducir', 'corridas', '09_autogravedad', 'sg__long_a0_1e-2.par')
TFIN, DT, SALIDA = 10000, 0.05, 200
RVENT = (3.0, 15.0)
# Banda del isócrono desnudo para J^2 exp(-J^2/0.1^2), umbral del 1 % (eta.py).
OM_MIN, OM_MAX, DJ_SOP = 0.0513, 0.0705, 0.270
TAU1 = 2*np.pi/(OM_MAX - OM_MIN)

CORRIDAS = [
    ('B0', 0.01, 25, False, 'no self-gravity: free mixing'),
    ('B1', 0.001, 25, True, r'$a_0=10^{-3}$'),
    ('B2', 0.01, 25, True, r'$a_0=10^{-2}$ (the largest in 09\_autogravedad)'),
    ('B3', 0.03, 25, True, r'$a_0=0.03$'),
    ('B4', 0.1, 25, True, r'$a_0=0.1$'),
    ('B5', 0.3, 25, True, r'$a_0=0.3$'),
    ('B4q', 0.1, 100, True, r'$a_0=0.1$, $N_{pc}=100$'),
]
# Segunda tanda (2026-10-19): la transición en masa entre a0 = 0.01 y 0.1, y la convergencia
# del residuo de B2 (|h1| ~ 3e-3, 40 veces la mezcla libre) y de un caso de la transición.
CORRIDAS += [
    ('B015', 0.015, 25, True, r'$a_0=0.015$'),
    ('B02', 0.02, 25, True, r'$a_0=0.02$'),
    ('B025', 0.025, 25, True, r'$a_0=0.025$'),
    ('B04', 0.04, 25, True, r'$a_0=0.04$'),
    ('B05', 0.05, 25, True, r'$a_0=0.05$'),
    ('B07', 0.07, 25, True, r'$a_0=0.07$'),
    ('B2q', 0.01, 100, True, r'$a_0=10^{-2}$, $N_{pc}=100$'),
    ('B2j', 0.01, 25, True, r'$a_0=10^{-2}$, $N_{rc}=800$'),
    ('B02q', 0.02, 100, True, r'$a_0=0.02$, $N_{pc}=100$'),
]
CAMBIOS = {'B2j': {'Nrc': '800'}}              # cambios en el .par fuera de los de preparar()
# Tercera tanda (2026-10-19): si la transición cae en el mismo mu cuando cambia el ancho del bulto
# en ángulo (sp = s_Q) o en acción (sr = s_J), y el mismo dato en el potencial de Kepler (el de
# Chaturvedi y Luk), con las partículas colocadas desde Python en sus variables de ángulo-acción.
for sq, npc, masas in ((0.2, 50, (0.01, 0.02, 0.04, 0.08)), (0.4, 25, (0.01, 0.02, 0.04, 0.08))):
    for m in masas:
        n = f'Q{int(sq*10)}_{m:g}'
        CORRIDAS.append((n, m, npc, True, rf'$s_Q={sq:g}$, $a_0={m:g}$'))
        CAMBIOS[n] = {'sp': str(sq)}
for sj, nrc, masas in ((0.05, 400, (0.005, 0.01, 0.02)), (0.2, 800, (0.02, 0.04, 0.08))):
    for m in masas:
        n = f'J{int(round(sj*100)):02d}_{m:g}'
        CORRIDAS.append((n, m, 25, True, rf'$s_J={sj:g}$, $a_0={m:g}$'))
        CAMBIOS[n] = {'sr': str(sj), 'Nrc': str(nrc)}
        if sj < 0.1:                           # banda más estrecha, tau_1 = 589: el doble de tiempo
            CAMBIOS[n]['Nt'] = str(2*int(round(TFIN/DT)))
FONDO = {}                                     # 'puntual': masa puntual (BGtype = sphere, r > 1)
for m in (0.01, 0.02, 0.04, 0.08):
    n = f'K_{m:g}'
    CORRIDAS.append((n, m, 25, True, rf'Kepler, $a_0={m:g}$'))
    CAMBIOS[n] = {'BGtype': 'sphere', 'state': 'checkpoint', 'checkpointfile': f'demo_eta/bulto/ic/{n}.dat'}
    FONDO[n] = 'puntual'
info = {c[0]: dict(a0=c[1], npc=c[2], sg=c[3], desc=c[4]) for c in CORRIDAS}


def parametro(n, k, defecto):
    return type(defecto)(CAMBIOS.get(n, {}).get(k, defecto))


def banda(n):
    """(Omega_min, Omega_max, ancho del soporte en J, tau_1) del fondo desnudo para el bulto n.
    Con s_J = 0.1 en el isócrono son las constantes de siempre."""
    sr, fondo = parametro(n, 'sr', 0.1), FONDO.get(n, 'isocrono')
    if sr == 0.1 and fondo == 'isocrono':
        return OM_MIN, OM_MAX, DJ_SOP, TAU1
    from scipy.optimize import brentq
    g = lambda x: x*np.exp(1 - x) - 0.01                  # (J/s_J)^2 e^{1-(J/s_J)^2} = 0.01
    xlo, xhi = brentq(g, 1e-9, 1.0), brentq(g, 1.0, 50.0)
    jlo, jhi = sr*np.sqrt(xlo), sr*np.sqrt(xhi)
    c = 0.5*(2.0 + np.sqrt(8.0)) if fondo == 'isocrono' else 2.0   # Omega = (J + c)^-3, L = 2
    omin, omax = (jhi + c)**-3, (jlo + c)**-3
    return omin, omax, jhi - jlo, 2*np.pi/(omax - omin)


def dato_kepler(n):
    """Estado inicial del bulto en el potencial de Kepler (GM = 1, L = 2): la misma F0(Q, J) que
    aa_quad, en una retícula de puntos medios de las variables de ángulo-acción de Kepler,
    J en [0, 6 s_J] y Q en [0, 2 pi); anomalía media Q = eta - e sin(eta)."""
    nrc, npc = parametro(n, 'Nrc', 400), info[n]['npc']
    sr, sp, L = parametro(n, 'sr', 0.1), parametro(n, 'sp', 0.1), 2.0
    J = (np.arange(nrc) + 0.5)*6*sr/nrc
    Q = (np.arange(npc) + 0.5)*2*np.pi/npc
    JJ, QQ = np.meshgrid(J, Q, indexing='ij')
    a = (JJ + L)**2
    e = np.sqrt(np.clip(1 - L**2/a, 0, None))
    eta = QQ.copy()
    for _ in range(60):
        eta = eta - (eta - e*np.sin(eta) - QQ)/(1 - e*np.cos(eta))
    r = a*(1 - e*np.cos(eta))
    p = np.sqrt(a)*e*np.sin(eta)/r
    F = np.exp(-np.sin(0.5*QQ)**2/sp**2)*np.exp(-JJ**2/sr**2)*JJ**2
    return r.ravel(), p.ravel(), F.ravel(), QQ.ravel(), JJ.ravel()


def ruta(*p):
    return os.path.join(BASE, *p)


def preparar():
    os.makedirs(BASE, exist_ok=True)
    lineas0 = open(PLANTILLA).read().split('\n')
    for n, a0, npc, sg, _ in CORRIDAS:
        v = {'courant': '1.0', 'Nt': str(int(round(TFIN/DT))), 'Npc': str(npc),
             'time_output': str(SALIDA), 'spatial_output': str(SALIDA), 'field_output': str(SALIDA),
             'directory': f'demo_eta/bulto/{n}', 'a0': str(a0),
             'autointeraction': '.true.' if sg else '.false.'}
        v.update(CAMBIOS.get(n, {}))
        out = [f'# Bulto {n}: el dato de 09_autogravedad (aa_quad, gauss, sr = sp = 0.1),',
               f'# a0 = {a0}, Npc = {npc}, autogravedad = {sg}. Generado por bulto.py.']
        if n in CAMBIOS:
            out.append('# Cambios: ' + ', '.join(f'{k} = {x}' for k, x in CAMBIOS[n].items()) + '.')
        for l in lineas0:
            if l.startswith('#') or '=' not in l:
                continue
            k = l.split('=')[0].strip()
            out.append(f'{k:<16} = {v[k]}' if k in v else l)
        claves = {l.split('=')[0].strip() for l in lineas0 if '=' in l and not l.startswith('#')}
        out += [f'{k:<16} = {x}' for k, x in CAMBIOS.get(n, {}).items() if k not in claves]
        open(os.path.join(PARDIR, f'bulto__{n}.par'), 'w').write('\n'.join(out) + '\n')
    for n in FONDO:
        f = os.path.join(EXE, CAMBIOS[n]['checkpointfile'])
        if not os.path.exists(f):
            os.makedirs(os.path.dirname(f), exist_ok=True)
            r, p, F, _, _ = dato_kepler(n)
            np.savetxt(f, np.column_stack([r, p, F]), fmt='%.17e')
    print('preparado:', len(CORRIDAS))


def correr():
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    for n, *_ in CORRIDAS:
        if os.path.exists(ruta(f'{n}.ok')):
            continue
        t0 = time.time()
        r = subprocess.run(['./VP_PIC', os.path.join(PARDIR, f'bulto__{n}.par')], cwd=EXE,
                           stdout=open(ruta(f'{n}.log'), 'w'), stderr=subprocess.STDOUT, env=env)
        print(f'{"OK" if r.returncode == 0 else "FALLO":5} {n} ({time.time()-t0:.0f} s)', flush=True)
        if r.returncode == 0:
            open(ruta(f'{n}.ok'), 'w').write('')
            metadatos(n, os.path.join(PARDIR, f'bulto__{n}.par'), ruta(f'{n}.meta'))


def _mapear(tarea):
    import h5py
    from aa_numerico import MapaAA, phi_iso
    h5, claves, sg, fondo = tarea
    f = h5py.File(h5, 'r')
    rg = f['grid']['r'][:]
    desnudo = MapaAA(L=2.0, fondo=fondo)
    # potencial de fondo del código: isócrono, o BGtype = sphere (esfera homogénea de radio 1)
    phib = phi_iso(rg) if fondo == 'isocrono' else np.where(np.abs(rg) < 1, 0.5*(rg**2 - 3), -1/np.maximum(np.abs(rg), 1e-12))
    sal = []
    for k in claves:
        g = f[k]
        r, p = g['r_part'][:], g['p_part'][:]
        ps = g['potential'][:] - phib if sg else np.zeros_like(rg)
        if fondo != 'isocrono':
            ps = np.where(rg < 1, ps[np.argmin(np.abs(rg - 1))], ps)   # las órbitas no bajan de r = 1
        Q, J, E = (MapaAA(rg, ps, L=2.0, fondo=fondo) if sg else desnudo)(r, p)
        Qd, Jd = (desnudo(r, p)[:2] if sg else (Q, J))
        sal.append((g.attrs['time'], r, p, Q, J, ps, Qd, Jd))
    f.close()
    return sal


def analizar(solo=None, procesos=4):
    import h5py
    from multiprocessing import Pool
    for n, *_ in CORRIDAS:
        if (solo and n not in solo) or os.path.exists(ruta(n, 'bulto.npz')):
            continue
        t0 = time.time()
        h5 = ruta(n, 'vlasov_output.h5')
        f = h5py.File(h5, 'r')
        pasos = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
        rg, w = f['grid']['r'][:], f[pasos[0]]['f'][:]
        f.close()
        grupos = [(h5, pasos[i:i+25], info[n]['sg'], FONDO.get(n, 'isocrono')) for i in range(0, len(pasos), 25)]
        with Pool(procesos) as pool:
            res = [x for b in pool.map(_mapear, grupos) for x in b]
        t = np.array([x[0] for x in res])
        R, P, Q, J, ps, Qd, Jd = (np.array([x[i] for x in res]) for i in range(1, 8))
        hk = np.stack([np.sum(w*np.exp(-1j*m*Q), axis=1) for m in range(6)], axis=1)
        hk = hk/hk[:, :1].real
        hkd = np.stack([np.sum(w*np.exp(-1j*m*Qd), axis=1) for m in range(6)], axis=1)
        hkd = hkd/hkd[:, :1].real
        np.savez(ruta(n, 'bulto.npz'), t=t, hk=hk, hk_desnudo=hkd, phi_self=ps, r=rg, w=w,
                 R=R.astype(np.float32), P=P.astype(np.float32), Q=Q.astype(np.float32),
                 J=J.astype(np.float32))
        print(f'  analizado {n}: {len(t)} instantáneas ({time.time()-t0:.0f} s)', flush=True)


def cargar(n):
    return np.load(ruta(n, 'bulto.npz'), mmap_mode='r')


def potencial(n):
    """Phi_self(r, t) en la ventana radial, su media en la segunda mitad y la
    fluctuación alrededor de ella, relativa a max|media|."""
    d = cargar(n)
    t, rg, ps = d['t'], d['r'], np.asarray(d['phi_self'])
    z = (rg >= RVENT[0]) & (rg <= RVENT[1])
    media = ps[len(t)//2:].mean(0)
    esc = np.abs(media[z]).max() if np.abs(media[z]).max() > 0 else 1.0
    return t, rg, ps, media, norma((ps - media)[:, z])/esc, esc


def escalas(n):
    """omega_b del bulto a partir de la oscilación de su potencial en la primera
    vuelta radial, y mu = omega_b/Delta Omega con la banda del isócrono."""
    t, rg, ps, media, fl, esc = potencial(n)
    omin, omax, djs, _ = banda(n)
    v = t <= 2*np.pi/omin
    z = (rg >= RVENT[0]) & (rg <= RVENT[1])
    A = 0.5*(ps[v][:, z].max(0) - ps[v][:, z].min(0)).max()
    wb = np.sqrt((omax - omin)/djs*A)
    return A, wb, wb/(omax - omin)


def resumen():
    with open(ruta('resumen.txt'), 'w') as fo:
        fo.write(f'tau_1 = {TAU1:.0f} (banda del isócrono desnudo [{OM_MIN}, {OM_MAX}])\n')
        for n in info:
            d = cargar(n); t = d['t']; hk = np.asarray(d['hk'])
            ln = f'{n:4} a0={info[n]["a0"]:<6g} Npc={info[n]["npc"]:<4}'
            if info[n]['sg']:
                A, wb, mu = escalas(n)
                ln += f' |dPhi|_1a vuelta={A:.2e} omega_b={wb:.2e} mu={mu:.2f}'
            fo.write(ln + '\n')
            for a, b in ((0, 1), (1, 3), (3, 10), (10, 20), (20, 31)):
                v = (t >= a*TAU1) & (t <= b*TAU1)
                fo.write(f'   [{a:>2},{b:>2}] tau1: ' + '  '.join(
                    f'|h{k}| {np.abs(hk[v, k]).mean():.3e}' for k in (1, 2, 3, 4)) + '\n')
            Jf = np.asarray(d['J'][-1]); J0 = np.asarray(d['J'][0]); w = np.asarray(d['w'])
            ok = np.isfinite(Jf)
            fo.write(f'   ligadas al final: {w[ok].sum()/w.sum():.4f} de la masa;'
                     f' <|J-J0|> = {np.average(np.abs(Jf-J0)[ok], weights=w[ok]):.3e}\n')
    print(open(ruta('resumen.txt')).read())


def transicion(solo=None):
    """Tabla de la transición: mu, |h1| tardío (potencial del instante y del isócrono desnudo) y
    fracción de la masa con J < 0.06 al final (el núcleo, con el mismo corte que la demo).
    Escribe exe/demo_eta/bulto/transicion.txt."""
    filas = []
    for n in info:
        if (solo and n not in solo) or not os.path.exists(ruta(n, 'bulto.npz')):
            continue
        d = cargar(n); t = d['t']; hk = np.asarray(d['hk']); hkd = np.asarray(d['hk_desnudo'])
        w = np.asarray(d['w']); Jf = np.asarray(d['J'][-1]); ok = np.isfinite(Jf)
        tau = banda(n)[3]
        v = (t >= 10*tau) & (t <= 20*tau); u = t >= 20*tau
        mu = escalas(n)[2] if info[n]['sg'] else 0.0
        nucleo = w[ok & (Jf < 0.06)].sum()/w.sum()
        filas.append((info[n]['a0'] if info[n]['sg'] else 0.0, n, mu, np.abs(hk[v, 1]).mean(), np.abs(hk[u, 1]).mean(),
                      np.abs(hkd[v, 1]).mean(), nucleo, info[n]['npc']))
    filas.sort()
    with open(ruta('transicion.txt'), 'w') as fo:
        fo.write('corrida   a0       Npc   mu     |h1| 10-20 tau1   |h1| desde 20 tau1   ídem desnudo 10-20   masa con J<0.06   tau1\n')
        for a0, n, mu, h1, h1b, h1d, nuc, npc in filas:
            fo.write(f'{n:8}  {a0:<7g}  {npc:<4}  {mu:5.2f}  {h1:14.3e}  {h1b:18.3e}  {h1d:18.3e}  {nuc:12.3f}   {banda(n)[3]:6.0f}\n')
    print(open(ruta('transicion.txt')).read())


def colores(w):
    return np.log10(np.maximum(w/w.max(), 1e-4))


def figuras():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9})
    for n in info:
        d = cargar(n); t = d['t']; c = colores(np.asarray(d['w'])); orden = np.argsort(c)
        cols = [int(np.argmin(np.abs(t - x))) for x in (0, 0.5*TAU1, 3*TAU1, t[-1])]
        fig, ax = plt.subplots(2, 4, figsize=(11, 5.4), constrained_layout=True)
        Jtop = max(0.62, float(np.nanmax(np.asarray(d['J'])[cols])) + 0.01)
        for k, m in enumerate(cols):
            R, P = np.asarray(d['R'][m])[orden], np.asarray(d['P'][m])[orden]
            Q, J = np.asarray(d['Q'][m])[orden], np.asarray(d['J'][m])[orden]
            ax[0, k].scatter(R, P, c=c[orden], s=0.8, cmap='magma_r', vmin=-4, vmax=0, lw=0, rasterized=True)
            sc = ax[1, k].scatter(Q, J, c=c[orden], s=0.8, cmap='magma_r', vmin=-4, vmax=0, lw=0,
                                  rasterized=True)
            ax[0, k].set_title(f'$t={t[m]:.0f}$ ($t/\\tau_1={t[m]/TAU1:.1f}$)')
            ax[0, k].set_xlim(0, 16); ax[0, k].set_ylim(-0.45, 0.45); ax[0, k].set_xlabel('$r$')
            ax[1, k].set_xlim(0, 2*np.pi); ax[1, k].set_ylim(0, Jtop); ax[1, k].set_xlabel('$Q$')
            ax[1, k].set_xticks([0, np.pi, 2*np.pi], ['0', r'$\pi$', r'$2\pi$'])
            if k:
                ax[0, k].set_yticklabels([]); ax[1, k].set_yticklabels([])
        ax[0, 0].set_ylabel('$p_r$'); ax[1, 0].set_ylabel('$J$ (instantaneous potential)')
        fig.colorbar(sc, ax=ax, label=r'$\log_{10}$ of the relative weight', shrink=0.8)
        fig.suptitle(f'{n}: ' + info[n]['desc'].replace('\\_', '_'), fontsize=10)
        fig.savefig(os.path.join(FIGDIR, f'bulto_fase_{n}.jpg'), dpi=130, pil_kwargs={'quality': 88})
        plt.close(fig)
    # Series: |h_1|, |h_2| y la fluctuación del potencial para toda la serie en masa.
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), constrained_layout=True)
    for i, n in enumerate(['B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B4q']):
        d = cargar(n); t = d['t']; hk = np.asarray(d['hk'])
        ls = ':' if n == 'B4q' else ('--' if n == 'B0' else '-')
        col = 'k' if n == 'B0' else ('C3' if n == 'B4q' else f'C{i-1}')
        lab = 'no self-gravity' if n == 'B0' else f'{n}: $a_0={info[n]["a0"]:g}$' + (', $N_{pc}=100$' if n == 'B4q' else '')
        for j, k in enumerate((1, 2)):
            ax[j].semilogy(t/TAU1, envolvente(np.abs(hk[:, k]), t), ls=ls, color=col, lw=1.0, label=lab)
        if info[n]['sg']:
            tt, _, _, _, fl, _ = potencial(n)
            ax[2].semilogy(tt/TAU1, envolvente(fl, tt), ls=ls, color=col, lw=1.0, label=lab)
    ax[0].set_ylabel('$|h_1|$ (envelope)'); ax[1].set_ylabel('$|h_2|$ (envelope)')
    ax[2].set_ylabel(r'$\|\Phi_{\rm self}-\overline{\Phi}_{\rm self}\|/\max|\overline{\Phi}_{\rm self}|$')
    for a in ax:
        a.grid(alpha=0.3); a.set_xlabel(r'$t/\tau_1$'); a.set_xlim(0, TFIN/TAU1)
    ax[0].legend(fontsize=7, loc='lower left')
    fig.savefig(os.path.join(FIGDIR, 'bulto_series.pdf'))
    plt.close(fig)
    resumen()


def _video(n):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter
    plt.rcParams.update({'font.size': 9})
    d = cargar(n); t = d['t']; hk = np.asarray(d['hk']); c = colores(np.asarray(d['w'])); orden = np.argsort(c)
    tt, rg, ps, media, fl, esc = potencial(n)
    z = (rg >= RVENT[0]) & (rg <= RVENT[1])
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    g = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.38, wspace=0.2,
                         left=0.06, right=0.98, top=0.86, bottom=0.08)
    a_rp, a_qj = fig.add_subplot(g[0, 0]), fig.add_subplot(g[0, 1])
    a_pr, a_ts = fig.add_subplot(g[1, 0]), fig.add_subplot(g[1, 1])
    kw = dict(c=c[orden], s=1.5, cmap='magma_r', vmin=-4, vmax=0, lw=0)
    sc1 = a_rp.scatter(np.asarray(d['R'][0])[orden], np.asarray(d['P'][0])[orden], **kw)
    sc2 = a_qj.scatter(np.asarray(d['Q'][0])[orden], np.asarray(d['J'][0])[orden], **kw)
    a_rp.set_xlim(0, 16); a_rp.set_ylim(-0.45, 0.45); a_rp.set_xlabel('$r$'); a_rp.set_ylabel('$p_r$')
    a_rp.set_title(r'phase space $(r,p_r)$')
    Jtop = max(0.62, float(np.nanmax(np.asarray(d['J'][::20]))) + 0.01)
    a_qj.set_xlim(0, 2*np.pi); a_qj.set_ylim(0, Jtop); a_qj.set_xlabel('$Q$'); a_qj.set_ylabel('$J$')
    a_qj.set_xticks([0, np.pi/2, np.pi, 1.5*np.pi, 2*np.pi], ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    a_qj.set_title(r'angle-action $(Q,J)$ in the instantaneous potential')
    lp, = a_pr.plot(rg[z], ps[0][z], lw=1.3, color='C0', label=r'$\Phi_{\rm self}(r,t)$')
    a_pr.plot(rg[z], media[z], color='k', ls='--', lw=0.9, label='late-time mean')
    lo, hi = ps[:, z].min(), ps[:, z].max()
    a_pr.set_ylim(lo - 0.05*(hi - lo) - 1e-12, hi + 0.05*(hi - lo) + 1e-12); a_pr.set_xlim(*RVENT)
    a_pr.set_xlabel('$r$'); a_pr.set_title('self-potential'); a_pr.legend(fontsize=8, loc='lower right')
    a_pr.grid(alpha=0.3)
    for k in (1, 2, 3):
        a_ts.semilogy(t/TAU1, np.abs(hk[:, k]), lw=0.7, label=f'$|h_{k}|$')
    marca = a_ts.axvline(0, color='C3', lw=0.9)
    a_ts.set_xlim(0, t[-1]/TAU1); a_ts.set_xlabel(r'$t/\tau_1$'); a_ts.grid(alpha=0.3)
    a_ts.set_title(r'harmonics of the blob, $h_k=\sum f e^{-ikQ}/\sum f$'); a_ts.legend(fontsize=8, loc='lower left')
    fig.text(0.06, 0.955, f'{n}: blob of 09_autogravedad, ' + info[n]['desc'].replace('\\_', '_'),
             fontsize=12, ha='left', va='center')
    reloj = fig.text(0.98, 0.955, '', fontsize=12, ha='right', va='center')
    fig.text(0.06, 0.918, r'colour: $\log_{10}$ of the particle weight (dark = where the mass is)',
             fontsize=9, ha='left', va='center', color='0.3')
    os.makedirs(VIDDIR, exist_ok=True)
    w = FFMpegWriter(fps=30, bitrate=3000, codec='libx264', extra_args=['-pix_fmt', 'yuv420p'])
    with w.saving(fig, os.path.join(VIDDIR, f'bulto_{n}.mp4'), dpi=100):
        for m in range(len(t)):
            sc1.set_offsets(np.c_[np.asarray(d['R'][m])[orden], np.asarray(d['P'][m])[orden]])
            sc2.set_offsets(np.c_[np.asarray(d['Q'][m])[orden], np.asarray(d['J'][m])[orden]])
            lp.set_ydata(ps[m][z]); marca.set_xdata([t[m]/TAU1]*2)
            reloj.set_text(f'$t={t[m]:.0f}$  ($t/\\tau_1={t[m]/TAU1:.1f}$)')
            w.grab_frame()
    plt.close(fig)
    return n


def videos(solo=None, procesos=4):
    from multiprocessing import Pool
    pend = [n for n in info if (not solo or n in solo)
            and not os.path.exists(os.path.join(VIDDIR, f'bulto_{n}.mp4'))]
    with Pool(procesos) as pool:
        for n in pool.imap_unordered(_video, pend):
            print(f'  video {n}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['preparar', 'correr', 'analizar', 'figuras', 'videos', 'transicion'])
    ap.add_argument('--solo', nargs='*')
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    print(f'=== {a.paso}', flush=True)
    globals()[a.paso](a.solo) if a.paso in ('analizar', 'videos', 'transicion') else globals()[a.paso]()
