"""Demo de la batería eta: ocho corridas PIC con la Maxwelliana de Wilson, sus
referencias con eps = 0, las soluciones lineales, el análisis, las figuras y un
video por corrida.

Pasos (en orden; cada uno reutiliza lo que ya exista):
    python3 demo_eta.py preparar   estados iniciales (equilibrio.py) y .par
    python3 demo_eta.py lineal     soluciones lineales hasta el t_fin de cada caso
    python3 demo_eta.py correr     corridas PIC, una detrás de otra, 4 hilos
    python3 demo_eta.py analizar   h_k y dPhi en el marco del equilibrio, menos eps = 0
    python3 demo_eta.py figuras    figuras para el documento
    python3 demo_eta.py videos     un video por corrida

Las corridas y los análisis se escriben en exe/demo_eta/, las figuras en
docs/demo_eta/figuras/ y los videos en reproducir/videos/demo_eta/. Las corridas usan courant = 1 (dt = 0.05),
dr = 0.1, r in [0, 20], yoshida4, Nrc x Npc = 400 x 25 y una instantánea cada
200 pasos (10 unidades de tiempo).
"""
import os, re, sys, subprocess, time, argparse, warnings, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
warnings.filterwarnings('ignore', category=RuntimeWarning)

RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
EXE = os.path.join(RAIZ, 'exe')
BASE = os.path.join(EXE, 'demo_eta')
PARDIR = os.path.join(RAIZ, 'reproducir', 'corridas', '12_demo_eta')
FIGDIR = os.path.join(RAIZ, 'docs', 'demo_eta', 'figuras')
VIDDIR = os.path.join(RAIZ, 'reproducir', 'videos', 'demo_eta')
PLANTILLA = os.path.join(RAIZ, 'reproducir', 'corridas', '11_landau', 'landau__L_a1e-2_n400_e0.1.par')
JT, W0, NRC, NPC, DT, SALIDA = 0.138, 3.0, 400, 25, 0.05, 200
J1 = 0.05                                       # centro y ancho de B(J)

# Equilibrios (los del bloque L): a0, exponente del borde g, eta, x y gamma lineales.
CASOS = {
    'A1':  dict(a0=0.0065, g=2.0, eta=0.10, tau1=555),
    'A3':  dict(a0=0.042,  g=2.0, eta=0.59, tau1=460),
    'A4':  dict(a0=0.075,  g=2.0, eta=1.00, tau1=395),
    'L5':  dict(a0=0.097,  g=2.0, eta=1.24, tau1=360),
    'L6':  dict(a0=0.170,  g=2.0, eta=2.00, tau1=276),
    'G1a': dict(a0=0.0069, g=1.0, eta=0.10, tau1=557),
}
# Corridas: nombre, caso, eps, t_fin, descripción. Las Z son las referencias eps = 0.
CORRIDAS = [
    ('D2', 'A3', 0.065, 4600, r'$\eta=0.6$, en la banda, lineal ($\nu=0.3$)'),
    ('Z_A3', 'A3', 0.0, 4600, 'referencia de D2'),
    ('D1', 'A1', 1.0, 11100, r'$\eta=0.1$, amortiguamiento rápido, lineal ($\nu\approx0.2$)'),
    ('Z_A1', 'A1', 0.0, 11100, 'referencia de D1'),
    ('D3', 'L5', 0.1, 7200, r'$\eta=1.24$, modo discreto junto al borde, amplitud chica'),
    ('D7', 'L5', 1.0, 7200, r'$\eta=1.24$, el mismo modo a amplitud grande'),
    ('Z_L5', 'L5', 0.0, 7200, 'referencia de D3 y D7'),
    ('D4', 'L6', 0.1, 5520, r'$\eta=2$, modo discreto separado de la banda'),
    ('Z_L6', 'L6', 0.0, 5520, 'referencia de D4'),
    ('D5', 'A4', 0.075, 17200, r'$\eta=1$, en el borde, no lineal ($\nu=3$)'),
    ('D6', 'A4', 0.84, 5200, r'$\eta=1$, en el borde, muy no lineal ($\nu=10$)'),
    ('Z_A4', 'A4', 0.0, 17200, 'referencia de D5 y D6'),
    ('D8', 'G1a', 1.0, 11140, r'$\eta=0.1$, borde abrupto (King, $g=1$)'),
    ('Z_G1a', 'G1a', 0.0, 11140, 'referencia de D8'),
]
DEMO = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8']
REF = {'D1': 'Z_A1', 'D2': 'Z_A3', 'D3': 'Z_L5', 'D7': 'Z_L5', 'D4': 'Z_L6',
       'D5': 'Z_A4', 'D6': 'Z_A4', 'D8': 'Z_G1a'}
RMED = np.linspace(3, 15, 121)
info = {c[0]: dict(caso=c[1], eps=c[2], tfin=c[3], desc=c[4]) for c in CORRIDAS}


def ruta(*p):
    return os.path.join(BASE, *p)


# ------------------------------------------------------------------ preparar
def preparar():
    os.makedirs(ruta('ic'), exist_ok=True); os.makedirs(PARDIR, exist_ok=True)
    plantilla = open(PLANTILLA).read().split('\n')
    for nombre, caso, eps, tfin, desc in CORRIDAS:
        c = CASOS[caso]
        dat = ruta('ic', f'{nombre}.dat')
        if not os.path.exists(dat):
            cmd = [sys.executable, os.path.join(AQUI, 'equilibrio.py'), '--forma', 'maxwell',
                   '--jt', str(JT), '--k', str(c['g']), '--w0', str(W0), '--a0', str(c['a0']),
                   '--eps', str(eps), '--pert', 'suave3', '--nrc', str(NRC), '--npc', str(NPC),
                   '--salida', dat]
            t0 = time.time()
            subprocess.run(cmd, check=True, stdout=open(ruta('ic', f'{nombre}.log'), 'w'),
                           stderr=subprocess.STDOUT)
            print(f'  estado inicial {nombre} ({time.time()-t0:.0f} s)', flush=True)
        valores = {'courant': '1.0', 'Nt': str(int(round(tfin/DT))), 'time_output': str(SALIDA),
                   'spatial_output': str(SALIDA), 'field_output': str(SALIDA),
                   'directory': f'demo_eta/{nombre}', 'a0': str(c['a0']),
                   'checkpointfile': f'demo_eta/ic/{nombre}.dat', 'j1': f'{J1:.6f}',
                   'sj1': f'{J1:.6f}', 'state': 'checkpoint'}
        lineas = [f'# Demo eta, corrida {nombre}: {caso}, eps = {eps}, t_fin = {tfin}.',
                  '# Generado por reproducir/scripts/demo_eta.py a partir de 11_landau.']
        for l in plantilla:
            if l.startswith('#') or '=' not in l:
                continue
            clave = l.split('=')[0].strip()
            lineas.append(f'{clave:<16} = {valores[clave]}' if clave in valores else l)
        open(os.path.join(PARDIR, f'demo__{nombre}.par'), 'w').write('\n'.join(lineas) + '\n')
    print('preparado:', len(CORRIDAS), 'corridas')


# ------------------------------------------------------------------ lineal
def lineal():
    from lineal import resolver
    os.makedirs(ruta('lineal'), exist_ok=True)
    for caso in CASOS:
        sal = ruta('lineal', f'{caso}.npz')
        if os.path.exists(sal):
            continue
        dcor = [n for n, c, e, t, _ in CORRIDAS if c == caso and e > 0][0]
        tmax = max(t for n, c, e, t, _ in CORRIDAS if c == caso)
        t0 = time.time()
        t, h1, h2, rmed, dphi = resolver(ruta('ic', f'{dcor}_equilibrio.npz'), 1600, 32, 0.5,
                                         tmax, verboso=False, j1=J1, sj1=J1)
        np.savez(sal, t=t, h1=h1, h2=h2, r=rmed, dphi=dphi)
        print(f'  lineal {caso} hasta t = {tmax} ({time.time()-t0:.0f} s)', flush=True)


# ------------------------------------------------------------------ correr
def correr():
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    for nombre, *_ in CORRIDAS:
        destino = ruta(nombre)
        if os.path.exists(os.path.join(destino, 'vlasov_output.h5')) and \
           os.path.exists(ruta(f'{nombre}.ok')):
            continue
        par = os.path.join(PARDIR, f'demo__{nombre}.par')
        t0 = time.time()
        r = subprocess.run(['./VP_PIC', par], cwd=EXE, stdout=open(ruta(f'{nombre}.log'), 'w'),
                           stderr=subprocess.STDOUT, env=env)
        estado = 'OK' if r.returncode == 0 else 'FALLO'
        print(f'{estado:5} {nombre} ({time.time()-t0:.0f} s)', flush=True)
        if r.returncode == 0:
            open(ruta(f'{nombre}.ok'), 'w').write('')


# ------------------------------------------------------------------ analizar
def _mapear(tarea):
    """(Q, J) en el mapa fijo del equilibrio y dPhi en la malla, para un grupo de
    instantáneas. Cada instantánea se mapea una sola vez: cuesta ~0.5 s."""
    import h5py
    from aa_numerico import MapaAA, phi_iso
    h5, eqf, claves = tarea
    eq = np.load(eqf)
    mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']))
    f = h5py.File(h5, 'r')
    rg = f['grid']['r'][:]
    fondo = phi_iso(rg) + np.interp(rg, eq['r'], eq['phi_self'])
    sal = []
    for k in claves:
        g = f[k]
        r, p = g['r_part'][:], g['p_part'][:]
        Q, J, _ = mapa(r, p)
        sal.append((g.attrs['time'], r, p, Q, J, g['potential'][:] - fondo))
    f.close()
    return sal


def analizar(solo=None, procesos=4):
    """Para cada corrida: landau.npz (h_k, dPhi, deriva de J; lo mismo que
    landau_analisis.py) y fase.npz (todas las instantáneas en (r,p) y (Q,J))."""
    import h5py
    from multiprocessing import Pool
    from landau_analisis import pesos_prueba
    for nombre, *_ in CORRIDAS:
        d = ruta(nombre)
        if (solo and nombre not in solo) or os.path.exists(os.path.join(d, 'fase.npz')):
            continue
        t0 = time.time()
        eqf = ruta('ic', f'{nombre}_equilibrio.npz')
        h5 = os.path.join(d, 'vlasov_output.h5')
        f = h5py.File(h5, 'r')
        pasos = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
        rg, w = f['grid']['r'][:], f[pasos[0]]['f'][:]
        f.close()
        grupos = [(h5, eqf, pasos[i:i+25]) for i in range(0, len(pasos), 25)]
        with Pool(procesos) as pool:
            res = [x for bloque in pool.map(_mapear, grupos) for x in bloque]
        t = np.array([x[0] for x in res])
        R, P, Q, J, dphi = (np.array([x[i] for x in res]) for i in range(1, 6))
        j1, sj1 = pesos_prueba(d)
        pesos = w*np.exp(-(J - j1)**2/sj1**2)*J**2
        hk = np.stack([np.sum(pesos*np.exp(-1j*m*Q), axis=1) for m in range(5)], axis=1)
        hk = hk/hk[0, 0].real
        eq = np.load(eqf)
        zona = (rg >= 3) & (rg <= 15)
        escala = np.max(np.abs(np.interp(rg, eq['r'], eq['phi_self'])[zona]))
        derivaJ = np.average(np.abs(J - J[0]), weights=w, axis=1)
        np.savez(os.path.join(d, 'landau.npz'), t=t, hk=hk, dphi=dphi, r=rg,
                 derivaJ=derivaJ, escala=escala)
        np.savez(os.path.join(d, 'fase.npz'), t=t, r=R.astype(np.float32), p=P.astype(np.float32),
                 Q=Q.astype(np.float32), J=J.astype(np.float32), w=w)
        print(f'  analizado {nombre}: {len(t)} instantáneas ({time.time()-t0:.0f} s)', flush=True)


def serie(nombre):
    """Señal de la corrida menos su referencia eps = 0, dividida por eps, en los
    radios RMED; y la solución lineal del caso en los mismos tiempos."""
    d = np.load(ruta(nombre, 'landau.npz'))
    t, hk, dphi, rg = d['t'], d['hk'], d['dphi'], d['r']
    eps = info[nombre]['eps']
    z = np.load(ruta(REF[nombre], 'landau.npz'))
    n = min(len(t), len(z['t']))
    t, h1 = t[:n], (hk[:n, 1] - z['hk'][:n, 1])/eps
    dp = (dphi[:n] - z['dphi'][:n])/eps
    dp = np.array([np.interp(RMED, rg, fila) for fila in dp])
    lin = np.load(ruta('lineal', f'{info[nombre]["caso"]}.npz'))
    idx = np.clip(np.searchsorted(lin['t'], t - 1e-9), 0, len(lin['t']) - 1)
    return dict(t=t, h1=h1, dphi=dp, h1_lin=lin['h1'][idx], dphi_lin=lin['dphi'][idx])


def norma(x):
    return np.sqrt(np.mean(x**2, axis=-1))


# ------------------------------------------------------------------ figuras
def colores(nombre, w):
    """Color de cada partícula: su perturbación de peso, F_eq s cos Q0, normalizada
    (en las referencias eps = 0, el ángulo inicial). Las partículas están en el
    orden de los nodos, (i-1)*Npc + j."""
    i, jq = np.divmod(np.arange(NRC*NPC), NPC)
    Jn, Qn = (i + 0.5)*JT/NRC, (jq + 0.5)*2*np.pi/NPC
    eps = info[nombre]['eps']
    if eps == 0:
        return np.cos(Qn), r'ángulo inicial $\cos Q_0$ (una etiqueta)'
    s = (Jn/JT)**1.5
    df = w - w/(1 + eps*s*np.cos(Qn))
    return df/np.abs(df).max(), r'perturbación de peso $\delta f$ (normalizada)'


def limites(fz):
    rr, pp = fz['r'], fz['p']
    return ((rr.min() - 0.2, rr.max() + 0.2), (1.08*pp.min(), 1.08*pp.max()))


def instantes(nombre, t):
    """Índices de t = 0, tau_1/2, 2 tau_1 y t_fin."""
    tau1 = CASOS[info[nombre]['caso']]['tau1']
    return [int(np.argmin(np.abs(t - x))) for x in (0, 0.5*tau1, 2*tau1, t[-1])]


def figura_fase(nombre, plt):
    fz = np.load(ruta(nombre, 'fase.npz'))
    c, etiqueta = colores(nombre, fz['w'])
    xl, yl = limites(fz)
    jtop = max(1.02*JT, float(fz['J'].max()) + 0.002)
    cols = instantes(nombre, fz['t'])
    fig, ax = plt.subplots(2, 4, figsize=(11, 5.4), constrained_layout=True)
    for k, m in enumerate(cols):
        ax[0, k].scatter(fz['r'][m], fz['p'][m], c=c, s=0.8, cmap='RdBu_r', vmin=-1, vmax=1,
                         lw=0, rasterized=True)
        ax[1, k].scatter(fz['Q'][m], fz['J'][m], c=c, s=0.8, cmap='RdBu_r', vmin=-1, vmax=1,
                         lw=0, rasterized=True)
        ax[0, k].set_title(f'$t={fz["t"][m]:.0f}$')
        ax[0, k].set_xlim(*xl); ax[0, k].set_ylim(*yl); ax[0, k].set_xlabel('$r$')
        ax[1, k].set_xlim(0, 2*np.pi); ax[1, k].set_ylim(0, jtop); ax[1, k].set_xlabel('$Q$')
        ax[1, k].axhline(JT, color='0.5', lw=0.5, ls=':')
        ax[1, k].set_xticks([0, np.pi, 2*np.pi], ['0', r'$\pi$', r'$2\pi$'])
        if k:
            ax[0, k].set_yticklabels([]); ax[1, k].set_yticklabels([])
    ax[0, 0].set_ylabel('$p_r$'); ax[1, 0].set_ylabel('$J$')
    i = info[nombre]
    fig.suptitle(f'{nombre}: {i["caso"]} ($\\eta={CASOS[i["caso"]]["eta"]:g}$), '
                 f'$\\varepsilon={i["eps"]:g}$.  Color: {etiqueta}', fontsize=10)
    fig.savefig(os.path.join(FIGDIR, f'fase_{nombre}.jpg'), dpi=130, pil_kwargs={'quality': 88})
    plt.close(fig)


def figura_islas(plt, lista=('D3', 'D5', 'D6', 'D7')):
    """Ampliación del borde de la banda en (Q, J) al final, coloreada por la acción
    inicial: una fila que deja de ser horizontal ha cambiado de acción."""
    fig, ax = plt.subplots(1, len(lista), figsize=(11, 3.1), constrained_layout=True, sharey=True)
    jtop = max(float(np.asarray(np.load(ruta(n, 'fase.npz'), mmap_mode='r')['J'][-1]).max()) for n in lista)
    for a, n in zip(ax, lista):
        fz = np.load(ruta(n, 'fase.npz'), mmap_mode='r')
        m = len(fz['t']) - 1
        J0, Q, J = (np.asarray(fz[k][i]) for k, i in (('J', 0), ('Q', m), ('J', m)))
        sel = J0 > 0.09
        sc = a.scatter(Q[sel], J[sel], c=J0[sel], s=1.2, cmap='viridis', vmin=0.09, vmax=JT, lw=0,
                       rasterized=True)
        a.set_title(f'{n} ($\\varepsilon={info[n]["eps"]:g}$), $t={fz["t"][m]:.0f}$')
        a.set_xlim(0, 2*np.pi); a.set_ylim(0.07, jtop + 0.003); a.set_xlabel('$Q$')
        a.axhline(JT, color='0.4', lw=0.6, ls=':')
        a.set_xticks([0, np.pi, 2*np.pi], ['0', r'$\pi$', r'$2\pi$'])
    ax[0].set_ylabel('$J$')
    fig.colorbar(sc, ax=ax, label='$J$ inicial', shrink=0.9)
    fig.savefig(os.path.join(FIGDIR, 'islas.jpg'), dpi=150, pil_kwargs={'quality': 90})
    plt.close(fig)


def figuras():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    os.makedirs(FIGDIR, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9})
    for nombre in info:
        figura_fase(nombre, plt)
    figura_islas(plt)
    # La transición a amplitud baja: una fila por corrida, PIC contra teoría lineal.
    lista = ['D1', 'D2', 'D3', 'D4']
    fig, ax = plt.subplots(len(lista), 2, figsize=(10, 8.4), constrained_layout=True)
    for k, nombre in enumerate(lista):
        s = serie(nombre)
        tau1 = CASOS[info[nombre]['caso']]['tau1']
        x = s['t']/tau1
        n0, h0 = norma(s['dphi_lin'][0]), np.abs(s['h1_lin'][0])
        ax[k, 0].semilogy(x, envolvente(norma(s['dphi']), s['t'])/n0, lw=1.0, color='C0', label='PIC')
        ax[k, 0].semilogy(x, envolvente(norma(s['dphi_lin']), s['t'])/n0, color='k', ls='--', lw=0.9,
                          label='lineal')
        ax[k, 1].semilogy(x, np.abs(s['h1'])/h0, lw=0.7, color='C0', label='PIC')
        ax[k, 1].semilogy(x, np.abs(s['h1_lin'])/h0, color='k', ls='--', lw=0.7, label='lineal')
        ax[k, 0].set_ylabel(r'envolvente de $\|\delta\Phi\|$')
        ax[k, 1].set_ylabel(r'$|h_1|/|h_1(0)|$')
        ax[k, 0].set_title(f'{nombre}: ' + info[nombre]['desc'], loc='left')
        for a in ax[k]:
            a.grid(alpha=0.3); a.set_xlim(0, x[-1])
    ax[0, 1].legend(fontsize=8, loc='upper right')
    ax[-1, 0].set_xlabel(r'$t/\tau_1$'); ax[-1, 1].set_xlabel(r'$t/\tau_1$')
    fig.savefig(os.path.join(FIGDIR, 'serie_transicion.pdf'))
    plt.close(fig)
    # Pares: el mismo equilibrio a dos amplitudes, o dos bordes con el mismo eta.
    pares = {'nolineal': ['D5', 'D6'], 'discreto': ['D3', 'D7'], 'borde': ['D1', 'D8']}
    for clave, lista in pares.items():
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.3), constrained_layout=True)
        for k, nombre in enumerate(lista):
            s = serie(nombre)
            tau1 = CASOS[info[nombre]['caso']]['tau1']
            x, col = s['t']/tau1, f'C{k}'
            n0, h0 = norma(s['dphi_lin'][0]), np.abs(s['h1_lin'][0])
            e = f'{nombre} ($\\varepsilon={info[nombre]["eps"]:g}$)' if clave != 'borde' else \
                f'{nombre} ($g={CASOS[info[nombre]["caso"]]["g"]:g}$)'
            ax[0].semilogy(x, envolvente(norma(s['dphi']), s['t'])/n0, lw=1.0, color=col, label=e)
            ax[1].semilogy(x, np.abs(s['h1'])/h0, lw=0.7, color=col, label=e)
            lin_igual = clave != 'borde'
            if not lin_igual or k == 0:
                c = 'k' if lin_igual else col
                ax[0].semilogy(x, envolvente(norma(s['dphi_lin']), s['t'])/n0, color=c, ls='--', lw=0.9,
                               label='lineal' if lin_igual else f'lineal {info[nombre]["caso"]}')
                ax[1].semilogy(x, np.abs(s['h1_lin'])/h0, color=c, ls='--', lw=0.8)
            eP = envolvente(norma(s['dphi']), s['t']); eL = envolvente(norma(s['dphi_lin']), s['t'])
            ax[2].semilogy(x, eP/eL, lw=1.0, color=col, label=nombre)
            xmax = max(xmax, x[-1]) if k else x[-1]
        ax[0].set_ylabel(r'envolvente de $\|\delta\Phi\|$'); ax[1].set_ylabel(r'$|h_1|/|h_1(0)|$')
        ax[2].set_ylabel('envolvente PIC / lineal'); ax[2].axhline(1, color='k', lw=0.6)
        ax[0].legend(fontsize=7, loc='lower left')
        for a in ax:
            a.grid(alpha=0.3); a.set_xlim(0, xmax); a.set_xlabel(r'$t/\tau_1$')
        fig.savefig(os.path.join(FIGDIR, f'serie_{clave}.pdf'))
        plt.close(fig)
    resumen()


def envolvente(x, t, ancho=50.0):
    """Máximo de x en |t' - t| <= ancho: una vuelta radial (~90) cabe en la ventana."""
    n = int(round(ancho/(t[1] - t[0])))
    return np.array([x[max(0, k - n):k + n + 1].max() for k in range(len(x))])


def resumen():
    """Tabla de la demo: PIC contra teoría lineal por ventanas de tau_1.

    R      = envolvente de ||dPhi|| PIC / envolvente lineal;
    est    = parte estática (media temporal) del residuo PIC - lineal;
    osc    = parte oscilante del residuo;
    t_nl   = primer t en que la envolvente del residuo supera el 10 % de la lineal;
    polo   = matrix pencil (M = 3) de h_1, PIC y lineal."""
    from landau_cola import ajustar
    filas = []
    with open(ruta('resumen.txt'), 'w') as fo:
        for nombre in DEMO:
            s = serie(nombre)
            t, p, L = s['t'], s['dphi'], s['dphi_lin']
            tau1 = CASOS[info[nombre]['caso']]['tau1']
            n0 = norma(L[0])
            eP, eL = envolvente(norma(p), t), envolvente(norma(L), t)
            eR = envolvente(norma(p - L), t)
            nl = eR > 0.1*eL
            tnl = t[np.argmax(nl)] if nl.any() else np.nan
            fo.write(f'{nombre}  t_nl = {tnl:.0f} ({tnl/tau1:.1f} tau1)   t_fin = {t[-1]:.0f} '
                     f'({t[-1]/tau1:.1f} tau1)\n')
            for a, b in ((0, 1), (1, 3), (3, 5), (5, 10), (10, 20), (20, 30), (30, 44)):
                v = (t >= a*tau1) & (t <= min(b*tau1, t[-1]))
                if v.sum() < 20 or a*tau1 >= t[-1] - 100:
                    continue
                r = p[v] - L[v]
                fo.write(f'   [{a:>2},{b:>2}] tau1: PIC {eP[v].max()/n0:.3e}  lin {eL[v].max()/n0:.3e}'
                         f'  R = {eP[v].max()/eL[v].max():.3f}   est = {norma(r.mean(0))/n0:.1e}'
                         f'  osc = {norma(r - r.mean(0)).mean()/n0:.1e}\n')
            for a, b in ((1, 3), (3, 6), (5, 9), (10, 20), (20, 44)):
                if b*tau1 > t[-1] + 50:
                    continue
                e = [ajustar(t, x, a*tau1, b*tau1)['pencil M=3'] for x in (s['h1'], s['h1_lin'])]
                fo.write(f'   polo h1 [{a},{b}] tau1: PIC w = {e[0][0]:.5f} g = {e[0][1]:+.2e}'
                         f' | lin w = {e[1][0]:.5f} g = {e[1][1]:+.2e}\n')
            filas.append((nombre, tnl))
    print(open(ruta('resumen.txt')).read())
    return filas


# ------------------------------------------------------------------ videos
def _video(nombre):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter
    plt.rcParams.update({'font.size': 9})
    i = info[nombre]; caso, eps = i['caso'], i['eps']
    sal = os.path.join(VIDDIR, f'{nombre}.mp4')
    os.makedirs(VIDDIR, exist_ok=True)
    fz = np.load(ruta(nombre, 'fase.npz'))
    c, etiqueta = colores(nombre, fz['w'])
    la = np.load(ruta(nombre, 'landau.npz'))
    if eps > 0:
        s = serie(nombre); tt, dp, dpl = s['t'], s['dphi'], s['dphi_lin']
        etq = r'$[\delta\Phi(\varepsilon)-\delta\Phi(0)]/\varepsilon$'
    else:
        tt, dpl = la['t'], None
        dp = np.array([np.interp(RMED, la['r'], fila) for fila in la['dphi']])
        etq = r'$\delta\Phi$ de la corrida sin perturbación'
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    g = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.38, wspace=0.2,
                         left=0.06, right=0.98, top=0.86, bottom=0.08)
    a_rp, a_qj = fig.add_subplot(g[0, 0]), fig.add_subplot(g[0, 1])
    a_pr, a_ts = fig.add_subplot(g[1, 0]), fig.add_subplot(g[1, 1])
    sc1 = a_rp.scatter(fz['r'][0], fz['p'][0], c=c, s=1.5, cmap='RdBu_r', vmin=-1, vmax=1, lw=0)
    sc2 = a_qj.scatter(fz['Q'][0], fz['J'][0], c=fz['J'][0], s=1.5, cmap='viridis', vmin=0, vmax=JT,
                       lw=0)
    xl, yl = limites(fz)
    a_rp.set_xlim(*xl); a_rp.set_ylim(*yl)
    a_rp.set_xlabel('$r$'); a_rp.set_ylabel('$p_r$'); a_rp.set_title(r'espacio fase $(r,p_r)$')
    a_qj.set_xlim(0, 2*np.pi); a_qj.set_ylim(0, max(1.02*JT, float(fz['J'].max()) + 0.002))
    a_qj.axhline(JT, color='0.5', lw=0.6, ls=':')
    a_qj.set_xlabel('$Q$'); a_qj.set_ylabel('$J$')
    a_qj.set_xticks([0, np.pi/2, np.pi, 1.5*np.pi, 2*np.pi], ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    a_qj.set_title(r'ángulo-acción $(Q,J)$ del equilibrio; color: $J$ inicial')
    l_pic, = a_pr.plot(RMED, dp[0], lw=1.3, color='C0', label='PIC')
    pico = np.abs(dp).max(axis=1)
    if dpl is not None:
        l_lin, = a_pr.plot(RMED, dpl[0], color='k', ls='--', lw=1.0, label='teoría lineal')
        pico = np.maximum(pico, np.abs(dpl).max(axis=1))
    # El eje vertical sigue a la envolvente, suavizada sobre ~1.5 órbitas (+-8 instantáneas).
    env = np.array([pico[max(0, n - 8):n + 9].max() for n in range(len(pico))])
    a_pr.set_ylim(-1.2*env[0], 1.2*env[0]); a_pr.set_xlim(RMED[0], RMED[-1]); a_pr.set_xlabel('$r$')
    a_pr.set_title(etq + ' (eje reescalado)'); a_pr.legend(fontsize=8, loc='upper right')
    a_pr.grid(alpha=0.3)
    nrm = norma(dp)
    if dpl is not None:
        a_ts.semilogy(tt, nrm/nrm[0], lw=0.8, color='C0', label='PIC')
        a_ts.semilogy(tt, norma(dpl)/norma(dpl[0]), color='k', ls='--', lw=0.8, label='teoría lineal')
        a_ts.set_ylabel(r'$\|\delta\Phi\|/\|\delta\Phi(0)\|$')
    else:
        a_ts.semilogy(tt, nrm/float(la['escala']), lw=0.8, color='C0', label='PIC')
        a_ts.set_ylabel(r'$\|\delta\Phi\|/\max|\Phi_{\rm self,eq}|$')
    marca = a_ts.axvline(0, color='C3', lw=0.9)
    a_ts.set_xlim(0, tt[-1]); a_ts.set_xlabel('$t$')
    a_ts.set_title(r'amplitud de $\delta\Phi$ (rms en $3\leq r\leq15$)'); a_ts.grid(alpha=0.3)
    a_ts.legend(fontsize=8, loc='upper right')
    desc = re.sub(r'^\$\\eta=[0-9.]+\$,? ?', '', i['desc'])
    fig.text(0.06, 0.955, f'{nombre}: equilibrio {caso} ($\\eta={CASOS[caso]["eta"]:g}$), '
             f'$\\varepsilon={eps:g}$.  {desc}', fontsize=12, ha='left', va='center')
    reloj = fig.text(0.98, 0.955, '', fontsize=12, ha='right', va='center')
    fig.text(0.06, 0.918, f'color en $(r,p_r)$: {etiqueta}; en $(Q,J)$: acción inicial (una fila que se '
             'deforma cambió de acción; la línea punteada es el borde $J_t$)', fontsize=9, ha='left',
             va='center', color='0.3')
    tau1 = CASOS[caso]['tau1']
    w = FFMpegWriter(fps=30, bitrate=3000, codec='libx264', extra_args=['-pix_fmt', 'yuv420p'])
    with w.saving(fig, sal, dpi=100):
        for m in range(len(fz['t'])):
            tm = fz['t'][m]
            sc1.set_offsets(np.c_[fz['r'][m], fz['p'][m]])
            sc2.set_offsets(np.c_[fz['Q'][m], fz['J'][m]])
            n = min(int(np.searchsorted(tt, tm - 1e-9)), len(tt) - 1)
            l_pic.set_ydata(dp[n])
            a_pr.set_ylim(-1.2*env[n], 1.2*env[n])
            if dpl is not None:
                l_lin.set_ydata(dpl[n])
            marca.set_xdata([tm, tm])
            reloj.set_text(f'$t={tm:.0f}$  ($t/\\tau_1={tm/tau1:.1f}$)')
            w.grab_frame()
    plt.close(fig)
    return nombre


def videos(solo=None, procesos=4):
    from multiprocessing import Pool
    os.makedirs(VIDDIR, exist_ok=True)
    pend = [n for n in info if (not solo or n in solo)
            and not os.path.exists(os.path.join(VIDDIR, f'{n}.mp4'))]
    with Pool(procesos) as pool:
        for n in pool.imap_unordered(_video, pend):
            print(f'  video {n}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['preparar', 'lineal', 'correr', 'analizar', 'figuras', 'videos', 'todo'])
    ap.add_argument('--solo', nargs='*')
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    pasos = ['preparar', 'lineal', 'correr', 'analizar', 'figuras', 'videos'] if a.paso == 'todo' else [a.paso]
    for p in pasos:
        print(f'=== {p}', flush=True)
        globals()[p](a.solo) if p in ('analizar', 'videos') else globals()[p]()
