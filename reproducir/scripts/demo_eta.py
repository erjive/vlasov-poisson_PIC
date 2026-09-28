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
    'M5':  dict(a0=0.5,    g=2.0, eta=5.11, tau1=125),
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
    ('D9', 'M5', 0.1, 5000, r'$\eta=5.1$, masa comparable a la del fondo, amplitud chica'),
    ('D10', 'M5', 1.0, 5000, r'$\eta=5.1$, masa comparable a la del fondo, amplitud grande'),
    ('Z_M5', 'M5', 0.0, 5000, 'referencia de D9 y D10'),
]
# Controles pedidos por la revisión (auditoria_deepseekpro): barrido en eps del modo
# de L5, y convergencia en N (400x25 -> 800x50) y en dt (courant 1 -> 0.5) de D3, D5 y
# D6, cada una con su referencia eps = 0 a la misma resolución. No tienen video.
EXTRA = [
    ('D3e003', 'L5', 0.003, 7200, r'D3 con $\varepsilon=0.003$'),
    ('D3e01', 'L5', 0.01, 7200, r'D3 con $\varepsilon=0.01$'),
    ('D3e03', 'L5', 0.03, 7200, r'D3 con $\varepsilon=0.03$'),
    ('D3N', 'L5', 0.1, 7200, r'D3 con $N=800\times50$'),
    ('Z_L5N', 'L5', 0.0, 7200, 'referencia de D3N'),
    ('D3dt', 'L5', 0.1, 7200, r'D3 con $\Delta t/2$'),
    ('Z_L5dt', 'L5', 0.0, 7200, 'referencia de D3dt'),
    ('D6N', 'A4', 0.84, 5200, r'D6 con $N=800\times50$'),
    ('D6dt', 'A4', 0.84, 5200, r'D6 con $\Delta t/2$'),
    ('D5dt', 'A4', 0.075, 17200, r'D5 con $\Delta t/2$'),
    ('Z_A4dt', 'A4', 0.0, 17200, 'referencia de D5dt y D6dt'),
    ('D5N', 'A4', 0.075, 17200, r'D5 con $N=800\times50$'),
    ('Z_A4N', 'A4', 0.0, 17200, 'referencia de D5N y D6N'),
    ('D6L', 'A4', 0.84, 17200, r'D6 hasta $t=17200$, para varios periodos de rebote'),
]
NUM = {n: ((800, 50, 1.0) if n.endswith('N') else (400, 25, 0.5) if n.endswith('dt') else (400, 25, 1.0))
       for n, *_ in EXTRA}
CORRIDAS = CORRIDAS + EXTRA
DEMO = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10']
REF = {'D1': 'Z_A1', 'D2': 'Z_A3', 'D3': 'Z_L5', 'D7': 'Z_L5', 'D4': 'Z_L6',
       'D5': 'Z_A4', 'D6': 'Z_A4', 'D8': 'Z_G1a', 'D9': 'Z_M5', 'D10': 'Z_M5',
       'D3e003': 'Z_L5', 'D3e01': 'Z_L5', 'D3e03': 'Z_L5', 'D3N': 'Z_L5N', 'D3dt': 'Z_L5dt',
       'D5N': 'Z_A4N', 'D6N': 'Z_A4N', 'D5dt': 'Z_A4dt', 'D6dt': 'Z_A4dt', 'D6L': 'Z_A4'}
RMED = np.linspace(3, 15, 121)
info = {c[0]: dict(caso=c[1], eps=c[2], tfin=c[3], desc=c[4],
                   nrc=NUM.get(c[0], (NRC, NPC, 1.0))[0], npc=NUM.get(c[0], (NRC, NPC, 1.0))[1],
                   courant=NUM.get(c[0], (NRC, NPC, 1.0))[2], extra=c[0] in NUM) for c in CORRIDAS}


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
                   '--eps', str(eps), '--pert', 'suave3', '--nrc', str(info[nombre]['nrc']),
                   '--npc', str(info[nombre]['npc']),
                   '--salida', dat]
            t0 = time.time()
            subprocess.run(cmd, check=True, stdout=open(ruta('ic', f'{nombre}.log'), 'w'),
                           stderr=subprocess.STDOUT)
            print(f'  estado inicial {nombre} ({time.time()-t0:.0f} s)', flush=True)
        cou = info[nombre]['courant']; sal = int(round(SALIDA/cou))   # una instantánea cada 10
        valores = {'courant': str(cou), 'Nt': str(int(round(tfin/(DT*cou)))), 'time_output': str(sal),
                   'spatial_output': str(sal), 'field_output': str(sal),
                   'Nrc': str(info[nombre]['nrc']), 'Npc': str(info[nombre]['npc']),
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
def metadatos(nombre, par, salida):
    """Versión exacta de lo que produjo la corrida: sha256 del ejecutable y del .par,
    commit del repositorio (y si src/ tenía cambios sin commit) y fecha."""
    import hashlib, datetime
    sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
    git = lambda *a: subprocess.run(['git', *a], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
    ahora = datetime.datetime.now().isoformat(timespec='minutes')
    src = git('log', '-1', '--format=%h', '--', 'src/') + ('+cambios' if git('status', '--porcelain', 'src/') else '')
    with open(os.path.join(PARDIR, 'METADATOS.txt'), 'a') as fo:           # copia versionada
        fo.write(f'{nombre:8} {ahora}  par {sha(par)[:16]}  VP_PIC {sha(os.path.join(EXE, "VP_PIC"))[:16]}'
                 f'  src {src}  HEAD {git("rev-parse", "--short", "HEAD")}\n')
    open(salida, 'w').write(
        f'corrida    {nombre}\nfecha      {datetime.datetime.now().isoformat(timespec="seconds")}\n'
        f'VP_PIC     {sha(os.path.join(EXE, "VP_PIC"))}\npar        {sha(par)}  {os.path.relpath(par, RAIZ)}\n'
        f'HEAD       {git("rev-parse", "--short", "HEAD")}\n'
        f'src/       {git("log", "-1", "--format=%h", "--", "src/")}'
        f'{" (con cambios sin commit)" if git("status", "--porcelain", "src/") else ""}\n')


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
            metadatos(nombre, par, ruta(f'{nombre}.meta'))


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


def cercano(tref, t):
    """Índice de la muestra de tref más cercana a cada t. (Con searchsorted, el redondeo
    acumulado en tref hacía tomar la muestra siguiente, 2 unidades de tiempo después,
    en ~60 % de las instantáneas.)"""
    k = np.clip(np.searchsorted(tref, t), 1, len(tref) - 1)
    return k - ((t - tref[k - 1]) < (tref[k] - t))


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
    idx = cercano(lin['t'], t)
    return dict(t=t, h1=h1, dphi=dp, h1_lin=lin['h1'][idx], dphi_lin=lin['dphi'][idx])


def norma(x):
    return np.sqrt(np.mean(x**2, axis=-1))


# ------------------------------------------------------------------ figuras
def colores(nombre, w):
    """Color de cada partícula: su perturbación de peso, F_eq s cos Q0, normalizada
    (en las referencias eps = 0, el ángulo inicial). Las partículas están en el
    orden de los nodos, (i-1)*Npc + j."""
    nrc, npc = info[nombre]['nrc'], info[nombre]['npc']
    i, jq = np.divmod(np.arange(nrc*npc), npc)
    Jn, Qn = (i + 0.5)*JT/nrc, (jq + 0.5)*2*np.pi/npc
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
        if not info[nombre]['extra']:
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
    pares = {'nolineal': ['D5', 'D6'], 'discreto': ['D3', 'D7'], 'borde': ['D1', 'D8'],
             'masa': ['D9', 'D10']}
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


def polo_con_error(t, h, lo, hi):
    """Polo (omega, gamma) de h en [lo, hi] con matrix pencil M = 3, y su
    incertidumbre: la mayor desviación entre M = 2, 3 y tres recortes de la
    ventana (completa, sin el primer 10 %, sin el último 10 %). Es una cota del
    error del ajuste, no del error numérico de la corrida."""
    from landau_cola import ajustar
    w0, g0 = ajustar(t, h, lo, hi)['pencil M=3']
    d = 0.1*(hi - lo)
    est = [ajustar(t, h, a, b)[f'pencil M={M}'] for M in (2, 3)
           for a, b in ((lo, hi), (lo + d, hi), (lo, hi - d))]
    return w0, max(abs(w - w0) for w, _ in est), g0, max(abs(g - g0) for _, g in est)


def suavizar(x, sigma=0.5):
    """Suavizado gaussiano en r (sobre RMED) de perfiles x[..., r]."""
    K = np.exp(-0.5*((RMED[:, None] - RMED[None, :])/sigma)**2); K /= K.sum(1, keepdims=True)
    return x @ K.T


def metricas(n, ventanas):
    """Envolventes, polo de h1 con error, partes lisa y fina y piso de ruido de una corrida."""
    s = serie(n); t, p, L = s['t'], s['dphi'], s['dphi_lin']
    tau = CASOS[info[n]['caso']]['tau1']; n0 = norma(L[0])
    eP, eL = envolvente(norma(p), t)/n0, envolvente(norma(L), t)/n0
    k5 = np.argmin(abs(t - 5*tau)); perf = L[k5]/np.linalg.norm(L[k5])
    z = np.load(ruta(REF[n], 'landau.npz')); zz = np.array([np.interp(RMED, z['r'], f) for f in z['dphi'][:len(t)]])
    piso = norma(zz - zz.mean(0)).mean()/info[n]['eps']/n0          # fluctuación de la referencia / eps
    out = dict(t=t, eP=eP, eL=eL, tau=tau, piso=piso)
    for a, b in ventanas:
        v = (t >= a*tau) & (t <= min(b*tau, t[-1]))
        e = polo_con_error(t, s['h1'], a*tau, min(b*tau, t[-1]))
        el = polo_con_error(t, s['h1_lin'], a*tau, min(b*tau, t[-1]))
        ps = suavizar(p[v])
        out[(a, b)] = dict(R=eP[v].max()/eL[v].max(), w=e[0], dw=e[1], g=e[2], dg=e[3], wl=el[0], gl=el[2], dgl=el[3],
                           liso=norma(ps).max()/norma(suavizar(L[v])).max(), fino=norma(p[v] - ps).max()/n0,
                           proy=np.abs(p[v] @ perf).max()/np.abs(L[v] @ perf).max())
    return out


def perdida(m, a=5, b=20):
    """Pérdida relativa de la envolvente PIC entre a y b tau_1, descontada la de la lineal."""
    t, tau = m['t'], m['tau']
    ka, kb = np.argmin(abs(t - a*tau)), np.argmin(abs(t - b*tau))
    return 1 - (m['eP'][kb]/m['eP'][ka])/(m['eL'][kb]/m['eL'][ka])


def borde(n, jt=JT):
    """Acción máxima al final y cambio máximo de acción de las filas del borde."""
    fz = np.load(ruta(n, 'fase.npz'), mmap_mode='r')
    J0, Jf = np.asarray(fz['J'][0]), np.asarray(fz['J'][-1])
    return float(Jf.max()), float(np.abs(Jf - J0)[J0 > 0.12].max())


def controles():
    """Barrido en eps de D3, convergencia en N y dt de D3, D5 y D6, y la isla de D6L.
    Números en exe/demo_eta/controles.txt; figuras controles_*.pdf."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9})
    fo = open(ruta('controles.txt'), 'w')
    def w(x=''):
        fo.write(x + '\n'); print(x)
    # --- barrido en eps
    w('== Barrido en eps de D3 (L5, eta = 1.24): omega_b(eps=1) = 4.76e-3, distancia al borde 2.5e-4')
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4), constrained_layout=True)
    filas = []
    for n in ('D3e003', 'D3e01', 'D3e03', 'D3'):
        m = metricas(n, [(10, 20)]); e = info[n]['eps']; q = m[(10, 20)]
        Jm, dJ = borde(n); cr = 4.76e-3*np.sqrt(e)/2.5e-4
        filas.append((e, perdida(m), q['g'], q['dg'], m['piso'], cr))
        w(f'{n:7} eps={e:<6g} wb/dist={cr:5.2f}  perdida 5->20 tau1 = {perdida(m):+.3f}  gamma[10,20] = {q["g"]:+.2e} +- {q["dg"]:.0e}'
          f' (lin {q["gl"]:+.1e})  piso/eps = {m["piso"]:.1e}  Jmax = {Jm:.4f}  dJ borde = {dJ:.4f}')
        ax[0].semilogy(m['t']/m['tau'], m['eP'], lw=1.0, label=f'$\\varepsilon={e:g}$')
    ax[0].semilogy(m['t']/m['tau'], m['eL'], 'k--', lw=0.9, label='lineal')
    ax[0].set_xlabel(r'$t/\tau_1$'); ax[0].set_ylabel(r'envolvente de $\|\delta\Phi_\varepsilon\|$'); ax[0].legend(fontsize=7)
    f = np.array(filas)
    ax[1].errorbar(f[:, 0], 1e5*f[:, 2], yerr=1e5*f[:, 3], fmt='o-', color='C3', capsize=3)
    ax[1].axhline(0, color='k', lw=0.6); ax[1].set_xscale('log')
    ax[1].set_xlabel(r'$\varepsilon$'); ax[1].set_ylabel(r'$\gamma$ de $h_1$ en $10$--$20\,\tau_1$ [$10^{-5}$]')
    sec = ax[1].secondary_xaxis('top', functions=(lambda x: 4.76e-3*np.sqrt(np.maximum(x, 1e-12))/2.5e-4,
                                                  lambda y: (np.maximum(y, 1e-12)*2.5e-4/4.76e-3)**2))
    sec.set_xscale('log'); sec.set_xticks([1, 2, 3, 6], ['1', '2', '3', '6']); sec.minorticks_off()
    sec.set_xlabel(r'$\omega_b/(\Omega_{\min}-\omega)$')
    for a in ax: a.grid(alpha=0.3)
    fig.savefig(os.path.join(FIGDIR, 'controles_eps.pdf')); plt.close(fig)
    # --- convergencia
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), constrained_layout=True)
    for k, (base, ven) in enumerate((('D3', [(10, 20)]), ('D5', [(10, 20), (20, 30), (30, 43.5)]), ('D6', [(5, 10), (10, 13)]))):
        w(f'== Convergencia de {base}')
        for n, et in ((base, 'base 400x25, dt'), (base + 'N', 'N x4 (800x50)'), (base + 'dt', 'dt/2')):
            m = metricas(n, ven)
            txt = '  '.join(f'[{a},{b}] R={m[(a, b)]["R"]:.3f} liso={m[(a, b)]["liso"]:.2f} proy={m[(a, b)]["proy"]:.2f}'
                            f' fino={m[(a, b)]["fino"]:.1e} w={m[(a, b)]["w"]:.5f}+-{m[(a, b)]["dw"]:.0e}'
                            f' g={m[(a, b)]["g"]:+.1e}+-{m[(a, b)]["dg"]:.0e}' for a, b in ven)
            extra = f'  perdida 5->20 = {perdida(m):+.3f}' if base == 'D3' else ''
            Jm, dJ = borde(n)
            w(f'  {n:6} {et:16} {txt}{extra}  Jmax={Jm:.4f} dJborde={dJ:.4f}')
            ax[k].semilogy(m['t']/m['tau'], m['eP'], lw=1.0, label=et)
        ax[k].semilogy(m['t']/m['tau'], m['eL'], 'k--', lw=0.9, label='lineal')
        ax[k].set_title(base); ax[k].set_xlabel(r'$t/\tau_1$'); ax[k].grid(alpha=0.3)
    ax[0].set_ylabel(r'envolvente de $\|\delta\Phi_\varepsilon\|$'); ax[0].legend(fontsize=7)
    fig.savefig(os.path.join(FIGDIR, 'controles_convergencia.pdf')); plt.close(fig)
    fo.close()


def isla(n='D6L', ventana=(10, 43)):
    """La isla de D6L contra el péndulo. Amplitud A(t) = |dPhi_1(J_r, t)| a lo largo de la
    órbita resonante, |Omega'(J_r)| del equilibrio, y las partículas: atrapadas si su
    fase psi = Q - omega t libra (recorre menos de 2 pi) en la ventana; frecuencia de
    rebote medida con la FFT de J(t) de las atrapadas."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from aa_numerico import MapaAA
    from equilibrio import invertir
    fo = open(ruta('isla.txt'), 'w')
    def w(x=''):
        fo.write(x + '\n'); print(x)
    s = serie(n); t = s['t']; tau = CASOS[info[n]['caso']]['tau1']
    om, dom, _, _ = polo_con_error(t, s['h1'], 5*tau, 15*tau)
    eq = np.load(ruta('ic', f'{n}_equilibrio.npz')); Jt, Et = eq['J_t'], eq['E_t']
    Om = np.gradient(Et, Jt); dOm = np.gradient(Om, Jt)
    o = np.argsort(Om); Jr = float(np.interp(om, Om[o], Jt[o])); dOr = abs(float(np.interp(Jr, Jt, dOm)))
    mapa = MapaAA(eq['r'], eq['phi_self'], L=2.0)
    Qo = (np.arange(64) + 0.5)*2*np.pi/64
    rr, _ = invertir(mapa, Et, Jt, Qo, np.full(64, Jr))
    la = np.load(ruta(n, 'landau.npz')); z = np.load(ruta(REF[n], 'landau.npz'))
    k = min(len(la['t']), len(z['t'])); dp = la['dphi'][:k] - z['dphi'][:k]
    A = np.array([abs(np.mean(np.interp(rr, la['r'], f)*np.exp(-1j*Qo))) for f in dp])
    fz = np.load(ruta(n, 'fase.npz'), mmap_mode='r'); tf = fz['t']
    v = (tf >= ventana[0]*tau) & (tf <= ventana[1]*tau)
    Q, J = np.asarray(fz['Q'])[v].astype(float), np.asarray(fz['J'])[v].astype(float)
    wgt = np.asarray(fz['w'])
    psi = np.unwrap(Q - om*tf[v][:, None], axis=0)
    atr = (psi.max(0) - psi.min(0)) < 2*np.pi
    Av = A[(la['t'][:k] >= ventana[0]*tau) & (la['t'][:k] <= ventana[1]*tau)].mean()
    w(f'{n}: omega del modo = {om:.5f} +- {dom:.0e}; J_r = {Jr:.4f}; |Omega\'(J_r)| = {dOr:.3f}')
    w(f'  A = |dPhi_1(J_r)|: t=0 {A[0]:.2e}; media en {ventana} tau1 {Av:.2e}; min {A.min():.2e} max {A[len(A)//4:].max():.2e}')
    w(f'  péndulo con A media: semiancho {2*np.sqrt(Av/dOr):.4f}, omega_b {np.sqrt(Av*dOr):.2e} (periodo {2*np.pi/np.sqrt(Av*dOr):.0f})')
    w('  oscilación de la amplitud (envolvente lenta, lóbulos principales): mínimos en 6.6 y 30.2 tau1,'
      f' máximos en 18.0 y 40.1 tau1: periodo ~{0.5*(23.6 + 22.1)*tau:.0f} ({0.5*(23.6 + 22.1):.1f} tau1)')
    w(f'  atrapadas (psi libra en la ventana): {wgt[atr].sum()/wgt.sum():.3f} de la masa, {atr.sum()} partículas;'
      f' J0 de ellas en [{np.asarray(fz["J"][0])[atr].min() if atr.any() else 0:.3f}, {np.asarray(fz["J"][0])[atr].max() if atr.any() else 0:.3f}]')
    if atr.any():
        semi = 0.5*(J[:, atr].max(0) - J[:, atr].min(0))
        dt = tf[1] - tf[0]
        Jc = J[:, atr] - J[:, atr].mean(0)
        esp = np.abs(np.fft.rfft(Jc, axis=0))**2; fr = 2*np.pi*np.fft.rfftfreq(Jc.shape[0], dt)
        wb = fr[1:][np.argmax(esp[1:], axis=0)]
        w(f'  semiancho medido: máx {semi.max():.4f}, mediana {np.median(semi):.4f}')
        w(f'  omega_b medida (FFT de J de las atrapadas): mediana {np.median(wb):.2e}, cuartiles [{np.percentile(wb, 25):.2e}, {np.percentile(wb, 75):.2e}]')
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), constrained_layout=True)
    n0 = norma(s['dphi_lin'][0])
    ax[0].semilogy(t/tau, envolvente(norma(s['dphi']), t)/n0, lw=1.0, label=n)
    ax[0].semilogy(t/tau, envolvente(norma(s['dphi_lin']), t)/n0, 'k--', lw=0.9, label='lineal')
    ax[0].set_xlabel(r'$t/\tau_1$'); ax[0].set_ylabel(r'envolvente de $\|\delta\Phi_\varepsilon\|$'); ax[0].legend(fontsize=7)
    tt = la['t'][:k]
    # Envolvente lenta (ventana de +-300, unas tres vueltas) y periodo de la oscilación de amplitud
    el = envolvente(norma(s['dphi']), t, ancho=300.0)/n0
    Tb = 2*np.pi/np.sqrt(Av*dOr)
    ax[1].plot(t/tau, el, color='C0', lw=1.0, label=r'envolvente lenta de $\|\delta\Phi_\varepsilon\|$')
    for x0 in (6.6, 30.2):
        ax[1].axvline(x0, color='C3', lw=0.6, ls=':')
    ax[1].annotate('', xy=(6.6 + Tb/tau, 0.12), xytext=(6.6, 0.12), arrowprops=dict(arrowstyle='<->', color='C1'))
    ax[1].text(6.6 + 0.5*Tb/tau, 0.13, f'$2\\pi/\\omega_b$ del péndulo = {Tb/tau:.0f}' + r'$\,\tau_1$', color='C1',
               ha='center', fontsize=7)
    ax[1].set_ylim(0.1, 0.85); ax[1].set_xlabel(r'$t/\tau_1$'); ax[1].legend(fontsize=7, loc='upper right')
    m = -1; ps = np.mod(Q[m] - 0*om, 2*np.pi)
    ax[2].scatter(np.mod(psi[m], 2*np.pi)[~atr], J[m][~atr], s=0.6, c='0.7', lw=0, label='libres', rasterized=True)
    ax[2].scatter(np.mod(psi[m], 2*np.pi)[atr], J[m][atr], s=1.5, c='C3', lw=0, label='atrapadas', rasterized=True)
    ax[2].set_xlabel(r'$\psi=Q-\omega t$'); ax[2].set_ylabel('$J$'); ax[2].set_ylim(0.06, 0.16); ax[2].legend(fontsize=7, markerscale=5)
    ax[2].set_xlim(0, 2*np.pi)
    for a_ in ax: a_.grid(alpha=0.3)
    fig.savefig(os.path.join(FIGDIR, 'controles_isla.pdf')); plt.close(fig)
    fo.close()


def energia():
    """Conservación de la energía total (atributo total_energy del HDF5) en las
    corridas de la demo y del bulto: tabla en exe/demo_eta/energia.txt y figura
    docs/demo_eta/figuras/energia.pdf."""
    import h5py, matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import bulto
    corr = [(n, ruta(n)) for n in info if os.path.exists(ruta(f'{n}.ok'))] + \
           [(n, bulto.ruta(n)) for n in bulto.info if os.path.exists(bulto.ruta(f'{n}.ok'))]
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), constrained_layout=True, sharey=True)
    with open(ruta('energia.txt'), 'w') as fo:
        fo.write(f"{'corrida':8} {'t_fin':>7} {'E0':>14} {'max|dE/E|':>10} {'final':>10}\n")
        for n, d in corr:
            f = h5py.File(os.path.join(d, 'vlasov_output.h5'), 'r')
            ks = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
            E = np.array([f[k].attrs['total_energy'] for k in ks])
            t = np.array([f[k].attrs['time'] for k in ks])
            dE = (E - E[0])/abs(E[0])
            fo.write(f'{n:8} {t[-1]:7.0f} {E[0]:+14.6e} {np.abs(dE).max():10.1e} {dE[-1]:+10.1e}\n')
            a = ax[1] if n.startswith('B') else ax[0]
            a.semilogy(t, np.maximum(np.abs(dE), 1e-13), lw=0.7, label=n)
    ax[0].set_title('demo $\\eta$'); ax[1].set_title('bulto')
    ax[0].set_ylabel(r'$|E(t)-E(0)|/|E(0)|$')
    for a in ax:
        a.set_xlabel('$t$'); a.grid(alpha=0.3); a.legend(fontsize=6, ncol=3, loc='lower right')
    fig.savefig(os.path.join(FIGDIR, 'energia.pdf')); plt.close(fig)
    print(open(ruta('energia.txt')).read())


def resolucion():
    """Curva de resolución del ajuste: señales sintéticas con gamma conocida a la
    frecuencia del modo de D3, muestreadas cada 10 y ajustadas en 10-20 tau_1 con
    polo_con_error; limpias y con un batido del 5 % (0.0735, gamma = 1e-4) más ruido
    complejo de 1e-3."""
    rng = np.random.default_rng(1)
    tau, w = 360, 0.07250
    t = np.arange(0, 7200.1, 10.0)
    with open(ruta('resolucion.txt'), 'w') as fo:
        fo.write('gamma_real  gamma_limpio  err  gamma_ruido  err  |dw|_ruido\n')
        for g in (0.0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
            h = np.exp(-(g + 1j*w)*t)
            hr = h + 0.05*np.exp(-(1e-4 + 1j*0.0735)*t) + 1e-3*(rng.normal(size=t.size) + 1j*rng.normal(size=t.size))
            a = polo_con_error(t, h, 10*tau, 20*tau); b = polo_con_error(t, hr, 10*tau, 20*tau)
            fo.write(f'{g:.0e}  {a[2]:+.2e}  {a[3]:.0e}  {b[2]:+.2e}  {b[3]:.0e}  {abs(b[0]-w):.0e}\n')
    print(open(ruta('resolucion.txt')).read())


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
            for a, b in ((1, 3), (3, 6), (5, 9), (10, 20), (20, 40)):
                if b*tau1 > t[-1] + 50:
                    continue
                e = [polo_con_error(t, x, a*tau1, b*tau1) for x in (s['h1'], s['h1_lin'])]
                v = (t >= a*tau1) & (t <= b*tau1); m = np.flatnonzero(v); m1, m2 = m[:len(m)//2], m[len(m)//2:]
                Rs = [eP[k].max()/eL[k].max() for k in (m, m1, m2)]
                fo.write(f'   polo h1 [{a},{b}] tau1: PIC w = {e[0][0]:.5f}+-{e[0][1]:.0e} g = {e[0][2]:+.2e}+-{e[0][3]:.0e}'
                         f' | lin w = {e[1][0]:.5f}+-{e[1][1]:.0e} g = {e[1][2]:+.2e}+-{e[1][3]:.0e}'
                         f' | R = {Rs[0]:.3f}+-{abs(Rs[1]-Rs[2])/2:.3f}\n')
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
            n = int(cercano(tt, np.array([tm]))[0])
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
    pend = [n for n in info if not info[n]['extra'] and (not solo or n in solo)
            and not os.path.exists(os.path.join(VIDDIR, f'{n}.mp4'))]
    with Pool(procesos) as pool:
        for n in pool.imap_unordered(_video, pend):
            print(f'  video {n}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['preparar', 'lineal', 'correr', 'analizar', 'figuras', 'videos', 'resolucion', 'energia', 'controles', 'isla',
                                     'todo'])
    ap.add_argument('--solo', nargs='*')
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    pasos = ['preparar', 'lineal', 'correr', 'analizar', 'figuras', 'videos'] if a.paso == 'todo' else [a.paso]
    for p in pasos:
        print(f'=== {p}', flush=True)
        globals()[p](a.solo) if p in ('analizar', 'videos') else globals()[p]()
