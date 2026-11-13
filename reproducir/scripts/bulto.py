"""El bulto: el dato inicial de 09_autogravedad a varias masas.

F0 = exp(-sin^2(Q/2)/sp^2) J^2 exp(-J^2/sr^2), sp = sr = 0.1, en la malla regular
(Q, J) del isócrono desnudo (state = aa_quad, dftype = gauss: lo genera el propio
código). No es un equilibrio más una perturbación: toda la masa está concentrada
en una fase orbital, con armónicos hasta k ~ 1/sp. La pregunta es si, al crecer
la masa, la autogravedad del bulto lo mantiene unido contra la cizalla de la
mezcla de fases.

Corridas: B0 sin autogravedad (la mezcla libre), B1..B5 con a0 = 1e-3 ... 0.3 y
Npc = 25 como el original, y B4q (a0 = 0.1, Npc = 100) para la resolución en Q.
Las tandas posteriores (transición en masa, otras anchuras, Kepler y las pruebas
numéricas en la transición) están comentadas junto a CORRIDAS.

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

def recorte_reticula(nj, nq, umbral, sj=0.1, sq=0.1):
    """Filas J_a y columnas Q_b de la retícula de puntos medios con nj filas en [0, 6 s_J] y nq
    columnas en [0, 2 pi) cuyo factor del dato inicial, relativo a su máximo, no baja de umbral."""
    J = (np.arange(nj) + 0.5)*6*sj/nj
    Q = (np.arange(nq) + 0.5)*2*np.pi/nq
    fj = (J/sj)**2*np.exp(1 - (J/sj)**2)
    fq = np.exp(-np.sin(0.5*Q)**2/sq**2)
    return J[fj >= umbral], Q[fq >= umbral]


# Cuarta tanda (2026-11-11): pruebas numéricas a ambos lados de la transición, B02 (mu = 0.43) y
# B025 (mu = 0.48), con un solo cambio cada vez respecto a la corrida base: doble de filas (j),
# mitad del paso radial con el mismo paso de tiempo (r) y mitad del paso de tiempo (t). Además la
# retícula fina (h): filas con dJ = 6 s_J/3200 y columnas con dQ = 2 pi/100, sin las partículas cuyo
# peso es menor que 1e-7 del máximo en alguno de los dos factores (quedan 2379 x 26, colocadas desde
# Python). El número efectivo de partículas, (sum w)^2/sum w^2, es 243 con 400 x 25, 892 con
# 400 x 100 y unas 7000 con la retícula fina; la recurrencia de las filas se retrasa ocho veces.
RETICULA = {}                                  # corrida: (filas y columnas de la retícula completa, umbral)
for base, m in (('B02', 0.02), ('B025', 0.025)):
    CORRIDAS.append((base + 'j', m, 25, True, rf'$a_0={m:g}$, $N_{{rc}}=800$'))
    CAMBIOS[base + 'j'] = {'Nrc': '800'}
    CORRIDAS.append((base + 'r', m, 25, True, rf'$a_0={m:g}$, $\Delta r=0.05$'))
    CAMBIOS[base + 'r'] = {'dr': '0.05', 'courant': '2.0'}
    CORRIDAS.append((base + 't', m, 25, True, rf'$a_0={m:g}$, $\Delta t=0.025$'))
    CAMBIOS[base + 't'] = {'courant': '0.5', 'Nt': str(2*int(round(TFIN/DT))), 'time_output': str(2*SALIDA),
                           'spatial_output': str(2*SALIDA), 'field_output': str(2*SALIDA)}
    n = base + 'h'
    RETICULA[n] = (3200, 100, 1e-7)
    filas, columnas = recorte_reticula(*RETICULA[n])
    CORRIDAS.append((n, m, len(columnas), True, rf'$a_0={m:g}$, fine lattice ($3200\times 100$, truncated)'))
    CAMBIOS[n] = {'Nrc': str(len(filas)), 'state': 'checkpoint', 'checkpointfile': f'demo_eta/bulto/ic/{n}.dat'}
# Quinta tanda (2026-11-11): masas intermedias de la familia de referencia entre B02 (mu = 0.43) y
# B025 (mu = 0.48), para distinguir un umbral de un cruce continuo. Retícula de 400 filas y columnas
# con dQ = 2 pi/100, recortada como en la cuarta tanda (quedan 299 x 26): con a0 = 0.02 la retícula
# de 400 x 100 da lo mismo que la fina en tres cifras, y recortada cuesta lo que una corrida base.
for m in (0.021, 0.022, 0.023, 0.024):
    n = f'B{int(round(m*1000)):03d}'
    RETICULA[n] = (400, 100, 1e-7)
    filas, columnas = recorte_reticula(*RETICULA[n])
    CORRIDAS.append((n, m, len(columnas), True, rf'$a_0={m:g}$, lattice $400\times 100$, truncated'))
    CAMBIOS[n] = {'Nrc': str(len(filas)), 'state': 'checkpoint', 'checkpointfile': f'demo_eta/bulto/ic/{n}.dat'}
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


def dato_isocrono(n, reticula=None):
    """Estado inicial del bulto en el isócrono (G = M = b = 1, L = 2) en la retícula recortada de la
    corrida n (RETICULA), o en la retícula (nj, nq, umbral) que se pase: la misma F0(Q, J) y la misma
    inversión que el estado aa_quad del código (Binney y Tremaine, sección 3.5.2):
    s = 2 + c (1 - e cos eta), r^2 = s (s - 2), Q = eta - e c/(c + 1) sin eta, c = (J + c_L)^2 - 1."""
    nj, nq, umbral = reticula or RETICULA[n]
    sr, sp, L = parametro(n, 'sr', 0.1), parametro(n, 'sp', 0.1), 2.0
    J, Q = recorte_reticula(nj, nq, umbral, sr, sp)
    JJ, QQ = np.meshgrid(J, Q, indexing='ij')
    E = -0.5/(JJ + 0.5*(L + np.sqrt(L*L + 4)))**2
    c = -0.5/E - 1.0
    e = np.sqrt(np.clip(1 - (L*L/c)*(1 + 1/c), 0, None))
    k = e*c/(c + 1)
    eta = QQ.copy()
    for _ in range(60):
        eta = eta - (eta - k*np.sin(eta) - QQ)/(1 - k*np.cos(eta))
    s_ = 2 + c*(1 - e*np.cos(eta))
    r = np.sqrt(s_*(s_ - 2))
    p = np.sqrt(np.clip(2*(E + 1/(1 + np.sqrt(1 + r*r)) - 0.5*L*L/(r*r)), 0, None))*np.where(np.sin(eta) >= 0, 1.0, -1.0)
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
    for n in RETICULA:
        f = os.path.join(EXE, CAMBIOS[n]['checkpointfile'])
        if not os.path.exists(f):
            os.makedirs(os.path.dirname(f), exist_ok=True)
            r, p, F, _, _ = dato_isocrono(n)
            np.savetxt(f, np.column_stack([r, p, F]), fmt='%.17e')
    print('preparado:', len(CORRIDAS))


def correr(solo=None):
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    for n, *_ in CORRIDAS:
        if (solo and n not in solo) or os.path.exists(ruta(f'{n}.ok')):
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
        # Las sumas h_k se hacen en doble precisión con cada instantánea y las partículas se guardan
        # en precisión simple al vuelo: con 6e4 partículas las mil instantáneas en doble precisión
        # no caben en la memoria.
        t, R, P, Q, J, ps, hk, hkd = ([] for _ in range(8))
        with Pool(procesos) as pool:
            for b in pool.imap(_mapear, grupos):
                for ti, r, p, q, j, pot, qd, jd in b:
                    t.append(ti); ps.append(pot)
                    hk.append([np.sum(w*np.exp(-1j*m*q)) for m in range(6)])
                    hkd.append([np.sum(w*np.exp(-1j*m*qd)) for m in range(6)])
                    R.append(r.astype(np.float32)); P.append(p.astype(np.float32))
                    Q.append(q.astype(np.float32)); J.append(j.astype(np.float32))
        t, ps, hk, hkd = np.array(t), np.array(ps), np.array(hk), np.array(hkd)
        hk = hk/hk[:, :1].real
        hkd = hkd/hkd[:, :1].real
        np.savez(ruta(n, 'bulto.npz'), t=t, hk=hk, hk_desnudo=hkd, phi_self=ps, r=rg, w=w,
                 R=np.array(R), P=np.array(P), Q=np.array(Q), J=np.array(J))
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


# Familias de corridas para la Sección 7 de las notas (docs/articulo): nombre, corridas con la
# resolución por omisión y corridas de convergencia (N_Q = 100 o N_J = 800).
FAMILIAS = [
    (r'$\sigma_Q=\sigma_J=0.1$', ['B0', 'B1', 'B2', 'B015', 'B02', 'B025', 'B3', 'B04', 'B05', 'B07', 'B4', 'B5'],
     ['B2q', 'B2j', 'B02q', 'B021', 'B022', 'B023', 'B024', 'B025h', 'B4q']),
    (r'$\sigma_J=0.05$', ['J05_0.005', 'J05_0.01', 'J05_0.02'], []),
    (r'$\sigma_J=0.2$', ['J20_0.02', 'J20_0.04', 'J20_0.08'], []),
    (r'$\sigma_Q=0.2$', ['Q2_0.01', 'Q2_0.02', 'Q2_0.04', 'Q2_0.08'], []),
    (r'$\sigma_Q=0.4$', ['Q4_0.01', 'Q4_0.02', 'Q4_0.04', 'Q4_0.08'], []),
    ('Kepler', ['K_0.01', 'K_0.02', 'K_0.04', 'K_0.08'], []),
]


def h1_ventana(n, a=10, b=20):
    """Media de |h_1| (mapa del potencial del instante) en a <= t/tau_1 <= b, con el tau_1 de la corrida."""
    d = cargar(n); t = d['t']; hk = np.asarray(d['hk'])
    tau = banda(n)[3]
    v = (t >= a*tau) & (t <= b*tau)
    return np.abs(hk[v, 1]).mean()


def campo_particulas(n):
    """Desviación del potencial que genera el gas respecto a su promedio en la segunda mitad,
    D(t) = ||Phi_gas(t) - <Phi_gas>|| / max|<Phi_gas>| en 3 <= r <= 15, con Phi_gas calculado desde
    las partículas (potencial exacto de cascarones, -sum_j m_j/max(r, r_j)). Sin autogravedad (B0)
    es el potencial que el gas generaría sin actuar sobre sí mismo. Coincide con el del código a
    1.5e-3 en las corridas con autogravedad. Se guarda en <n>/campo.npz."""
    sal = ruta(n, 'campo.npz')
    if not os.path.exists(sal):
        d = cargar(n); t = d['t']; w = np.asarray(d['w'], float); R = np.asarray(d['R'], float)
        r = np.linspace(RVENT[0], RVENT[1], 121)
        P = np.empty((len(t), len(r)))
        for m in range(len(t)):
            o = np.argsort(R[m]); Rs, ws = R[m][o], w[o]
            mc = np.concatenate([[0], np.cumsum(ws)])                       # masa con r_j < r
            pc = np.concatenate([np.cumsum((ws/Rs)[::-1])[::-1], [0]])      # sum_{r_j >= r} m_j/r_j
            k = np.searchsorted(Rs, r)
            P[m] = -(mc[k]/r + pc[k])
        media = P[t >= t[-1]/2].mean(0)
        np.savez(sal, t=t, D=np.sqrt(np.mean((P - media)**2, axis=1))/np.abs(media).max())
    d = np.load(sal)
    return d['t'], d['D']


# Estado tardío junto a la órbita circular (Sección 6 de las notas). Franjas de J/sigma_J del perfil.
XPERFIL = (0, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0)
VENTANAS = ((1, 2), (2, 4), (4, 7), (7, 10), (10, 15), (15, 20), (20, 1e9))


def nucleo(n):
    """Estado tardío de la corrida n junto a la órbita circular, con las acciones del potencial del
    instante. Se guarda en <n>/nucleo.npz:
      perfil   densidad promediada en ángulo, dm/(2 pi dJ), en las franjas XPERFIL sigma_J, media en
               t >= 20 tau_1, en unidades del máximo Fbar_max de la densidad promediada en ángulo del
               dato inicial (el mayor valor que puede alcanzar un gas que se mezcla sin cambiar de
               acción); un trozo de la nube sin mezclar puede llegar a f_max/Fbar_max = 1/(e^-z I_0(z));
      rho_c    lo mismo en J < 0.025 sigma_J, en cada instante;
      m02      masa con J < 0.2 sigma_J en cada instante (fracción del total);
      mnuc     masa dentro de una elipse del espacio de fases (r, p_r) de área 2 pi 0.05 sigma_J
               centrada en el máximo de densidad, media en t >= 20 tau_1; rho_p es la densidad media
               en esa elipse en unidades de Fbar_max (un gas mezclado da ~1, con ruido 0.8-1.2);
      rc, dr   radio medio y ancho rms, en cada instante, de las partículas que están en esa elipse
               en la última instantánea: la oscilación del núcleo."""
    sal = ruta(n, 'nucleo.npz')
    if os.path.exists(sal):
        return np.load(sal)
    d = np.load(ruta(n, 'bulto.npz'))
    t, w = d['t'], d['w'].astype(float)
    J, R, P = d['J'], d['R'], d['P']
    W = w.sum()
    sj, tau = parametro(n, 'sr', 0.1), banda(n)[3]
    n0max = np.exp(-1.0)/(sj*np.sqrt(np.pi)/4)            # máximo de la masa por unidad de J del dato inicial
    tarde = np.where(t >= 20*tau)[0]
    x = np.array(XPERFIL)*sj
    H = np.mean([np.histogram(J[m], bins=x, weights=w)[0] for m in tarde[::5]], axis=0)/W
    jc = 0.025*sj
    rho_c = np.array([w[J[m] < jc].sum() for m in range(len(t))])/W/jc/n0max
    m02 = np.array([w[J[m] < 0.2*sj].sum() for m in range(len(t))])/W
    # Elipse de acción 0.05 sigma_J alrededor del máximo de densidad en (r, p_r); kappa es la
    # frecuencia de la órbita circular del potencial externo.
    kappa = (0.5*(2.0 + np.sqrt(8.0)) if FONDO.get(n, 'isocrono') == 'isocrono' else 2.0)**-3
    ar, ap = np.sqrt(2*0.05*sj/kappa), np.sqrt(2*0.05*sj*kappa)
    br, bp = np.arange(2.0, 20.0 + ar, 0.5*ar), np.arange(-0.5, 0.5 + ap, 0.5*ap)

    def elipse(m):
        r, p = R[m].astype(float), P[m].astype(float)
        h = np.histogram2d(r, p, bins=(br, bp), weights=w)[0]
        k = h.copy()
        k[1:] += h[:-1]; k[:-1] += h[1:]                   # suma con las celdas vecinas en r
        k2 = k.copy()
        k2[:, 1:] += k[:, :-1]; k2[:, :-1] += k[:, 1:]     # y en p
        i, j = np.unravel_index(np.argmax(k2), k2.shape)
        r0, p0 = 0.5*(br[i] + br[i + 1]), 0.5*(bp[j] + bp[j + 1])
        for _ in range(8):
            sel = ((r - r0)/ar)**2 + ((p - p0)/ap)**2 < 1
            if w[sel].sum() == 0:
                break
            r0, p0 = np.average(r[sel], weights=w[sel]), np.average(p[sel], weights=w[sel])
        return ((r - r0)/ar)**2 + ((p - p0)/ap)**2 < 1, (r0, p0)

    mnuc = np.mean([w[elipse(m)[0]].sum() for m in tarde[::5]])/W
    sel, centro = elipse(len(t) - 1)
    ms = w[sel]
    rc = (R[:, sel].astype(float) @ ms)/ms.sum()
    dr = np.sqrt((((R[:, sel].astype(float) - rc[:, None])**2) @ ms)/ms.sum())
    np.savez(sal, t=t, x=np.array(XPERFIL), perfil=H/np.diff(x)/n0max, rho_c=rho_c, m02=m02, mnuc=mnuc, rc=rc, dr=dr,
             msel=ms.sum()/W, centro=np.array(centro), rho_p=mnuc/(0.05*sj*n0max))
    return np.load(sal)


def amplitudes(n):
    """Semiamplitud de la oscilación radial del núcleo en las ventanas VENTANAS (en tau_1)."""
    d = nucleo(n); t, rc = d['t'], d['rc']; tau = banda(n)[3]
    return [0.5*np.ptp(rc[(t >= a*tau) & (t <= b*tau)]) for a, b in VENTANAS]


def amortiguamiento(n, a=2.0, b=14.0):
    """Tiempo de caída a 1/e (en tau_1) de la oscilación radial del núcleo, de un ajuste exponencial
    de su envolvente en a <= t/tau_1 <= b, y el cociente entre la envolvente en b y en a."""
    d = nucleo(n); t = d['t']; tau = banda(n)[3]
    k = int(round(45/(t[1] - t[0])))
    media = np.array([d['rc'][max(0, i - k):i + k + 1].mean() for i in range(len(t))])
    env = envolvente(np.abs(d['rc'] - media), t)
    v = (t >= a*tau) & (t <= b*tau)
    p = np.polyfit(t[v]/tau, np.log(env[v]), 1)
    return -1/p[0], env[np.argmin(np.abs(t - b*tau))]/env[np.argmin(np.abs(t - a*tau))]


def estado(solo=None):
    """Tabla del estado tardío de todas las corridas: exe/demo_eta/bulto/estado.txt."""
    from scipy.special import ive
    filas = []
    for n in info:
        if (solo and n not in solo) or not os.path.exists(ruta(n, 'bulto.npz')):
            continue
        d = nucleo(n); t = d['t']; tau = banda(n)[3]
        u = t >= 20*tau
        mu = escalas(n)[2] if info[n]['sg'] else 0.0
        filas.append((FONDO.get(n, 'isocrono'), parametro(n, 'sp', 0.1), parametro(n, 'sr', 0.1), mu, n, info[n]['npc'],
                      h1_ventana(n), d['rho_c'][u].mean(), d['m02'][u].mean(), float(d['mnuc']), float(d['rho_p']), amplitudes(n),
                      d['dr'][u].mean(), 1/ive(0, 0.5/parametro(n, 'sp', 0.1)**2)))
    with open(ruta('estado.txt'), 'w') as fo:
        fo.write('Estado tardío (t >= 20 tau_1) de cada corrida. rho_c: densidad promediada en ángulo en J < 0.025 sigma_J,\n'
                 'en unidades de Fbar_max (un gas mezclado sin cambio de acción no pasa de 1; f_max/Fbar_max en la última\n'
                 'columna). m02: masa con J < 0.2 sigma_J. mnuc: masa en la elipse de acción 0.05 sigma_J alrededor del\n'
                 'máximo de densidad, y rho_p la densidad media en ella en unidades de Fbar_max. amp: semiamplitud de la oscilación radial del núcleo en las ventanas\n'
                 + ', '.join(f'[{a}, {b if b < 1e8 else "fin"}]' for a, b in VENTANAS) + ' tau_1. dr: ancho rms del núcleo.\n')
        fo.write(f'{"corrida":10} {"mu":>5} {"N_Q":>4} {"|h1|":>9} {"rho_c":>6} {"m02":>6} {"mnuc":>6} {"rho_p":>6}  '
                 + ' '.join(f'{"amp":>5}' for _ in VENTANAS) + f' {"dr":>5} {"fmax/Fb":>7}\n')
        for fondo, sq, sj, mu, n, npc, h1, rc, m02, mn, rp, amp, dr, cont in sorted(filas):
            fo.write(f'{n:10} {mu:5.2f} {npc:4d} {h1:9.2e} {rc:6.2f} {m02:6.3f} {mn:6.3f} {rp:6.2f}  '
                     + ' '.join(f'{x:5.2f}' for x in amp) + f' {dr:5.2f} {cont:7.1f}\n')
        fo.write('\nOscilación del núcleo en las corridas que se relajan con núcleo: tiempo de caída a 1/e de su envolvente\n'
                 '(ajuste en 2 <= t/tau_1 <= 14) y cociente entre la envolvente en 14 y en 2 tau_1.\n')
        for n in ('B02', 'B02q', 'Q2_0.02', 'J20_0.04', 'B025', 'B3', 'B04'):
            if os.path.exists(ruta(n, 'bulto.npz')) and (not solo or n in solo):
                te, q = amortiguamiento(n)
                fo.write(f'{n:10} mu = {escalas(n)[2]:4.2f}   tiempo = {te:6.2f} tau_1   cociente = {q:.3f}\n')
    print(open(ruta('estado.txt')).read())


def articulo():
    """Figuras de la Sección 6 de las notas, en docs/articulo/figuras/: bulto_fase.pdf (espacio de
    fases de B2q, B02q y B4q), bulto_series.pdf (|h_1| y fluctuación del potencial de la familia de
    referencia), bulto_estado.pdf (perfil de densidad junto a la órbita circular y oscilación del
    núcleo) y bulto_transicion.pdf (|h_1| tardío y densidad máxima frente a mu, todas las familias)."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    destino = os.path.join(RAIZ, 'docs', 'articulo', 'figuras')
    plt.rcParams.update({'font.size': 9})
    # 1. Espacio de fases (r, p_r) de tres nubes con N_Q = 100: se mezcla (B2q), deja un núcleo en
    #    reposo (B02q) y deja un núcleo que oscila (B4q). La última columna amplía el entorno del
    #    máximo de densidad en la última instantánea.
    fig = plt.figure(figsize=(7.2, 4.9), layout='constrained')
    gs = fig.add_gridspec(3, 5)
    for i, n in enumerate(('B2q', 'B02q', 'B4q')):
        d = np.load(ruta(n, 'bulto.npz')); t = d['t']; c = colores(d['w']); orden = np.argsort(c)
        R, P = d['R'], d['P']
        mu = escalas(n)[2]
        for k, x in enumerate((0, 0.5*TAU1, 3*TAU1, t[-1])):
            a = fig.add_subplot(gs[i, k])
            m = int(np.argmin(np.abs(t - x)))
            sc = a.scatter(R[m][orden], P[m][orden], c=c[orden], s=0.3, cmap='magma_r', vmin=-4, vmax=0, lw=0,
                           rasterized=True)
            a.set_xlim(2, 16); a.set_ylim(-0.3, 0.3)
            if i == 0:
                a.set_title(f'$t/\\tau_1={t[m]/TAU1:.1f}$', fontsize=8)
            if i == 2:
                a.set_xlabel('$r$')
            else:
                a.set_xticklabels([])
            if k:
                a.set_yticklabels([])
            else:
                a.set_ylabel('$p_r$')
                a.text(0.97, 0.95, f'$\\mu={mu:.2f}$', transform=a.transAxes, fontsize=7, va='top', ha='right')
        nu = nucleo(n)
        if nu['rho_p'] > 1.5:                              # hay núcleo: se amplía su entorno
            r0, p0 = nu['centro']
        else:                                              # no lo hay: el entorno de la órbita circular
            r0, p0 = float(R[-1][np.argmin(d['J'][-1])]), 0.0
        a = fig.add_subplot(gs[i, 4])
        a.scatter(R[-1][orden], P[-1][orden], c=c[orden], s=1.5, cmap='magma_r', vmin=-4, vmax=0, lw=0, rasterized=True)
        a.set_xlim(r0 - 1.2, r0 + 1.2); a.set_ylim(p0 - 0.08, p0 + 0.08)
        a.tick_params(labelsize=6); a.yaxis.tick_right()
        if i == 0:
            a.set_title(f'$t/\\tau_1={t[-1]/TAU1:.1f}$, zoom', fontsize=8)
        if i == 2:
            a.set_xlabel('$r$')
    fig.colorbar(sc, ax=fig.axes, label='$\\log_{10}$ of the relative mass', shrink=0.7, pad=0.06)
    fig.savefig(os.path.join(destino, 'bulto_fase.pdf'), dpi=220); plt.close(fig)
    # 2. |h_1| y fluctuación del potencial propio, familia de referencia.
    sel = ['B1', 'B2', 'B02', 'B025', 'B3', 'B4', 'B5']
    mus = {n: escalas(n)[2] for n in sel}
    cmap = plt.get_cmap('viridis')
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4), layout='constrained')
    d = cargar('B0'); t = d['t']
    ax[0].semilogy(t/TAU1, envolvente(np.abs(np.asarray(d['hk'])[:, 1]), t), 'k--', lw=1.0, label='no self-gravity')
    tc, Dc = campo_particulas('B0')
    ax[1].semilogy(tc/TAU1, envolvente(Dc, tc), 'k--', lw=1.0)
    for n in sel:
        col = cmap(0.85*mus[n]/1.16)
        d = cargar(n); t = d['t']
        ax[0].semilogy(t/TAU1, envolvente(np.abs(np.asarray(d['hk'])[:, 1]), t), color=col, lw=1.0,
                       label=f'$\\mu={mus[n]:.2f}$')
        tc, Dc = campo_particulas(n)
        ax[1].semilogy(tc/TAU1, envolvente(Dc, tc), color=col, lw=1.0)
        q = {'B2': 'B2q', 'B02': 'B02q', 'B4': 'B4q'}.get(n)
        if q:
            d = cargar(q); t = d['t']
            ax[0].semilogy(t/TAU1, envolvente(np.abs(np.asarray(d['hk'])[:, 1]), t), color=col, lw=1.0, ls=':')
            tc, Dc = campo_particulas(q)
            ax[1].semilogy(tc/TAU1, envolvente(Dc, tc), color=col, lw=1.0, ls=':')
    ax[0].set_ylabel('$|h_1|/h_0$'); ax[0].set_ylim(1e-5, 2)
    ax[1].set_ylabel('$\\|\\Phi_{gas}-\\langle\\Phi_{gas}\\rangle\\|/\\max|\\langle\\Phi_{gas}\\rangle|$')
    ax[1].set_ylim(1e-5, 0.3)
    for a in ax:
        a.set_xlabel('$t/\\tau_1$'); a.set_xlim(0, TFIN/TAU1); a.grid(alpha=0.3)
    ley = fig.legend(*ax[0].get_legend_handles_labels(), loc='outside upper center', ncol=8, fontsize=7, frameon=False,
                     handlelength=1.6, columnspacing=1.0)
    for l in ley.get_lines():
        l.set_linewidth(2.0)
    for a, l in zip(ax, 'ab'):
        a.text(0.97, 0.95, f'({l})', transform=a.transAxes, ha='right', va='top')
    fig.savefig(os.path.join(destino, 'bulto_series.pdf')); plt.close(fig)
    # 3. Estado tardío: perfil de densidad promediada en ángulo junto a la órbita circular, y
    #    oscilación radial del núcleo.
    from scipy.special import ive
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2), layout='constrained')
    xs = np.geomspace(1e-2, 3, 200)
    ax[0].loglog(xs*0.1, xs**2*np.exp(1 - xs**2), color='0.6', lw=2.5, label='initial')
    for n in ('B0', 'B2q', 'B015', 'B02q', 'B3', 'B05'):
        d = nucleo(n); mu = escalas(n)[2] if info[n]['sg'] else 0.0
        xc = 0.5*(d['x'][1:] + d['x'][:-1])*parametro(n, 'sr', 0.1)
        ax[0].loglog(xc, d['perfil'], marker='o', ms=2.5, lw=1.0, color='k' if n == 'B0' else cmap(0.85*mu/1.16),
                     ls='--' if n == 'B0' else '-', label='no self-gravity' if n == 'B0' else f'$\\mu={mu:.2f}$')
    ax[0].axhline(1, color='k', lw=0.6, ls=':'); ax[0].axhline(1/ive(0, 50.0), color='k', lw=0.6, ls=':')
    ax[0].text(0.28, 1.15, '$\\hat{F}_{max}$', fontsize=7, ha='right'); ax[0].text(0.28, 20.5, '$F_{max}$', fontsize=7, ha='right')
    ax[0].set_xlim(5e-4, 0.3); ax[0].set_ylim(1e-3, 40)
    ax[0].set_xlabel('$J$'); ax[0].set_ylabel('angle-averaged DF$/\\hat{F}_{max}$')
    ax[0].legend(fontsize=6.5, frameon=False, loc='lower right')
    # La serie de masas intermedias de la quinta tanda (N_Q = 100; 0.025 con la retícula fina), con un tono de azul
    # por mu del más claro al más oscuro, y dos corridas de mu mayor en gris.
    cruce = ('B02q', 'B021', 'B022', 'B023', 'B024', 'B025h')
    azules = plt.get_cmap('Blues')(np.linspace(0.45, 1.0, len(cruce)))
    for n, col, ls in [(n, c, '-') for n, c in zip(cruce, azules)] + [('B3', '0.55', '--'), ('B4q', '0.3', '--')]:
        d = nucleo(n); t = d['t']; mu = escalas(n)[2]
        k = int(round(45/(t[1] - t[0])))
        media = np.array([d['rc'][max(0, i - k):i + k + 1].mean() for i in range(len(t))])
        ax[1].semilogy(t/TAU1, envolvente(np.abs(d['rc'] - media), t), color=col, ls=ls, lw=1.0, label=f'$\\mu={mu:.2f}$')
    ax[1].set_xlim(0, TFIN/TAU1); ax[1].set_ylim(1.5e-3, 2)
    ax[1].set_xlabel('$t/\\tau_1$'); ax[1].set_ylabel('amplitude $a_c$ of the oscillation of the core')
    ax[1].legend(fontsize=6, frameon=False, loc='lower left', ncol=2)
    for a, l in zip(ax, 'ab'):
        a.grid(alpha=0.3); a.set_title(f'({l})', loc='left', fontsize=9)
    fig.savefig(os.path.join(destino, 'bulto_estado.pdf')); plt.close(fig)
    # 4. Transición: |h_1| medio en 10 <= t/tau_1 <= 20 y densidad en torno al máximo, frente a mu.
    estilos = [('o', 'k'), ('v', '#0072B2'), ('^', '#009E73'), ('s', '#E69F00'), ('D', '#CC79A7'), ('*', '#D55E00')]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2), layout='constrained')
    filas = []
    for (nombre, base, conv), (mk, col) in zip(FAMILIAS, estilos):
        x = [escalas(n)[2] if info[n]['sg'] else 0.0 for n in base]
        xc = [escalas(n)[2] for n in conv]
        for a, fun in zip(ax, (h1_ventana, lambda n: float(nucleo(n)['rho_p']))):
            a.semilogy(x, [fun(n) for n in base], mk, color=col, ms=6 if mk == '*' else 4.5, label=nombre, ls='-', lw=0.6)
            if conv:
                a.semilogy(xc, [fun(n) for n in conv], mk, color=col, mfc='none', ms=6.5, ls='none',
                           label='finer lattices')
        filas += [(n, m, h1_ventana(n)) for n, m in zip(base + conv, x + xc)]
    ax[0].set_ylabel('$|h_1|/h_0$ for $10\\leq t/\\tau_1\\leq 20$'); ax[0].set_ylim(3e-5, 2)
    ax[1].set_ylabel('$m_c/m_{mix}$ for $t\\geq 20\\tau_1$'); ax[1].set_ylim(0.5, 30)
    ax[1].axhline(1, color='k', lw=0.6, ls=':')
    for a, l in zip(ax, 'ab'):
        a.set_xlabel('$\\mu$'); a.set_xlim(-0.03, 1.2); a.grid(alpha=0.3)
        a.set_title(f'({l})', loc='left', fontsize=9)
    ax[0].legend(fontsize=6.5, loc='lower right', frameon=False)
    fig.savefig(os.path.join(destino, 'bulto_transicion.pdf')); plt.close(fig)
    with open(ruta('articulo.txt'), 'w') as fo:
        fo.write('corrida      mu     |h1| en 10-20 tau_1\n')
        for n, x, y in sorted(filas, key=lambda r: (r[1], r[0])):
            fo.write(f'{n:10} {x:5.2f}  {y:.3e}\n')
        fo.write('\nRelajación del potencial que genera el gas (campo_particulas): primera vez que la envolvente\n'
                 '(ventana de +-50) de D baja de 1e-2, 1e-3 y 3e-4, en unidades de tau_1, y media de D en 10-20 tau_1.\n')
        for n in ['B0', 'B1', 'B2', 'B2q', 'B015', 'B02', 'B02q', 'B025', 'B3', 'B04', 'B05', 'B07', 'B4', 'B4q', 'B5']:
            tc, Dc = campo_particulas(n)
            k = int(round(50/(tc[1] - tc[0])))
            env = np.array([Dc[max(0, i - k):i + k + 1].max() for i in range(len(Dc))])
            cruces = [f'{tc[np.argmax(env < x)]/TAU1:6.2f}' if (env < x).any() else f'{"--":>6}' for x in (1e-2, 1e-3, 3e-4)]
            v = (tc >= 10*TAU1) & (tc <= 20*TAU1)
            mu = escalas(n)[2] if info[n]['sg'] else 0.0
            fo.write(f'{n:6} mu = {mu:4.2f}  ' + ' '.join(cruces) + f'   D en 10-20 tau_1: {Dc[v].mean():.2e}\n')
    print('figuras en', destino)


def medidas_articulo(n):
    """Las medidas de la corrida n que usan el texto y las tablas de la Sección 6 del artículo: mu, |h_1| en
    10-20 tau_1, masa con J < 0.2 s_J y m_c (t >= 20 tau_1), m_c/m_mix, semiamplitud del radio del núcleo en
    [1, 2] y desde 20 tau_1 (amplitudes), D en 10-20 tau_1, y de la envolvente a_c de la oscilación del núcleo
    su primer mínimo local en 2-30 tau_1 (el menor valor en +-2 tau_1; y cuándo) y su rango en t >= 25 tau_1. Con
    mu >= 0.47 la envolvente tiene un segundo mínimo de valor parecido, así que el mínimo global salta de uno a
    otro entre retículas; el primero no."""
    d = nucleo(n); t = d['t']; tau = banda(n)[3]
    k = int(round(45/(t[1] - t[0])))
    media = np.array([d['rc'][max(0, i - k):i + k + 1].mean() for i in range(len(t))])
    env = envolvente(np.abs(d['rc'] - media), t)
    i = None
    for a in np.arange(2, 30, 0.5):                 # primer mínimo local: el menor valor en +-2 tau_1
        z = (t >= a*tau) & (t < (a + 2)*tau)
        k = int(np.argmin(np.where(z, env, np.inf)))
        if env[k] <= env[(t >= t[k] - 2*tau) & (t <= t[k] + 2*tau)].min():
            i = k
            break
    if i is None:                                   # sin mínimo local: decae hasta el final; el menor valor
        i = int(np.argmin(np.where((t >= 2*tau) & (t <= 30*tau), env, np.inf)))
    tc, Dc = campo_particulas(n)
    amp = amplitudes(n)
    return dict(mu=escalas(n)[2] if info[n]['sg'] else 0.0, h1=h1_ventana(n), m02=d['m02'][t >= 20*tau].mean(),
                mc=float(d['mnuc']), rho_p=float(d['rho_p']), a12=amp[0], a20=amp[-1],
                D=Dc[(tc >= 10*tau) & (tc <= 20*tau)].mean(), acmin=env[i], tmin=t[i]/tau,
                ac_lo=env[t >= 25*tau].min(), ac_hi=env[t >= 25*tau].max(), t=t/tau, env=env, rc=d['rc'] - media)


def _sci(x):
    if x >= 0.0995:
        return f'${x:.2f}$'
    e = int(np.floor(np.log10(x))); a = x/10.0**e
    if round(a, 1) >= 10:
        a, e = a/10, e + 1
    return rf'${a:.1f}\times 10^{{{e}}}$'


def _macro(nombre, archivo, filas):
    """Las filas de una tabla como la macro \\<nombre> en docs/articulo/<archivo>: main.tex la carga en el
    preámbulo (un \\input dentro de tabular deja una fila vacía)."""
    cuerpo = '%\n'.join(filas)
    open(os.path.join(RAIZ, 'docs', 'articulo', archivo), 'w').write(
        '% Escrito por reproducir/scripts/bulto.py tablas; no editar a mano.\n\\newcommand{\\' + nombre + '}{%\n' + cuerpo + '%\n}\n')


# Filas de las tablas de la Sección 6: estado tardío de la familia de referencia (la mejor retícula de cada
# masa), la serie fina entre mu = 0.43 y 0.48, y la ubicación del cambio en cada familia.
ESTADO = [('B0', 0), ('B1', 0.001), ('B2q', 0.01), ('B015', 0.015), ('B02q', 0.02), ('B025h', 0.025), ('B3', 0.03),
          ('B04', 0.04), ('B05', 0.05), ('B07', 0.07), ('B4q', 0.1), ('B5', 0.3)]
CRUCE = [('B02q', 0.02), ('B021', 0.021), ('B022', 0.022), ('B023', 0.023), ('B024', 0.024), ('B025h', 0.025)]
CRUCE_FAMILIAS = [('isochrone', r'$\sigma_Q = \sigma_J = 0.1$', 'B02q', 'B025h'),
                  ('', r'$\sigma_J = 0.05$', 'J05_0.005', 'J05_0.01'), ('', r'$\sigma_J = 0.2$', 'J20_0.04', 'J20_0.08'),
                  ('', r'$\sigma_Q = 0.2$', 'Q2_0.02', 'Q2_0.04'), ('', r'$\sigma_Q = 0.4$', 'Q4_0.04', 'Q4_0.08'),
                  ('Kepler', r'$\sigma_Q = \sigma_J = 0.1$', 'K_0.02', 'K_0.04')]


def tablas():
    """Escribe las tablas de la Sección 6 del artículo (docs/articulo/tab_estado.tex, tab_cruce.tex y
    tab_familias.tex) e imprime las cifras del texto que no están en ellas: las pruebas de resolución de la
    cuarta tanda a mu = 0.43 y 0.48 y el salto de fase de la oscilación del núcleo en su mínimo."""
    from scipy.signal import hilbert
    f2 = lambda x: f'{x:.2f}'
    filas = []
    for n, m in ESTADO:
        M = medidas_articulo(n)
        nuc = M['rho_p'] >= 1.5                    # hay un núcleo: m_c por encima de la cota de la mezcla
        filas.append(f'${m:g}$ & ${M["mu"]:.2f}$ & {_sci(M["h1"])} & ${M["m02"]:.3f}$ & ${M["mc"]:.3f}$ & ${M["rho_p"]:.1f}$ & '
                     + (f'${f2(M["a12"])}$ & ${f2(M["a20"])}$' if nuc else '-- & --') + ' \\\\')
    _macro('TabEstado', 'tab_estado.tex', filas)
    filas = []
    for n, m in CRUCE:
        M = medidas_articulo(n)
        filas.append(f'${m:g}$ & ${M["mu"]:.2f}$ & {_sci(M["h1"])} & {_sci(M["D"])} & ${M["mc"]:.3f}$ & '
                     f'${M["acmin"]:.3f}$ & ${M["tmin"]:.0f}$ & ${M["ac_lo"]:.2f}$--${M["ac_hi"]:.2f}$ \\\\')
    _macro('TabCruce', 'tab_cruce.tex', filas)
    filas = []
    for pot, anchos, a, b in CRUCE_FAMILIAS:
        A, B_ = medidas_articulo(a), medidas_articulo(b)
        filas.append(f'{pot} & {anchos} & ${A["mu"]:.2f}$ & {_sci(A["h1"])} & ${A["rho_p"]:.1f}$ & '
                     f'${B_["mu"]:.2f}$ & {_sci(B_["h1"])} & ${B_["rho_p"]:.1f}$ \\\\')
    _macro('TabFamilias', 'tab_familias.tex', filas)
    print('tablas escritas en docs/articulo/tab_*.tex')
    print('\nPruebas de resolución (cuarta tanda). a_c: mínimo de la envolvente en 2-30 tau_1 (cuándo) y rango en t >= 25 tau_1.')
    for n in ('B02', 'B02q', 'B02j', 'B02t', 'B02r', 'B02h', 'B025', 'B025j', 'B025t', 'B025r', 'B025h'):
        if os.path.exists(ruta(n, 'nucleo.npz')) or os.path.exists(ruta(n, 'bulto.npz')):
            M = medidas_articulo(n)
            print(f'  {n:6} mu {M["mu"]:.3f}  |h1| {M["h1"]:.3e}  m_c {M["mc"]:.3f}  D {M["D"]:.2e}  a_c mín {M["acmin"]:.4f} '
                  f'({M["tmin"]:.1f})  tarde {M["ac_lo"]:.3f}-{M["ac_hi"]:.3f}')
    print('\nFase de la oscilación del núcleo antes y después del mínimo (señal analítica de r_c - media).')
    for n in ('B02q', 'B02h', 'B02j'):
        M = medidas_articulo(n); t, x = M['t'], M['rc']
        fase = np.unwrap(np.angle(hilbert(x)))
        antes = (t >= M['tmin'] - 8) & (t <= M['tmin'] - 2); despues = (t >= M['tmin'] + 2) & (t <= M['tmin'] + 8)
        p1, p2 = np.polyfit(t[antes], fase[antes], 1), np.polyfit(t[despues], fase[despues], 1)
        salto = np.angle(np.exp(1j*(np.polyval(p2, M['tmin']) - np.polyval(p1, M['tmin']))))
        print(f'  {n:5} mínimo en {M["tmin"]:.1f} tau_1; frecuencia antes {p1[0]/TAU1:.5f} y después {p2[0]/TAU1:.5f};'
              f' salto de fase {salto/np.pi:+.2f} pi')


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
    ap.add_argument('paso', choices=['preparar', 'correr', 'analizar', 'figuras', 'videos', 'transicion', 'estado', 'articulo',
                                     'tablas'])
    ap.add_argument('--solo', nargs='*')
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    print(f'=== {a.paso}', flush=True)
    globals()[a.paso](a.solo) if a.paso in ('correr', 'analizar', 'videos', 'transicion', 'estado') else globals()[a.paso]()
