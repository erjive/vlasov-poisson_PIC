"""El escenario de Hadžić, Rein, Schrecker y Straub (ARMA 249, 45, 2025): |L| fijo, masa
puntual de masa 1 en el centro y politropos en la energía F = A (E_t - E)^k, con el borde
E_t = E(J_t). El teorema dice que, con masa pequeña, las perturbaciones se amortiguan si
k > 1 y no si 1/2 < k <= 1. Aquí se estudia con la teoría lineal y con el código PIC.

Pasos (cada uno reutiliza lo que ya exista):
    python3 hadzic.py lineal     mapa de lambda_edge(k, a0) y soluciones lineales en el tiempo
    python3 hadzic.py preparar   estados iniciales (equilibrio.py) y .par de las corridas
    python3 hadzic.py correr     corridas PIC, una detrás de otra, 4 hilos
    python3 hadzic.py serie [-j N] <corrida> ...   serie.npz de esas corridas, con N procesos
                                 (para el cluster: reproducir/cluster/serie.slurm)
    python3 hadzic.py analizar   dPhi y h_1 en el mapa del equilibrio, menos eps = 0
    python3 hadzic.py figuras    figuras del informe (docs/hadzic/figuras/)
    python3 hadzic.py videos     un video por corrida perturbada (exe/hadzic/videos/)

Parámetros: L0 = 2 (su L = 4), J_t = 0.7, es decir E_t = -0.0686 sin masa propia, dentro
de su condición de un solo hueco (-0.079 < E_t < 0; Omega_max/Omega_min = 2.46). En el
código, BGtype = "sphere" es una bola uniforme de masa 1 y radio 1, es decir una masa
puntual para r > 1, y la cáscara está en 2.4 < r < 12.2. Salidas en exe/hadzic/.
"""
import os, sys, subprocess, time, argparse
# analizar y videos reparten las corridas en cuatro procesos: un hilo de BLAS en cada uno. La
# variable tiene que estar antes de importar numpy; puesta después, cada proceso usaba ~1.6
# núcleos y los cuatro, unos 16 hilos en los 4 núcleos físicos.
if len(sys.argv) > 1 and sys.argv[1] in ('serie', 'analizar', 'videos'):
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[v] = '1'
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
EXE = os.path.join(RAIZ, 'exe')
BASE = os.path.join(EXE, 'hadzic')
PARDIR = os.path.join(RAIZ, 'reproducir', 'corridas', '13_hadzic')
PLANTILLA = os.path.join(RAIZ, 'reproducir', 'corridas', '11_landau', 'landau__L_a1e-2_n400_e0.1.par')
JT, L0 = 0.7, 2.0
J1, SJ1 = 0.35, 0.20                    # función de prueba de h_1: B = J^2 exp(-(J-J1)^2/SJ1^2)
# courant = 2 (el de 11_landau): dt = courant dr/pmax = 0.1. Con a0 = 0.01, una instantánea cada
# 200 pasos (20 unidades de tiempo) hasta t_fin = 8000, unos 95 tau_1. Con a0 = 1 la banda sube a
# 0.15-0.33 y tau_1 baja a ~37: una instantánea cada 50 pasos (5 unidades, para no submuestrear
# la frecuencia) hasta t_fin = 4000.
NRC, NPC, COURANT, DT = 400, 25, 2.0, 0.1
TFIN, SALIDA = {0.01: 8000.0, 1.0: 4000.0}, {0.01: 200, 1.0: 50}
# Radios de la cáscara para ||dPhi|| (la solución lineal da 3-15): con a0 = 1 la cáscara se
# contrae a 1.8 < r < 7.5.
ZONA = {0.01: (3.0, 12.0), 1.0: (3.0, 7.5)}
# Mapa lineal: exponentes del borde y masas.
K_MAPA = [0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
A0_MAPA = [0.01, 0.03, 0.1, 0.3, 0.6, 1.0]
# Corridas PIC: nombre, k, a0, eps y la referencia eps = 0 (None en las referencias).
#  serie 1, a0 = 0.01 y eps = 0.1: el régimen del teorema;
#  serie 2, a0 = 1 y eps = 0.1: el modo de k <= 1 bien ligado, y el umbral (k = 1.25 y 1.5
#           tienen modo, k = 2 no);
#  serie 3, a0 = 0.01 y eps = 0.5, con las referencias de la serie 1: piso de ruido más bajo.
CORRIDAS = []
for k in [0.75, 1.0, 1.5, 2.0]:
    CORRIDAS += [(f'D_k{k:g}', k, 0.01, 0.1, f'Z_k{k:g}'), (f'Z_k{k:g}', k, 0.01, 0.0, None)]
for k in [0.75, 1.0, 1.25, 1.5, 2.0]:
    CORRIDAS += [(f'D_k{k:g}_a1', k, 1.0, 0.1, f'Z_k{k:g}_a1'), (f'Z_k{k:g}_a1', k, 1.0, 0.0, None)]
for k in [0.75, 1.0, 1.5, 2.0]:
    CORRIDAS += [(f'D5_k{k:g}', k, 0.01, 0.5, f'Z_k{k:g}')]
#  serie 4, a0 = 1, k = 1.25 y 1.5: ¿la pérdida de amplitud del modo cercano al borde es física
#           o numérica? Cada corrida cambia una sola cosa respecto de D_k{k}_a1: eps = 0.03 o 0.3
#           (con la misma referencia), dr/2 con el mismo dt (courant = 4, mismo dato inicial) o
#           cuatro veces más partículas (800 x 50); estas dos, con su propia referencia.
CAMBIOS = {}                            # cambios en el .par y dato inicial ('ic') de la serie 4
for k in [1.25, 1.5]:
    b = f'k{k:g}_a1'
    CORRIDAS += [(f'De03_{b}', k, 1.0, 0.03, f'Z_{b}'), (f'De3_{b}', k, 1.0, 0.3, f'Z_{b}')]
    for s, c in (('dr', {'dr': '0.05', 'courant': '4.0'}), ('N', {'Nrc': '800', 'Npc': '50'})):
        CORRIDAS += [(f'D{s}_{b}', k, 1.0, 0.1, f'Z{s}_{b}'), (f'Z{s}_{b}', k, 1.0, 0.0, None)]
        for x in 'DZ':
            CAMBIOS[f'{x}{s}_{b}'] = dict(c, ic=f'{x}_{b}') if s == 'dr' else c
#  serie 5, piloto con diez veces más partículas: k = 1.25, a0 = 1, eps = 0.03, 1280 x 80
#           partículas (la proporción N_J/N_Q = 16 de 400 x 25 y de 800 x 50) y dt = 0.05
#           (courant = 1), con su propia referencia. Con eps = 0.1 la respuesta deja de ser
#           lineal desde t ~ 300, y con eps = 0.03 y 10^4 partículas el ruido la tapa desde
#           t ~ 1000: ¿con menos ruido se sigue la respuesta lineal hasta t_fin?
for k in [1.25]:
    b = f'k{k:g}_a1'
    CORRIDAS += [(f'DP_{b}', k, 1.0, 0.03, f'ZP_{b}'), (f'ZP_{b}', k, 1.0, 0.0, None)]
    for x in 'DZ':
        CAMBIOS[f'{x}P_{b}'] = {'Nrc': '1280', 'Npc': '80', 'courant': '1.0'}


def ic_de(nombre):
    """Nombre del dato inicial de la corrida (el suyo, salvo las de dr/2)."""
    return CAMBIOS.get(nombre, {}).get('ic', nombre)


def nombre_lin(k, a0):
    """Solución lineal en el tiempo del equilibrio (k, a0)."""
    return f'lin_k{k:g}' + ('' if a0 == 0.01 else f'_a{a0:g}')


def ruta(*p):
    return os.path.join(BASE, *p)


def equilibrio(salida, k, a0, eps, nrc=NRC, npc=NPC):
    """Genera, si falta, el estado inicial y el equilibrio con equilibrio.py."""
    if os.path.exists(salida):
        return
    cmd = [sys.executable, os.path.join(AQUI, 'equilibrio.py'), '--fondo', 'puntual', '--forma',
           'polE', '--jt', str(JT), '--k', str(k), '--a0', str(a0), '--l0', str(L0), '--eps',
           str(eps), '--pert', 'suave3', '--nrc', str(nrc), '--npc', str(npc), '--salida', salida]
    subprocess.run(cmd, check=True, stdout=open(os.path.splitext(salida)[0] + '.log', 'w'),
                   stderr=subprocess.STDOUT)


# ------------------------------------------------------------------ lineal
def lineal():
    """lambda_edge(k, a0), el modo discreto si lo hay y, para k <= 1, lambda cerca del
    borde; después, la solución lineal en el tiempo de los casos de las corridas PIC."""
    from lambda_borde import Lazo
    from lineal import resolver
    os.makedirs(ruta('lineal'), exist_ok=True)
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    if os.path.exists(ruta('lineal', 'lambda.txt')):        # el mapa ya está hecho
        print(open(ruta('lineal', 'lambda.txt')).read())
        K_MAPA_, A0_MAPA_ = [], []
    else:
        K_MAPA_, A0_MAPA_ = K_MAPA, A0_MAPA
    w(f'Masa puntual, polE, J_t = {JT}, L0 = {L0}. delta en anchos de banda bajo Omega_min.')
    w(f'{"k":>5} {"a0":>6} {"iter":>5} {"Omega_min":>10} {"Omax/Omin":>9} {"lambda_edge":>12} '
      f'{"l(1e-3)":>8} {"l(1e-6)":>8} {"omega_d":>9} {"x_d":>8}')
    for k in K_MAPA_:
        for a0 in A0_MAPA_:
            npz = ruta('lineal', f'P_k{k:g}_a{a0:g}_equilibrio.npz')
            try:
                equilibrio(ruta('lineal', f'P_k{k:g}_a{a0:g}.dat'), k, a0, 0.0)
            except subprocess.CalledProcessError:
                w(f'{k:5g} {a0:6g}  el equilibrio no converge'); continue
            log = open(ruta('lineal', f'P_k{k:g}_a{a0:g}.log')).read()
            it = log.count('iteración')
            z = Lazo(npz)
            an = z.Om_max - z.Om_min
            lb, l3, l6 = (z.lam(z.Om_min - d*an) for d in (0.0, 1e-3, 1e-6))
            wd = z.omega_modo()
            txt = f'{wd:9.5f} {(wd - z.Om_min)/an:+8.4f}' if wd is not None else f'{"--":>9} {"--":>8}'
            w(f'{k:5g} {a0:6g} {it:5d} {z.Om_min:10.5f} {z.Om_max/z.Om_min:9.3f} {lb:12.4f} '
              f'{l3:8.4f} {l6:8.4f} {txt}')
    if K_MAPA_:
        open(ruta('lineal', 'lambda.txt'), 'w').write('\n'.join(lineas) + '\n')
    tabla_eta()
    for k, a0 in sorted({(k, a0) for _, k, a0, eps, ref in CORRIDAS if ref}):
        sal = ruta('lineal', nombre_lin(k, a0) + '.npz')
        if os.path.exists(sal):
            continue
        t0 = time.time()
        t, h1, h2, rmed, dphi = resolver(ruta('lineal', f'P_k{k:g}_a{a0:g}_equilibrio.npz'),
                                         1600, 32, 0.5, TFIN[a0], verboso=False, j1=J1, sj1=SJ1)
        np.savez(sal, t=t, h1=h1, h2=h2, r=rmed, dphi=dphi)
        print(f'  lineal k = {k:g}, a0 = {a0:g} hasta t = {TFIN[a0]:g} ({time.time()-t0:.0f} s)', flush=True)


def tabla_eta():
    """eta = dOmega/(ancho/media) de cada equilibrio del mapa, como en eta.py: medias pesadas
    por la masa con la tabla E(J) del equilibrio, y la del fondo desnudo con eta.banda(0)."""
    from eta import banda
    from equilibrio import F_E, borde_E
    lineas = ['eta de los equilibrios del mapa (masa puntual, polE, J_t = 0.7).',
              f'{"k":>5} {"a0":>6} {"eta":>7} {"dOmega":>7} {"ancho/media":>11}']
    for k in K_MAPA:
        b0 = banda(0.0, None, JT, L0, forma='polE', jt=JT, k=k, fondo='puntual')
        for a0 in A0_MAPA:
            d = np.load(ruta('lineal', f'P_k{k:g}_a{a0:g}_equilibrio.npz'))
            J = np.linspace(0.0, JT, 4000)
            E = np.interp(J, d['J_t'], d['E_t']); Om = np.gradient(E, J)
            F = F_E('polE', E, *borde_E('polE', d['E_t'], d['J_t'], JT, None), k)
            media = np.sum(F*Om)/np.sum(F)
            dO = (media - b0['Om_media'])/b0['Om_media']; rel = (Om.max() - Om.min())/media
            lineas.append(f'{k:5g} {a0:6g} {dO/rel:7.3f} {dO:7.3f} {rel:11.3f}')
    open(ruta('lineal', 'eta.txt'), 'w').write('\n'.join(lineas) + '\n')
    print('\n'.join(lineas))


# ------------------------------------------------------------------ preparar
def escribir_par(nombre, k, a0, eps, destino=PARDIR):
    """.par de la corrida a partir de la plantilla de 11_landau, con los cambios de CAMBIOS."""
    c = CAMBIOS.get(nombre, {})
    # dt = courant dr/pmax (pmax = 2): 0.1 en todas las corridas salvo las de la serie 5. Las
    # instantáneas salen en los mismos tiempos con cualquier dt, cada SALIDA[a0] DT.
    dt = float(c.get('courant', COURANT))*float(c.get('dr', 0.1))/2.0
    cada = SALIDA[a0]*DT/dt
    assert abs(cada - round(cada)) < 1e-9 and abs(TFIN[a0]/dt - round(TFIN[a0]/dt)) < 1e-6
    sal = str(int(round(cada)))
    valores = {'courant': str(COURANT), 'Nt': str(int(round(TFIN[a0]/dt))), 'time_output': sal,
               'spatial_output': sal, 'field_output': sal,
               'Nrc': str(NRC), 'Npc': str(NPC), 'Lfix': str(L0),
               'directory': f'hadzic/{nombre}', 'a0': str(a0), 'BGtype': 'sphere',
               'checkpointfile': f'hadzic/ic/{ic_de(nombre)}.dat', 'j1': f'{J1:.6f}',
               'sj1': f'{SJ1:.6f}', 'state': 'checkpoint'}
    valores.update({x: v for x, v in c.items() if x != 'ic'})
    lineas = [f'# Escenario de Hadžić, corrida {nombre}: polE con k = {k:g}, J_t = {JT}, '
              f'a0 = {a0:g}, eps = {eps}, t_fin = {TFIN[a0]:g}.',
              '# Masa puntual: BGtype = sphere (masa 1, radio 1). Generado por '
              'reproducir/scripts/hadzic.py a partir de 11_landau.']
    if c:
        lineas.append('# Cambios respecto de la corrida base: '
                      + ', '.join(f'{x} = {v}' for x, v in c.items() if x != 'ic')
                      + (f'; dato inicial de {c["ic"]}.' if 'ic' in c else '.'))
    for l in open(PLANTILLA).read().split('\n'):
        if l.startswith('#') or '=' not in l:
            continue
        clave = l.split('=')[0].strip()
        lineas.append(f'{clave:<16} = {valores[clave]}' if clave in valores else l)
    open(os.path.join(destino, f'hadzic__{nombre}.par'), 'w').write('\n'.join(lineas) + '\n')


def preparar():
    os.makedirs(ruta('ic'), exist_ok=True); os.makedirs(PARDIR, exist_ok=True)
    for nombre, k, a0, eps, ref in CORRIDAS:
        if os.path.exists(ruta(f'{nombre}.ok')):       # corrida hecha: su .par queda como se usó
            continue
        c = CAMBIOS.get(nombre, {})
        dat = ruta('ic', f'{ic_de(nombre)}.dat')
        if not os.path.exists(dat):
            t0 = time.time()
            equilibrio(dat, k, a0, eps, int(c.get('Nrc', NRC)), int(c.get('Npc', NPC)))
            print(f'  estado inicial {ic_de(nombre)} ({time.time()-t0:.0f} s)', flush=True)
        escribir_par(nombre, k, a0, eps)
    print('preparado:', len(CORRIDAS), 'corridas')


# ------------------------------------------------------------------ correr
def metadatos(nombre, par, salida):
    """sha256 del ejecutable y del .par, commit del repositorio y fecha (como en demo_eta)."""
    import hashlib, datetime
    sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
    git = lambda *a: subprocess.run(['git', *a], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
    src = git('log', '-1', '--format=%h', '--', 'src/') + ('+cambios' if git('status', '--porcelain', 'src/') else '')
    with open(os.path.join(PARDIR, 'METADATOS.txt'), 'a') as fo:
        fo.write(f'{nombre:8} {datetime.datetime.now().isoformat(timespec="minutes")}  par {sha(par)[:16]}'
                 f'  VP_PIC {sha(os.path.join(EXE, "VP_PIC"))[:16]}  src {src}  HEAD {git("rev-parse", "--short", "HEAD")}\n')
    open(salida, 'w').write(f'corrida    {nombre}\nVP_PIC     {sha(os.path.join(EXE, "VP_PIC"))}\n'
                            f'par        {sha(par)}  {os.path.relpath(par, RAIZ)}\n'
                            f'HEAD       {git("rev-parse", "--short", "HEAD")}\nsrc/       {src}\n')


def correr():
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    for nombre, *_ in CORRIDAS:
        if os.path.exists(ruta(f'{nombre}.ok')):
            continue
        par = os.path.join(PARDIR, f'hadzic__{nombre}.par')
        t0 = time.time()
        r = subprocess.run(['./VP_PIC', par], cwd=EXE, stdout=open(ruta(f'{nombre}.log'), 'w'),
                           stderr=subprocess.STDOUT, env=env)
        print(f'{"OK" if r.returncode == 0 else "FALLO":5} {nombre} ({time.time()-t0:.0f} s)', flush=True)
        if r.returncode == 0:
            open(ruta(f'{nombre}.ok'), 'w').write('')
            metadatos(nombre, par, ruta(f'{nombre}.meta'))


# ------------------------------------------------------------------ analizar
def _tramo(arg):
    """t, h_1 sin normalizar, dPhi y energía de un tramo de instantáneas de una corrida."""
    import h5py
    from aa_numerico import MapaAA
    nombre, primero, pasos = arg
    eq = np.load(ruta('ic', f'{ic_de(nombre)}_equilibrio.npz'))
    mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']), fondo='puntual')
    t, h1, dphi, E = [], [], [], []
    with h5py.File(ruta(nombre, 'vlasov_output.h5'), 'r') as f:
        rg = f['grid']['r'][:]
        fondo = -1.0/rg + np.interp(rg, eq['r'], eq['phi_self'])
        w = f[primero]['f'][:]
        for c in pasos:
            g = f[c]
            Q, J, _ = mapa(g['r_part'][:], g['p_part'][:])
            B = J**2*np.exp(-(J - J1)**2/SJ1**2)
            t.append(g.attrs['time']); E.append(g.attrs['total_energy'])
            h1.append(np.sum(w*B*np.exp(-1j*Q)))
            dphi.append(g['potential'][:] - fondo)
    return t, h1, dphi, E


def serie(nombre, procesos=1):
    """t, h_1 (en el mapa del equilibrio), dPhi(r, t) en la malla, y la energía total;
    se guarda en exe/hadzic/<nombre>/serie.npz. Con procesos > 1 las instantáneas se
    reparten en tramos, uno por proceso; el resultado es el mismo."""
    import h5py
    from aa_numerico import MapaAA
    sal = ruta(nombre, 'serie.npz')
    if os.path.exists(sal):
        return np.load(sal)
    with h5py.File(ruta(nombre, 'vlasov_output.h5'), 'r') as f:
        pasos = sorted([c for c in f if c.startswith('step_')], key=lambda c: int(c.split('_')[1]))
        rg = f['grid']['r'][:]
        w = f[pasos[0]]['f'][:]
    n = max(1, min(procesos, len(pasos)))
    tramos = [(nombre, pasos[0], [str(c) for c in tr]) for tr in np.array_split(np.array(pasos), n)]
    if n > 1:
        from multiprocessing import Pool
        with Pool(n) as pool:
            res = pool.map(_tramo, tramos)
    else:
        res = [_tramo(tramos[0])]
    t, h1, dphi, E = (sum((list(r[i]) for r in res), []) for i in range(4))
    eq = np.load(ruta('ic', f'{ic_de(nombre)}_equilibrio.npz'))
    mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']), fondo='puntual')
    Q0, J0, _ = mapa(*np.loadtxt(ruta('ic', f'{ic_de(nombre)}.dat'), usecols=(0, 1), unpack=True))
    norma = np.sum(w*J0**2*np.exp(-(J0 - J1)**2/SJ1**2))
    np.savez(sal, t=np.array(t), h1=np.array(h1)/norma, dphi=np.array(dphi), r=rg, E=np.array(E))
    return np.load(sal)


def series_de(nombres, procesos=4):
    """Paso serie: serie.npz de las corridas pedidas, una detrás de otra."""
    for nombre in nombres:
        if not os.path.exists(ruta(nombre, 'vlasov_output.h5')):
            print(f'{nombre}: no hay vlasov_output.h5', flush=True)
            continue
        t0 = time.time()
        serie(nombre, procesos)
        print(f'{nombre}: serie.npz ({time.time()-t0:.0f} s)', flush=True)


def _serie(nombre):
    return nombre, dict(serie(nombre))


def frecuencia(t, x, a):
    """Frecuencia del máximo del periodograma (ventana de Hann) de x en t >= a."""
    s = t >= a
    om = np.linspace(0.005, 0.6, 11901)
    return om[np.argmax(np.abs(np.exp(-1j*np.outer(om, t[s])) @ (x[s]*np.hanning(s.sum()))))]


def modos():
    """omega_d de lambda.txt para cada (k, a0); None si no hay modo."""
    res = {}
    for l in open(ruta('lineal', 'lambda.txt')).readlines()[2:]:
        p = l.split()
        if len(p) == 10:
            res[(float(p[0]), float(p[1]))] = None if p[8] == '--' else float(p[8])
    return res


def analizar(procesos=4):
    from multiprocessing import Pool
    hechas = [n for n, *_ in CORRIDAS if os.path.exists(ruta(f'{n}.ok'))]
    with Pool(procesos) as pool:                 # una corrida por proceso
        series = dict(pool.map(_serie, hechas))
    wd = modos()
    lineas, polos = [], {}
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Perturbación = (D - Z)/eps. Envolventes: máximo en ventanas de 100 de |h_1| y de ||dPhi_eps|| '
      '(rms en la zona de la cáscara), PIC y lineal; pendiente de ln(envolvente) frente a ln t en las '
      'cuatro últimas ventanas; frecuencia tardía (segunda mitad), PIC y lineal; omega_d de lambda.')
    for nombre, k, a0, eps, ref in CORRIDAS:
        if not ref or nombre not in series or ref not in series:
            continue
        d, z = series[nombre], series[ref]
        lin = np.load(ruta('lineal', nombre_lin(k, a0) + '.npz'))
        t, tl = d['t'], lin['t']
        h = (d['h1'] - z['h1'])/eps
        zr = ZONA[a0]
        zona = (d['r'] >= zr[0]) & (d['r'] <= zr[1])
        nphi = np.sqrt(np.mean(((d['dphi'] - z['dphi'])/eps)[:, zona]**2, axis=1))
        zl = (lin['r'] >= zr[0]) & (lin['r'] <= zr[1])
        nlin = np.sqrt(np.mean(lin['dphi'][:, zl]**2, axis=1))
        T = TFIN[a0]
        ventanas = (0, T/16, T/8, T/4, T/2, T - 100)
        env = lambda x, tt, a: np.max(np.abs(x[(tt >= a) & (tt < a + 100)]))
        pend = lambda x, tt: np.polyfit(np.log(np.array(ventanas[2:]) + 50),
                                        np.log([env(x, tt, a) for a in ventanas[2:]]), 1)[0]
        dE = max(np.max(np.abs(d['E']/d['E'][0] - 1)), np.max(np.abs(z['E']/z['E'][0] - 1)))
        eq = np.load(ruta('ic', f'{ic_de(nombre)}_equilibrio.npz'))
        om_min = float(np.gradient(eq['E_t'], eq['J_t'])[np.searchsorted(eq['J_t'], JT)])
        m = wd.get((k, a0))
        w(f'{nombre}: k = {k:g}, a0 = {a0:g}, eps = {eps:g}   max|dE/E| = {dE:.1e}   Omega_min = {om_min:.5f}'
          f'   omega_d = {"--" if m is None else f"{m:.5f}"}')
        w(f'   pendientes: |h_1| PIC {pend(h, t):+.2f}, lineal {pend(lin["h1"], tl):+.2f};  ||dPhi|| PIC '
          f'{pend(nphi, t):+.2f}, lineal {pend(nlin, tl):+.2f};  frecuencia tardía de h_1: PIC '
          f'{frecuencia(t, h, T/2):.5f}, lineal {frecuencia(tl, lin["h1"], T/2):.5f}')
        for a in ventanas:
            w(f'   t {a:5.0f}-{a + 100:5.0f}: |h_1| PIC {env(h, t, a):.2e}  lineal {env(lin["h1"], tl, a):.2e}'
              f'   ||dPhi|| PIC {env(nphi, t, a):.2e}  lineal {env(nlin, tl, a):.2e}')
        if a0 == 1.0:                 # un polo (matrix pencil, M = 3) antes de que domine el ruido
            from landau_cola import ajustar
            for lo, hi in ((300, 1000), (300, 1500)) if k < 2 else ((100, 400), (100, 600)):
                wp, gp = ajustar(t, h, lo, hi)['pencil M=3']
                wl, gl = ajustar(tl, lin['h1'], lo, hi)['pencil M=3']
                polos[nombre, lo, hi] = (wp, gp)
                w(f'   polo en [{lo}, {hi}]: PIC omega = {wp:.5f}, gamma = {gp:+.1e};  lineal omega = '
                  f'{wl:.5f}, gamma = {gl:+.1e}')
    # Series 4 y 5: el polo de cada variante junto al de la corrida base, y el ruido de la
    # referencia.
    w('\nSeries 4 y 5 (a0 = 1): polo de (D - Z)/eps en [300, 1000] y [300, 1500], y |h_1| de la '
      'referencia Z (máximo en ventanas de 100).')
    w(f'{"corrida":>15} {"k":>5} {"eps":>5} {"cambio":>36} {"omega":>8} {"gamma":>8} {"omega":>8} '
      f'{"gamma":>8} {"|Z| t=1000":>10} {"|Z| t=3900":>10}')
    datos = dict((n, (e, r)) for n, _, _, e, r in CORRIDAS)
    for k in [1.25, 1.5]:
        b = f'k{k:g}_a1'
        for nombre in (f'D_{b}', f'De03_{b}', f'De3_{b}', f'Ddr_{b}', f'DN_{b}', f'DP_{b}'):
            if (nombre, 300, 1500) not in polos:
                continue
            eps, ref = datos[nombre]
            c = ', '.join(f'{x} = {v}' for x, v in CAMBIOS.get(nombre, {}).items() if x != 'ic') or '--'
            zt, zh = series[ref]['t'], np.abs(series[ref]['h1'])
            zv = [zh[(zt >= a) & (zt < a + 100)].max() for a in (1000, 3900)]
            (w1, g1), (w2, g2) = polos[nombre, 300, 1000], polos[nombre, 300, 1500]
            w(f'{nombre:>15} {k:5g} {eps:5g} {c:>36} {w1:8.5f} {g1:+8.1e} {w2:8.5f} {g2:+8.1e} '
              f'{zv[0]:10.1e} {zv[1]:10.1e}')
    open(ruta('resumen.txt'), 'w').write('\n'.join(lineas) + '\n')


def envolvente(t, x, ancho=100.0):
    """Máximo de |x| en |t' - t| <= ancho/2."""
    x = np.abs(x)
    return np.array([x[(t >= a - ancho/2) & (t <= a + ancho/2)].max() for a in t])


def figuras():
    """Figuras del informe en docs/hadzic/figuras/."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    destino = os.path.join(RAIZ, 'docs', 'hadzic', 'figuras')
    os.makedirs(destino, exist_ok=True)
    plt.rcParams.update({'font.size': 9})
    # 1. lambda_edge(a0) para cada k; con k <= 1, lambda a 1e-6 anchos de banda del borde.
    tabla = {}
    for l in open(ruta('lineal', 'lambda.txt')).readlines()[2:]:
        q = l.split()
        if len(q) == 10:
            tabla[(float(q[0]), float(q[1]))] = (float(q[5]), float(q[7]))
    fig, ax = plt.subplots(figsize=(5.2, 3.6), constrained_layout=True)
    for k, c in zip(K_MAPA, ['C3', 'C1', 'C2', 'C0', 'C4', 'C5']):
        a = np.array(A0_MAPA)
        if k <= 1:
            ax.loglog(a, [tabla[(k, x)][1] for x in a], '--o', color=c, ms=3,
                      label=f'$k={k:g}$: $\\lambda$ at $\\delta=10^{{-6}}$ (diverges at the edge)')
        else:
            ax.loglog(a, [tabla[(k, x)][0] for x in a], '-o', color=c, ms=3, label=f'$k={k:g}$')
    ax.axhline(1, color='k', lw=0.8)
    ax.set_xlabel('$a_0$'); ax.set_ylabel('$\\lambda_{\\rm edge}$')
    ax.legend(fontsize=7, frameon=False); ax.grid(alpha=0.3, which='both')
    fig.savefig(os.path.join(destino, 'lambda_mapa.pdf')); plt.close(fig)
    # 2. a0 = 0.01: |h_1| lineal y PIC con eps = 0.1 y 0.5.
    fig, axs = plt.subplots(2, 2, figsize=(7, 5), constrained_layout=True, sharex=True, sharey=True)
    for ax, k in zip(axs.flat, [0.75, 1.0, 1.5, 2.0]):
        lin = np.load(ruta('lineal', nombre_lin(k, 0.01) + '.npz'))
        ax.loglog(lin['t'][1:], envolvente(lin['t'], lin['h1'])[1:], 'k', lw=1.2, label='linear')
        z = np.load(ruta(f'Z_k{k:g}', 'serie.npz'))
        for nombre, eps, c in ((f'D_k{k:g}', 0.1, 'C0'), (f'D5_k{k:g}', 0.5, 'C1')):
            d = np.load(ruta(nombre, 'serie.npz'))
            ax.loglog(d['t'][1:], envolvente(d['t'], (d['h1'] - z['h1'])/eps)[1:], color=c, lw=0.9,
                      label=f'PIC, $\\varepsilon={eps:g}$')
        ax.set_title(f'$k={k:g}$, $a_0=0.01$'); ax.grid(alpha=0.3, which='both')
        ax.set_ylim(1e-8, 1)
    for ax in axs[1]: ax.set_xlabel('$t$')
    for ax in axs[:, 0]: ax.set_ylabel('$|h_1|$ (envelope)')
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.savefig(os.path.join(destino, 'pic_a001.pdf')); plt.close(fig)
    # 3. a0 = 1: |h_1| lineal y PIC, y la referencia sola (ruido de la corrida sin perturbar).
    fig, axs = plt.subplots(2, 3, figsize=(8, 5), constrained_layout=True, sharex=True, sharey=True)
    for ax, k in zip(axs.flat, [0.75, 1.0, 1.25, 1.5, 2.0]):
        lin = np.load(ruta('lineal', nombre_lin(k, 1.0) + '.npz'))
        d, z = np.load(ruta(f'D_k{k:g}_a1', 'serie.npz')), np.load(ruta(f'Z_k{k:g}_a1', 'serie.npz'))
        ax.semilogy(lin['t'], envolvente(lin['t'], lin['h1']), 'k', lw=1.2, label='linear')
        ax.semilogy(d['t'], envolvente(d['t'], (d['h1'] - z['h1'])/0.1), 'C0', lw=0.9, label='PIC, $(D-Z)/\\varepsilon$')
        ax.semilogy(z['t'], envolvente(z['t'], z['h1']/0.1), color='0.6', lw=0.8, label='$Z/\\varepsilon$ (noise)')
        ax.set_title(f'$k={k:g}$, $a_0=1$'); ax.grid(alpha=0.3); ax.set_ylim(1e-6, 1)
    axs.flat[-1].axis('off')
    for ax in axs[1]: ax.set_xlabel('$t$')
    for ax in axs[:, 0]: ax.set_ylabel('$|h_1|$ (envelope)')
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.savefig(os.path.join(destino, 'pic_a1.pdf')); plt.close(fig)
    # 4. Serie 4, k = 1.25 y 1.5 con a0 = 1: arriba, la amplitud (eps = 0.03, 0.1, 0.3); abajo, la
    #    malla y las partículas (dr/2, 4N), con el ruido de las referencias Z de la base y de 4N.
    fig, axs = plt.subplots(2, 2, figsize=(7, 5), constrained_layout=True, sharex=True, sharey=True)
    for col, k in enumerate([1.25, 1.5]):
        b = f'k{k:g}_a1'
        lin = np.load(ruta('lineal', nombre_lin(k, 1.0) + '.npz'))
        filas = [((f'De03_{b}', f'Z_{b}', 0.03, 'C1', '$\\varepsilon=0.03$'),
                  (f'D_{b}', f'Z_{b}', 0.1, 'C0', '$\\varepsilon=0.1$ (base)'),
                  (f'De3_{b}', f'Z_{b}', 0.3, 'C2', '$\\varepsilon=0.3$')),
                 ((f'D_{b}', f'Z_{b}', 0.1, 'C0', 'base'),
                  (f'Ddr_{b}', f'Zdr_{b}', 0.1, 'C3', '$\\Delta r/2$'),
                  (f'DN_{b}', f'ZN_{b}', 0.1, 'C4', '$4N$'))]
        for fila, curvas in enumerate(filas):
            ax = axs[fila, col]
            ax.semilogy(lin['t'], envolvente(lin['t'], lin['h1']), 'k', lw=1.4, label='linear')
            for nombre, ref, eps, c, et in curvas:
                d, z = np.load(ruta(nombre, 'serie.npz')), np.load(ruta(ref, 'serie.npz'))
                ax.semilogy(d['t'], envolvente(d['t'], (d['h1'] - z['h1'])/eps), color=c, lw=1.0, label=et)
            if fila == 1:
                for ref, ls, et in ((f'Z_{b}', '-', '$Z/\\varepsilon$, base'), (f'ZN_{b}', '--', '$Z/\\varepsilon$, $4N$')):
                    z = np.load(ruta(ref, 'serie.npz'))
                    ax.semilogy(z['t'], envolvente(z['t'], z['h1']/0.1), color='0.55', ls=ls, lw=0.8, label=et)
            ax.set_title(f'$k={k:g}$, $a_0=1$'); ax.grid(alpha=0.3); ax.set_xlim(0, 4000)
            ax.set_ylim(5e-3, 0.4)
    for ax in axs[1]: ax.set_xlabel('$t$')
    for ax in axs[:, 0]: ax.set_ylabel('$|h_1|$ (envelope)')
    for ax, donde, nc in ((axs[0, 0], 'lower right', 1), (axs[1, 0], 'upper center', 3)):
        for l in ax.legend(fontsize=7, frameon=False, loc=donde, ncol=nc).get_lines():
            l.set_linewidth(1.6)
    fig.savefig(os.path.join(destino, 'pic_barrido.pdf')); plt.close(fig)
    print('figuras en', destino)


# ------------------------------------------------------------------ videos
# La perturbación se dibuja en ángulo-acción, a partir de sus armónicos en Q. En (r, p_r) la
# diferencia de histogramas de D y Z está dominada por la red inicial, que se estira en
# espirales finas: al separarse un poco las partículas de D y Z, la diferencia de celda a celda
# es mucho mayor que la señal. En (Q, J) cada anillo de la red tiene 25 partículas igualmente
# espaciadas en Q que giran juntas, y sus armónicos k < 12 casi no tienen ruido (como h_1).
KMAX, NJC, SJ = 6, 100, 0.007            # armónicos, centros en J y ancho del núcleo en J


def armonicos(nombre):
    """f_k(J_c, t) = (dQ dJ/2pi) sum_j F_j e^{-ikQ_j} K(J_j - J_c)/N(J_c), k = 0..KMAX, con K una
    gaussiana de ancho SJ y N su integral en [0, J_max] (corrige el borde J = 0), de modo que
    f(Q, J) = f_0 + 2 Re sum_{k>=1} f_k e^{ikQ}. Se guarda en exe/hadzic/<nombre>/armonicos.npz."""
    import h5py
    from scipy.special import erf
    from aa_numerico import MapaAA
    sal = ruta(nombre, 'armonicos.npz')
    if os.path.exists(sal):
        return np.load(sal)
    eq = np.load(ruta('ic', f'{ic_de(nombre)}_equilibrio.npz'))
    mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']), fondo='puntual')
    jmx = float(eq['J_max'])
    dQdJ = 2*np.pi/int(eq['npc'])*jmx/int(eq['nrc'])
    Jc = (np.arange(NJC) + 0.5)*JT/NJC
    norma = 0.5*(erf((jmx - Jc)/(np.sqrt(2)*SJ)) + erf(Jc/(np.sqrt(2)*SJ)))
    kk = np.arange(KMAX + 1)
    f = h5py.File(ruta(nombre, 'vlasov_output.h5'), 'r')
    pasos = sorted([c for c in f if c.startswith('step_')], key=lambda c: int(c.split('_')[1]))
    w = f[pasos[0]]['f'][:]
    t, fk = [], np.empty((len(pasos), KMAX + 1, NJC), np.complex64)
    for i, c in enumerate(pasos):
        g = f[c]
        Q, J, _ = mapa(g['r_part'][:], g['p_part'][:])
        K = np.exp(-0.5*((J[:, None] - Jc[None, :])/SJ)**2)/(np.sqrt(2*np.pi)*SJ)
        fk[i] = dQdJ/(2*np.pi)*((np.exp(-1j*np.outer(kk, Q))*w[None, :]) @ K)/norma[None, :]
        t.append(g.attrs['time'])
    f.close()
    np.savez(sal, t=np.array(t), Jc=Jc, fk=fk, Fmax=w.max())
    return np.load(sal)


def _armonicos(nombre):
    t0 = time.time()
    armonicos(nombre)
    return nombre, f'armónicos en {time.time()-t0:.0f} s'


def borde_rp(eq):
    """Órbita del borde del soporte, E = E(J_t), en (r, p_r >= 0), y el radio de la circular."""
    r = np.linspace(1.0, 20.0, 8000)                        # masa puntual para r > 1
    phi = -1.0/r + np.interp(r, eq['r'], eq['phi_self']) + L0**2/(2*r**2)
    Eb = np.interp(JT, eq['J_t'], eq['E_t'])
    d = Eb > phi
    return r[d], np.sqrt(2*(Eb - phi[d])), r[np.argmin(phi)]


class Fotogramas:
    """Fotogramas del video de una corrida D: sus partículas en (r, p_r) coloreadas por F, el
    armónico k = 1 de la perturbación, dF_1(Q, J) = 2 Re[(f_D - f_Z)_1(J) e^{iQ}]/eps, y |h_1(t)|.
    F y dF van en unidades del máximo de F_eq de los pesos del código (que reescala los del dato
    inicial por un factor constante, 0.7997 en a0 = 1).

    Solo k = 1: el dato inicial es puro k = 1 y es lo que mide h_1. Con 10^4 partículas, en
    D_k2_a1 a t = 1500 los armónicos 0 y 2 tienen un rms en J de 0.95 y 1.1 veces el máximo
    inicial de k = 1, casi todo ruido, y sumarlos lo dobla."""

    def __init__(self, nombre, k, a0, eps, ref, nq=128):
        import h5py
        self.nombre, self.k, self.a0, self.eps, self.ref = nombre, k, a0, eps, ref
        eq = np.load(ruta('ic', f'{ic_de(nombre)}_equilibrio.npz'))
        self.rb, self.pb, self.rc = borde_rp(eq)
        self.r0, self.r1, self.p1 = self.rb.min() - 0.3, self.rb.max() + 0.3, 1.12*self.pb.max()
        self.fD = h5py.File(ruta(nombre, 'vlasov_output.h5'), 'r')
        self.pasos = sorted([c for c in self.fD if c.startswith('step_')], key=lambda c: int(c.split('_')[1]))
        self.w = self.fD[self.pasos[0]]['f'][:]
        aD, aZ = armonicos(nombre), armonicos(ref)
        self.Fmax, self.Jc = float(aZ['Fmax']), aZ['Jc']
        self.dfk = (aD['fk'] - aZ['fk'])/eps/self.Fmax            # (t, k, J)
        self.Q = (np.arange(nq) + 0.5)*2*np.pi/nq
        self.base = 2*np.exp(1j*self.Q)                            # k = 1
        d, z = np.load(ruta(nombre, 'serie.npz')), np.load(ruta(ref, 'serie.npz'))
        self.t, self.h, self.hz = d['t'], np.abs((d['h1'] - z['h1'])/eps), np.abs(z['h1'])/eps
        lin = np.load(ruta('lineal', nombre_lin(k, a0) + '.npz'))
        self.tl, self.hl = lin['t'], np.abs(lin['h1'])

    def dF(self, i):
        """dF_1(J, Q) en el fotograma i."""
        return np.real(self.dfk[i, 1][:, None]*self.base[None, :])

    def figura(self, lim):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10})
        fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
        gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1], hspace=0.34, wspace=0.26,
                              left=0.06, right=0.945, top=0.87, bottom=0.08)
        a1, a2, a3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
        g = self.fD[self.pasos[0]]
        self.orden = np.argsort(self.w)                         # los de F grande, encima
        s = 0.6 if len(self.w) <= 10000 else 0.15
        self.pts = a1.scatter(g['r_part'][:][self.orden], g['p_part'][:][self.orden],
                              c=self.w[self.orden]/self.Fmax, s=s, cmap='viridis', vmin=0, vmax=1,
                              linewidths=0, rasterized=True)
        fig.colorbar(self.pts, ax=a1, fraction=0.04, pad=0.02, label='$F/F_{\\max}$')
        for sg in (1, -1):
            a1.plot(self.rb, sg*self.pb, color='0.4', lw=0.7)
        a1.plot([self.rc], [0], '+', color='0.3', ms=6)
        a1.set_xlim(self.r0, self.r1); a1.set_ylim(-self.p1, self.p1)
        a1.set_xlabel('$r$'); a1.set_ylabel('$p_r$')
        a1.set_title(f'particles of run D ({len(self.w)}); grey: edge of the support', fontsize=10)
        self.img = a2.imshow(self.dF(0), origin='lower', aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim,
                             extent=(0, 2*np.pi, 0, JT), interpolation='bilinear')
        fig.colorbar(self.img, ax=a2, fraction=0.04, pad=0.02, label='$\\delta f_1/(\\varepsilon F_{\\max})$')
        a2.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
        a2.set_xticklabels(['0\n(pericentre)', '$\\pi/2$', '$\\pi$\n(apocentre)', '$3\\pi/2$', '$2\\pi$'])
        a2.set_xlabel('angle $Q$'); a2.set_ylabel('action $J$')
        a2.set_title(f'harmonic $k=1$ of the perturbation in angle–action variables '
                     f'(smoothed over $\\Delta J={SJ}$)', fontsize=10)
        a3.semilogy(self.tl, self.hl, 'k', lw=1.2, label='linear theory')
        a3.semilogy(self.t, self.hz, color='0.65', lw=0.8, label='reference alone, $|h_{1,Z}|/\\varepsilon$')
        a3.semilogy(self.t, self.h, 'C0', lw=0.9, label='PIC, $|h_{1,D}-h_{1,Z}|/\\varepsilon$')
        self.cursor = a3.axvline(0, color='C3', lw=1)
        a3.set_xlim(0, self.t[-1])
        a3.set_ylim(max(1e-7, 0.3*min(self.h[1:].min(), self.hl[1:].min())), 1)
        a3.set_xlabel('$t$'); a3.set_ylabel('$|h_1|$'); a3.grid(alpha=0.3)
        a3.legend(fontsize=8, frameon=False, loc='upper right', ncol=3)
        wd = modos().get((self.k, self.a0))
        teoria = ('no discrete mode, the perturbation damps' if wd is None else
                  f'discrete mode below the band at $\\omega_d={wd:.5f}$')
        c = ', '.join(f'{x} = {v}' for x, v in CAMBIOS.get(self.nombre, {}).items() if x != 'ic')
        fig.suptitle(f'Hadžić setting, run {self.nombre}:  edge exponent $k={self.k:g}$,  shell mass '
                     f'$a_0={self.a0:g}$,  $\\varepsilon={self.eps:g}$' + (f'  ({c})' if c else ''), fontsize=12)
        fig.text(0.5, 0.915, f'linear theory: {teoria}', ha='center', fontsize=10, color='0.25')
        self.reloj = fig.text(0.06, 0.915, '', fontsize=11, family='monospace')
        return fig

    def poner(self, i):
        g = self.fD[self.pasos[i]]
        self.pts.set_offsets(np.c_[g['r_part'][:][self.orden], g['p_part'][:][self.orden]])
        self.img.set_data(self.dF(i))
        t = g.attrs['time']
        self.cursor.set_xdata([t, t])
        self.reloj.set_text(f't = {t:6.0f}')


def video(caso, fps=25):
    """Video de una corrida D en exe/hadzic/videos/<nombre>.mp4 (H.264, 1280 x 720, 25 fps)."""
    from matplotlib.animation import FFMpegWriter
    nombre = caso[0]
    sal = ruta('videos', f'{nombre}.mp4')
    if os.path.exists(sal):
        return nombre, 'ya estaba'
    t0 = time.time()
    v = Fotogramas(*caso)
    lim = 0.6*np.abs(v.dF(0)).max()                 # escala fija: lo que se amortigua, se apaga
    fig = v.figura(lim)
    escritor = FFMpegWriter(fps=fps, codec='libx264', bitrate=-1,
                            extra_args=['-pix_fmt', 'yuv420p', '-crf', '20', '-preset', 'slow'])
    with escritor.saving(fig, sal + '.tmp.mp4', dpi=100):
        for i in range(len(v.pasos)):
            v.poner(i); escritor.grab_frame()
    os.replace(sal + '.tmp.mp4', sal)
    return nombre, f'{len(v.pasos)} fotogramas en {time.time()-t0:.0f} s'


def _video(caso):
    try:
        return video(caso)
    except Exception as e:
        return caso[0], f'FALLO: {e!r}'


def videos(procesos=4):
    """Un video por corrida D hecha, con su referencia: primero los armónicos de todas las
    corridas (en caché) y después los videos, cuatro procesos a la vez."""
    from multiprocessing import Pool
    os.makedirs(ruta('videos'), exist_ok=True)
    casos = [c for c in CORRIDAS if c[4] and os.path.exists(ruta(f'{c[0]}.ok'))
             and os.path.exists(ruta(f'{c[4]}.ok'))]
    nombres = sorted({n for c in casos for n in (c[0], c[4])},
                     key=lambda n: -os.path.getsize(ruta(n, 'vlasov_output.h5')))   # las grandes primero
    with Pool(procesos) as pool:
        for nombre, msg in pool.imap_unordered(_armonicos, nombres):
            print(f'{nombre:16} {msg}', flush=True)
        for nombre, msg in pool.imap_unordered(_video, casos):
            print(f'{nombre:16} {msg}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['lineal', 'preparar', 'correr', 'serie', 'analizar', 'figuras', 'videos'])
    ap.add_argument('nombres', nargs='*', help='corridas del paso serie')
    ap.add_argument('-j', '--procesos', type=int, default=4, help='procesos del paso serie')
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    if a.paso == 'serie':
        series_de(a.nombres, a.procesos)
    else:
        globals()[a.paso]()
