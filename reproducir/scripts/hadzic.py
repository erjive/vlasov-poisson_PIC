"""El escenario de Hadžić, Rein, Schrecker y Straub (ARMA 249, 45, 2025): |L| fijo, masa
puntual de masa 1 en el centro y politropos en la energía F = A (E_t - E)^k, con el borde
E_t = E(J_t). El teorema dice que, con masa pequeña, las perturbaciones se amortiguan si
k > 1 y no si 1/2 < k <= 1. Aquí se estudia con la teoría lineal y con el código PIC.

Pasos (cada uno reutiliza lo que ya exista):
    python3 hadzic.py lineal     mapa de lambda_edge(k, a0) y soluciones lineales en el tiempo
    python3 hadzic.py preparar   estados iniciales (equilibrio.py) y .par de las corridas
    python3 hadzic.py correr     corridas PIC, una detrás de otra, 4 hilos
    python3 hadzic.py analizar   dPhi y h_1 en el mapa del equilibrio, menos eps = 0

Parámetros: L0 = 2 (su L = 4), J_t = 0.7, es decir E_t = -0.0686 sin masa propia, dentro
de su condición de un solo hueco (-0.079 < E_t < 0; Omega_max/Omega_min = 2.46). En el
código, BGtype = "sphere" es una bola uniforme de masa 1 y radio 1, es decir una masa
puntual para r > 1, y la cáscara está en 2.4 < r < 12.2. Salidas en exe/hadzic/.
"""
import os, sys, subprocess, time, argparse, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
EXE = os.path.join(RAIZ, 'exe')
BASE = os.path.join(EXE, 'hadzic')
PARDIR = os.path.join(RAIZ, 'reproducir', 'corridas', '13_hadzic')
PLANTILLA = os.path.join(RAIZ, 'reproducir', 'corridas', '11_landau', 'landau__L_a1e-2_n400_e0.1.par')
JT, L0 = 0.7, 2.0
J1, SJ1 = 0.35, 0.20                    # función de prueba de h_1: B = J^2 exp(-(J-J1)^2/SJ1^2)
# courant = 2 (el de 11_landau): dt = courant dr/pmax = 0.1; una instantánea cada 200 pasos,
# es decir cada 20 unidades de tiempo; t_fin = 8000, unos 95 tau_1.
NRC, NPC, COURANT, DT, SALIDA, TFIN = 400, 25, 2.0, 0.1, 200, 8000.0
ZONA = (3.0, 12.0)                      # radios de la cáscara para ||dPhi|| (la solución lineal da 3-15)
# Mapa lineal: exponentes del borde y masas.
K_MAPA = [0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
A0_MAPA = [0.01, 0.03, 0.1, 0.3, 0.6, 1.0]
# Corridas PIC: a0 = 0.01; D con eps = 0.1 y su referencia Z con eps = 0.
K_PIC = [0.75, 1.0, 1.5, 2.0]
A0_PIC, EPS = 0.01, 0.1
CORRIDAS = [(f'{p}_k{k:g}', k, e) for k in K_PIC for p, e in (('D', EPS), ('Z', 0.0))]


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
    w(f'Masa puntual, polE, J_t = {JT}, L0 = {L0}. delta en anchos de banda bajo Omega_min.')
    w(f'{"k":>5} {"a0":>6} {"iter":>5} {"Omega_min":>10} {"Omax/Omin":>9} {"lambda_edge":>12} '
      f'{"l(1e-3)":>8} {"l(1e-6)":>8} {"omega_d":>9} {"x_d":>8}')
    for k in K_MAPA:
        for a0 in A0_MAPA:
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
    open(ruta('lineal', 'lambda.txt'), 'w').write('\n'.join(lineas) + '\n')
    for k in K_PIC:
        sal = ruta('lineal', f'lin_k{k:g}.npz')
        if os.path.exists(sal):
            continue
        t0 = time.time()
        t, h1, h2, rmed, dphi = resolver(ruta('lineal', f'P_k{k:g}_a{A0_PIC:g}_equilibrio.npz'),
                                         1600, 32, 0.5, TFIN, verboso=False, j1=J1, sj1=SJ1)
        np.savez(sal, t=t, h1=h1, h2=h2, r=rmed, dphi=dphi)
        print(f'  lineal k = {k:g} hasta t = {TFIN:g} ({time.time()-t0:.0f} s)', flush=True)


# ------------------------------------------------------------------ preparar
def preparar():
    os.makedirs(ruta('ic'), exist_ok=True); os.makedirs(PARDIR, exist_ok=True)
    plantilla = open(PLANTILLA).read().split('\n')
    for nombre, k, eps in CORRIDAS:
        dat = ruta('ic', f'{nombre}.dat')
        if not os.path.exists(dat):
            t0 = time.time()
            equilibrio(dat, k, A0_PIC, eps)
            print(f'  estado inicial {nombre} ({time.time()-t0:.0f} s)', flush=True)
        valores = {'courant': str(COURANT), 'Nt': str(int(round(TFIN/DT))), 'time_output': str(SALIDA),
                   'spatial_output': str(SALIDA), 'field_output': str(SALIDA),
                   'Nrc': str(NRC), 'Npc': str(NPC), 'Lfix': str(L0),
                   'directory': f'hadzic/{nombre}', 'a0': str(A0_PIC), 'BGtype': 'sphere',
                   'checkpointfile': f'hadzic/ic/{nombre}.dat', 'j1': f'{J1:.6f}',
                   'sj1': f'{SJ1:.6f}', 'state': 'checkpoint'}
        lineas = [f'# Escenario de Hadžić, corrida {nombre}: polE con k = {k:g}, J_t = {JT}, '
                  f'a0 = {A0_PIC}, eps = {eps}, t_fin = {TFIN:g}.',
                  '# Masa puntual: BGtype = sphere (masa 1, radio 1). Generado por '
                  'reproducir/scripts/hadzic.py a partir de 11_landau.']
        for l in plantilla:
            if l.startswith('#') or '=' not in l:
                continue
            clave = l.split('=')[0].strip()
            lineas.append(f'{clave:<16} = {valores[clave]}' if clave in valores else l)
        open(os.path.join(PARDIR, f'hadzic__{nombre}.par'), 'w').write('\n'.join(lineas) + '\n')
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
def serie(nombre):
    """t, h_1 (en el mapa del equilibrio), dPhi(r, t) en la malla, y la energía total;
    se guarda en exe/hadzic/<nombre>/serie.npz."""
    import h5py
    from aa_numerico import MapaAA
    sal = ruta(nombre, 'serie.npz')
    if os.path.exists(sal):
        return np.load(sal)
    eq = np.load(ruta('ic', f'{nombre}_equilibrio.npz'))
    mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']), fondo='puntual')
    f = h5py.File(ruta(nombre, 'vlasov_output.h5'), 'r')
    pasos = sorted([c for c in f if c.startswith('step_')], key=lambda c: int(c.split('_')[1]))
    rg = f['grid']['r'][:]
    fondo = -1.0/rg + np.interp(rg, eq['r'], eq['phi_self'])
    w = f[pasos[0]]['f'][:]
    t, h1, dphi, E = [], [], [], []
    for c in pasos:
        g = f[c]
        Q, J, _ = mapa(g['r_part'][:], g['p_part'][:])
        B = J**2*np.exp(-(J - J1)**2/SJ1**2)
        t.append(g.attrs['time']); E.append(g.attrs['total_energy'])
        h1.append(np.sum(w*B*np.exp(-1j*Q)))
        dphi.append(g['potential'][:] - fondo)
    f.close()
    Q0, J0, _ = mapa(*np.loadtxt(ruta('ic', f'{nombre}.dat'), usecols=(0, 1), unpack=True))
    norma = np.sum(w*J0**2*np.exp(-(J0 - J1)**2/SJ1**2))
    np.savez(sal, t=np.array(t), h1=np.array(h1)/norma, dphi=np.array(dphi), r=rg, E=np.array(E))
    return np.load(sal)


def _serie(nombre):
    return nombre, dict(serie(nombre))


def analizar(procesos=4):
    from multiprocessing import Pool
    os.environ['OMP_NUM_THREADS'] = '1'
    with Pool(procesos) as pool:                 # una corrida por proceso
        series = dict(pool.map(_serie, [n for n, *_ in CORRIDAS]))
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w(f'Corridas PIC (a0 = {A0_PIC}, eps = {EPS}): envolventes (máximo en ventanas de 100) de |h_1| y '
      f'de ||dPhi_eps|| (rms en {ZONA[0]} <= r <= {ZONA[1]}), PIC y lineal; pendiente de ln(envolvente) '
      f'frente a ln t entre t = 1000 y 8000; error de la energía.')
    ventanas = (0, 500, 1000, 2000, 4000, 7900)
    for k in K_PIC:
        d, z = series[f'D_k{k:g}'], series[f'Z_k{k:g}']
        lin = np.load(ruta('lineal', f'lin_k{k:g}.npz'))
        t, tl = d['t'], lin['t']
        h = (d['h1'] - z['h1'])/EPS
        zona = (d['r'] >= ZONA[0]) & (d['r'] <= ZONA[1])
        nphi = np.sqrt(np.mean(((d['dphi'] - z['dphi'])/EPS)[:, zona]**2, axis=1))
        zl = (lin['r'] >= ZONA[0]) & (lin['r'] <= ZONA[1])
        nlin = np.sqrt(np.mean(lin['dphi'][:, zl]**2, axis=1))
        env = lambda x, tt, a: np.max(np.abs(x[(tt >= a) & (tt < a + 100)]))
        pend = lambda x, tt: np.polyfit(np.log(np.array(ventanas[2:]) + 50),
                                        np.log([env(x, tt, a) for a in ventanas[2:]]), 1)[0]
        dE = max(np.max(np.abs(d['E']/d['E'][0] - 1)), np.max(np.abs(z['E']/z['E'][0] - 1)))
        w(f'k = {k:g}   max|dE/E| = {dE:.1e}   pendientes: |h_1| PIC {pend(h, t):+.2f}, lineal '
          f'{pend(lin["h1"], tl):+.2f};  ||dPhi|| PIC {pend(nphi, t):+.2f}, lineal {pend(nlin, tl):+.2f}')
        for a in ventanas:
            w(f'   t {a:4d}-{a + 100:4d}: |h_1| PIC {env(h, t, a):.2e}  lineal {env(lin["h1"], tl, a):.2e}'
              f'   ||dPhi|| PIC {env(nphi, t, a):.2e}  lineal {env(nlin, tl, a):.2e}')
    open(ruta('resumen.txt'), 'w').write('\n'.join(lineas) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['lineal', 'preparar', 'correr', 'analizar'])
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    globals()[a.paso]()
