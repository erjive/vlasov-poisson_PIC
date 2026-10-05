"""Corridas de validación del artículo (docs/articulo, Sección 5), con el código auditado.

Dos bloques, todo con L0 = 2 y G = M = b = 1:

  libre  Mezcla sin autogravedad en el isócrono frente a la solución exacta, con el pulso del
         borrador anterior, F0 = exp(-sin^2(Q/2)/sQ^2) J^2 exp(-J^2/sJ^2), sQ = sJ = 0.1, en una
         retícula de N_J x N_Q nodos (state = aa_quad). La función de prueba es la propia F0.
         Barridos en N_J, en N_Q y en el paso de tiempo (courant = 4, 2, 1, 0.5 y el avance
         exacto, integrator = analytic), y una corrida larga hasta t = 10^4.
  paso   Paso de tiempo con autogravedad en el caso más exigente del informe de Hadžić: masa
         puntual, a0 = 1, k = 1 (modo discreto). Las corridas D (eps = 0.1) y Z (eps = 0) con
         courant = 4, 1 y 0.5, con el mismo dato inicial que las de courant = 2 de exe/hadzic.

En el código dt = courant dr/pmax con pmax = 2, así que courant = 2 es dt = 0.1.

    python3 validacion.py correr     las corridas, una detrás de otra, con 4 hilos
    python3 validacion.py analizar   tablas: exe/validacion/resumen.txt
    python3 validacion.py armonico   chi_1 del bloque paso en el mapa del equilibrio (4 procesos,
                                     unos minutos): exe/validacion/armonico.txt
    python3 validacion.py lineal     convergencia del solucionador lineal (dt/2, 2 N_Q, 2 N_J) en el
                                     politropo k = 1.25, a0 = 1: exe/validacion/lineal.txt (5 min)
    python3 validacion.py residuo    parte estática de h_1 con el mapa del isócrono y con el del
                                     potencial total (corridas de exe/sg): exe/validacion/residuo.txt
    python3 validacion.py respuesta  corridas de masa pequeña de exe/hadzic frente a la solución
                                     lineal: exe/validacion/respuesta.txt
    python3 validacion.py referencia corridas de referencia (eps = 0) de exe/hadzic con tres números
                                     de partículas y dr/2: exe/validacion/referencia.txt
    python3 validacion.py figuras    figuras del artículo (docs/articulo/figuras/)

Salidas en exe/validacion/. Las corridas hechas no se repiten.

Los demás números de la Sección 5 salen de otros pasos:
  - tres amplitudes con las mismas 1.024e5 partículas (kappa): hadzic.py analizar,
    exe/hadzic/resumen.txt;
  - frecuencia de la solución lineal frente a omega_d: hadzic.py analizar y
    exe/hadzic/lineal/lambda.txt;
  - ganancia con dispersión en el límite de L fijo: lambda_L.py validar;
  - escala con dr del error del campo (modelo de Wilson, a0 = 0.5): informe de la demo,
    docs/demo_eta, corridas E_M5*;
  - mezcla libre, estado estacionario con los dos mapas y colapso frío del código con L:
    VlasovPoisson_PIC_sp (reproducir/corridas/08_verificacion y 11_equilibrio, verificacion I2).
"""
import os, sys, subprocess, time
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
EXE = os.path.join(RAIZ, 'exe')
BASE = os.path.join(EXE, 'validacion')
BASE_LIBRE = os.path.join(RAIZ, 'reproducir', 'corridas', 'base', 'base_dftest.par')
BASE_PASO = os.path.join(EXE, 'hadzic', 'D_k1_a1', 'params_usados.par')
SJ, SQ = 0.1, 0.1                      # anchos del pulso y de la función de prueba
TFIN_LIBRE, CADA = 2000.0, 2.0         # h_k cada 2 unidades de tiempo
TFIN_PASO, CADA_PASO, EPS = 2000.0, 5.0, 0.1
PULSO = dict(dftype='gauss', sr=SJ, sp=SQ, state='aa_quad', j1=0.0, sj1=SJ, sq1=SQ)


def tiempo(courant, tfin, cada):
    """Nt y la cadencia de salida para dt = courant 0.1/2, con los mismos tiempos de salida."""
    dt = courant*0.1/2.0
    return dict(courant=courant, Nt=int(round(tfin/dt)), spatial_output=int(round(cada/dt)),
                time_output=int(round(cada/dt)))


CORRIDAS = {}                           # nombre -> (archivo base, cambios)
for nj in (25, 50, 100, 200, 400, 800):
    CORRIDAS[f'libre/nj{nj}'] = (BASE_LIBRE, dict(PULSO, Nrc=nj, Npc=64, **tiempo(2.0, TFIN_LIBRE, CADA)))
for nq in (16, 32):
    CORRIDAS[f'libre/nq{nq}'] = (BASE_LIBRE, dict(PULSO, Nrc=400, Npc=nq, **tiempo(2.0, TFIN_LIBRE, CADA)))
for c in (4.0, 1.0, 0.5):
    CORRIDAS[f'libre/c{c:g}'] = (BASE_LIBRE, dict(PULSO, Nrc=400, Npc=64, **tiempo(c, TFIN_LIBRE, CADA)))
CORRIDAS['libre/exacto'] = (BASE_LIBRE, dict(PULSO, Nrc=400, Npc=64, integrator='analytic',
                                             **tiempo(2.0, TFIN_LIBRE, CADA)))
CORRIDAS['libre/largo'] = (BASE_LIBRE, dict(PULSO, Nrc=800, Npc=64, field_output=1000,
                                            **tiempo(2.0, 1.0e4, CADA)))
for c in (4.0, 1.0, 0.5):
    for x in 'DZ':
        cambios = dict(tiempo(c, TFIN_PASO, CADA_PASO), checkpointfile=f'hadzic/ic/{x}_k1_a1.dat')
        cambios['field_output'] = cambios['spatial_output']
        CORRIDAS[f'paso/{x}_c{c:g}'] = (BASE_PASO, cambios)


def ruta(*p):
    return os.path.join(BASE, *p)


def correr():
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    for nombre, (base, cambios) in CORRIDAS.items():
        if os.path.exists(ruta(nombre + '.ok')):
            continue
        os.makedirs(os.path.dirname(ruta(nombre)), exist_ok=True)
        args = [f'{k}={v}' for k, v in dict(cambios, directory=f'validacion/{nombre}').items()]
        t0 = time.time()
        r = subprocess.run(['./VP_PIC', base] + args, cwd=EXE, stdout=open(ruta(nombre + '.log'), 'w'),
                           stderr=subprocess.STDOUT, env=env)
        print(f'{"OK" if r.returncode == 0 else "FALLO":5} {nombre} ({time.time()-t0:.0f} s)', flush=True)
        if r.returncode == 0:
            open(ruta(nombre + '.ok'), 'w').write('')


# ------------------------------------------------------------------ analizar
def coef():
    """a_n/a_0, n = 0..4, de la parte angular de la función de prueba. El código escribe
    a_n h_n, el término n de la suma del observable; aquí se da h_n, como en el artículo."""
    from exact import ak_test
    return np.array([ak_test(SQ, k)/ak_test(SQ, 0) for k in range(5)])


def hk(nombre):
    """t y h_n, n = 0..4, con el peso B(J) de la primera función de prueba (hk1_complex.tl)."""
    a = np.loadtxt(ruta(nombre, 'hk1_complex.tl'))
    return a[:, 0], (a[:, 1::2] + 1j*a[:, 2::2])/coef()


def pasos(nombre, base=BASE):
    import h5py
    f = h5py.File(os.path.join(base, nombre, 'vlasov_output.h5'), 'r')
    st = sorted([s for s in f if s.startswith('step_')], key=lambda s: int(s.split('_')[1]))
    return f, st, np.array([f[s].attrs['time'] for s in st])


def analizar_libre(w):
    from exact import make_hk
    hex_ = make_hk('gauss', 1.0e-4, 0.0, SJ, SQ, nJ=200001)
    t, h = hk('libre/nj400')
    ex = np.array([hex_(k, t) for k in range(5)]).T/coef()
    h0 = abs(ex[0, 0])
    w('Mezcla libre, pulso gaussiano. h_k con el peso B(J) = J^2 exp(-J^2/sJ^2), sin el coeficiente a_k de '
      'la función de prueba.\nError = max_t |h_k - h_k exacto| / h_0 en 0 <= t <= 2000.')
    w(f'h_0 exacto = {h0:.6e};  |h_k(0)|/h_0, k = 1..4: ' + ' '.join(f'{abs(ex[0, k])/h0:.4f}' for k in range(1, 5)))
    w(f'{"corrida":>14} {"N_J":>5} {"N_Q":>4} {"dt":>6} ' + ' '.join(f'{"k=" + str(k):>9}' for k in range(5)))
    for nombre, (_, c) in CORRIDAS.items():
        if not nombre.startswith('libre/') or nombre == 'libre/largo' or not os.path.exists(ruta(nombre + '.ok')):
            continue
        t1, h1 = hk(nombre)
        assert np.allclose(t1, t)
        err = np.max(np.abs(h1 - ex), axis=0)/h0
        dt = 'exacto' if c.get('integrator') == 'analytic' else f'{c["courant"]*0.05:g}'
        w(f'{nombre[6:]:>14} {c["Nrc"]:5d} {c["Npc"]:4d} {dt:>6} ' + ' '.join(f'{e:9.1e}' for e in err))
    # Error de la dinámica: frente al avance exacto de los mismos nodos.
    if os.path.exists(ruta('libre/exacto.ok')):
        _, ha = hk('libre/exacto')
        w('\nError del integrador: max_t |h_k - h_k con avance exacto| / h_0, mismos nodos (400 x 64).')
        w(f'{"dt":>6} ' + ' '.join(f'{"k=" + str(k):>9}' for k in range(5)) + '   cociente k = 1')
        previo = None
        for n, dt in (('libre/c4', 0.2), ('libre/nj400', 0.1), ('libre/c1', 0.05), ('libre/c0.5', 0.025)):
            if not os.path.exists(ruta(n + '.ok')):
                continue
            err = np.max(np.abs(hk(n)[1] - ha), axis=0)/h0
            w(f'{dt:6g} ' + ' '.join(f'{e:9.1e}' for e in err)
              + ('' if previo is None else f'   {previo/err[1]:6.1f}'))
            previo = err[1]
    if os.path.exists(ruta('libre/largo.ok')):
        tl, hl = hk('libre/largo')
        # 20001 nodos en J bastan hasta t = 10^4 (0.09 rad por nodo con k = 4) y ahorran tiempo.
        hex_ = make_hk('gauss', 1.0e-4, 0.0, SJ, SQ, nJ=20001)
        exl = np.array([hex_(k, tl) for k in range(5)]).T/coef()
        w('\nCorrida larga (800 x 64, dt = 0.1): |h_k|/h_0 del código y exacto, y el error máximo.')
        for a, b in ((0, 500), (500, 2000), (2000, 5000), (5000, 10000)):
            v = (tl >= a) & (tl <= b)
            w(f'   t en [{a:5d}, {b:5d}]: ' + '  '.join(
                f'k={k}: {np.max(np.abs(hl[v, k]))/h0:.1e} / {np.max(np.abs(exl[v, k]))/h0:.1e} / '
                f'{np.max(np.abs(hl[v, k] - exl[v, k]))/h0:.1e}' for k in (1, 2, 4)))
        f, st, ts = pasos('libre/largo')
        E = np.array([f[s].attrs['total_energy'] for s in st])
        w(f'   h_0: max|h_0(t)/h_0(0) - 1| = {np.max(np.abs(hl[:, 0]/hl[0, 0] - 1)):.1e};  '
          f'energía: max|E/E(0) - 1| = {np.max(np.abs(E/E[0] - 1)):.1e}')


def datos_paso():
    """Por cada courant con la corrida D hecha: tiempos, potencial de D en la zona de la cáscara,
    energía y archivos, hasta TFIN_PASO; y la respuesta (D - Z)/eps si también está Z. Las de
    courant = 2 son las de exe/hadzic."""
    HAD = os.path.join(EXE, 'hadzic')
    casos = [(4.0, BASE, 'paso/{}_c4'), (2.0, HAD, '{}_k1_a1'), (1.0, BASE, 'paso/{}_c1'), (0.5, BASE, 'paso/{}_c0.5')]
    hecha = lambda base, nombre: os.path.exists(os.path.join(base, nombre + '.ok'))
    datos = {}
    for c, base, patron in casos:
        if not hecha(base, patron.format('D')):
            continue
        fd, sd, td = pasos(patron.format('D'), base)
        v = td <= TFIN_PASO + 1e-6
        sd, td = [s for s, ok in zip(sd, v) if ok], td[v]
        rg = fd['grid']['r'][:]
        zona = (rg >= 3.0) & (rg <= 7.5)
        datos[c] = dict(t=td, f=fd, s=sd, pot=np.array([fd[a]['potential'][:] for a in sd])[:, zona],
                        E=np.array([fd[a].attrs['total_energy'] for a in sd]),
                        pmax=max(np.abs(fd[a]['p_part'][:]).max() for a in sd[::20]))
        if hecha(base, patron.format('Z')):
            fz, sz, tz = pasos(patron.format('Z'), base)
            sz = [s for s, ok in zip(sz, v) if ok]
            datos[c]['res'] = (datos[c]['pot'] - np.array([fz[b]['potential'][:] for b in sz])[:, zona])/EPS
    return datos


def serie_paso(nombre):
    """t y h_1 de una corrida del bloque paso en el mapa del equilibrio, como hadzic.serie (misma
    función de prueba, J1 = 0.35 y ancho 0.2); se guarda en la carpeta de la corrida."""
    import h5py
    from aa_numerico import MapaAA
    sal = ruta(nombre, 'serie.npz')
    if os.path.exists(sal):
        return np.load(sal)
    J1, SJ1 = 0.35, 0.20
    eq = np.load(os.path.join(EXE, 'hadzic', 'ic', 'D_k1_a1_equilibrio.npz'))
    mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']), fondo='puntual')
    f, st, ts = pasos(nombre)
    peso = f[st[0]]['f'][:]
    B = lambda J: J**2*np.exp(-(J - J1)**2/SJ1**2)
    h1 = []
    for a in st:
        Q, J, _ = mapa(f[a]['r_part'][:], f[a]['p_part'][:])
        h1.append(np.sum(peso*B(J)*np.exp(-1j*Q)))
    Q0, J0, _ = mapa(f[st[0]]['r_part'][:], f[st[0]]['p_part'][:])
    np.savez(sal, t=ts, h1=np.array(h1)/np.sum(peso*B(J0)))
    return np.load(sal)


def chi_paso(procesos=4):
    """chi_1 = (h_1 de D - h_1 de Z)/eps para cada courant disponible; la de courant = 2, de
    exe/hadzic/*/serie.npz."""
    from multiprocessing import Pool
    nuevos = [f'paso/{x}_c{c:g}' for c in (4.0, 1.0, 0.5) for x in 'DZ'
              if os.path.exists(ruta(f'paso/D_c{c:g}.ok')) and os.path.exists(ruta(f'paso/Z_c{c:g}.ok'))]
    with Pool(procesos) as pool:
        pool.map(serie_paso, nuevos)
    chi = {}
    for c in (4.0, 1.0, 0.5):
        if os.path.exists(ruta(f'paso/D_c{c:g}.ok')) and os.path.exists(ruta(f'paso/Z_c{c:g}.ok')):
            d, z = serie_paso(f'paso/D_c{c:g}'), serie_paso(f'paso/Z_c{c:g}')
            chi[c] = (d['t'], (d['h1'] - z['h1'])/EPS)
    d, z = (np.load(os.path.join(EXE, 'hadzic', f'{x}_k1_a1', 'serie.npz')) for x in 'DZ')
    m = d['t'] <= TFIN_PASO + 1e-6
    chi[2.0] = (d['t'][m], ((d['h1'] - z['h1'])/EPS)[m])
    return chi


VENTANAS = ((0, 250), (250, 500), (500, 750), (750, 1000), (1000, 1250), (1250, 1500), (1500, 2000))
NORMA = lambda x: np.sqrt(np.mean(x**2, axis=1))


def analizar_paso(w):
    """Error del paso de tiempo con autogravedad: la corrida D frente a la de dt más pequeño, en
    unidades de la respuesta, y la respuesta (D - Z)/eps entre las parejas completas."""
    datos = datos_paso()
    if len(datos) < 2:
        return
    from hadzic import frecuencia
    ref = min(datos)
    d0 = datos[ref]
    con_res = sorted(c for c in datos if 'res' in datos[c])
    escala = datos[con_res[0]]['res']              # respuesta con el dt menor que tiene pareja
    en = lambda x, a, b: np.max(NORMA(x[(d0['t'] >= a) & (d0['t'] <= b)]))
    w(f'\nPaso de tiempo con autogravedad: masa puntual, a0 = 1, k = 1, eps = 0.1, 400 x 25 partículas, '
      f't <= {TFIN_PASO:g}.\nNormas: rms en 3 <= r <= 7.5. Referencia: courant = {ref:g} (dt = {ref*0.05:g}).')
    w('\nCorrida D: ||Phi_D - Phi_D(ref)|| / (eps max ||dPhi_res||), máximo por ventana: el error del paso '
      'en unidades de la respuesta.')
    w(f'{"courant":>8} {"dt":>6} {"|p| dt/dr":>10} {"max|dE/E|":>10} ' + ' '.join(f'{f"[{a},{b}]":>12}' for a, b in VENTANAS))
    for c in sorted(datos, reverse=True):
        d = datos[c]
        assert np.allclose(d['t'], d0['t'])
        fila = '' if c == ref else ' '.join(f'{en(d["pot"] - d0["pot"], a, b)/(EPS*en(escala, a, b)):12.1e}'
                                             for a, b in VENTANAS)
        w(f'{c:8g} {c*0.05:6g} {d["pmax"]*c*0.05/0.1:10.3f} {np.max(np.abs(d["E"]/d["E"][0] - 1)):10.1e} {fila}')
    w('\nDiferencia máxima de radio entre partículas homólogas de D, frente a la referencia:')
    tiempos = (250, 500, 750, 1000, 1500, 2000)
    w(f'{"courant":>8} ' + ' '.join(f'{"t=" + str(x):>10}' for x in tiempos))
    for c in sorted(datos, reverse=True):
        if c == ref:
            continue
        d = datos[c]
        fila = []
        for tm in tiempos:
            i = np.argmin(abs(d['t'] - tm))
            fila.append(np.max(np.abs(d['f'][d['s'][i]]['r_part'][:] - d0['f'][d0['s'][i]]['r_part'][:])))
        w(f'{c:8g} ' + ' '.join(f'{x:10.1e}' for x in fila))
    if len(con_res) >= 2:
        r0 = datos[con_res[0]]['res']
        w(f'\nRespuesta: ||dPhi_res - dPhi_res(courant = {con_res[0]:g})|| / max ||dPhi_res|| por ventana, y '
          f'frecuencia del máximo del periodograma en t >= 500.')
        w(f'{"courant":>8} ' + ' '.join(f'{f"[{a},{b}]":>12}' for a, b in VENTANAS) + f' {"frecuencia":>11}')
        for c in sorted(con_res, reverse=True):
            d = datos[c]
            k = np.argmax(np.var(d['res'], axis=0))
            x = d['res'][:, k] - d['res'][d['t'] >= 500.0, k].mean()      # sin la parte estática
            fila = ' '.join(f'{"":>12}' if c == con_res[0] else f'{en(d["res"] - r0, a, b)/en(r0, a, b):12.1e}'
                            for a, b in VENTANAS)
            w(f'{c:8g} {fila} {frecuencia(d["t"], x, 500.0):11.5f}')
    # La referencia sin perturbar no depende del paso: su deriva es la misma con todos.
    HAD = os.path.join(EXE, 'hadzic')
    w('\nCorrida Z (eps = 0): ||Phi(t) - Phi(0)||, máximo por ventana.')
    w(f'{"courant":>8} ' + ' '.join(f'{f"[{a},{b}]":>12}' for a, b in VENTANAS))
    for c, base, nombre in ((4.0, BASE, 'paso/Z_c4'), (2.0, HAD, 'Z_k1_a1'), (1.0, BASE, 'paso/Z_c1'), (0.5, BASE, 'paso/Z_c0.5')):
        if not os.path.exists(os.path.join(base, nombre + '.ok')):
            continue
        f, st, ts = pasos(nombre, base)
        v = ts <= TFIN_PASO + 1e-6
        rg = f['grid']['r'][:]
        zona = (rg >= 3.0) & (rg <= 7.5)
        pot = np.array([f[a]['potential'][:] for a, ok in zip(st, v) if ok])[:, zona]
        w(f'{c:8g} ' + ' '.join(f'{en(pot - pot[0], a, b):12.2e}' for a, b in VENTANAS))


def armonico():
    """chi_1 = (h_1 de D - h_1 de Z)/eps en el mapa del equilibrio para cada courant con pareja
    completa, frente a la de paso más pequeño y a la teoría lineal."""
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    chi = chi_paso()
    from landau_cola import ajustar
    lin = np.load(os.path.join(EXE, 'hadzic', 'lineal', 'lin_k1_a1.npz'))
    ref = min(chi)
    t0, x0 = chi[ref]
    w('chi_1 = (h_1 de D - h_1 de Z)/eps en el mapa del equilibrio. Por ventanas: max|chi_1| de la teoría '
      f'lineal y de courant = {ref:g},\ny max|chi_1 - chi_1(courant = {ref:g})| / max|chi_1| de los demás. '
      'Polo: matrix pencil K = 3 en [300, 1000].')
    w(f'{"":>8} ' + ' '.join(f'{f"[{a},{b}]":>12}' for a, b in VENTANAS) + f' {"omega":>9} {"gamma":>9}')
    wl, gl = ajustar(lin['t'], lin['h1'], 300, 1000)['pencil M=3']
    w(f'{"lineal":>8} ' + ' '.join(f'{np.max(np.abs(lin["h1"][(lin["t"] >= a) & (lin["t"] <= b)])):12.3e}'
                                    for a, b in VENTANAS) + f' {wl:9.5f} {gl:+9.1e}')
    for c in sorted(chi):
        t, x = chi[c]
        assert np.allclose(t, t0)
        wp, gp = ajustar(t, x, 300, 1000)['pencil M=3']
        v = lambda y, a, b: np.max(np.abs(y[(t >= a) & (t <= b)]))
        fila = [v(x, a, b) if c == ref else v(x - x0, a, b)/v(x0, a, b) for a, b in VENTANAS]
        w(f'{"c=" + format(c, "g"):>8} ' + ' '.join(f'{y:12.3e}' if c == ref else f'{y:12.1e}' for y in fila)
          + f' {wp:9.5f} {gp:+9.1e}')
    open(ruta('armonico.txt'), 'w').write('\n'.join(lineas) + '\n')


# ------------------------------------------------------------------ residuo estático
def residuo():
    """Parte estática de h_1 del pulso con autogravedad (reproducir/corridas/09_autogravedad,
    salidas en exe/sg) según el mapa ángulo-acción: el del isócrono, que usa el código, y el del
    potencial total promediado en t >= 2000 (aa_meseta.py). La parte estática es el promedio
    de h_1/h_0 en 2000 <= t <= 20000; el resto, su rms en la misma ventana."""
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    SG = os.path.join(EXE, 'sg')
    rms = lambda x: np.sqrt(np.mean(np.abs(x)**2))
    v = lambda t: (t >= 2000.0) & (t <= 20000.0)
    w('Pulso gaussiano con autogravedad, isócrono, L0 = 2, dt = 0.1, t <= 20000. h_1/h_0 en 2000 <= t <= 20000.')
    w('Con el mapa del isócrono (h_1 del código):')
    w(f'{"corrida":>23} {"N_J x N_Q":>10} {"M_gas/M":>8} {"|<h_1>|/h_0":>12} {"rms del resto":>14}')
    for n in ('long20k_nosg', 'long_a0_1e-4', 'long_a0_1e-3', 'long_a0_1e-2', 'long20k_nrc800', 'long20k_npc50',
              'long20k_a0_1e-2_nrc800'):
        a = np.loadtxt(os.path.join(SG, n, 'hk1_complex.tl'))
        t, z = a[:, 0], (a[:, 3] + 1j*a[:, 4])/(a[0, 1] + 1j*a[0, 2])
        par = dict((x.split('=')[0].strip(), x.split('=')[1].split('#')[0].strip())
                   for x in open(os.path.join(SG, n, 'params_usados.par')) if '=' in x and not x.startswith('#'))
        masa = '0' if par.get('autointeraction', '.true.') == '.false.' else par['a0']
        z = z[v(t)]
        w(f'{n:>23} {par["Nrc"] + " x " + par["Npc"]:>10} {masa:>8} {abs(z.mean()):12.3e} {rms(z - z.mean()):14.3e}')
    d = np.load(os.path.join(SG, 'long20k_snap', 'aa_meseta.npz'))
    w('\nM_gas/M = 1e-3, 400 x 25 (long20k_snap, 501 instantáneas), según el mapa:')
    w(f'{"mapa":>23} {"|<h_1>|/h_0":>12} {"rms del resto":>14}')
    for k, nombre in (('iso', 'isócrono'), ('promedio', 'potencial promediado'), ('instante', 'potencial instantáneo')):
        z = d[f'h1_{k}'][v(d['t'])]
        w(f'{nombre:>23} {abs(z.mean()):12.3e} {rms(z - z.mean()):14.3e}')
    open(ruta('residuo.txt'), 'w').write('\n'.join(lineas) + '\n')


# ------------------------------------------------------------------ corridas de referencia
REFERENCIAS = (('Z_k1.25_a1', 'H_k1.25'), ('Zdr_k1.25_a1', 'H_k1.25_dr'), ('ZN_k1.25_a1', 'H_k1.25_N'),
               ('ZP_k1.25_a1', 'H_k1.25_P'))


def referencia():
    """Corridas de referencia (eps = 0) del politropo k = 1.25, a0 = 1 (exe/hadzic): |h_1|/h_0, que
    se anula en un estado estacionario (mayor valor en [t, t + 100]); el desplazamiento rms de la
    acción de las partículas (exe/relajacion, de relajacion_N.py); y la energía."""
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    T, TJ = (0, 500, 1000, 2000, 3900), (1000, 2000, 4000)
    w('Masa puntual, politropo k = 1.25, a0 = 1, eps = 0. |h_1|/h_0 (1e-3), mayor valor en [t, t + 100]; '
      'desplazamiento rms de la acción sobre J_max (1e-3); max|E/E(0) - 1|.')
    w(f'{"corrida":>13} {"N":>7} {"dr":>5} {"dt":>5} ' + ' '.join(f'{f"h1({a})":>9}' for a in T) + '  '
      + ' '.join(f'{f"dJ({a})":>9}' for a in TJ) + f' {"energía":>9}')
    for corrida, rel in REFERENCIAS:
        z = np.load(os.path.join(EXE, 'hadzic', corrida, 'serie.npz'))
        d = np.load(os.path.join(EXE, 'relajacion', rel + '.npz'))
        par = dict((x.split('=')[0].strip(), x.split('=')[1].split('#')[0].strip())
                   for x in open(os.path.join(EXE, 'hadzic', corrida, 'params_usados.par'))
                   if '=' in x and not x.startswith('#'))
        dr, dt = float(par['dr']), float(par['courant'])*float(par['dr'])/float(par['pmax'])
        h = np.abs(z['h1'])
        w(f'{corrida:>13} {int(d["N"]):7d} {dr:5g} {dt:5g} '
          + ' '.join(f'{1e3*h[(z["t"] >= a) & (z["t"] < a + 100)].max():9.2f}' for a in T) + '  '
          + ' '.join(f'{1e3*d["d2"][np.argmin(np.abs(d["t"] - a))]:9.2f}' for a in TJ)
          + f' {np.max(np.abs(z["E"]/z["E"][0] - 1)):9.1e}')
    eq = np.load(os.path.join(EXE, 'hadzic', 'ic', 'Z_k1.25_a1_equilibrio.npz'))
    J = np.linspace(0.0, float(eq['J_borde']), 2001)
    dom = np.abs(np.gradient(np.gradient(eq['E_t'], eq['J_t']), eq['J_t']))
    dmax = np.max(np.interp(J, eq['J_t'], dom))
    w(f'\nmax |dOmega/dJ| en 0 <= J <= J_max = {dmax:.3f}: recurrencia de la retícula de 400 filas en '
      f't = 2 pi/(|dOmega/dJ| dJ) = {2*np.pi/(dmax*float(eq["J_borde"])/400):.0f} (n = 1).')
    open(ruta('referencia.txt'), 'w').write('\n'.join(lineas) + '\n')


# ------------------------------------------------------------------ respuesta frente a la lineal
def respuesta():
    """Las corridas de masa pequeña (a0 = 0.01, 10^4 partículas) frente a la solución lineal, con
    el cociente complejo de amplitudes de hadzic.ganancia (ventana de Hann de ancho 800) y la
    diferencia punto a punto en los tiempos comunes. La comparación con a0 = 1 y tres amplitudes
    está en exe/hadzic/resumen.txt (hadzic.py analizar)."""
    import hadzic as H
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    T0, V = (500, 1000, 1500, 2000), ((0, 500), (500, 1000), (1000, 2000))
    w('Masa puntual, a0 = 0.01, 400 x 25 partículas, dt = 0.1. kappa(t0) = <chi_1, chi_lin>/<chi_lin, chi_lin> '
      'con ventana de Hann de ancho 800;\nmax|chi_1 - chi_lin|/max|chi_lin| por ventana; |chi_lin| al final '
      'de cada ventana sobre su valor inicial; y max|chi_1 - chi_lin|/|chi_lin(0)| por ventana (abs).')
    w(f'{"k":>5} {"eps":>4} ' + ' '.join(f'{f"|kappa|({a})":>12} {"arg":>7}' for a in T0) + '  '
      + ' '.join(f'{f"[{a},{b}]":>11}' for a, b in V) + '  ' + ' '.join(f'{f"lin({b})":>9}' for a, b in V)
      + '  ' + ' '.join(f'{f"abs[{a},{b}]":>15}' for a, b in V))
    for k in (0.75, 1.0, 1.5, 2.0):
        lin = np.load(H.ruta('lineal', H.nombre_lin(k, 0.01) + '.npz'))
        tl, xl = lin['t'], lin['h1']
        z = np.load(H.ruta(f'Z_k{k:g}', 'serie.npz'))
        t = z['t']
        _, i, j = np.intersect1d(np.round(t, 6), np.round(tl, 6), return_indices=True)
        tc, cl = t[i], xl[j]
        mx = lambda y, a, b: np.max(np.abs(y[(tc >= a) & (tc <= b)]))
        cola = ' '.join(f'{mx(cl, b - 100, b)/abs(cl[0]):9.1e}' for a, b in V)
        for nombre, eps in ((f'D_k{k:g}', 0.1), (f'D5_k{k:g}', 0.5)):
            x = (np.load(H.ruta(nombre, 'serie.npz'))['h1'] - z['h1'])[i]/eps
            g = [H.ganancia(tc, x, cl, a, 800.0) for a in T0]
            w(f'{k:5g} {eps:4g} ' + ' '.join(f'{abs(v):12.4f} {np.angle(v):+7.4f}' for v in g) + '  '
              + ' '.join(f'{mx(x - cl, a, b)/mx(cl, a, b):11.1e}' for a, b in V) + '  ' + cola
              + '  ' + ' '.join(f'{mx(x - cl, a, b)/abs(cl[0]):15.1e}' for a, b in V))
    open(ruta('respuesta.txt'), 'w').write('\n'.join(lineas) + '\n')


# ------------------------------------------------------------------ solucionador lineal
EQ_LIN = os.path.join(EXE, 'hadzic', 'lineal', 'P_k1.25_a1_equilibrio.npz')
BASE_LIN = (1600, 32, 0.5)                       # N_J, N_Q y dt de las soluciones del artículo
VARIANTES_LIN = [(1600, 32, 0.25), (1600, 64, 0.5), (3200, 32, 0.5)]
TFIN_LIN = 4000.0


def _lineal(arg):
    """h_1 de la solución lineal del politropo k = 1.25, a0 = 1 con otra resolución."""
    from lineal import resolver
    nj, nq, dt = arg
    sal = ruta('lineal', f'lin_nj{nj}_nq{nq}_dt{dt:g}.npz')
    if not os.path.exists(sal):
        t, h1, _, _, _ = resolver(EQ_LIN, nj, nq, dt, TFIN_LIN, verboso=False, j1=0.35, sj1=0.20)
        np.savez(sal, t=t, h1=h1)
    return arg


def lineal():
    """Convergencia del solucionador lineal: la solución con dt/2, con 2 N_Q y con 2 N_J frente
    a la del artículo (exe/hadzic/lineal/lin_k1.25_a1.npz)."""
    from multiprocessing import Pool
    from landau_cola import ajustar
    from hadzic import ganancia
    os.makedirs(ruta('lineal'), exist_ok=True)
    with Pool(len(VARIANTES_LIN)) as pool:
        pool.map(_lineal, VARIANTES_LIN)
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    base = np.load(os.path.join(EXE, 'hadzic', 'lineal', 'lin_k1.25_a1.npz'))
    t, x0 = base['t'], base['h1']
    V = ((0, 250), (250, 1000), (1000, 2000), (2000, 4000))
    mx = lambda y, a, b: np.max(np.abs(y[(t >= a) & (t <= b)]))
    polo = lambda x: ajustar(t, x, 300, 1500)['pencil M=3']
    w('Solucionador lineal, politropo k = 1.25, a0 = 1, masa puntual. Base: N_J x N_Q = '
      f'{BASE_LIN[0]} x {BASE_LIN[1]}, dt = {BASE_LIN[2]}.')
    w('max|chi_1 - chi_1 base| / max|chi_1 base| por ventana; polo (matrix pencil, K = 3) en [300, 1500]; '
      'media de |chi_1| en [3000, 4000];\ncociente de amplitudes con la base (mínimo y máximo de |kappa(t0)|, '
      '300 <= t0 <= 3900) y frecuencia menos la de la base (pendiente de -arg kappa).')
    w(f'{"N_J":>5} {"N_Q":>4} {"dt":>5} ' + ' '.join(f'{f"[{a},{b}]":>12}' for a, b in V)
      + f' {"omega":>9} {"gamma":>9} {"|chi_1|":>9} {"amplitud":>15} {"frecuencia":>11}')
    t0 = np.arange(300.0, 3901.0, 100.0)
    for nj, nq, dt in [BASE_LIN] + VARIANTES_LIN:
        if (nj, nq, dt) == BASE_LIN:
            x, fila, cola = x0, ' '.join(f'{"":>12}' for _ in V), ''
        else:
            d = np.load(ruta('lineal', f'lin_nj{nj}_nq{nq}_dt{dt:g}.npz'))
            assert np.allclose(d['t'], t)
            x = d['h1']
            fila = ' '.join(f'{mx(x - x0, a, b)/mx(x0, a, b):12.1e}' for a, b in V)
            g = np.array([ganancia(t, x, x0, a) for a in t0])
            cola = f' {abs(g).min():7.4f}-{abs(g).max():6.4f} {-np.polyfit(t0, np.unwrap(np.angle(g)), 1)[0]:+11.1e}'
        om, ga = polo(x)
        w(f'{nj:5d} {nq:4d} {dt:5g} {fila} {om:9.5f} {ga:+9.1e} {np.mean(np.abs(x[t >= 3000])):9.5f}{cola}')
    open(ruta('lineal.txt'), 'w').write('\n'.join(lineas) + '\n')


def figuras():
    """Figuras de la Sección 5 del artículo, en docs/articulo/figuras/."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from exact import make_hk
    destino = os.path.join(RAIZ, 'docs', 'articulo', 'figuras')
    os.makedirs(destino, exist_ok=True)
    plt.rcParams.update({'font.size': 9})
    # 1. Mezcla libre: |h_n|/h_0 de la corrida larga frente a la solución exacta, y el error del
    #    integrador frente al paso de tiempo.
    tl, hl = hk('libre/largo')
    hex_ = make_hk('gauss', 1.0e-4, 0.0, SJ, SQ, nJ=20001)
    te = np.unique(np.round(np.geomspace(10, 1.0e4, 60)/CADA)*CADA)
    h0 = abs(hex_(0, np.array([0.0]))[0])
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.1), constrained_layout=True)
    for k, c in zip((1, 2, 3, 4), ('C0', 'C1', 'C2', 'C3')):
        axs[0].loglog(tl[1:], np.abs(hl[1:, k])/h0, color=c, lw=1.0, label=f'$n={k}$')
        axs[0].loglog(te, np.abs(hex_(k, te))/coef()[k]/h0, 'o', color=c, ms=3, mfc='none', mew=0.7)
    axs[0].loglog([], [], 'ko', ms=3, mfc='none', mew=0.7, label='exact')
    axs[0].set_xlim(10, 1.0e4); axs[0].set_ylim(1e-11, 3)
    axs[0].set_xlabel('$t$'); axs[0].set_ylabel('$|h_n|/h_0$'); axs[0].grid(alpha=0.3)
    axs[0].legend(fontsize=7, frameon=False, loc='lower left')
    _, ha = hk('libre/exacto')
    dts, errs = [], []
    for nombre, dt in (('libre/c4', 0.2), ('libre/nj400', 0.1), ('libre/c1', 0.05), ('libre/c0.5', 0.025)):
        dts.append(dt); errs.append(np.max(np.abs(hk(nombre)[1] - ha), axis=0)/h0)
    dts, errs = np.array(dts), np.array(errs)
    for k, c in zip((1, 4), ('C0', 'C3')):
        axs[1].loglog(dts, errs[:, k], 'o-', color=c, ms=4, lw=1.0, label=f'$n={k}$')
    axs[1].loglog(dts, errs[1, 1]*(dts/dts[1])**4, 'k--', lw=0.8, label='$\\propto\\Delta t^4$')
    axs[1].set_xticks(dts); axs[1].set_xticks([], minor=True)
    axs[1].set_xticklabels([f'{x:g}' for x in dts])
    axs[1].set_xlabel('$\\Delta t$'); axs[1].set_ylabel('error of $h_n/h_0$'); axs[1].grid(alpha=0.3)
    axs[1].legend(fontsize=7, frameon=False)
    fig.savefig(os.path.join(destino, 'validacion_libre.pdf')); plt.close(fig)
    # 2. El pulso en el espacio fase, en (r, p_r) y en (Q, J), a cuatro tiempos de la corrida larga.
    #    Cada partícula es un punto con el color de su peso; solo las de peso > 1e-3 del máximo.
    #    Más tarde las franjas son más finas que la separación entre filas de la retícula y el
    #    dibujo, no la dinámica, deja de resolverlas.
    from df0 import rp_to_QJ
    f, st, ts = pasos('libre/largo')
    fig, axs = plt.subplots(2, 4, figsize=(7, 3.6), constrained_layout=True, sharey='row')
    for col, tt in enumerate((0.0, 100.0, 300.0, 1000.0)):
        g = f[st[np.argmin(abs(ts - tt))]]
        r, p, peso = g['r_part'][:], g['p_part'][:], g['f'][:]
        m = peso > 1e-3*peso.max()
        r, p, peso = r[m], p[m], peso[m]/peso.max()
        o = np.argsort(peso)
        Q, J = rp_to_QJ(r, p)
        for ax, x, y in ((axs[0, col], r, p), (axs[1, col], Q, J)):
            sc = ax.scatter(x[o], y[o], c=peso[o], s=0.6, lw=0, cmap='viridis', vmin=0, vmax=1, rasterized=True)
        axs[0, col].set_title(f'$t={tt:g}$', fontsize=9)
        axs[0, col].set_xlabel('$r$'); axs[1, col].set_xlabel('$Q$')
        axs[1, col].set_xlim(0, 2*np.pi); axs[1, col].set_xticks([0, np.pi, 2*np.pi])
        axs[1, col].set_xticklabels(['0', '$\\pi$', '$2\\pi$'])
    axs[0, 0].set_ylabel('$p_r$'); axs[1, 0].set_ylabel('$J$')
    fig.colorbar(sc, ax=axs, shrink=0.6, pad=0.01, label='$\\mathcal{F}_0/\\max\\mathcal{F}_0$')
    fig.savefig(os.path.join(destino, 'validacion_pulso.pdf'), dpi=300); plt.close(fig)
    # 3. Paso de tiempo con autogravedad: diferencia de la respuesta del potencial con la de la
    #    corrida de paso más pequeño.
    datos = {c: d for c, d in datos_paso().items() if 'res' in d}
    if len(datos) >= 2:
        ref = min(datos); d0 = datos[ref]
        fig, ax = plt.subplots(figsize=(4.2, 3.1), constrained_layout=True)
        for c, col in zip(sorted(datos, reverse=True), ('C3', 'C1', 'C0', 'C2')):
            if c != ref:
                ax.semilogy(d0['t'][1:], (NORMA(datos[c]['res'] - d0['res'])/np.max(NORMA(d0['res'])))[1:],
                            color=col, lw=1.0, label=f'$\\Delta t={c*0.05:g}$')
        ax.set_xlabel('$t$'); ax.set_ylabel(f'difference with $\\Delta t={ref*0.05:g}$'); ax.grid(alpha=0.3)
        ax.set_xlim(0, TFIN_PASO); ax.legend(fontsize=7, frameon=False, loc='lower right')
        fig.savefig(os.path.join(destino, 'validacion_paso.pdf')); plt.close(fig)
    # 4. Corridas de referencia (eps = 0) del politropo k = 1.25, a0 = 1 (exe/hadzic): |h_1|/h_0, que
    #    se anula en un estado estacionario, y el desplazamiento rms de la acción de las partículas
    #    (relajacion_N.py), con tres números de partículas y con dr/2.
    from hadzic import envolvente
    HAD, REL = os.path.join(EXE, 'hadzic'), os.path.join(EXE, 'relajacion')
    estilos = (('$N=10^4$', 'C0', '-'), ('$N=10^4$, $\\Delta r=0.05$', 'C0', '--'),
               ('$N=4\\times10^4$', 'C1', '-'), ('$N=1.024\\times10^5$', 'C2', '-'))
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.0), constrained_layout=True)
    for (corrida, rel), (etiqueta, c, ls) in zip(REFERENCIAS, estilos):
        z = np.load(os.path.join(HAD, corrida, 'serie.npz'))
        axs[0].semilogy(z['t'], envolvente(z['t'], z['h1']), color=c, ls=ls, lw=1.0, label=etiqueta)
        d = np.load(os.path.join(REL, rel + '.npz'))
        axs[1].semilogy(d['t'][1:], d['d2'][1:], color=c, ls=ls, lw=1.0)
    axs[0].set_ylabel('$|h_1^{(0)}|/h_0$'); axs[1].set_ylabel('$\\Delta J_{rms}/J_{max}$')
    axs[0].set_ylim(2e-4, 3e-2); axs[1].set_ylim(1e-4, 1e-1)
    for ax in axs:
        ax.set_xlabel('$t$'); ax.set_xlim(0, 4000); ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7, frameon=False, loc='upper left')
    fig.savefig(os.path.join(destino, 'validacion_referencia.pdf')); plt.close(fig)
    print('figuras en', destino)


def analizar():
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    analizar_libre(w)
    analizar_paso(w)
    open(ruta('resumen.txt'), 'w').write('\n'.join(lineas) + '\n')


if __name__ == '__main__':
    pasos_ = {'correr': correr, 'analizar': analizar, 'armonico': armonico, 'lineal': lineal, 'residuo': residuo,
              'respuesta': respuesta, 'referencia': referencia, 'figuras': figuras}
    if len(sys.argv) != 2 or sys.argv[1] not in pasos_:
        sys.exit(__doc__)
    pasos_[sys.argv[1]]()
