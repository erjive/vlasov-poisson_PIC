"""Resultados con L fijo para el artículo (docs/articulo, Sección 6).

Reúne, con las mismas definiciones, las dos familias de estados estacionarios:

  W   maxwellianas rebajadas (Wilson, g = 2, W0 = 3) en el isócrono, L0 = 2, con J_max = 0.138
      (y 0.083 y 0.276): equilibrios de la batería (exe/eta_lineal) y corridas de la demo de
      eta (demo_eta.py, exe/demo_eta);
  P   politropos (E_0 - E)^k con masa puntual, J_max = 0.7, L0 = 2: equilibrios y corridas del
      escenario de Hadžić (hadzic.py, exe/hadzic).

Pasos, en este orden (cada uno reutiliza lo que ya exista en exe/resultados/):
    python3 resultados.py estados    banda, lambda(0), lambda_edge, segundo autovalor y modo de
                                     cada estado estacionario (estados.txt; 12 min)
    python3 resultados.py umbrales   masa en la que lambda_edge = 1, por familia; calcula los
                                     equilibrios intermedios (umbrales.txt; 8 min)
    python3 resultados.py ganancia   lambda(omega) cerca del borde y leyes del borde, politropos
                                     con a0 = 1 (ganancia.txt; 1 min)
    python3 resultados.py respuesta  respuestas sin modo (polos PIC y lineal, mezcla libre,
                                     colas algebraicas con y sin autogravedad) y frecuencias de
                                     los modos con tres métodos (respuesta.txt; la primera vez
                                     resuelve 23 problemas lineales, unos 50 min)
    python3 resultados.py modos      el péndulo de las órbitas del borde en cada modo: omega_b,
                                     eps_c, J_r (modos.txt; 3 min)
    python3 resultados.py finita     corridas con modo: s, mayor acción frente a la del péndulo
                                     y cociente de amplitudes kappa (finita.txt)
    python3 resultados.py figuras    figuras de la Sección 6 (docs/articulo/figuras/)

No lanza simulaciones. Necesita los pasos "rebote" y "orbitas" de hadzic.py y las corridas de
demo_eta.py y hadzic.py ya analizadas. De dónde sale cada número de la Sección 6:
    6.A y 6.B  tablas de estados y de umbrales, lambda(0), segundo autovalor,
               leyes del borde y distancia del modo al borde             estados, umbrales, ganancia
    6.C        tiempos de caída con y sin autogravedad, polos, colas
               t^-(k+1) y t^-k, potencial acumulado en el borde          respuesta
    6.D        mesetas y frecuencias de los modos                        respuesta
    6.E        péndulo, eps_c, mayor acción y pérdida de amplitud,
               saturación del amortiguamiento                            modos, finita, respuesta
"""
import os, sys, subprocess
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')                # los pasos reparten casos en procesos
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
EXE = os.path.join(RAIZ, 'exe')
BASE = os.path.join(EXE, 'resultados')
ETA = os.path.join(EXE, 'eta_lineal')
import hadzic as H
import demo_eta as D

JT_W, JT_P = D.JT, H.JT
# Familia W: los equilibrios de la batería, por grupos.
W_GRUPOS = {'ref': ['L1', 'A1', 'L2', 'A2', 'L3', 'A3', 'L4', 'A4', 'L5', 'A5', 'L6', 'M5'],   # J_max = 0.138
            'ancho': ['A6', 'A7'], 'estrecho': ['A8', 'A9'],                      # J_max = 0.276 y 0.083
            'borde': ['G3', 'G15', 'G1', 'G075'],                                 # g = 3, 1.5, 1, 0.75; masa ~0.075
            'debil': ['G075a', 'G075b', 'G1a', 'G1b']}                            # g <= 1 con masa pequeña
# Modos de la familia W: caso -> (corrida con su equilibrio, ventanas del ajuste de c_+).
MODOS_W = {'L5': ('D3', ((3000, 5000), (5000, 7200))),
           'L6': ('D4', ((2000, 3500), (3500, 5520))),
           'M5': ('D9', ((1000, 3000), (3000, 5000)))}
# Corridas con modo discreto. W: corrida de la demo. P: (corrida, tiempo hasta el que su
# referencia conserva la retícula; después ni la mayor acción ni kappa miden el atrapamiento).
# CONTROLES: las que repiten otra con más partículas, otro paso u otra malla; no van a la figura.
FINITA_W = ['D3e003', 'D3e01', 'D3e03', 'D3', 'D3N', 'D3dt', 'D7', 'D4', 'D9', 'D10']
FINITA_P = [('DPe01_k1.25_a1', 4000), ('DP_k1.25_a1', 4000), ('DPe1_k1.25_a1', 4000), ('D_k1.25_a1', 1000),
            ('DN_k1.25_a1', 2000), ('De3_k1.25_a1', 1000), ('D_k1_a1', 1000), ('D_k0.75_a1', 1000)]
CONTROLES = {'D3N', 'D3dt', 'D_k1.25_a1', 'DN_k1.25_a1'}


def ruta(*p):
    return os.path.join(BASE, *p)


def escribir(nombre, lineas):
    os.makedirs(BASE, exist_ok=True)
    open(ruta(nombre), 'w').write('\n'.join(lineas) + '\n')


def eq_w(caso):
    return os.path.join(ETA, f'{caso}_equilibrio.npz')


def eq_p(k, a0):
    return H.ruta('lineal', f'P_k{k:g}_a{a0:g}_equilibrio.npz')


# ------------------------------------------------------------------ estados estacionarios
def _estado(arg):
    """Un estado estacionario: banda, extensión radial, lambda(0), lambda a 1e-3, 1e-6, 1e-9 y
    1e-12 anchos de banda del borde, lambda_edge (infinito si g <= 1), el segundo autovalor en
    el borde (a 1e-12 del borde si g <= 1) y el modo."""
    from lambda_borde import Lazo
    fam, etiqueta, npz = arg
    d = np.load(npz)
    z = Lazo(npz)
    g, an = float(d['k_borde']), z.Om_max - z.Om_min
    def dos(w):
        h = np.sqrt(np.abs(z.R(w)))
        return np.linalg.eigvalsh(h[:, None]*(-z.M)*h[None, :])[-2:]
    l3, l6, l9, l12 = (z.lam(z.Om_min - x*an) for x in (1e-3, 1e-6, 1e-9, 1e-12))
    e = dos(z.Om_min if g > 1 else z.Om_min - 1e-12*an)
    dd = z.delta_modo()
    r = d['r'][d['rho'] > 1e-10*d['rho'].max()]
    return dict(fam=fam, caso=etiqueta, a0=z.a0, g=g, jt=float(d['J_borde']), om_min=z.Om_min, om_max=z.Om_max,
                rin=r.min(), rout=r.max(), phi=float(np.abs(d['phi_self']).max()), lam0=z.lam(0.0), l3=l3, l6=l6,
                l9=l9, l12=l12, lame=float(e[1]) if g > 1 else np.inf, l2=float(e[0]),
                dd=np.nan if dd is None else dd, wd=np.nan if dd is None else z.Om_min - dd*an)


def estados_datos():
    """Los estados estacionarios de las dos familias (se calculan una vez; estados.npz)."""
    sal = ruta('estados.npz')
    if not os.path.exists(sal):
        from multiprocessing import Pool
        os.makedirs(BASE, exist_ok=True)
        casos = [('W' + g, c, eq_w(c)) for g, cs in W_GRUPOS.items() for c in cs]
        casos += [('P', f'k{k:g}_a{a0:g}', eq_p(k, a0)) for k in H.K_MAPA for a0 in H.A0_MAPA]
        with Pool(4) as pool:
            res = pool.map(_estado, casos, chunksize=1)
        np.savez(sal, **{k: np.array([r[k] for r in res]) for k in res[0]})
    d = np.load(sal)
    return [{k: (str(d[k][i]) if d[k].dtype.kind == 'U' else float(d[k][i])) for k in d.files}
            for i in range(len(d['caso']))]


def extrapolar_dd(e):
    """delta_d con la ley del borde cuando el modo queda a menos de 1e-12 del borde (g <= 1):
    lambda = a + b delta^(g-1) (g < 1) o a + b ln(1/delta) (g = 1), por los valores en 1e-9 y
    1e-12. Devuelve log10(delta_d)."""
    g = e['g']
    if g < 1:
        b = (e['l12'] - e['l9'])/(1e-12**(g - 1) - 1e-9**(g - 1))
        a = e['l12'] - b*1e-12**(g - 1)
        return np.log10((1 - a)/b)/(g - 1)
    b = (e['l12'] - e['l9'])/np.log(1e3)
    return -(np.log(1e12) + (1 - e['l12'])/b)/np.log(10)


def estados():
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Estados estacionarios. l(x): lambda a x anchos de banda bajo el borde; lam_edge = inf si g <= 1; l2: segundo '
      'autovalor en el borde (a 1e-12 del borde si g <= 1);')
    w('delta_d = (Omega_min - omega_d)/(Omega_max - Omega_min); entre paréntesis, log10(delta_d) extrapolado con la '
      'ley del borde.')
    w(f'{"fam":>9} {"caso":>10} {"masa":>7} {"g":>5} {"J_max":>6} {"r_in":>6} {"r_out":>6} {"Omega_min":>9} '
      f'{"Omega_max":>9} {"lam(0)":>7} {"l(1e-3)":>8} {"l(1e-6)":>8} {"l(1e-9)":>8} {"l(1e-12)":>8} {"lam_edge":>9} {"l2":>7} '
      f'{"omega_d":>9} {"delta_d":>9}')
    es = estados_datos()
    for e in es:
        if np.isnan(e['dd']):
            modo = f'{"--":>9} ' + (f'({extrapolar_dd(e):7.1f})' if e['g'] <= 1 else f'{"--":>9}')
        else:
            modo = f'{e["wd"]:9.6f} {e["dd"]:9.2e}'
        w(f'{e["fam"]:>9} {e["caso"]:>10} {e["a0"]:7g} {e["g"]:5g} {e["jt"]:6.3f} {e["rin"]:6.2f} {e["rout"]:6.2f} '
          f'{e["om_min"]:9.6f} {e["om_max"]:9.6f} {e["lam0"]:7.4f} {e["l3"]:8.4f} {e["l6"]:8.4f} {e["l9"]:8.4f} '
          f'{e["l12"]:8.4f} {e["lame"]:9.4f} {e["l2"]:7.4f} {modo}')
    w(f'\nmáximo de lambda(0): {max(e["lam0"] for e in es):.4f};  máximo del segundo autovalor: '
      f'{max(e["l2"] for e in es):.4f}')
    escribir('estados.txt', lineas)


# ------------------------------------------------------------------ masa umbral
# Familias: (etiqueta, familia, parámetro, intervalo de masas en el que se busca lambda_edge = 1).
UMBRALES = [('W, J_max = 0.138', 'W', 0.138, (0.075, 0.097)), ('W, J_max = 0.083', 'W', 0.083, (0.0109, 0.044)),
            ('W, J_max = 0.276', 'W', 0.276, (0.148, 1.0)), ('P, k = 1.25', 'P', 1.25, (0.3, 0.6)),
            ('P, k = 1.5', 'P', 1.5, (0.6, 1.0)), ('P, k = 2', 'P', 2.0, (1.0, 3.0)), ('P, k = 3', 'P', 3.0, (1.0, 5.0))]


def lam_borde(fam, par, a0):
    """lambda_edge y la banda del estado de masa a0 de una familia; genera el equilibrio si
    falta (exe/resultados/eq, con una retícula mínima: solo se usa el equilibrio)."""
    from lambda_borde import Lazo
    os.makedirs(ruta('eq'), exist_ok=True)
    base = ruta('eq', f'{fam}_{par:g}_a{a0:.6g}')
    if not os.path.exists(base + '_equilibrio.npz'):
        forma = (['--forma', 'maxwell', '--jt', str(par), '--k', '2', '--w0', str(D.W0)] if fam == 'W' else
                 ['--fondo', 'puntual', '--forma', 'polE', '--jt', str(JT_P), '--k', str(par), '--l0', str(H.L0)])
        subprocess.run([sys.executable, os.path.join(AQUI, 'equilibrio.py')] + forma
                       + ['--a0', f'{a0:.6g}', '--eps', '0', '--pert', 'suave3', '--nrc', '40', '--npc', '8',
                          '--maxit', '400', '--salida', base + '.dat'], check=True, stdout=open(base + '.log', 'w'),
                       stderr=subprocess.STDOUT)
    z = Lazo(base + '_equilibrio.npz')
    return z.lam(z.Om_min), z.Om_min, z.Om_max


def _umbral(arg):
    """Masa con lambda_edge = 1 por bisección-secante en ln(masa)."""
    from scipy.optimize import brentq
    etiqueta, fam, par, (lo, hi) = arg
    vistos = {}
    def f(x):
        a0 = float(np.exp(x))
        vistos[a0] = lam_borde(fam, par, a0)
        return vistos[a0][0] - 1
    try:
        if f(np.log(hi)) < 0:
            return etiqueta, None, vistos
        x = brentq(f, np.log(lo), np.log(hi), xtol=2e-4, rtol=1e-8)
    except subprocess.CalledProcessError:
        return etiqueta, None, vistos
    return etiqueta, float(np.exp(x)), vistos


def umbrales():
    from multiprocessing import Pool
    with Pool(4) as pool:
        res = pool.map(_umbral, UMBRALES, chunksize=1)
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Masa umbral, lambda_edge = 1, de cada familia, y los estados que se evaluaron al buscarla.')
    for etiqueta, a0, vistos in res:
        w(f'{etiqueta}: ' + (f'umbral = {a0:.4f}' if a0 else f'lambda_edge < 1 hasta la masa {max(vistos):g}'
                             if vistos else 'el equilibrio no converge'))
        for m in sorted(vistos):
            w(f'      masa {m:9.5f}   lambda_edge = {vistos[m][0]:8.5f}   Omega_min = {vistos[m][1]:.5f}   '
              f'Omega_max = {vistos[m][2]:.5f}')
    escribir('umbrales.txt', lineas)


# ------------------------------------------------------------------ ganancia cerca del borde
DELTAS = np.logspace(-9, np.log10(0.5), 53)


def _ganancia(k):
    from lambda_borde import Lazo
    z = Lazo(eq_p(k, 1.0))
    an = z.Om_max - z.Om_min
    return k, z.Om_min, z.Om_max, np.array([z.lam(z.Om_min - x*an) for x in DELTAS]), z.lam(z.Om_min) if k > 1 else np.inf


def ganancia_datos():
    sal = ruta('ganancia.npz')
    if not os.path.exists(sal):
        from multiprocessing import Pool
        os.makedirs(BASE, exist_ok=True)
        with Pool(4) as pool:
            res = pool.map(_ganancia, H.K_MAPA, chunksize=1)
        np.savez(sal, k=np.array([r[0] for r in res]), om_min=np.array([r[1] for r in res]),
                 om_max=np.array([r[2] for r in res]), lam=np.array([r[3] for r in res]),
                 lame=np.array([r[4] for r in res]), delta=DELTAS)
    return np.load(sal)


def ganancia():
    """lambda(omega) de los politropos con a0 = 1 a la distancia delta (en anchos de banda) bajo
    el borde, y el exponente local de la ley del borde."""
    d = ganancia_datos()
    lineas = ['lambda a la distancia delta (anchos de banda) bajo el borde; politropos, masa puntual, a0 = 1.',
              f'{"delta":>9} ' + ' '.join(f'{f"k={k:g}":>9}' for k in d['k'])]
    for i in range(0, len(d['delta']), 4):
        lineas.append(f'{d["delta"][i]:9.1e} ' + ' '.join(f'{d["lam"][j, i]:9.4f}' for j in range(len(d['k']))))
    lineas.append(f'{"borde":>9} ' + ' '.join(f'{x:9.4f}' for x in d['lame']))
    lineas.append('\nLey del borde entre delta = 1e-8 y 1e-5 (pendiente local, mínimo y máximo). k < 1: d ln(lambda)/d ln(delta);')
    lineas.append('k = 1: d lambda/d ln(1/delta); k > 1: d ln(lambda_edge - lambda)/d ln(delta), a comparar con min(k - 1, 1).')
    v = (d['delta'] >= 1e-8) & (d['delta'] <= 1e-5)
    x = np.log(d['delta'][v])
    for j, k in enumerate(d['k']):
        y = d['lam'][j, v]
        p = np.gradient(np.log(y), x) if k < 1 else -np.gradient(y, x) if k == 1 else np.gradient(np.log(d['lame'][j] - y), x)
        lineas.append(f'   k = {k:g}: {p.min():.4f} a {p.max():.4f}')
    escribir('ganancia.txt', lineas)
    print('\n'.join(lineas))


# ------------------------------------------------------------------ respuesta PIC y lineal
def serie_w(corrida, comun=True):
    """t, chi_1 = (h_1 de D - h_1 de Z)/eps de una corrida de la demo y la solución lineal de su
    caso; con comun, las dos en los tiempos comunes. Si no, (t, chi), (t_lin, chi_lin)."""
    caso, eps = D.info[corrida]['caso'], D.info[corrida]['eps']
    d, z = np.load(D.ruta(corrida, 'landau.npz')), np.load(D.ruta(D.REF[corrida], 'landau.npz'))
    n = min(len(d['t']), len(z['t']))
    t, x = d['t'][:n], (d['hk'][:n, 1] - z['hk'][:n, 1])/eps
    lin = np.load(D.ruta('lineal', caso + '.npz'))
    if not comun:
        return (t, x), (lin['t'], lin['h1'])
    _, i, j = np.intersect1d(np.round(t, 6), np.round(lin['t'], 6), return_indices=True)
    return t[i], x[i], lin['h1'][j]


def datos_p(nombre):
    return next((k, a0, e, r) for n, k, a0, e, r in H.CORRIDAS if n == nombre)


def serie_p(nombre, comun=True):
    """Lo mismo para una corrida de hadzic.py (serie.npz y lineal/lin_k*.npz)."""
    k, a0, eps, ref = datos_p(nombre)
    d, z = np.load(H.ruta(nombre, 'serie.npz')), np.load(H.ruta(ref, 'serie.npz'))
    t, x = d['t'], (d['h1'] - z['h1'])/eps
    lin = np.load(H.ruta('lineal', H.nombre_lin(k, a0) + '.npz'))
    if not comun:
        return (t, x), (lin['t'], lin['h1'])
    _, i, j = np.intersect1d(np.round(t, 6), np.round(lin['t'], 6), return_indices=True)
    return t[i], x[i], lin['h1'][j]


def particulas_p(nombre):
    c = H.CAMBIOS.get(nombre, {})
    return int(c.get('Nrc', H.NRC))*int(c.get('Npc', H.NPC))


def kappas(t, x, xl, t0, ancho=200.0):
    return np.array([H.ganancia(t, x, xl, a, ancho) for a in t0])


# Respuestas sin modo: familia, caso, corrida y ventanas del ajuste del polo.
AMORTIGUADAS = [('W', 'A1', 'D1', ((555, 1665),)), ('W', 'A3', 'D2', ((460, 1380), (1380, 2760))),
                ('W', 'A4', 'D5', ((1000, 2500), (2000, 4000), (4000, 7900))),
                ('P', 'k = 2', 'D_k2_a1', ((100, 400), (100, 600)))]
# Modos: familia, caso, ventana del ajuste de la frecuencia en la solución lineal, y corridas
# (de menor a mayor amplitud) con la suya, dentro del tiempo en que la referencia es fiable.
FRECUENCIAS = [('W', 'L5', (3600, 7200), [('D3e01', (3600, 7200)), ('D3e03', (3600, 7200)), ('D3', (3600, 7200))]),
               ('W', 'L6', (2760, 5520), [('D4', (2760, 5520))]),
               ('W', 'M5', (2500, 5000), [('D9', (2500, 5000)), ('D10', (2500, 5000))]),
               ('P', 0.75, (1000, 4000), [('D_k0.75_a1', (300, 1000))]),
               ('P', 1.0, (1000, 4000), [('D_k1_a1', (300, 1000))]),
               ('P', 1.25, (1000, 4000), [('DPe01_k1.25_a1', (500, 3900)), ('DP_k1.25_a1', (500, 3900)),
                                          ('DPe1_k1.25_a1', (500, 3900))]),
               ('P', 1.5, (1000, 4000), [])]


def cola_libre(k, a0=0.01):
    """Coeficiente C de la cola de la mezcla libre, |chi_1| ~ C t^-(k+1) (Ec. Tail del artículo
    con el estado estacionario autoconsistente): C = Gamma(k+1) alpha Omega_min^k B(J_max)
    / (2 |Omega'|^(k+1) int F_eq B dJ), para la perturbación F_eq s(J) cos Q con s(J_max) = 1."""
    from math import gamma
    from scipy.interpolate import CubicSpline
    from equilibrio import F_E
    d = np.load(eq_p(k, a0))
    sp = CubicSpline(d['J_t'], d['E_t'])
    x, w = np.polynomial.legendre.leggauss(400)
    s = 0.5*(x + 1)
    J, wJ = JT_P*(1 - s**2), 0.5*w*2*JT_P*s
    B = lambda j: j**2*np.exp(-(j - H.J1)**2/H.SJ1**2)
    F = F_E('polE', sp(J), float(sp(JT_P)), None, k)                 # sin la constante alpha
    om, dom = float(sp(JT_P, 1)), abs(float(sp(JT_P, 2)))
    return gamma(k + 1)*om**k*B(JT_P)/(2*dom**(k + 1)*np.sum(wJ*F*B(J))), om


# Colas: soluciones lineales de los politropos de masa 0.01 hasta t = 8000, con y sin la
# autogravedad de la perturbación (libre) y con 1600 y 3200 filas, para separar el efecto de
# la autogravedad en la cola del error de la retícula; la del politropo k = 2 con masa 1 con
# 3200 filas; y la solución lineal larga del modo pegado al borde (k = 1.5, masa 1). De las
# soluciones con autogravedad se guarda también el primer armónico del potencial de la
# perturbación en el borde, phi_1(t, J_max), extrapolado de las dos últimas filas.
K_COLAS = (0.75, 1.0, 1.5, 2.0)
# Las soluciones largas (masa 1, 1600 filas) llegan a t = 20000, pero solo se usan hasta 14000:
# con masa 1, max |Omega'| ~ 0.3, y el desfase entre filas vecinas del armónico n = 2 llega a pi
# en t ~ 12000-15000; después la suma sobre la retícula deja de representar la integral (en la
# de k = 2 la envolvente vuelve a crecer desde t ~ 14000).
T_LARGO = 20000.0


def _cola(arg):
    from lineal import resolver
    k, a0, nj, libre, tmax = arg
    res = resolver(eq_p(k, a0), nj, 32, 0.5, tmax, libre=libre, verboso=False, j1=H.J1, sj1=H.SJ1, phi1=not libre)
    if libre:
        return res[1], None
    J, om, ph = res[5]
    return res[1], (ph[:, -1] + (ph[:, -1] - ph[:, -2])*(JT_P - J[-1])/(J[-1] - J[-2])).real


def colas_datos():
    sal = ruta('colas.npz')
    if not os.path.exists(sal):
        from multiprocessing import Pool
        os.makedirs(BASE, exist_ok=True)
        casos = [(k, 0.01, nj, libre, 8000.0) for nj in (1600, 3200) for libre in (True, False) for k in K_COLAS]
        casos += [(2.0, 1.0, 3200, False, 8000.0), (1.5, 1.0, 1600, False, T_LARGO)]
        with Pool(4) as pool:
            res = pool.map(_cola, casos, chunksize=1)
        d = {}
        for (k, a0, nj, libre, _), (h, ph) in zip(casos[:-1], res[:-1]):
            n = f'k{k:g}_{nj}_{"libre" if libre else "sg"}' + ('' if a0 == 0.01 else f'_a{a0:g}')
            d[n] = h
            if ph is not None:
                d[n + '_phi'] = ph
        np.savez(sal, t=np.arange(len(res[0][0]))*2.0, largo=res[-1][0], **d)
    return np.load(sal)


# Otro estado sin modo y con autogravedad fuerte, con k no entero, para el exponente de la cola
# (3200 filas, hasta t = 8000). Las soluciones lineales tienen un piso numérico de ~3e-7
# |chi_1(0)|: la cola de k = 3 con masa 1 queda por debajo desde t ~ 3000 y no sirve, y la de
# k = 2 con masa 0.01 lo alcanza hacia t = 6000.
COLAS_EXTRA = [(1.5, 0.6)]


def colas_extra():
    sal = ruta('colas_extra.npz')
    if not os.path.exists(sal):
        from multiprocessing import Pool
        casos = [(k, a0, 3200, False, 8000.0) for k, a0 in COLAS_EXTRA]
        with Pool(len(casos)) as pool:
            res = pool.map(_cola, casos, chunksize=1)
        d = {}
        for (k, a0, *_), (h, ph) in zip(casos, res):
            d[f'k{k:g}_3200_sg_a{a0:g}'] = h
            d[f'k{k:g}_3200_sg_a{a0:g}_phi'] = ph
        np.savez(sal, t=np.arange(len(res[0][0]))*2.0, **d)
    return np.load(sal)


def phi_borde(c, nombre, k, a0):
    """Phi(t) = int_0^t phi_1(t', J_max) e^{i Omega_min t'} dt' de una solución de colas_datos, y
    |Omega'(J_max)| y Omega_min del estado."""
    from scipy.interpolate import CubicSpline
    eq = np.load(eq_p(k, a0))
    sp = CubicSpline(eq['J_t'], eq['E_t'])
    om_min, dom = float(sp(JT_P, 1)), abs(float(sp(JT_P, 2)))
    f = c[nombre + '_phi']*np.exp(1j*om_min*c['t'])
    return np.concatenate([[0], np.cumsum(0.5*(f[1:] + f[:-1])*(c['t'][1] - c['t'][0]))]), dom, om_min


def _libre_w(arg):
    from lineal import resolver
    caso, corrida = arg
    return resolver(D.ruta('ic', f'{corrida}_equilibrio.npz'), 1600, 32, 0.5, 6000.0, libre=True, verboso=False,
                    j1=D.J1, sj1=D.J1)[:2]


def libres_w():
    """Mezcla libre (sin la autogravedad de la perturbación) de los estados de la familia W bajo
    el umbral, para compararla con la respuesta completa (libres_w.npz)."""
    sal = ruta('libres_w.npz')
    casos = [(c, n) for _, c, n, _ in AMORTIGUADAS if c in D.CASOS]
    if not os.path.exists(sal):
        from multiprocessing import Pool
        with Pool(len(casos)) as pool:
            res = pool.map(_libre_w, casos)
        np.savez(sal, t=res[0][0], **{c: h for (c, _), (_, h) in zip(casos, res)})
    return np.load(sal)


def largo_k2():
    """Solución lineal del politropo k = 2 con masa 1 hasta T_LARGO, para la cola (largo_k2.npz)."""
    sal = ruta('largo_k2.npz')
    if not os.path.exists(sal):
        from lineal import resolver
        t, h1 = resolver(eq_p(2.0, 1.0), 1600, 32, 0.5, T_LARGO, verboso=False, j1=H.J1, sj1=H.SJ1)[:2]
        np.savez(sal, t=t, h1=h1)
    return np.load(sal)


def respuesta():
    """Respuestas sin modo (polos y colas) y frecuencias de los modos."""
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    est = {e['caso']: e for e in estados_datos()}
    w('Respuesta amortiguada: polo (omega, gamma) de chi_1 con matrix pencil (K = 3; incertidumbre con K = 2, 3 y la '
      'ventana recortada un 10 %), PIC y lineal en la misma ventana, y |kappa| (ventana de Hann de 200) en ella.')
    for fam, caso, corrida, ventanas in AMORTIGUADAS:
        t, x, xl = serie_w(corrida) if fam == 'W' else serie_p(corrida)
        eps = D.info[corrida]['eps'] if fam == 'W' else datos_p(corrida)[2]
        e = est[caso] if fam == 'W' else est['k2_a1']
        w(f'  {fam} {caso:>6} {corrida:>8} eps = {eps:<5g} lam_edge = {e["lame"]:.4f}  Omega_min = {e["om_min"]:.5f}  '
          f'Omega_max = {e["om_max"]:.5f}  |chi_1(0)| = {abs(xl[0]):.4f}')
        for lo, hi in ventanas:
            p, pl = D.polo_con_error(t, x, lo, hi), D.polo_con_error(t, xl, lo, hi)
            kap = np.abs(kappas(t, x, xl, np.arange(lo + 100, hi - 99, 50.0)))
            w(f'     [{lo:5d}, {hi:5d}]  PIC    omega = {p[0]:.5f} +- {p[1]:.0e}   gamma = {p[2]:.3e} +- {p[3]:.0e}')
            w(f'                     lineal omega = {pl[0]:.5f} +- {pl[1]:.0e}   gamma = {pl[2]:.3e} +- {pl[3]:.0e}   '
              f'|kappa| entre {kap.min():.3f} y {kap.max():.3f};  (omega - Omega_min)/ancho = '
              f'{(pl[0] - e["om_min"])/(e["om_max"] - e["om_min"]):+.3f};  |chi_1| lineal al final / inicial = '
              f'{abs(xl[np.argmin(np.abs(t - hi))])/abs(xl[0]):.1e}')

    lw = libres_w()
    w('\nTiempo en que |chi_1| baja de 0.1, 0.01 y 0.001 de su valor inicial: mezcla libre (sin la autogravedad de la '
      'perturbación) y solución lineal completa.')
    for caso in lw.files[1:]:
        lin = np.load(D.ruta('lineal', caso + '.npz'))
        fila = []
        for t, h in ((lw['t'], lw[caso]), (lin['t'], lin['h1'])):
            x = np.abs(h)/abs(h[0])
            fila.append(', '.join(f'{t[np.argmax(x < f)]:.0f}' if np.any(x < f) else '--' for f in (0.1, 0.01, 0.001)))
        w(f'   {caso}: libre {fila[0]};  completa {fila[1]}')
    lin = np.load(D.ruta('lineal', 'A4.npz'))
    env = H.envolvente(lin['t'], lin['h1'], 200.0)
    w('   A4, solución lineal: envolvente/|chi_1(0)| y tasa local de decaimiento en t = 4000, 8000, 12000, 14000, 16000: '
      + ', '.join(f'{env[i]/abs(lin["h1"][0]):.1e} ({-np.log(env[j]/env[i])/(lin["t"][j] - lin["t"][i]):.1e})'
                  for i, j in ((np.argmin(np.abs(lin['t'] - a)), np.argmin(np.abs(lin['t'] - 1.1*a)))
                               for a in (4000, 8000, 12000, 14000, 16000))))
    l4 = np.load(os.path.join(ETA, 'L4_lineal.npz'))
    for lo, hi in ((1300, 2600), (2000, 3900), (1300, 3900)):
        p = D.polo_con_error(l4['t'], l4['h1'], lo, hi)
        e = est['L4']
        w(f'   L4 (masa 0.058, lam_edge = {e["lame"]:.4f}), solución lineal de la batería, [{lo}, {hi}]: omega = {p[0]:.5f} '
          f'+- {p[1]:.0e}  gamma = {p[2]:.3e} +- {p[3]:.0e}  (omega - Omega_min)/ancho = '
          f'{(p[0] - e["om_min"])/(e["om_max"] - e["om_min"]):+.3f}')
    w('   A4 con eps = 0.075: |kappa| en t0 = 2000, 4000, 5000, 6000, 7000 con 400 x 25, con 800 x 50 y con dt/2:')
    for corrida in ('D5', 'D5N', 'D5dt'):
        t, x, xl = serie_w(corrida)
        w(f'      {corrida:>5}: ' + ', '.join(f'{abs(H.ganancia(t, x, xl, a)):.3f}' for a in (2000, 4000, 5000, 6000, 7000)))
    w('   A4 con eps = 0.84 (saturación): envolvente de |chi_1| / |chi_1(0)| (máximo en una ventana de 600) en sus '
      'extremos, y en t = 3000, 4000, 5000 con 400 x 25, 800 x 50 y dt/2 (ventana de 100):')
    from scipy.signal import argrelextrema
    (t, x), (tl, xl) = serie_w('D6L', comun=False)
    e2 = H.envolvente(t, x, 600.0)/abs(xl[0])
    for nombre, f in (('mínimos', np.less_equal), ('máximos', np.greater_equal)):
        idx = [i for i in argrelextrema(e2, f, order=60)[0] if 500 < t[i] < t[-1] - 300]
        idx = [i for n, i in enumerate(idx) if n == 0 or t[i] - t[idx[n - 1]] > 1000]
        w(f'      {nombre}: ' + ', '.join(f't = {t[i]:.0f}: {e2[i]:.3f}' for i in idx))
    for corrida in ('D6L', 'D6N', 'D6dt'):
        (t, x), _ = serie_w(corrida, comun=False)
        e1 = H.envolvente(t, x, 100.0)/abs(xl[0])
        w(f'      {corrida:>5}: ' + ', '.join(f'{e1[np.argmin(np.abs(t - a))]:.4f}' for a in (3000, 4000, 5000)))

    c = colas_datos()
    tc = c['t']
    env = lambda x, ancho=150.0: H.envolvente(tc, x, ancho)
    t10 = tc
    en = lambda y, T: ' '.join(f'{y[np.argmin(np.abs(t10 - a))]:6.3f}' for a in T)
    w('\nColas algebraicas, politropos con masa 0.01 (soluciones lineales hasta t = 8000). Exponente: pendiente de '
      'ln(envolvente de |chi_1|) frente a ln t en [2000, 4000], con 3200 filas (y con 1600); envolvente: máximo en una '
      'ventana de 150.')
    w('libre/cola: envolvente de la mezcla libre (sin la autogravedad de la perturbación) entre C t^-(k+1), con 3200 y '
      'con 1600 filas; sg/libre: cociente de las envolventes con y sin autogravedad.')
    T = (1000, 2000, 4000, 6000, 8000)
    w(f'{"k":>5} {"-(k+1)":>7} {"exp. sg":>16} {"exp. libre":>16} {"C":>10} | {"libre/cola 3200:":>17} '
      + ' '.join(f'{a:>6}' for a in T) + f' | {"1600:":>6} ' + ' '.join(f'{a:>6}' for a in T)
      + f' | {"sg/libre 3200:":>14} ' + ' '.join(f'{a:>6}' for a in T) + f' | {"1600:":>6} ' + ' '.join(f'{a:>6}' for a in T))
    envs = {}
    for k in K_COLAS:
        e = {n: env(c[f'k{k:g}_{n}']) for n in ('1600_libre', '3200_libre', '1600_sg', '3200_sg')}
        envs[k] = e
        v = (t10 >= 2000) & (t10 <= 4000)
        pend = {n: np.polyfit(np.log(t10[v]), np.log(e[n][v]), 1)[0] for n in e}
        C, om = cola_libre(k)
        cola = C*np.maximum(t10, 1.0)**-(k + 1)
        w(f'{k:5g} {-(k + 1):7.2f} {pend["3200_sg"]:7.3f} ({pend["1600_sg"]:6.3f}) {pend["3200_libre"]:7.3f} '
          f'({pend["1600_libre"]:6.3f}) {C:10.3e} | {"":>17} {en(e["3200_libre"]/cola, T)} | {"":>6} '
          f'{en(e["1600_libre"]/cola, T)} | {"":>14} {en(e["3200_sg"]/e["3200_libre"], T)} | {"":>6} '
          f'{en(e["1600_sg"]/e["1600_libre"], T)}')
    w('\nCola con autogravedad a primer orden en la masa: chi_1 ~ chi_1 libre x (1 + 2 |Omega\'| Phi(t) t), con Phi(t) = '
      'int_0^t phi_1(t\', J_max) e^{i Omega_min t\'} dt\' (3200 filas). t_x = 1/(2 |Omega\'| |Phi|): cruce a t^-k.')
    w('pred: C t^-(k+1) |1 + 2 |Omega\'| Phi t|; env/pred: la envolvente de la solución lineal entre esa predicción.')
    T2 = (2000, 4000, 8000)
    w(f'{"k":>5} {"masa":>5} {"|Omega\'|":>8} | {"|Phi|:":>7} ' + ' '.join(f'{a:>9}' for a in T2) + f' | {"arg Phi:":>8} '
      + ' '.join(f'{a:>6}' for a in T2) + f' | {"2|Om\'||Phi|t:":>13} ' + ' '.join(f'{a:>7}' for a in T2)
      + f' | {"t_x(8000)":>9} | {"env/pred:":>9} ' + ' '.join(f'{a:>6}' for a in T2) + f' | {"env t^k:":>8} '
      + ' '.join(f'{a:>9}' for a in T2))
    for k, a0, nombre in [(k, 0.01, f'k{k:g}_3200_sg') for k in K_COLAS] + [(2.0, 1.0, 'k2_3200_sg_a1')]:
        Phi, dom, om_min = phi_borde(c, nombre, k, a0)
        C, _ = cola_libre(k, a0)
        e = H.envolvente(t10, c[nombre], 150.0)
        pred = C*np.maximum(t10, 1.0)**-(k + 1)*np.abs(1 + 2*dom*Phi*t10)
        i = [int(np.argmin(np.abs(t10 - a))) for a in T2]
        w(f'{k:5g} {a0:5g} {dom:8.4f} | {"":>7} ' + ' '.join(f'{abs(Phi[j]):9.3e}' for j in i) + f' | {"":>8} '
          + ' '.join(f'{np.angle(Phi[j]):+6.3f}' for j in i) + f' | {"":>13} '
          + ' '.join(f'{2*dom*abs(Phi[j])*t10[j]:7.3f}' for j in i) + f' | {1/(2*dom*abs(Phi[-1])):9.2e} | {"":>9} '
          + ' '.join(f'{e[j]/pred[j]:6.3f}' for j in i) + f' | {"":>8} ' + ' '.join(f'{e[j]*t10[j]**k:9.3e}' for j in i))
    e = H.envolvente(t10, c['k2_3200_sg_a1'], 150.0)
    w('Politropo k = 2 con masa 1, 3200 filas: exponente de la envolvente por intervalos: '
      + ', '.join(f'[{a}, {b}]: {np.polyfit(np.log(t10[v]), np.log(e[v]), 1)[0]:.2f}'
                  for a, b, v in ((a, b, (t10 >= a) & (t10 <= b)) for a, b in ((1500, 2500), (2000, 4000), (4000, 8000))))
      + ';  envolvente/|chi_1(0)| en t = 2000, 8000: '
      + ', '.join(f'{e[np.argmin(np.abs(t10 - a))]/abs(c["k2_3200_sg_a1"][0]):.2e}' for a in (2000, 8000)))
    cx = colas_extra()
    for k, a0 in COLAS_EXTRA:
        nombre = f'k{k:g}_3200_sg_a{a0:g}'
        e = H.envolvente(cx['t'], cx[nombre], 150.0)
        Phi, dom, om_min = phi_borde(cx, nombre, k, a0)
        C, _ = cola_libre(k, a0)
        w(f'Politropo k = {k:g} con masa {a0:g} (lam_edge = {est[f"k{k:g}_a{a0:g}"]["lame"]:.4f}), 3200 filas: exponente de la '
          'envolvente por intervalos: '
          + ', '.join(f'[{a}, {b}]: {np.polyfit(np.log(cx["t"][v]), np.log(e[v]), 1)[0]:.2f}'
                      for a, b, v in ((a, b, (cx['t'] >= a) & (cx['t'] <= b)) for a, b in
                                      ((1500, 2500), (2000, 4000), (4000, 8000))))
          + f';  env/(C t^-(k+1)) en t = 2000, 8000: '
          + ', '.join(f'{e[np.argmin(np.abs(cx["t"] - a))]/(C*a**-(k + 1)):.3g}' for a in (2000, 8000))
          + f';  |Phi| = {abs(Phi[-1]):.3e}, t_x = {1/(2*dom*abs(Phi[-1])):.3g};  env/pred en 4000: '
          + f'{e[np.argmin(np.abs(cx["t"] - 4000))]/(C*4000.0**-(k + 1)*abs(1 + 2*dom*Phi[np.argmin(np.abs(cx["t"] - 4000))]*4000)):.3g}')
    w('PIC (eps = 0.1, 400 x 25) frente a la lineal de 1600 filas: |kappa| en 300 <= t0 <= 1000 y la mayor '
      'diferencia en t <= 2000, en unidades de |chi_1(0)|; y |chi_1| lineal en t = 1000 y 2000 entre |chi_1(0)|.')
    for k in K_COLAS:
        t, x, xl = serie_p(f'D_k{k:g}')
        kap = np.abs(kappas(t, x, xl, np.arange(400.0, 901.0, 50.0)))
        dif = np.abs(x - xl)/abs(xl[0])
        w(f'   k = {k:g}: |kappa| entre {kap.min():.3f} y {kap.max():.3f};  diferencia máxima {dif[t <= 2000].max():.1e};  '
          + ', '.join(f'{H.envolvente(t, xl, 150.0)[np.argmin(np.abs(t - a))]/abs(xl[0]):.1e}' for a in (1000, 2000)))

    lk = largo_k2()
    env = H.envolvente(lk['t'], lk['h1'], 150.0)
    w(f'\nPolitropo k = 2 con masa 1 (lam_edge = {est["k2_a1"]["lame"]:.4f}): exponente de la envolvente de la solución '
      'lineal por intervalos: '
      + ', '.join(f'[{a}, {b}]: {np.polyfit(np.log(lk["t"][v]), np.log(env[v]), 1)[0]:.2f}'
                  for a, b, v in ((a, b, (lk['t'] >= a) & (lk['t'] <= b)) for a, b in
                                  ((1500, 2500), (2000, 4000), (4000, 8000), (8000, 14000), (14000, 20000)))))

    w('\nMesetas de los modos: envolvente de |chi_1| lineal (ventana de 100) al final, entre |chi_1(0)|.')
    for fam, caso, vl, corridas in FRECUENCIAS:
        (_, _), (tl, xl) = serie_w(MODOS_W[caso][0], comun=False) if fam == 'W' else serie_p(f'D_k{caso:g}_a1', comun=False)
        env = H.envolvente(tl, xl, 100.0)
        w(f'   {fam} {caso}: |chi_1(0)| = {abs(xl[0]):.4f};  envolvente final / inicial = {env[-60]/abs(xl[0]):.3f};  '
          f'periodos hasta el final: {tl[-1]*est[caso if fam == "W" else f"k{caso:g}_a1"]["wd"]/(2*np.pi):.0f}')

    w('\nFrecuencia de los modos: omega_d de lambda = 1; matrix pencil (K = 3) de la solución lineal y de cada corrida, '
      'cada una en su ventana; el máximo del periodograma de la corrida; y omega - omega_lineal de la fase de kappa.')
    for fam, caso, vl, corridas in FRECUENCIAS:
        e = est[caso] if fam == 'W' else est[f'k{caso:g}_a1']
        w(f'  {fam} {caso if fam == "W" else f"k = {caso:g}":>8}: Omega_min = {e["om_min"]:.6f}  omega_d = {e["wd"]:.6f}  '
          f'delta_d = {e["dd"]:.4f}')
        (_, _), (tl, xl) = serie_w(MODOS_W[caso][0], comun=False) if fam == 'W' else serie_p(f'D_k{caso:g}_a1', comun=False)
        pl = D.polo_con_error(tl, xl, *vl)
        w(f'       lineal [{vl[0]}, {vl[1]}]:{"":>27} omega = {pl[0]:.6f} +- {pl[1]:.0e}   gamma = {pl[2]:+.1e}   '
          f'|chi_1| = {np.abs(xl[(tl >= vl[0]) & (tl <= vl[1])]).mean():.4f}')
        if fam == 'P' and caso == 1.5:
            tg = np.arange(len(c['largo']))*2.0
            eg = H.envolvente(tg, c['largo'], 100.0)
            w('       lineal larga, envolvente (ventana de 100) en t = 500, 1000, 2000, 4000, 6000, 8000, 10000, 14000: '
              + ', '.join(f'{eg[np.argmin(np.abs(tg - a))]:.4f}' for a in (500, 1000, 2000, 4000, 6000, 8000, 10000, 14000))
              + f';  entre {eg[(tg >= 4000) & (tg <= 14000)].min():.4f} y {eg[(tg >= 4000) & (tg <= 14000)].max():.4f} en '
              f'[4000, 14000];  |chi_1(0)| = {abs(c["largo"][0]):.4f}')
            for a, b in ((4000, 8000), (8000, 14000), (14000, 20000)):
                pl = D.polo_con_error(tg, c['largo'], a, b)
                w(f'       lineal hasta {T_LARGO:g} [{a}, {b}]:{"":>10} omega = {pl[0]:.6f} +- {pl[1]:.0e}   gamma = {pl[2]:+.1e}   '
                  f'|chi_1| = {np.abs(c["largo"][(tg >= a) & (tg <= b)]).mean():.4f}')
        for corrida, (lo, hi) in corridas:
            t, x, xl = serie_w(corrida) if fam == 'W' else serie_p(corrida)
            eps = D.info[corrida]['eps'] if fam == 'W' else datos_p(corrida)[2]
            part = D.info[corrida]['nrc']*D.info[corrida]['npc'] if fam == 'W' else particulas_p(corrida)
            p = D.polo_con_error(t, x, lo, hi)
            v = (t >= lo) & (t <= hi)
            om = np.linspace(p[0] - 2e-3, p[0] + 2e-3, 4001)
            per = om[np.argmax(np.abs(np.exp(-1j*np.outer(om, t[v])) @ (x[v]*np.hanning(v.sum()))))]
            t0 = np.arange(lo + 100, hi - 99, 50.0)
            pend = -np.polyfit(t0, np.unwrap(np.angle(kappas(t, x, xl, t0))), 1)[0]
            w(f'       {corrida:>15} eps = {eps:<5g} N = {part:6d} [{lo}, {hi}]: omega = {p[0]:.6f} +- {p[1]:.0e}   '
              f'gamma = {p[2]:+.1e}   periodograma {per:.6f}   omega - omega_lineal = {pend:+.1e}')
    escribir('respuesta.txt', lineas)


# ------------------------------------------------------------------ modos y péndulo
def _modo_w(caso):
    """c_+ en el borde de un modo de la familia W: la parte e^{-i omega_d t} del primer armónico
    del potencial de la perturbación en la fila del borde (lineal.resolver con phi1), ajustada
    en dos ventanas."""
    from lineal import resolver
    corrida, ventanas = MODOS_W[caso]
    e = next(x for x in estados_datos() if x['caso'] == caso)
    lin = np.load(D.ruta('lineal', caso + '.npz'))
    t, h1, _, _, _, (J, om, ph) = resolver(D.ruta('ic', f'{corrida}_equilibrio.npz'), 1600, 32, 0.5,
                                           float(lin['t'][-1]), verboso=False, j1=D.J1, sj1=D.J1, phi1=True)
    assert np.array_equal(h1, lin['h1']), 'la solución lineal no es la guardada'
    cs, res = [], []
    for a, b in ventanas:
        v = (t >= a) & (t <= b)
        M = np.stack([np.exp(-1j*e['wd']*t[v]), np.exp(1j*e['wd']*t[v])], axis=1)
        c = np.linalg.lstsq(M, ph[v, -1], rcond=None)[0]
        cs.append(c[0]); res.append(np.sqrt(np.sum(np.abs(ph[v, -1] - M @ c)**2)/np.sum(np.abs(ph[v, -1])**2)))
    return dict(caso=caso, masa=e['a0'], om_min=e['om_min'], om_max=e['om_max'], wd=e['wd'],
                dom=abs(np.gradient(om, J)[-1]), cmas=cs[0], cmas2=cs[1], res=res[0], res2=res[1])


def modos_w():
    """Los modos de la familia W (se calculan una vez; exe/resultados/modos_w.npz)."""
    sal = ruta('modos_w.npz')
    if not os.path.exists(sal):
        from multiprocessing import Pool
        estados_datos()
        with Pool(len(MODOS_W)) as pool:
            res = pool.map(_modo_w, list(MODOS_W))
        np.savez(sal, **{k: np.array([r[k] for r in res]) for k in res[0]})
    d = np.load(sal)
    return {str(c): {k: d[k][i] for k in d.files} for i, c in enumerate(d['caso'])}


def pendulo(familia, caso, eps):
    """Péndulo de las órbitas del borde en el modo de un estado estacionario, con amplitud eps:
    omega_b, s = 2 omega_b/(Omega_min - omega_d), J_r, semiancho de la separatriz, la mayor
    acción que alcanza una órbita del borde, la fase del punto O y omega_d. familia 'P':
    caso = k (hadzic.pendulo); familia 'W': caso = 'L5', 'L6' o 'M5'."""
    if familia == 'P':
        return H.pendulo(caso, eps)
    m = modos_w()[caso]
    hueco, dom = float(m['om_min'] - m['wd']), float(m['dom'])
    wb = float(np.sqrt(2*abs(m['cmas'])*dom*eps))
    d, dJs = hueco/dom, 2*wb/dom
    jmax = JT_W + d + dJs if dJs > d else JT_W + d - np.sqrt(d**2 - dJs**2)
    return dict(wb=wb, s=2*wb/hueco, Jr=JT_W + d, dJs=dJs, Jmax=jmax, fase_O=-float(np.angle(m['cmas'])),
                wd=float(m['wd']))


def modos():
    """Tabla del péndulo de los modos discretos de las dos familias."""
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Péndulo de las órbitas del borde. hueco = Omega_min - omega_d; c_+ en el borde (dos ventanas, con el residuo '
      'del ajuste); omega_b = sqrt(2 eps |c_+| |Omega\'|); eps_c = hueco^2/(8 |c_+| |Omega\'|); J_r = J_max + hueco/|Omega\'|.')
    w(f'{"fam":>3} {"caso":>7} {"masa":>6} {"Omega_min":>9} {"omega_d":>9} {"hueco":>9} {"|dOm/dJ|":>8} {"|c_+|":>9} '
      f'{"res":>6} {"|c_+| 2":>9} {"res":>6} {"w_b/sqrt(eps)":>13} {"eps_c":>8} {"J_r":>7} {"J_r/J_max":>9}')
    for caso, m in modos_w().items():
        hueco = float(m['om_min'] - m['wd'])
        wb1 = np.sqrt(2*abs(m['cmas'])*m['dom'])
        w(f'{"W":>3} {caso:>7} {float(m["masa"]):6g} {float(m["om_min"]):9.6f} {float(m["wd"]):9.6f} {hueco:9.2e} '
          f'{float(m["dom"]):8.4f} {abs(m["cmas"]):9.2e} {float(m["res"]):6.3f} {abs(m["cmas2"]):9.2e} '
          f'{float(m["res2"]):6.3f} {wb1:13.5f} {(hueco/(2*wb1))**2:8.4f} {JT_W + hueco/float(m["dom"]):7.4f} '
          f'{1 + hueco/float(m["dom"])/JT_W:9.4f}')
    rb = np.load(H.ruta('lineal', 'rebote.npz'))
    om_min = {e['g']: e['om_min'] for e in estados_datos() if e['fam'] == 'P' and e['a0'] == 1.0}
    for i, k in enumerate(rb['k']):
        hueco, wb1, dom = float(rb['hueco'][i]), float(rb['wb1'][i]), float(rb['dom'][i])
        w(f'{"P":>3} {f"k={k:g}":>7} {1:6g} {om_min[float(k)]:9.6f} {float(rb["wd"][i]):9.6f} {hueco:9.2e} {dom:8.4f} '
          f'{abs(rb["cmas"][i]):9.2e} {"":>6} {"":>9} {"":>6} {wb1:13.5f} {(hueco/(2*wb1))**2:8.4f} '
          f'{JT_P + hueco/dom:7.4f} {1 + hueco/dom/JT_P:9.4f}')
    escribir('modos.txt', lineas)


# ------------------------------------------------------------------ amplitud finita
def casos_finita():
    """Las corridas con modo de las dos familias: (familia, caso, corrida, eps, partículas, tiempo
    fiable, serie, tiempos y mayor acción de las partículas en cada uno, J_max del estado)."""
    lista = []
    for corrida in FINITA_W:
        caso, eps = D.info[corrida]['caso'], D.info[corrida]['eps']
        fz = np.load(D.ruta(corrida, 'fase.npz'))
        lista.append(('W', caso, corrida, eps, fz['J'].shape[1], float(fz['t'][-1]), serie_w(corrida), fz['t'],
                      fz['J'].max(axis=1).astype(float), JT_W))
    for nombre, tmax in FINITA_P:
        k, a0, eps, ref = datos_p(nombre)
        H._orbitas(nombre)
        o = np.load(H.ruta(nombre, 'orbitas.npz'))
        lista.append(('P', k, nombre, eps, particulas_p(nombre), tmax, serie_p(nombre), o['t'], o['J'].max(axis=1), JT_P))
    return lista


def finita():
    """Las corridas con modo discreto: s, la mayor acción de las partículas junto a la del
    péndulo, y el cociente de amplitudes |kappa| frente al tiempo en unidades del rebote, todo
    dentro del tiempo fiable t_f de cada corrida."""
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Corridas con modo discreto. J_+: mayor acción de las partículas hasta t_f y la del péndulo; dif: su diferencia '
      'en unidades de la excursión J_+ - J_max del péndulo.')
    w('kappa: cociente de amplitudes PIC/lineal (ventana de Hann de 200) en omega_b t0 = 1, 2, 3, 4 (-- si cae fuera de '
      't_f); mín y máx en 200 <= t0 <= t_f - 100; primer mínimo: el menor valor con omega_b t0 <= 6')
    w('y dónde (con * si queda a menos de 0.3/omega_b del final: sigue bajando). c: control de otra corrida.')
    w(f'{"fam":>3} {"caso":>6} {"corrida":>15} {"partíc":>7} {"eps":>6} {"s":>6} {"omega_b":>8} {"t_f":>5} {"w_b t_f":>7} '
      f'{"J_+ PIC":>8} {"péndulo":>8} {"dif":>6} ' + ' '.join(f'{f"k({x})":>6}' for x in (1, 2, 3, 4))
      + f' {"mín":>6} {"máx":>6} {"1er mín":>7} {"w_b t0":>7}')
    filas = []
    for fam, caso, corrida, eps, part, tf, (t, x, xl), tj, jm, jt in casos_finita():
        pd = pendulo(fam, caso, eps)
        t0 = np.arange(200.0, tf - 99.0, 10.0)
        kap = np.abs(kappas(t, x, xl, t0))
        en = [f'{abs(H.ganancia(t, x, xl, y/pd["wb"])):6.3f}' if 100 <= y/pd['wb'] <= tf - 100 else f'{"--":>6}'
              for y in (1, 2, 3, 4)]
        jmax = float(jm[tj <= tf + 1].max())
        v = pd['wb']*t0 <= 6
        i = int(np.argmin(kap[v]))
        marca = '*' if pd['wb']*(t0[-1] - t0[i]) < 0.3 else ' '
        etiqueta = caso if fam == 'W' else f'k={caso:g}'
        filas.append((fam, pd['s'], f'{fam:>3} {etiqueta:>6} {corrida:>15} {part:7d} {eps:6g} {pd["s"]:6.2f} '
                      f'{pd["wb"]:8.5f} {tf:5.0f} {pd["wb"]*tf:7.2f} {jmax:8.4f} {pd["Jmax"]:8.4f} '
                      f'{(jmax - pd["Jmax"])/(pd["Jmax"] - jt):+6.2f} ' + ' '.join(en)
                      + f' {kap.min():6.3f} {kap.max():6.3f} {kap[v][i]:7.3f} {pd["wb"]*t0[i]:6.2f}{marca}'
                      + (' c' if corrida in CONTROLES else '')))
    for fam in 'WP':
        for _, _, f in sorted(x for x in filas if x[0] == fam):
            w(f)
    escribir('finita.txt', lineas)


# ------------------------------------------------------------------ figuras
def figuras():
    """Figuras de la Sección 6 del artículo, en docs/articulo/figuras/."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    destino = os.path.join(RAIZ, 'docs', 'articulo', 'figuras')
    os.makedirs(destino, exist_ok=True)
    plt.rcParams.update({'font.size': 9})
    es = estados_datos()

    # 1. La ganancia: lambda_edge frente a la masa en las dos familias, y lambda(omega) junto al
    #    borde para los politropos de masa 1.
    fig, axs = plt.subplots(1, 3, figsize=(7, 2.55), constrained_layout=True)
    ax = axs[0]
    for grupo, jt, c, m in (('Westrecho', 0.083, 'C1', 's'), ('Wref', 0.138, 'C0', 'o'), ('Wancho', 0.276, 'C2', '^')):
        pts = sorted((e['a0'], e['lame']) for e in es if e['fam'] == grupo)
        nombre = f'W, J_max = {jt:g}'
        extra = sorted((float(l.split()[1]), float(l.split()[4])) for l in bloque_umbral(nombre))
        pts = sorted(set(pts) | set(extra))
        ax.loglog(*zip(*pts), '-', color=c, lw=0.9)
        ax.loglog(*zip(*[q for q in pts if q not in extra]), m, color=c, ms=3.5, label=f'$J_{{max}}={jt:g}$')
    ax.set_xlim(2e-3, 1.5); ax.set_ylim(0.03, 4)
    ax.set_title('lowered Maxwellians, $g=2$', fontsize=9)
    ax = axs[1]
    for k, c, m in ((1.25, 'C0', 'o'), (1.5, 'C1', 's'), (2.0, 'C2', '^'), (3.0, 'C3', 'v')):
        pts = sorted((e['a0'], e['lame']) for e in es if e['fam'] == 'P' and e['g'] == k)
        extra = sorted((float(l.split()[1]), float(l.split()[4])) for l in bloque_umbral(f'P, k = {k:g}'))
        todos = sorted(set(pts) | set(extra))
        ax.loglog(*zip(*todos), '-', color=c, lw=0.9)
        ax.loglog(*zip(*pts), m, color=c, ms=3.5, label=f'$k={k:g}$')
    ax.set_xlim(7e-3, 6); ax.set_ylim(0.01, 3)
    ax.set_title('polytropes', fontsize=9)
    for ax in axs[:2]:
        ax.axhline(1.0, color='k', lw=0.7, ls='--')
        ax.set_xlabel('$M_{gas}/M$'); ax.grid(alpha=0.3)
        ax.legend(fontsize=7, frameon=False, loc='lower right')
    axs[0].set_ylabel('$\\lambda_{edge}$')
    g = ganancia_datos()
    ax = axs[2]
    for j, (k, c) in enumerate(zip(g['k'], ('C4', 'C5', 'C0', 'C1', 'C2', 'C3'))):
        ax.loglog(g['delta'], g['lam'][j], color=c, lw=1.0, label=f'$k={k:g}$')
    ax.axhline(1.0, color='k', lw=0.7, ls='--')
    ax.set_xlim(1e-9, 0.5); ax.set_ylim(0.35, 300)
    ax.set_xlabel('$(\\Omega_{min}-\\omega)/(\\Omega_{max}-\\Omega_{min})$'); ax.set_ylabel('$\\lambda(\\omega)$')
    ax.set_title('polytropes, $M_{gas}/M=1$', fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=7, frameon=False, loc='upper right', ncol=2, columnspacing=1.0)
    for ax, letra in zip(axs, 'abc'):
        ax.text(0.04, 0.93, f'({letra})', transform=ax.transAxes, fontsize=9, va='top')
    fig.savefig(os.path.join(destino, 'resultados_ganancia.pdf')); plt.close(fig)

    # 2. Respuesta amortiguada: familia W bajo el umbral, con la mezcla libre; colas de los
    #    politropos de masa 0.01 (PIC hasta t = 2000, donde alcanza el ruido de la referencia);
    #    y la cola del politropo k = 2 con masa 1, que va como t^-k.
    fig, axs = plt.subplots(1, 3, figsize=(7, 2.75), constrained_layout=True)
    ax = axs[0]
    lw = libres_w()
    for (caso, corrida, col), masa in zip((('A1', 'D1', 'C0'), ('A3', 'D2', 'C1'), ('A4', 'D5', 'C2')), (0.0065, 0.042, 0.075)):
        (t, x), (tl, xl) = serie_w(corrida, comun=False)
        ax.semilogy(lw['t'], np.abs(lw[caso]), color=col, lw=0.8, ls=':')
        ax.semilogy(t, np.abs(x), color=col, lw=1.0, label=f'${masa:g}$')
        ax.semilogy(tl, np.abs(xl), 'k--', lw=0.6)
    ax.set_xlim(0, 9000); ax.set_ylim(3e-5, 0.4)
    ax.set_xlabel('$t$'); ax.set_ylabel('$|\\chi_1|$'); ax.grid(alpha=0.3)
    ax.set_title('lowered Maxwellians, $g=2$', fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc='upper right', title='$M_{gas}/M$', title_fontsize=7)
    ax = axs[1]
    c = colas_datos()
    for k, col in zip((0.75, 1.0, 1.5, 2.0), ('C4', 'C5', 'C1', 'C2')):
        (t, x), _ = serie_p(f'D_k{k:g}', comun=False)
        v = (t > 0) & (t <= 2000)
        ax.loglog(t[v], np.abs(x[v]), color=col, lw=1.0, label=f'$k={k:g}$')
        ax.loglog(c['t'][1:], np.abs(c[f'k{k:g}_3200_sg'][1:]), 'k--', lw=0.6)
        C, _ = cola_libre(k)
        tt = np.array([800.0, 8000.0])
        ax.loglog(tt, C*tt**-(k + 1), color=col, lw=0.8, ls=':')
    ax.set_xlim(20, 8000); ax.set_ylim(1e-8, 0.4)
    ax.set_xlabel('$t$'); ax.grid(alpha=0.3)
    ax.set_title('polytropes, $M_{gas}/M=0.01$', fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc='lower left')
    ax = axs[2]
    (t, x), _ = serie_p('D_k2_a1', comun=False)
    v = (t > 0) & (t <= 600)
    ax.loglog(t[v], H.envolvente(t[v], x[v], 40.0), color='C2', lw=1.0, label='PIC')
    hl = c['k2_3200_sg_a1']
    el = H.envolvente(c['t'], hl, 40.0)
    ax.loglog(c['t'][1:], el[1:], 'k--', lw=0.6, label='linearized')
    C, _ = cola_libre(2.0, 1.0)
    tt = np.array([1500.0, 8000.0])
    ax.loglog(tt, el[np.argmin(np.abs(c['t'] - 4000))]*(tt/4000)**-2.0*1.6, color='0.3', lw=0.8, ls='-.', label='$\\propto t^{-2}$')
    tt = np.array([300.0, 8000.0])
    ax.loglog(tt, C*tt**-3.0, color='C2', lw=0.8, ls=':', label='$C\\,t^{-3}$')
    ax.set_xlim(20, 8000); ax.set_ylim(1e-9, 0.4)
    ax.set_xlabel('$t$'); ax.grid(alpha=0.3)
    ax.set_title('polytrope, $k=2$, $M_{gas}/M=1$', fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc='lower left')
    for ax, letra in zip(axs, 'abc'):
        ax.text(0.5, 0.96, f'({letra})', transform=ax.transAxes, fontsize=9, va='top', ha='center')
    fig.savefig(os.path.join(destino, 'resultados_amortiguado.pdf')); plt.close(fig)

    # 3. Modos discretos: la envolvente de la respuesta no decae, y la respuesta misma al final
    #    de una corrida.
    fig, axs = plt.subplots(1, 3, figsize=(7, 2.7), constrained_layout=True)
    ax = axs[0]
    for (corrida, col), masa in zip((('D5', 'C2'), ('D3e01', 'C3'), ('D4', 'C4'), ('D9', 'C5')), (0.075, 0.097, 0.17, 0.5)):
        (t, x), (tl, xl) = serie_w(corrida, comun=False)
        ax.plot(t, H.envolvente(t, x, 100.0), color=col, lw=1.0, label=f'${masa:g}$')
        ax.plot(tl, H.envolvente(tl, xl, 100.0), 'k--', lw=0.6)
    ax.set_xlim(0, 7200); ax.set_ylim(0, 0.2)
    ax.set_title('lowered Maxwellians, $g=2$', fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc='upper left', title='$M_{gas}/M$', title_fontsize=7, ncol=2,
              columnspacing=1.0)
    ax = axs[1]
    for k, corrida, tf, col in ((0.75, 'D_k0.75_a1', 1000, 'C4'), (1.0, 'D_k1_a1', 1000, 'C5'), (1.25, 'DPe01_k1.25_a1', 4000, 'C0'),
                                (2.0, 'D_k2_a1', 600, 'C2')):
        (t, x), (tl, xl) = serie_p(corrida, comun=False)
        v = t <= tf
        ax.semilogy(t[v], H.envolvente(t[v], x[v], 50.0), color=col, lw=1.0, label=f'$k={k:g}$')
        ax.semilogy(tl, H.envolvente(tl, xl, 50.0), 'k--', lw=0.6)
    ax.set_xlim(0, 4000); ax.set_ylim(1e-5, 0.5)
    ax.set_title('polytropes, $M_{gas}/M=1$', fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc='center right')
    for ax in axs[:2]:
        ax.set_xlabel('$t$'); ax.set_ylabel('envelope of $|\\chi_1|$'); ax.grid(alpha=0.3)
    ax = axs[2]
    (t, x), (tl, xl) = serie_w('D4', comun=False)
    v, vl = t >= 5120, tl >= 5120
    ax.plot(tl[vl], xl[vl].real, 'k-', lw=0.7, label='linearized')
    ax.plot(t[v], x[v].real, 'o', color='C4', ms=2.5, mfc='none', mew=0.7, label='PIC')
    ax.set_xlim(5120, 5520); ax.set_ylim(-0.17, 0.21)
    ax.set_xlabel('$t$'); ax.set_ylabel('Re$\\,\\chi_1$'); ax.grid(alpha=0.3)
    ax.set_title('$M_{gas}/M=0.17$', fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc='upper left', ncol=2, columnspacing=1.0)
    for ax, letra in zip(axs, 'abc'):
        ax.text(0.97, 0.96 if letra != 'c' else 0.09, f'({letra})', transform=ax.transAxes, fontsize=9, va='top', ha='right')
    fig.savefig(os.path.join(destino, 'resultados_modos.pdf')); plt.close(fig)

    # 4. Amplitud finita: cociente de amplitudes frente al tiempo en unidades del rebote, las
    #    partículas del borde con la separatriz del péndulo, y la mayor acción frente a la del péndulo.
    casos = {q[2]: q for q in casos_finita()}
    fig, axs = plt.subplots(2, 2, figsize=(7, 5.3), constrained_layout=True)
    for ax, lista, titulo, sitio in ((axs[0, 0], (('DPe01_k1.25_a1', 'C2'), ('DP_k1.25_a1', 'C1'), ('DPe1_k1.25_a1', 'C0')),
                                      'polytrope, $k=1.25$, $M_{gas}/M=1$', 'upper right'),
                                     (axs[0, 1], (('D3e03', 'C2'), ('D3', 'C1'), ('D7', 'C0')),
                                      'lowered Maxwellian, $M_{gas}/M=0.097$', 'lower right')):
        for corrida, col in lista:
            fam, caso, _, eps, part, tf, (t, x, xl), tj, jm, jt = casos[corrida]
            pd = pendulo(fam, caso, eps)
            t0 = np.arange(150.0, tf - 99.0, 10.0)
            ax.plot(pd['wb']*t0, np.abs(kappas(t, x, xl, t0)), color=col, lw=1.0,
                    label=f'$\\varepsilon={eps:g}$ ($s={pd["s"]:.2f}$)')
        ax.set_xlabel('$\\omega_b t_0$'); ax.set_ylabel('$|\\kappa(t_0)|$'); ax.grid(alpha=0.3)
        ax.set_title(titulo, fontsize=9); ax.set_xlim(0, 14); ax.set_ylim(0.5, 1.25)
        ax.legend(fontsize=7, frameon=False, loc=sitio)
    ax = axs[1, 0]
    nombre, k, eps = 'DPe1_k1.25_a1', 1.25, 0.1
    o, pd = np.load(H.ruta(nombre, 'orbitas.npz')), pendulo('P', k, eps)
    n = int(np.argmin(np.abs(o['t'] - np.pi/pd['wb'])))            # medio periodo de rebote
    fase = (o['Q'][n] - pd['wd']*o['t'][n] - pd['fase_O'] + np.pi) % (2*np.pi) - np.pi
    ax.plot(fase, o['J'][n], '.', color='C0', ms=1.0, rasterized=True)
    ph = np.linspace(-np.pi, np.pi, 401)
    for sg in (1, -1):
        ax.plot(ph, pd['Jr'] + sg*pd['dJs']*np.abs(np.cos(ph/2)), 'k', lw=0.8)
    ax.axhline(JT_P, color='0.4', lw=0.7, ls='--'); ax.axhline(pd['Jr'], color='0.4', lw=0.7, ls=':')
    ax.set_xlim(-np.pi, np.pi); ax.set_ylim(0.62, 0.78)
    ax.set_xticks([-np.pi, 0, np.pi]); ax.set_xticklabels(['$-\\pi$', '0', '$\\pi$'])
    ax.set_xlabel('$\\vartheta$'); ax.set_ylabel('$J$')
    ax.set_title(f'$k=1.25$, $\\varepsilon={eps:g}$, $t={o["t"][n]:.0f}$', fontsize=9)
    ax = axs[1, 1]
    for fam, m, col, etiqueta in (('W', 'o', 'C3', 'lowered Maxwellians'), ('P', 's', 'C0', 'polytropes')):
        xs, ys = [], []
        for corrida, (f, caso, _, eps, part, tf, _, tj, jm, jt) in casos.items():
            if f == fam and corrida not in CONTROLES:
                pd = pendulo(fam, caso, eps)
                xs.append(pd['Jmax']/jt - 1); ys.append(float(jm[tj <= tf + 1].max())/jt - 1)
        ax.loglog(xs, ys, m, color=col, ms=4, mfc='none', mew=0.9, label=etiqueta)
    ax.loglog([1e-3, 1], [1e-3, 1], 'k--', lw=0.7)
    ax.set_xlim(3e-3, 0.5); ax.set_ylim(3e-3, 0.5)
    ax.set_xlabel('pendulum estimate of $J/J_{max}-1$'); ax.set_ylabel('largest $J/J_{max}-1$ of the particles')
    ax.grid(alpha=0.3); ax.legend(fontsize=7, frameon=False, loc='lower right')
    for ax, letra in zip(axs.flat, 'abcd'):
        ax.text(0.03, 0.97, f'({letra})', transform=ax.transAxes, fontsize=9, va='top')
    fig.savefig(os.path.join(destino, 'resultados_atrapamiento.pdf'), dpi=300); plt.close(fig)
    print('figuras en', destino)


def bloque_umbral(nombre):
    """Las líneas "masa ... lambda_edge = ..." de una familia en umbrales.txt."""
    lineas, dentro = [], False
    for l in open(ruta('umbrales.txt')):
        if l.startswith(('W,', 'P,')):
            dentro = l.startswith(nombre + ':')
        elif dentro and l.strip().startswith('masa'):
            lineas.append(l)
    return lineas


if __name__ == '__main__':
    pasos = {'estados': estados, 'umbrales': umbrales, 'ganancia': ganancia, 'respuesta': respuesta,
             'modos': modos, 'finita': finita, 'figuras': figuras}
    if len(sys.argv) != 2 or sys.argv[1] not in pasos:
        sys.exit(__doc__)
    os.makedirs(BASE, exist_ok=True)
    pasos[sys.argv[1]]()
