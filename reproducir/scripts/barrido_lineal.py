"""Barrido en eta con la teoría lineal: gamma(eta), omega(eta) y la posición del
modo respecto de la banda, sin partículas.

Para cada caso: equilibrio autoconsistente (equilibrio.py), banda y eta
(eta.banda), Vlasov-Poisson linealizado (lineal.resolver) y ajuste de la cola
de h_1 con matrix pencil (landau_cola.ajustar, el de 11_landau) en ventanas
[a, 2a] tau_1 desde a = 1, con tau_1 = 2 pi/(ancho de la banda). Se reporta el par de
ventanas consecutivas que más coinciden y su discrepancia: si ninguna coincide
(cola algebraica de la mezcla de fases, dos polos comparables), no hay polo
dominante y el caso se marca así. Entre los pares que coinciden se toma el más
tardío, que es el comportamiento asintótico: un modo discreto puede quedar tapado
al principio por una componente amortiguada del continuo.

El mismo ajuste se repite sobre dPhi(r, t), que no depende de la función de
prueba B(J) de h_1: la serie temporal dominante de su descomposición en valores
singulares en la ventana, con matrix pencil sobre la señal real. El polo es
propiedad del operador y debe coincidir en los dos observables; si no, la
discrepancia es del ajuste (componentes que compiten en la cola).

Además estima la frecuencia de rebote con eps = 1 a partir del dPhi lineal,
omega_b = sqrt(|dOmega/dJ| max_r|dPhi|), en t = 0 y al comienzo de la ventana
del polo: con eps cualquiera, omega_b escala como sqrt(eps) y el parámetro de
O'Neil nu = omega_b/gamma decide si la teoría lineal vale a tiempos largos.

La posición x = (omega - Omega_min)/(Omega_max - Omega_min) dice dónde está el
modo: 0 < x < 1 dentro de la banda (resonante, amortiguado); x < 0 por debajo
(discreto, gamma debe ser ~0).

Uso:  python3 barrido_lineal.py [--casos A4 E1 ...] [--dir exe/eta_lineal]
      [--nj 1600] [--nq 32] [--dt 0.5]
"""
import os, sys, argparse, subprocess, time, warnings, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eta import banda
from lineal import resolver
from landau_cola import ajustar, matrix_pencil

warnings.filterwarnings('ignore', category=RuntimeWarning)
AQUI = os.path.dirname(os.path.abspath(__file__))
# En tau_1, desde t = tau_1: antes domina el transitorio de la mezcla de fases,
# que no es un polo. Las tardías solo se usan si la corrida las cubre.
VENTANAS = ((1, 2), (2, 4), (3, 6), (5, 9), (10, 20), (20, 40), (30, 60))
PISO = 1e-9                                      # |h_1| mínimo en la ventana
NMAX = 1000       # muestras por ventana: matrix pencil hace una SVD de N/2 x N/2


def wilson(jt, a0, w0=3.0, l0=2.0, g=2.0):
    return dict(forma='maxwell', jt=jt, k=g, w0=w0, l0=l0, a0=a0, pert='suave3')


def hueca(jt, a0, l0=2.0):
    return dict(forma='politropo', jt=jt, k=2.0, m=2.0, l0=l0, a0=a0, pert='plana')


def gauss(sigma, a0, l0=2.0):
    return dict(forma='gauss', sigma=sigma, l0=l0, a0=a0, pert='plana')


# Las corridas del documento, y puntos intermedios del barrido en eta (L*).
CASOS = {
    'L1': wilson(0.138, 0.0030), 'A1': wilson(0.138, 0.0065), 'L2': wilson(0.138, 0.0130),
    'A2': wilson(0.138, 0.020), 'L3': wilson(0.138, 0.031), 'A3': wilson(0.138, 0.042),
    'L4': wilson(0.138, 0.058), 'A4': wilson(0.138, 0.075), 'L5': wilson(0.138, 0.097),
    'A5': wilson(0.138, 0.121), 'L6': wilson(0.138, 0.170),
    'A6': wilson(0.276, 0.044), 'A7': wilson(0.276, 0.148),
    'A8': wilson(0.083, 0.0109), 'A9': wilson(0.083, 0.044),
    'D1': wilson(0.094, 0.018, l0=1.0), 'D2': wilson(0.094, 0.069, l0=1.0),
    'E1': hueca(0.138, 0.075), 'E2': gauss(0.05, 0.069), 'E3': wilson(0.138, 0.073, w0=6.0),
    'H9': hueca(0.083, 0.041), 'H3': hueca(0.083, 0.10),
    # El exponente del borde (el k de Hadžić) a eta = 1; g = 2 es A4 y g = 1 es King.
    'G075': wilson(0.138, 0.0775, g=0.75), 'G1': wilson(0.138, 0.0771, g=1.0),
    'G15': wilson(0.138, 0.0762, g=1.5), 'G3': wilson(0.138, 0.0739, g=3.0),
    # Bordes abruptos a eta chico: ¿hay modo discreto con acoplamiento débil?
    'G075a': wilson(0.138, 0.0070, g=0.75), 'G075b': wilson(0.138, 0.0215, g=0.75),
    'G1a': wilson(0.138, 0.0069, g=1.0), 'G1b': wilson(0.138, 0.0212, g=1.0),
    'landau': gauss(0.10, 0.01),                 # validación: 11_landau
}


def polo_dphi(t, dphi, lo, hi, M=4):
    """Polo dominante de dPhi(r, t) en [lo, hi]: serie temporal principal (SVD) y
    matrix pencil sobre la señal real; los polos vienen en pares conjugados."""
    v = (t >= lo) & (t <= hi)
    U, S, _ = np.linalg.svd(dphi[v], full_matrices=False)
    y = U[:, 0]*S[0]
    sp, a = matrix_pencil(y.astype(complex), t[1] - t[0], M)
    pos = -sp.imag > 0
    if not pos.any():
        return None
    k = np.flatnonzero(pos)[np.argmax(np.abs(a[pos]))]
    return -sp[k].imag, -sp[k].real


def mejor_par(ajustes, ancho, tol=0.2):
    """El par de ventanas consecutivas MÁS TARDÍO que coincide (discrepancia < tol):
    es el comportamiento asintótico. Si ninguno coincide, el de menor discrepancia.
    Devuelve (discrepancia, i, omega, gamma)."""
    pares = []
    for i in range(len(ajustes) - 1):
        if ajustes[i] is None or ajustes[i+1] is None:
            continue
        (w1, g1), (w2, g2) = ajustes[i], ajustes[i+1]
        disc = max(abs(w1 - w2)/ancho, abs(g1 - g2)/max(0.5*abs(g1 + g2), 2e-3*ancho))
        pares.append((disc, i, 0.5*(w1 + w2), 0.5*(g1 + g2)))
    if not pares:
        return None
    buenos = [p for p in pares if p[0] < tol]
    return max(buenos, key=lambda p: p[1]) if buenos else min(pares)


def equilibrio(nombre, c, carpeta):
    """Genera (si falta) el equilibrio del caso con equilibrio.py y devuelve su npz."""
    dat = os.path.join(carpeta, f'{nombre}.dat')
    npz = os.path.join(carpeta, f'{nombre}_equilibrio.npz')
    if os.path.exists(npz):
        return npz
    cmd = [sys.executable, os.path.join(AQUI, 'equilibrio.py'), '--forma', c['forma'],
           '--a0', str(c['a0']), '--l0', str(c['l0']), '--eps', '0.1', '--pert', c['pert'],
           '--nrc', '400', '--npc', '25', '--salida', dat]
    if c['forma'] == 'gauss':
        cmd += ['--sigma', str(c['sigma'])]
        if c['sigma'] != 0.10:
            cmd += ['--jmax', str(6*c['sigma'])]
    else:
        cmd += ['--jt', str(c['jt']), '--k', str(c['k'])]
        cmd += ['--w0', str(c['w0'])] if c['forma'] == 'maxwell' else ['--m', str(c['m'])]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return npz


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--casos', nargs='+', default=list(CASOS))
    ap.add_argument('--dir', default=os.path.join(AQUI, '..', '..', 'exe', 'eta_lineal'))
    ap.add_argument('--nj', type=int, default=1600)
    ap.add_argument('--nq', type=int, default=32)
    ap.add_argument('--dt', type=float, default=0.5)
    ap.add_argument('--ntau', type=float, default=10.0, help='duración en tau_1')
    ap.add_argument('--pb', type=float, default=2.0, help='exponente de J en B(J) de h_1')
    ap.add_argument('--j1f', type=float, default=None,
                    help='centro y ancho de B(J) como fracción de J_t (por omisión 1/2.763)')
    arg = ap.parse_args()
    os.makedirs(arg.dir, exist_ok=True)
    cab = (f'{"caso":>7} {"a0":>7} {"eta":>6} {"banda Omega":>19} {"tau_1":>6} | '
           f'{"omega":>8} {"gamma":>9} {"x":>6} {"gamma/dOm":>9} {"disc.":>6} {"vent.":>9} | '
           f'{"omega_P":>8} {"gamma_P":>9} {"d_obs":>6} | {"wb(0)":>8} {"wb(tv)":>8} | {"seg":>4}')
    print(cab); print('-'*len(cab))
    for nombre in arg.casos:
        c = CASOS[nombre]
        t0 = time.time()
        npz = equilibrio(nombre, c, arg.dir)
        kw = {k: c[k] for k in ('forma', 'jt', 'k', 'm', 'w0') if k in c}
        salida = os.path.join(arg.dir, f'{nombre}_lineal.npz')
        if os.path.exists(salida) and 'eta' in np.load(salida).files:
            d = np.load(salida)
            eta, b = float(d['eta']), dict(Om_min=float(d['Om_min']), Om_max=float(d['Om_max']))
            j1 = float(d['j1'])
        elif c['forma'] == 'gauss':
            jmax = 6*c['sigma']
            b0 = banda(0.0, c['sigma'], jmax, c['l0']); b = banda(c['a0'], c['sigma'], jmax, c['l0'])
            j1 = c['sigma'] if c['sigma'] != 0.10 else 0.10
        else:
            b0 = banda(0.0, None, c['jt'], c['l0'], **kw); b = banda(c['a0'], None, c['jt'], c['l0'], **kw)
            j1 = c['jt']/2.763 if arg.j1f is None else c['jt']*arg.j1f
        if 'Om_media' in b:
            eta = abs(b['Om_media'] - b0['Om_media'])/b0['Om_media']/((b['Om_max'] - b['Om_min'])/b['Om_media'])
        ancho = b['Om_max'] - b['Om_min']
        tau1 = 2*np.pi/ancho
        if os.path.exists(salida):
            d = np.load(salida); t, h1, dphi = d['t'], d['h1'], d['dphi']
        else:
            t, h1, h2, rmed, dphi = resolver(npz, arg.nj, arg.nq, arg.dt, arg.ntau*tau1,
                                             verboso=False, j1=j1, sj1=j1, pb=arg.pb)
            np.savez(salida, t=t, h1=h1, h2=h2, r=rmed, dphi=dphi, nj=arg.nj, nq=arg.nq,
                     dt=arg.dt, j1=j1, sj1=j1, eta=eta, Om_min=b['Om_min'], Om_max=b['Om_max'])
        h1n = np.abs(h1)/np.abs(h1[0])
        pn = np.sqrt(np.mean(dphi**2, axis=1)); pn = pn/pn[0]
        ajustes, ajustes_P = [], []
        vent_caso = [(lo, hi) for lo, hi in VENTANAS if hi*tau1 <= t[-1] + 1e-9]
        for lo, hi in vent_caso:
            v = (t >= lo*tau1) & (t <= hi*tau1)
            ok = v.sum() >= 20
            m = max(1, int(np.ceil(v.sum()/NMAX)))     # submuestreo: omega*m*dt << pi
            ajustes.append(ajustar(t[::m], h1[::m], lo*tau1, hi*tau1)['pencil M=3']
                           if ok and h1n[v].min() >= PISO else None)
            ajustes_P.append(polo_dphi(t[::m], dphi[::m], lo*tau1, hi*tau1)
                             if ok and pn[v].min() >= PISO else None)
        mejor = mejor_par(ajustes, ancho)
        mejor_P = mejor_par(ajustes_P, ancho)
        # omega_b con eps = 1: sqrt(|dOmega/dJ| max_r |dPhi|), |dOmega/dJ| ~ ancho/J_soporte
        jsop = c.get('jt', 2.763*c.get('sigma', 0.1))
        wb = lambda n: np.sqrt(ancho/jsop*np.max(np.abs(dphi[n])))
        if mejor is None:
            print(f'{nombre:>7} {c["a0"]:7.4f} {eta:6.3f} [{b["Om_min"]:.5f},{b["Om_max"]:.5f}] {tau1:6.0f} | '
                  f'{"sin ventanas válidas":>52} | {wb(0):8.2e} {"":>8} | {time.time()-t0:4.0f}', flush=True)
            continue
        disc, i, w, g = mejor
        x = (w - b['Om_min'])/ancho
        ntv = np.argmin(np.abs(t - vent_caso[i][0]*tau1))
        vent = f'{vent_caso[i][0]:g}-{vent_caso[i+1][1]:g}'
        if mejor_P is None:
            wP = gP = dobs = np.nan
        else:
            _, _, wP, gP = mejor_P
            dobs = max(abs(w - wP)/ancho, abs(g - gP)/max(0.5*abs(g + gP), 2e-3*ancho))
        marca = '' if (disc < 0.2 and dobs < 0.2) else ' *'
        print(f'{nombre:>7} {c["a0"]:7.4f} {eta:6.3f} [{b["Om_min"]:.5f},{b["Om_max"]:.5f}] {tau1:6.0f} | '
              f'{w:8.5f} {g:9.2e} {x:6.2f} {g/ancho:9.3f} {disc:6.2f} {vent:>9} | '
              f'{wP:8.5f} {gP:9.2e} {dobs:6.2f} | '
              f'{wb(0):8.2e} {wb(ntv):8.2e} | {time.time()-t0:4.0f}{marca}', flush=True)
    print('\nx = (omega - Omega_min)/ancho. disc. = discrepancia entre las dos ventanas consecutivas que más '
          'coinciden (vent., en tau_1),\nmáximo de |d omega|/ancho y |d gamma|/gamma; * = mayor que 0.2: sin polo '
          'dominante. omega_P, gamma_P: el polo en dPhi; d_obs: su discrepancia con el de h_1.\n'
          'wb = omega_b con eps = 1 en t = 0 y al comienzo de la ventana de h_1.')


if __name__ == '__main__':
    main()
