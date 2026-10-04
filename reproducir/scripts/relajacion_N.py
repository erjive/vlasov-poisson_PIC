"""Relajación por número finito de partículas en las corridas de referencia (eps = 0).

Para cada corrida, con el mapa ángulo-acción del equilibrio, mide en una muestra de
instantáneas (pesos w de las partículas, J_t el borde del soporte):

  d1(t)   = sum w |J - J(0)| / sum w / J_t        desplazamiento medio de la acción
  d2(t)   = sqrt(sum w (J - J(0))^2 / sum w) / J_t
  med(t)  = sum w (J - J(0)) / sum w / J_t        desplazamiento con signo
  W1(t)   = int |C_t(J) - C_0(J)| dJ / J_t        distancia de Wasserstein entre las
            distribuciones acumuladas de masa en J: el cambio neto de F(J)
  dF(t)   = (1/2) sum_b |m_b(t) - m_b(0)| / M     lo mismo con 40 celdas en [0, J_t]

d1 y d2 miden cuánto se mueve cada partícula (la difusión, que no se bloquea); W1 y dF
miden cuánto cambia la distribución (el flujo neto, que es lo que el bloqueo cinético
anula: Roule, Fouvry y Pichon 2022; Fouvry y Roule 2023). Si las partículas solo se
intercambian, d1 crece y W1 no.

Las corridas son de arranque silencioso y con pesos distintos por partícula, no una
muestra de Poisson de masas iguales: la ley en N que salga de aquí describe el ruido de
este esquema, no la teoría cinética.

Uso:  python3 relajacion_N.py            (escribe exe/relajacion/resumen.txt)
"""
import os, sys
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'                       # cuatro procesos, un hilo cada uno
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
EXE = os.path.abspath(os.path.join(AQUI, '..', '..', 'exe'))
SAL = os.path.join(EXE, 'relajacion')
NMUESTRA = 41                                 # instantáneas por corrida
# nombre, directorio de la corrida, equilibrio, partículas, descripción
CORRIDAS = [
    ('H_k1.25',    'hadzic/Z_k1.25_a1',    'hadzic/ic/Z_k1.25_a1',  'masa puntual, k=1.25, a0=1'),
    ('H_k1.25_N',  'hadzic/ZN_k1.25_a1',   'hadzic/ic/ZN_k1.25_a1', 'idem, 4N'),
    ('H_k1.25_dr', 'hadzic/Zdr_k1.25_a1',  'hadzic/ic/Z_k1.25_a1',  'idem, dr/2'),
    ('H_k1.5',     'hadzic/Z_k1.5_a1',     'hadzic/ic/Z_k1.5_a1',   'masa puntual, k=1.5, a0=1'),
    ('H_k1.5_N',   'hadzic/ZN_k1.5_a1',    'hadzic/ic/ZN_k1.5_a1',  'idem, 4N'),
    ('H_k1.5_dr',  'hadzic/Zdr_k1.5_a1',   'hadzic/ic/Z_k1.5_a1',   'idem, dr/2'),
    ('H_k2',       'hadzic/Z_k2_a1',       'hadzic/ic/Z_k2_a1',     'masa puntual, k=2, a0=1 (sin modo)'),
    ('H_k2_m',     'hadzic/Z_k2',          'hadzic/ic/Z_k2',        'masa puntual, k=2, a0=0.01'),
    ('A4',         'demo_eta/Z_A4',        'demo_eta/ic/Z_A4',      'isócrono, A4 (a0=0.075, amortiguado)'),
    ('A4_N',       'demo_eta/Z_A4N',       'demo_eta/ic/Z_A4N',     'idem, 4N'),
    ('L5',         'demo_eta/Z_L5',        'demo_eta/ic/Z_L5',      'isócrono, L5 (modo junto al borde)'),
    ('L5_N',       'demo_eta/Z_L5N',       'demo_eta/ic/Z_L5N',     'idem, 4N'),
    ('M5',         'demo_eta/edo/E_M5',    'demo_eta/edo/ic/E_M5',  'isócrono, M5 (a0=0.5, modo discreto)'),
    ('M5_N',       'demo_eta/edo/E_M5N',   'demo_eta/edo/ic/E_M5N', 'idem, 4N'),
    ('M5_dr',      'demo_eta/edo/E_M5dr',  'demo_eta/edo/ic/E_M5',  'idem, dr/2'),
    ('M5_drN',     'demo_eta/edo/E_M5drN', 'demo_eta/edo/ic/E_M5N', 'idem, dr/2 y 4N'),
]


def medir(caso):
    import h5py
    from aa_numerico import MapaAA
    nombre, corrida, ic, _ = caso
    sal = os.path.join(SAL, nombre + '.npz')
    if os.path.exists(sal):
        return nombre
    eq = np.load(os.path.join(EXE, ic + '_equilibrio.npz'))
    fondo = str(eq['fondo']) if 'fondo' in eq.files else 'isocrono'
    l0 = float(eq['L0']) if 'L0' in eq.files else 2.0
    jt = float(eq['J_borde'])
    mapa = MapaAA(eq['r'], eq['phi_self'], L=l0, fondo=fondo)
    f = h5py.File(os.path.join(EXE, corrida, 'vlasov_output.h5'), 'r')
    pasos = sorted([c for c in f if c.startswith('step_')], key=lambda c: int(c.split('_')[1]))
    elegidos = [pasos[i] for i in np.unique(np.linspace(0, len(pasos) - 1, NMUESTRA).round().astype(int))]
    w = f[pasos[0]]['f'][:]; M = w.sum()
    malla = np.linspace(0.0, 1.2*jt, 2401)                    # para las acumuladas
    bordes = np.linspace(0.0, jt, 41)

    def acumulada(J):
        o = np.argsort(J)
        return np.interp(malla, J[o], np.cumsum(w[o])/M, left=0.0, right=1.0)

    def celdas(J):
        return np.histogram(J, bins=np.r_[bordes, np.inf], weights=w)[0]/M

    filas = []
    for c in elegidos:
        g = f[c]
        _, J, _ = mapa(g['r_part'][:], g['p_part'][:])
        if c == elegidos[0]:
            J0, C0, m0 = J.copy(), acumulada(J), celdas(J)
        d = J - J0
        filas.append((g.attrs['time'], np.sum(w*np.abs(d))/M/jt, np.sqrt(np.sum(w*d*d)/M)/jt, np.sum(w*d)/M/jt,
                      np.sum(np.abs(acumulada(J) - C0))*(malla[1] - malla[0])/jt, 0.5*np.sum(np.abs(celdas(J) - m0))))
    f.close()
    a = np.array(filas)
    np.savez(sal, t=a[:, 0], d1=a[:, 1], d2=a[:, 2], med=a[:, 3], W1=a[:, 4], dF=a[:, 5], N=len(w), jt=jt)
    return nombre


if __name__ == '__main__':
    from multiprocessing import Pool
    os.makedirs(SAL, exist_ok=True)
    with Pool(4) as pool:
        for n in pool.imap_unordered(medir, CORRIDAS):
            print('medida', n, flush=True)
    lineas = ['Corridas de referencia (eps = 0). Todo en unidades de J_t (1e-3). W1 y dF: cambio de F(J).',
              f'{"corrida":>11} {"N":>6} {"t":>6} {"d1":>8} {"d2":>8} {"|med|":>8} {"W1":>8} {"dF":>8}   descripción']
    for nombre, _, _, desc in CORRIDAS:
        z = np.load(os.path.join(SAL, nombre + '.npz')); t = z['t']
        for tt in (t[-1]/4, t[-1]/2, t[-1]):
            i = int(np.argmin(np.abs(t - tt)))
            lineas.append(f'{nombre:>11} {int(z["N"]):6d} {t[i]:6.0f} {1e3*z["d1"][i]:8.3f} {1e3*z["d2"][i]:8.3f} '
                          f'{1e3*abs(z["med"][i]):8.4f} {1e3*z["W1"][i]:8.3f} {1e3*z["dF"][i]:8.3f}   {desc if tt == t[-1]/4 else ""}')
    open(os.path.join(SAL, 'resumen.txt'), 'w').write('\n'.join(lineas) + '\n')
    print('\n'.join(lineas))
