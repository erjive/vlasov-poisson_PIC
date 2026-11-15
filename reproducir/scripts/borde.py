"""La frontera de la FD y la relajación: exploración con el modelo epicíclico (propuesta del 2026-11-15).

Pregunta: ¿decide el comportamiento de la FD en la frontera de su soporte si una nube autogravitante se relaja,
como decide el amortiguamiento de las perturbaciones de un estado estacionario (Hadžić, Rein, Schrecker y Straub:
FD ~ (E_0 - E)^k, amortiguado si k > 1, modo que no se amortigua si 1/2 < k <= 1)? Tres pruebas con el modelo de
reducido.py (clase Modelo, esquema de Strang), todas cortas y en el portátil:

  A  estado estacionario en Kepler (L_0 = 2), FD = (E(J_0) - E(J))^k en 0 < J < J_0, uniforme en el ángulo, con un
     desplazamiento rígido pequeño de todas las partículas en el plano de zeta (perturbación m = 1). Se sigue el
     centro de masa C(t) = sum m zeta / M: ¿decae, o queda oscilando con una frecuencia fuera de la banda?
     Con g = 0 (sin autogravedad) C decae por la sola mezcla: es la referencia.
  B  la nube lejos del equilibrio en Kepler (ancho en ángulo s_Q = 0.1) con la FD J^2 (E(J_0) - E(J))^k, la misma
     frontera exterior, frente a la nube gaussiana del artículo: ¿se relaja?
  C  la nube de referencia del artículo (isócrono, J^p exp(-J^2/s_J^2), s_Q = s_J = 0.1) con otro exponente p en
     la órbita circular: ¿cambian el núcleo y su oscilación?

    python3 borde.py correr [A|B|C ...]     corre lo que falte (en serie, cuatro hilos)
    python3 borde.py resumen                 tabla de resultados (exe/demo_eta/reducido/borde/resumen.txt)
"""
import os, sys, time
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import reducido as R

SAL = os.path.join(R.DATOS, 'borde')
L0 = 2.0
GUARDA = 5                                      # una instantánea de (X, Y) cada GUARDA salidas (50 unidades de tiempo)


def energia(J, pot):
    return -0.5/(J + R.c_de_L(L0, pot))**2


def reticula(nj, nq, perfil, pot, J0=None, k=None, p=2, sj=0.1, sq=None, umbral=1e-7):
    """Partículas (Q, J, peso) en la retícula de puntos medios: perfil 'hrss' = (E(J0) - E(J))^k, 'borde' = J^p (E(J0)
    - E(J))^k, ambos en 0 < J < J0; 'gauss' = J^p exp(-J^2/sj^2) en 0 < J < 6 sj. En el ángulo, uniforme (sq None) o
    exp(-sin^2(Q/2)/sq^2). Se omiten los puntos con peso menor que umbral por el mayor."""
    if perfil in ('hrss', 'borde'):
        J = (np.arange(nj) + 0.5)*J0/nj
        fJ = (energia(J0, pot) - energia(J, pot))**k*(J**p if perfil == 'borde' else 1.0)
    else:
        J = (np.arange(nj) + 0.5)*6*sj/nj
        fJ = J**p*np.exp(-(J/sj)**2)
    Q = (np.arange(nq) + 0.5)*2*np.pi/nq
    fQ = np.ones(nq) if sq is None else np.exp(-np.sin(0.5*Q)**2/sq**2)
    W = fJ[:, None]*fQ[None, :]
    QQ, JJ = np.meshgrid(Q, J)
    ok = W >= umbral*W.max()
    return QQ[ok], JJ[ok], W[ok]


# nombre: (prueba, perfil, potencial, M_gas/M, k, p, J0, s_Q, desplazamiento, g (None = el del modelo), N_J, N_Q)
CASOS = {}
for kk in (0.5, 1.0, 1.5, 2.0, 3.0):
    for mm in (0.01, 0.03):
        CASOS[f'A_k{kk:g}_M{mm:g}'] = ('A', 'hrss', 'kepler', mm, kk, 0, 0.2, None, 0.02, None, 160, 64)
for kk in (0.5, 2.0):
    CASOS[f'A_k{kk:g}_libre'] = ('A', 'hrss', 'kepler', 0.01, kk, 0, 0.2, None, 0.02, 0.0, 160, 64)
# Segunda tanda (2026-11-15): la prueba A con masas mayores, bordes intermedios y el doble de tiempo.
TFIN = {}
for mm in (0.03, 0.1):
    for kk in (0.5, 0.75, 1.0, 1.25):
        CASOS[f'A2_k{kk:g}_M{mm:g}'] = ('A', 'hrss', 'kepler', mm, kk, 0, 0.2, None, 0.02, None, 160, 64)
        TFIN[f'A2_k{kk:g}_M{mm:g}'] = 2.0e4
CASOS['A2_k2_M0.1'] = ('A', 'hrss', 'kepler', 0.1, 2.0, 0, 0.2, None, 0.02, None, 160, 64); TFIN['A2_k2_M0.1'] = 2.0e4
for mm in (0.01, 0.02):
    for kk in (0.5, 1.0, 2.0):
        CASOS[f'B_k{kk:g}_M{mm:g}'] = ('B', 'borde', 'kepler', mm, kk, 2, 0.25, 0.1, 0.0, None, 400, 100)
    CASOS[f'B_gauss_M{mm:g}'] = ('B', 'gauss', 'kepler', mm, None, 2, None, 0.1, 0.0, None, 400, 100)
for pp in (0, 1, 2, 4):
    CASOS[f'C_p{pp}_M0.03'] = ('C', 'gauss', 'isocrono', 0.03, None, pp, None, 0.1, 0.0, None, 400, 100)


def correr(nombre, dt=5.0, tfin=1.0e4):
    prueba, perfil, pot, masa, k, p, J0, sq, desp, g, nj, nq = CASOS[nombre]
    Q, J, w = reticula(nj, nq, perfil, pot, J0, k, p, sq=sq)
    M = R.Modelo(Q, J, w, masa, g, 'strang', True, pot, None, L0)
    M.x = M.x + desp
    nsal = int(round(tfin/R.SALIDA)); npas = int(round(R.SALIDA/dt))
    C = np.empty(nsal + 1, complex); h1 = np.empty(nsal + 1, complex)
    X, Y, tg = [], [], []
    W = M.m.sum()

    def mide(i):
        z = M.x + 1j*M.y
        C[i] = (M.m*z).sum()/W
        h1[i] = (M.m*np.conj(z)/(np.abs(z) + 1e-300)).sum()/W
        if i % GUARDA == 0:
            X.append(M.x.astype(np.float32)); Y.append(M.y.astype(np.float32)); tg.append(i*R.SALIDA)
    mide(0)
    inv0 = M.invariantes(); t0 = time.time()
    for i in range(1, nsal + 1):
        for _ in range(npas):
            M.paso(dt)
        mide(i)
    inv1 = M.invariantes()
    os.makedirs(SAL, exist_ok=True)
    np.savez(os.path.join(SAL, nombre + '.npz'), t=R.SALIDA*np.arange(nsal + 1), C=C, h1=h1, tg=np.array(tg), X=np.array(X),
             Y=np.array(Y), m=M.m, c=M.c, kappa=M.kappa, g=M.g, dH=(inv1[0] + inv1[1] - inv0[0] - inv0[1])/inv0[1] if inv0[1] else 0.0)
    print(f'  {nombre}: {M.n} partículas, {time.time() - t0:.0f} s', flush=True)


def banda(d):
    """Frecuencias de las partículas (marco de kappa_0) en el último tercio, de la pendiente del ángulo entre
    instantáneas guardadas: cuantiles 1 % y 99 % por masa, y la media."""
    tg, X, Y, m = d['tg'], d['X'].astype(float), d['Y'].astype(float), d['m']
    k0 = int(len(tg)*2/3)
    ang = np.arctan2(Y[k0:], X[k0:])
    dang = np.angle(np.exp(1j*np.diff(ang, axis=0)))
    om = dang.sum(0)/(tg[-1] - tg[k0])
    o = np.argsort(om); cm = np.cumsum(m[o])/m.sum()
    return om[o][np.searchsorted(cm, 0.01)], om[o][np.searchsorted(cm, 0.99)]


def frecuencia(t, z, desde):
    sel = t >= desde
    ph = np.unwrap(np.angle(z[sel]))
    return np.polyfit(t[sel], ph, 1)[0]


def resumen():
    lineas = ['Exploración de la frontera de la FD con el modelo epicíclico. Frecuencias en el marco que rota con kappa_0.']
    for prueba in 'ABC':
        nombres = [n for n in CASOS if CASOS[n][0] == prueba and os.path.exists(os.path.join(SAL, n + '.npz'))]
        if not nombres:
            continue
        lineas.append('')
        if prueba == 'A':
            lineas.append('A  estado estacionario (E(J0) - E(J))^k en Kepler, J0 = 0.2, con desplazamiento 0.02: |C| inicial, medio en 20-40 % y en '
                          'el último 20 % de la corrida (t = 1e4), frecuencia de C en el último 20 % y banda de las partículas (1 % - 99 %)')
        elif prueba == 'B':
            lineas.append('B  nube lejos del equilibrio en Kepler, s_Q = 0.1: |h_1| medio en 1000-2000, en 3000-5000 y en el último 20 %; |C| igual')
        else:
            lineas.append('C  nube de referencia (isócrono, J^p exp(-J^2/0.01)), M = 0.03: |h_1| tardío, m_c, |C| tardío y su frecuencia')
        for n in nombres:
            d = np.load(os.path.join(SAL, n + '.npz')); t = d['t']; C, h1 = d['C'], d['h1']
            fin = t >= 0.8*t[-1]; med = (t >= 0.2*t[-1]) & (t <= 0.4*t[-1])
            if prueba == 'A':
                lo, hi = banda(d); w = frecuencia(t, C, 0.8*t[-1])
                lineas.append(f'  {n:16} |C| {abs(C[0]):.4f} -> {np.abs(C[med]).mean():.2e} -> {np.abs(C[fin]).mean():.2e} '
                              f'(cociente {np.abs(C[fin]).mean()/abs(C[0]):.2e})  w_C = {w:+.5f}  banda [{lo:+.5f}, {hi:+.5f}]'
                              f'  {"FUERA" if (w < lo or w > hi) else "dentro"}  dH/U {float(d["dH"]):+.1e}')
            elif prueba == 'B':
                a = (t >= 1000) & (t <= 2000); b = (t >= 3000) & (t <= 5000)
                lineas.append(f'  {n:16} |h_1| {np.abs(h1[a]).mean():.2e}  {np.abs(h1[b]).mean():.2e}  {np.abs(h1[fin]).mean():.2e}'
                              f'   |C| {np.abs(C[a]).mean():.2e}  {np.abs(C[b]).mean():.2e}  {np.abs(C[fin]).mean():.2e}  dH/U {float(d["dH"]):+.1e}')
            else:
                tg = d['tg']; M = R.medidas_plano(tg, d['X'], d['Y'], d['m'], 0.1, R.tau1(), float(d['kappa']))
                lineas.append(f'  {n:16} |h_1| (t >= 20 tau_1) {np.abs(h1[t >= 20*R.tau1()]).mean():.2e}  m_c {M["mnuc"]:.3f}'
                              f'  |C| tardío {np.abs(C[fin]).mean():.3f}  w_C = {frecuencia(t, C, 0.8*t[-1]):+.5f}')
    texto = '\n'.join(lineas)
    open(os.path.join(SAL, 'resumen.txt'), 'w').write(texto + '\n')
    print(texto)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'correr':
        pruebas = sys.argv[2:] or ['A', 'B', 'C']
        for n in CASOS:
            if CASOS[n][0] in pruebas and not os.path.exists(os.path.join(SAL, n + '.npz')):
                correr(n, tfin=TFIN.get(n, 1.0e4))
    elif len(sys.argv) > 1 and sys.argv[1] == 'resumen':
        resumen()
    else:
        print(__doc__)
