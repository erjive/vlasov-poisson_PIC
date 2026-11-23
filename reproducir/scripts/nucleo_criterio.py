"""Comprobaciones del criterio del núcleo de la Sección V con las corridas ya hechas.

No lanza simulaciones ni escribe en las carpetas de las corridas. La salida queda en
exe/demo_eta/bulto/criterio/.

Paso nulo: el valor de m_c en un gas mezclado a acción fija. Sin autogravedad cada fila J_a de la
retícula gira rígidamente con Omega(J_a), y la solución es exacta. Se evalúa el estimador de m_c (la
elipse alrededor del máximo de densidad, como bulto.nucleo) sobre esa mezcla libre, con la retícula y en
los instantes de cada corrida, y se compara con el valor de la corrida, para J_* = 0.025, 0.05 y 0.1 sigma_J.
"""
import os, sys, struct, zipfile, argparse, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import bulto as B
import mu_convenciones as C

SAL = B.ruta('criterio')
XS = (0.025, 0.05, 0.1)                        # J_*/sigma_J
MMIX = 4/(np.e*np.sqrt(np.pi))                 # m_mix = (J_*/sigma_J) MMIX


def mapa(n, clave):
    """Arreglo clave de <n>/bulto.npz como memmap (el .npz no está comprimido): no carga las 1001
    instantáneas en la memoria."""
    ruta = B.ruta(n, 'bulto.npz')
    with zipfile.ZipFile(ruta) as z:
        inf = z.getinfo(clave + '.npy')
    with open(ruta, 'rb') as f:
        f.seek(inf.header_offset)
        h = f.read(30)
        nl, xl = struct.unpack('<HH', h[26:30])
        f.seek(inf.header_offset + 30 + nl + xl)
        v = np.lib.format.read_magic(f)
        lee = np.lib.format.read_array_header_1_0 if v == (1, 0) else np.lib.format.read_array_header_2_0
        forma, fortran, tipo = lee(f)
        pos = f.tell()
    return np.memmap(ruta, dtype=tipo, mode='r', offset=pos, shape=forma, order='F' if fortran else 'C')


def pot(n):
    return 'kepler' if B.FONDO.get(n) == 'puntual' else 'isocrono'


NUEVAS = ['B011', 'B012', 'B013', 'B014', 'B015q', 'B016', 'B017', 'B018', 'B019']    # tandas 6 y 7 (2026-11-22): 400 x 100


def corridas():
    return list(dict.fromkeys([x for _, a, b in B.FAMILIAS for x in a + b] + NUEVAS))


def reticula(n):
    """Filas J, columnas Q y pesos F(Q_b, J_a) del dato inicial de la corrida n."""
    sj, sq = B.parametro(n, 'sr', 0.1), B.parametro(n, 'sp', 0.1)
    if n in B.RETICULA:
        J, Q = B.recorte_reticula(*B.RETICULA[n], sj, sq)
    else:
        nj, nq = B.parametro(n, 'Nrc', 400), B.info[n]['npc']
        J, Q = (np.arange(nj) + 0.5)*6*sj/nj, (np.arange(nq) + 0.5)*2*np.pi/nq
    F = np.outer(J**2*np.exp(-J**2/sj**2), np.exp(-np.sin(0.5*Q)**2/sq**2))
    return J, Q, F


def fase(Q, J, fondo, L=2.0):
    """(r, p_r) de los puntos (Q, J) en el potencial externo, como bulto.dato_isocrono y dato_kepler."""
    if fondo == 'kepler':
        a = (J + L)**2
        e = np.sqrt(np.clip(1 - L**2/a, 0, None))
        eta = Q.copy()
        for _ in range(14):
            eta = eta - (eta - e*np.sin(eta) - Q)/(1 - e*np.cos(eta))
        r = a*(1 - e*np.cos(eta))
        return r, np.sqrt(a)*e*np.sin(eta)/r
    E = -0.5/(J + 0.5*(L + np.sqrt(L*L + 4)))**2
    c = -0.5/E - 1.0
    e = np.sqrt(np.clip(1 - (L*L/c)*(1 + 1/c), 0, None))
    k = e*c/(c + 1)
    eta = Q.copy()
    for _ in range(14):
        eta = eta - (eta - k*np.sin(eta) - Q)/(1 - k*np.cos(eta))
    s = 2 + c*(1 - e*np.cos(eta))
    r = np.sqrt(s*(s - 2))
    p = np.sqrt(np.clip(2*(E + 1/(1 + np.sqrt(1 + r*r)) - 0.5*L*L/(r*r)), 0, None))*np.where(np.sin(eta) >= 0, 1.0, -1.0)
    return r, p


def ejes(n, x):
    sj = B.parametro(n, 'sr', 0.1)
    kappa = (0.5*(2.0 + np.sqrt(8.0)) if pot(n) == 'isocrono' else 2.0)**-3
    return np.sqrt(2*x*sj/kappa), np.sqrt(2*x*sj*kappa)


def elipse(r, p, w, ar, ap):
    """El estimador de bulto.nucleo: máscara de la elipse alrededor del máximo de densidad y su centro."""
    br, bp = np.arange(2.0, 20.0 + ar, 0.5*ar), np.arange(-0.5, 0.5 + ap, 0.5*ap)
    h = np.histogram2d(r, p, bins=(br, bp), weights=w)[0]
    k = h.copy()
    k[1:] += h[:-1]; k[:-1] += h[1:]
    k2 = k.copy()
    k2[:, 1:] += k[:, :-1]; k2[:, :-1] += k[:, 1:]
    i, j = np.unravel_index(np.argmax(k2), k2.shape)
    r0, p0 = 0.5*(br[i] + br[i + 1]), 0.5*(bp[j] + bp[j + 1])
    for _ in range(8):
        sel = ((r - r0)/ar)**2 + ((p - p0)/ap)**2 < 1
        if w[sel].sum() == 0:
            break
        r0, p0 = np.average(r[sel], weights=w[sel]), np.average(p[sel], weights=w[sel])
    return ((r - r0)/ar)**2 + ((p - p0)/ap)**2 < 1, (r0, p0)


def mc_corrida(n):
    """m_c(t) para cada J_* en las instantáneas t >= 20 tau_1 (una de cada cinco, como bulto.nucleo)."""
    t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
    w = np.asarray(mapa(n, 'w'), float); W = w.sum()
    R, P = mapa(n, 'R'), mapa(n, 'P')
    marcos = np.where(t >= 20*tau)[0][::5]
    M = np.empty((len(marcos), len(XS)))
    for a, m in enumerate(marcos):
        r, p = np.asarray(R[m], float), np.asarray(P[m], float)
        for b, x in enumerate(XS):
            M[a, b] = w[elipse(r, p, w, *ejes(n, x))[0]].sum()/W
    return M


_NULO = {}


def mc_nulo(n):
    """m_c(t) de la mezcla libre exacta de la retícula de la corrida n, en los mismos instantes. Una fase
    aleatoria por fila no sirve de referencia: las filas vecinas tienen fases correlacionadas (la espiral), y
    al descorrelacionarlas el ruido de disparo sube m_c de 0.051 a 0.073 en la retícula de B0."""
    J, Q, F = reticula(n)
    t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
    marcos = np.where(t >= 20*tau)[0][::5]
    clave = (len(J), len(Q), J[0], J[-1], Q[0], B.parametro(n, 'sp', 0.1), pot(n), len(marcos), t[-1])
    if clave in _NULO:
        return _NULO[clave]
    c = 0.5*(2.0 + np.sqrt(8.0)) if pot(n) == 'isocrono' else 2.0
    Om = (J + c)**-3
    JJ = np.repeat(J, len(Q)); QQ = np.tile(Q, len(J)); w = F.ravel(); W = w.sum()
    M = np.empty((len(marcos), len(XS)))
    for a, m in enumerate(marcos):
        r, p = fase(QQ + np.repeat(Om*t[m], len(Q)), JJ, pot(n))
        for b, x in enumerate(XS):
            M[a, b] = w[elipse(r, p, w, *ejes(n, x))[0]].sum()/W
    _NULO[clave] = M
    return M


def nulo(solo=None):
    os.makedirs(SAL, exist_ok=True)
    sal = {}
    print('m_c de la corrida (media en t >= 20 tau_1) frente a la mezcla libre exacta de la misma retícula y los mismos instantes.')
    print('nulo: media +- desviación y percentil 99 entre instantáneas. z = (m_c - media)/desviación.')
    print(f'{"corrida":10} {"M":>6} {"mu":>5} {"N":>6} {"Nef":>6} | ' + ' | '.join(
        f'x={x:<5g} m_c  nulo         p99     z ' for x in XS))
    for n in corridas():
        if (solo and n not in solo) or not os.path.exists(B.ruta(n, 'bulto.npz')):
            continue
        w = np.asarray(mapa(n, 'w'), float)
        F = reticula(n)[2]
        if w.size != F.size:
            print(f'  {n}: la retícula no coincide ({w.size} partículas frente a {F.size})', flush=True)
            continue
        if not np.allclose(np.sort(w)/w.max(), np.sort(F.ravel())/F.max(), atol=1e-5, rtol=1e-3):
            print(f'  {n}: los pesos no coinciden con los de la retícula', flush=True)
        Mc, Mn = mc_corrida(n), mc_nulo(n)
        sj, sq = B.parametro(n, 'sr', 0.1), B.parametro(n, 'sp', 0.1)
        m0 = B.info[n]['a0'] if B.info[n]['sg'] else 0.0
        mu = C.coef_pasivo(None, sj, sq, pot(n))[1]*np.sqrt(m0)
        trozos = []
        for b, x in enumerate(XS):
            mc, med, de, q = Mc[:, b].mean(), Mn[:, b].mean(), Mn[:, b].std(), np.percentile(Mn[:, b], 99)
            trozos.append(f'{mc:6.3f} {med:6.3f}+-{de:5.3f} {q:6.3f} {(mc - med)/de:5.1f}')
        print(f'{n:10} {m0:6.3g} {mu:5.2f} {w.size:6d} {w.sum()**2/(w**2).sum():6.0f} | ' + ' | '.join(trozos), flush=True)
        sal[n] = dict(M=m0, mu=mu, Mc=Mc, Mn=Mn)
    if not solo:
        np.savez(os.path.join(SAL, 'nulo.npz'), corridas=np.array(list(sal)), xs=np.array(XS),
                 M=np.array([v['M'] for v in sal.values()]), mu=np.array([v['mu'] for v in sal.values()]),
                 mc=np.array([v['Mc'].mean(0) for v in sal.values()]), mc_de=np.array([v['Mc'].std(0) for v in sal.values()]),
                 nulo=np.array([v['Mn'].mean(0) for v in sal.values()]), nulo_de=np.array([v['Mn'].std(0) for v in sal.values()]),
                 nulo_p99=np.array([np.percentile(v['Mn'], 99, axis=0) for v in sal.values()]))


XPERFIL = (0.0125, 0.025, 0.05, 0.1, 0.2)
CASOS = ('B0', 'B2q', 'B015', 'B02', 'B02q', 'J20_0.02', 'J20_0.04', 'J20_0.08', 'Q2_0.01', 'Q2_0.02', 'Q4_0.02', 'Q4_0.04',
         'K_0.02', 'J05_0.005')


def perfil(solo=None):
    """Densidad media en J < x sigma_J, con las acciones del potencial del instante, en unidades de Fbar_max, la
    mayor densidad de un gas mezclado a acción fija. Media en t >= 20 tau_1 (una instantánea de cada cinco), con su
    desviación entre instantáneas y el número efectivo de partículas dentro del intervalo."""
    print('Densidad en J < x sigma_J en unidades de Fbar_max: media +- desviación entre instantáneas [N efectivo].')
    print(f'{"corrida":10} {"M":>6} ' + ' '.join(f'{"x=" + str(x):>21}' for x in XPERFIL))
    for n in (solo or CASOS):
        t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
        w = np.asarray(mapa(n, 'w'), float); W = w.sum()
        J, sj = mapa(n, 'J'), B.parametro(n, 'sr', 0.1)
        marcos = np.where(t >= 20*tau)[0][::5]
        D = np.empty((len(marcos), len(XPERFIL))); N = D.copy()
        for a, m in enumerate(marcos):
            j = np.asarray(J[m], float)
            for b, x in enumerate(XPERFIL):
                s = j < x*sj
                D[a, b] = w[s].sum()/W/(x*MMIX)
                N[a, b] = w[s].sum()**2/max((w[s]**2).sum(), 1e-300)
        m0 = B.info[n]['a0'] if B.info[n]['sg'] else 0.0
        print(f'{n:10} {m0:6.3g} ' + ' '.join(f'{D[:, b].mean():6.2f}+-{D[:, b].std():4.2f} [{N[:, b].mean():5.1f}]'
                                              for b in range(len(XPERFIL))), flush=True)


SERIE = ('B011', 'B012', 'B013', 'B014', 'B015', 'B015q', 'B016', 'B017', 'B018', 'B019', 'B02q', 'B021', 'B022', 'B023', 'B024', 'B025h', 'B3', 'B04', 'B05', 'B07', 'B4q', 'B5')


def serie(n):
    """Series temporales de la corrida n con la elipse de J_* = 0.05 sigma_J en cada instantánea: centro (r_*, p_*),
    masa m_c(t), contribución a h_1 del gas dentro y fuera de la elipse (en unidades de la masa total), y radio medio
    de las partículas elegidas en la última instantánea, en 20 tau_1 y en 10 tau_1. Se guarda en criterio/<n>.npz."""
    sal = os.path.join(SAL, f'{n}.npz')
    if os.path.exists(sal):
        return np.load(sal)
    t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
    w = np.asarray(mapa(n, 'w'), float); W = w.sum()
    R, P, Q = mapa(n, 'R'), mapa(n, 'P'), mapa(n, 'Q')
    ar, ap = ejes(n, 0.05)
    elegir = lambda m: elipse(np.asarray(R[m], float), np.asarray(P[m], float), w, ar, ap)[0]
    i10, i20 = int(np.argmin(np.abs(t - 10*tau))), int(np.argmin(np.abs(t - 20*tau)))
    sf, s20, s10 = elegir(len(t) - 1), elegir(i20), elegir(i10)
    c = np.empty((len(t), 2)); mc = np.empty(len(t)); hn = np.empty(len(t), complex); hh = hn.copy()
    rc = np.empty((len(t), 3))
    for m in range(len(t)):
        r, p = np.asarray(R[m], float), np.asarray(P[m], float)
        s, c[m] = elipse(r, p, w, ar, ap)
        mc[m] = w[s].sum()/W
        e = w*np.exp(-1j*np.nan_to_num(np.asarray(Q[m], float)))
        hn[m], hh[m] = e[s].sum()/W, e[~s].sum()/W
        rc[m] = [np.average(r[x], weights=w[x]) for x in (sf, s20, s10)]
    sol = lambda a, b: w[a & b].sum()/w[a].sum()
    np.savez(sal, t=t, tau=tau, centro=c, mc=mc, hn=hn, hh=hh, rc=rc, sol20=sol(sf, s20), sol10=sol(sf, s10))
    return np.load(sal)


def env_ac(x, t, tau):
    """Envolvente de la desviación de x respecto a su media en una vuelta radial, su primer mínimo local y su rango
    en t >= 25 tau_1 (la definición de a_c de bulto.medidas_articulo)."""
    k = int(round(45/(t[1] - t[0])))
    media = np.array([x[max(0, i - k):i + k + 1].mean() for i in range(len(t))])
    env = B.envolvente(np.abs(x - media), t)
    i = None
    for a in np.arange(2, 30, 0.5):
        z = (t >= a*tau) & (t < (a + 2)*tau)
        j = int(np.argmin(np.where(z, env, np.inf)))
        if env[j] <= env[(t >= t[j] - 2*tau) & (t <= t[j] + 2*tau)].min():
            i = j
            break
    if i is None:
        i = int(np.argmin(np.where((t >= 2*tau) & (t <= 30*tau), env, np.inf)))
    return env, env[i], t[i]/tau, env[t >= 25*tau].min(), env[t >= 25*tau].max()


def oscilacion(solo=None):
    os.makedirs(SAL, exist_ok=True)
    print('a_c: primer mínimo (valor, t/tau_1) y rango en t >= 25 tau_1, según qué partículas siguen al núcleo.')
    print('sol: fracción de la masa elegida al final que también estaba en la elipse en 20 y 10 tau_1.')
    print(f'{"corrida":7} {"M":>6} | {"final":^28} | {"elegidas en 20 tau_1":^28} | {"elegidas en 10 tau_1":^28} | '
          f'{"centro de la elipse":^28} | sol20 sol10')
    for n in (solo or SERIE):
        d = serie(n); t, tau = d['t'], float(d['tau'])
        tr = []
        for x in (d['rc'][:, 0], d['rc'][:, 1], d['rc'][:, 2], d['centro'][:, 0]):
            _, a, b, lo, hi = env_ac(x, t, tau)
            tr.append(f'{a:6.3f} {b:5.1f}  {lo:5.3f}-{hi:5.3f}    ')
        print(f'{n:7} {B.info[n]["a0"]:6.3g} | ' + ' | '.join(tr) + f' | {float(d["sol20"]):5.2f} {float(d["sol10"]):5.2f}', flush=True)


def halo(solo=None):
    os.makedirs(SAL, exist_ok=True)
    print('Medias en 10 <= t/tau_1 <= 20. h_1 = h_nucleo + h_halo, en unidades de la masa total; |h_1| guardado para comprobar.')
    print(f'{"corrida":7} {"M":>6} {"m_c":>6} {"|h1|":>9} {"guardado":>9} {"|h_nuc|":>9} {"|h_halo|":>9} {"|h_nuc|/m_c":>11} {"a_c/0.38":>9}')
    for n in (solo or SERIE):
        d = serie(n); t, tau = d['t'], float(d['tau'])
        v = (t >= 10*tau) & (t <= 20*tau)
        h = d['hn'] + d['hh']
        hk = np.asarray(mapa(n, 'hk')[:, 1])
        env = env_ac(d['rc'][:, 0], t, tau)[0]
        mc = d['mc'][v].mean()
        print(f'{n:7} {B.info[n]["a0"]:6.3g} {mc:6.3f} {np.abs(h[v]).mean():9.2e} {np.abs(hk[v]).mean():9.2e} '
              f'{np.abs(d["hn"][v]).mean():9.2e} {np.abs(d["hh"][v]).mean():9.2e} {np.abs(d["hn"][v]).mean()/mc:11.3f} '
              f'{env[v].mean()/0.38:9.3f}', flush=True)


def central(solo=None):
    """Para todas las corridas: densidad en J < 0.025 y 0.05 sigma_J en unidades de Fbar_max (media en t >= 20 tau_1) y
    fracción de la masa que está en ese intervalo al final y ya lo estaba en 20 tau_1. Después, el ajuste lineal de
    m_c menos su valor en la mezcla libre, frente a la masa, en las cinco corridas de 400 x 100 entre 0.02 y 0.024."""
    fil = {}
    print(f'{"corrida":10} {"M":>6} {"mu":>5} | densidad x=0.025  x=0.05 | sigue en J<x desde 20 tau_1: x=0.025  x=0.05 | masa J<0.05 s_J')
    nu = np.load(os.path.join(SAL, 'nulo.npz'))
    for n, m0, mu in zip(nu['corridas'], nu['M'], nu['mu']):
        n = str(n)
        if solo and n not in solo:
            continue
        t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
        w = np.asarray(mapa(n, 'w'), float); W = w.sum()
        J, sj = mapa(n, 'J'), B.parametro(n, 'sr', 0.1)
        marcos = np.where(t >= 20*tau)[0][::5]
        D = np.mean([[w[np.asarray(J[m]) < x*sj].sum()/W/(x*MMIX) for x in (0.025, 0.05)] for m in marcos], axis=0)
        jf, j20 = np.asarray(J[-1], float), np.asarray(J[int(np.argmin(np.abs(t - 20*tau)))], float)
        S = [w[(jf < x*sj) & (j20 < x*sj)].sum()/max(w[jf < x*sj].sum(), 1e-300) for x in (0.025, 0.05)]
        mas = w[jf < 0.05*sj].sum()/W
        fil[n] = (m0, mu, D[0], D[1], S[0], S[1], mas)
        print(f'{n:10} {m0:6.3g} {mu:5.2f} | {D[0]:15.2f} {D[1]:7.2f} | {S[0]:30.2f} {S[1]:7.2f} | {mas:6.3f}', flush=True)
    if not solo:
        np.savez(os.path.join(SAL, 'central.npz'), corridas=np.array(list(fil)), datos=np.array(list(fil.values())))
    k = [list(map(str, nu['corridas'])).index(x) for x in ('B02q', 'B021', 'B022', 'B023', 'B024')]
    M, ex = nu['M'][k], nu['mc'][k, 1] - nu['nulo'][k, 1]
    p, cov = np.polyfit(M, ex, 1, cov=True)
    M0 = -p[1]/p[0]
    dM = abs(M0)*np.sqrt(cov[0, 0]/p[0]**2 + cov[1, 1]/p[1]**2 - 2*cov[0, 1]/(p[0]*p[1]))
    print(f'\nm_c - mezcla libre = {p[0]:.2f} (M - M0), M0 = {M0:.4f} +- {dM:.4f} (ajuste de 5 corridas); residuo rms '
          f'{np.sqrt(np.mean((np.polyval(p, M) - ex)**2)):.4f}')


def centrado(solo=None):
    """Si la elipse de J_* = 0.05 sigma_J está sobre la órbita circular o sobre un trozo de la espiral. f_J: fracción de
    su masa con J < 0.1 sigma_J (acciones del potencial del instante), media y mínimo en t >= 20 tau_1. sol20 y sol10:
    fracción de la masa de la elipse final que ya estaba en la elipse en 20 y 10 tau_1."""
    os.makedirs(SAL, exist_ok=True)
    fil = {}
    print(f'{"corrida":10} {"M":>6} {"mu":>5} {"columnas":>8} | f_J media  mínimo | sol20 sol10')
    for n in corridas():
        if (solo and n not in solo) or not os.path.exists(B.ruta(n, 'bulto.npz')):
            continue
        t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
        w = np.asarray(mapa(n, 'w'), float)
        R, P, J = mapa(n, 'R'), mapa(n, 'P'), mapa(n, 'J')
        sj, sq = B.parametro(n, 'sr', 0.1), B.parametro(n, 'sp', 0.1)
        ar, ap = ejes(n, 0.05)
        elegir = lambda m: elipse(np.asarray(R[m], float), np.asarray(P[m], float), w, ar, ap)[0]
        marcos = np.where(t >= 20*tau)[0][::5]
        f = np.array([w[s & (np.asarray(J[m], float) < 0.1*sj)].sum()/max(w[s].sum(), 1e-300)
                      for m, s in ((m, elegir(m)) for m in marcos)])
        i10, i20 = int(np.argmin(np.abs(t - 10*tau))), int(np.argmin(np.abs(t - 20*tau)))
        sf, s20, s10 = elegir(len(t) - 1), elegir(i20), elegir(i10)
        so20, so10 = w[sf & s20].sum()/w[sf].sum(), w[sf & s10].sum()/w[sf].sum()
        m0 = B.info[n]['a0'] if B.info[n]['sg'] else 0.0
        mu = C.coef_pasivo(None, sj, sq, pot(n))[1]*np.sqrt(m0)
        nq = len(reticula(n)[1])
        fil[n] = (m0, mu, nq, f.mean(), f.min(), so20, so10)
        print(f'{n:10} {m0:6.3g} {mu:5.2f} {nq:8d} | {f.mean():9.2f} {f.min():7.2f} | {so20:5.2f} {so10:5.2f}', flush=True)
    if not solo:
        np.savez(os.path.join(SAL, 'centrado.npz'), corridas=np.array(list(fil)), datos=np.array(list(fil.values())))


BR_CAPA = np.arange(0.5, 60.0001, 0.02)        # radios para dM/dr
_CAPA = {}


def _max_dMdr(r, w, K):
    """Mayor valor de dM/dr (por unidad de masa) promediado en una ventana de K celdas, y el radio donde ocurre."""
    c = np.concatenate([[0.0], np.cumsum(np.histogram(r, BR_CAPA, weights=w)[0])])/w.sum()
    d = (c[K:] - c[:-K])/(K*0.02)
    i = int(np.argmax(d))
    return d[i], 0.5*(BR_CAPA[i] + BR_CAPA[i + K])


def _capa_libre(n, t, marcos, K):
    """El mismo máximo para la mezcla libre exacta de la retícula de n en los mismos instantes."""
    J, Q, F = reticula(n)
    clave = (len(J), len(Q), J[0], J[-1], Q[0], B.parametro(n, 'sp', 0.1), pot(n), len(marcos), t[-1], K)
    if clave not in _CAPA:
        c = 0.5*(2.0 + np.sqrt(8.0)) if pot(n) == 'isocrono' else 2.0
        Om = (J + c)**-3
        JJ = np.repeat(J, len(Q)); QQ = np.tile(Q, len(J)); w = F.ravel()
        _CAPA[clave] = np.mean([_max_dMdr(fase(QQ + np.repeat(Om*t[m], len(Q)), JJ, pot(n))[0], w, K) for m in marcos], axis=0)
    return _CAPA[clave]


def capa(solo=None):
    """La cáscara en el espacio físico: mayor valor instantáneo de dM/dr, promediado en una ventana de ancho a = (2 J_*/kappa)^(1/2)
    (el semieje radial de la elipse del núcleo) y luego en t >= 20 tau_1, frente al de la mezcla libre exacta de la misma retícula.
    Se guarda en criterio/capa.npz: M, Gamma, F_c/Fhat_max, dM/dr de la corrida y de la mezcla libre, su cociente, radio del máximo
    y radio de la órbita circular."""
    os.makedirs(SAL, exist_ok=True)
    fil = {}
    print(f'{"corrida":10} {"M":>6} {"Gamma":>6} {"F_c":>6} | dM/dr: {"PIC":>6} {"libre":>6} {"cociente":>8} | r del máximo: {"PIC":>5} {"libre":>5} {"R_0":>5}')
    for n in corridas():
        if (solo and n not in solo) or not os.path.exists(B.ruta(n, 'bulto.npz')):
            continue
        t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
        w = np.asarray(mapa(n, 'w'), float); R = mapa(n, 'R')
        marcos = np.where(t >= 20*tau)[0][::5]
        K = max(1, int(round(ejes(n, 0.05)[0]/0.02)))
        a = np.mean([_max_dMdr(np.asarray(R[m], float), w, K) for m in marcos], axis=0)
        b = _capa_libre(n, t, marcos, K)
        R0 = 5.742 if pot(n) == 'isocrono' else 4.0
        m0 = B.info[n]['a0'] if B.info[n]['sg'] else 0.0
        fil[n] = (m0, B.gamma(n), B.densidad_central(n), a[0], b[0], a[0]/b[0], a[1], b[1], R0)
        print(f'{n:10} {m0:6.3g} {B.gamma(n):6.3f} {B.densidad_central(n):6.2f} |        {a[0]:6.3f} {b[0]:6.3f} {a[0]/b[0]:8.2f} |               {a[1]:5.2f} {b[1]:5.2f} {R0:5.2f}', flush=True)
    if not solo:
        np.savez(os.path.join(SAL, 'capa.npz'), corridas=np.array(list(fil)), datos=np.array(list(fil.values())))


DESTINO = os.path.join(B.RAIZ, 'docs', 'articulo', 'core_checks')
ESTILOS = [('o', 'k'), ('v', '#0072B2'), ('^', '#009E73'), ('s', '#E69F00'), ('D', '#CC79A7'), ('*', '#D55E00')]


def clasificar():
    """Clase de cada corrida. grupo: al menos 40 % de la masa de la elipse final ya estaba en la elipse en 20 tau_1 (entre
    las 48 corridas no hay valores entre 0.12 y 0.47; en la mezcla libre el valor es 0.19). sin_resolver: sin lo anterior,
    pero con densidad mayor que Fbar_max en J < 0.025 sigma_J. sin_grupo: ninguna de las dos. Esta clase dice si la elipse
    se queda sobre un mismo gas, no si existe gas junto a J = 0: la densidad central crece sin salto desde 0.01."""
    nu, ce, cn = (np.load(os.path.join(SAL, f)) for f in ('nulo.npz', 'central.npz', 'centrado.npz'))
    nom = [str(x) for x in nu['corridas']]
    assert nom == [str(x) for x in ce['corridas']] == [str(x) for x in cn['corridas']]
    clase = {}
    for i, n in enumerate(nom):
        clase[n] = 'grupo' if cn['datos'][i, 5] >= 0.4 else 'sin_resolver' if ce['datos'][i, 2] > 1 else 'sin_grupo'
    return nom, clase, nu, ce


def figuras(solo=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Circle
    os.makedirs(DESTINO, exist_ok=True)
    nom, clase, nu, ce = clasificar()
    print({c: sum(v == c for v in clase.values()) for c in ('sin_grupo', 'sin_resolver', 'grupo')})
    print('sin resolver:', [n for n in nom if clase[n] == 'sin_resolver'])
    i = {n: k for k, n in enumerate(nom)}

    # 1. Los dos indicadores frente a mu, por familia. Relleno: núcleo; vacío: sin núcleo.
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3), layout='constrained')
    for (nombre, base, conv), (mk, col) in zip(B.FAMILIAS, ESTILOS):
        rs = [n for n in dict.fromkeys(base + conv + (NUEVAS if nombre == B.FAMILIAS[0][0] else [])) if n in i and B.info[n]['sg']]
        for a, y in zip(ax, (lambda n: nu['mc'][i[n], 1]/nu['nulo'][i[n], 1], lambda n: max(ce['datos'][i[n], 2], 5e-3))):
            for n in rs:
                a.plot(nu['mu'][i[n]], y(n), mk, ms=7 if mk == '*' else 5, mec=col, mew=0.9,
                       mfc=col if clase[n] == 'grupo' else 'none', ls='none')
            a.plot([], [], mk, color=col, ms=5, ls='none', label=nombre)
    ax[0].plot([], [], 'o', mec='0.3', mfc='0.3', ms=5, ls='none', label='persistent ellipse')
    ax[0].plot([], [], 'o', mec='0.3', mfc='none', ms=5, ls='none', label='non-persistent ellipse')
    ax[0].axhline(1, color='k', lw=0.7, ls=':')
    ax[1].axhline(1, color='k', lw=0.7, ls=':')
    ax[0].set_yscale('log'); ax[1].set_yscale('log')
    ax[0].set_ylabel('$m_c$ / value for exact free mixing'); ax[1].set_ylabel(r'density in $J<0.025\,\sigma_J$ / $\hat{F}_{max}$')
    ax[1].set_ylim(4e-3, 30)
    for a, l in zip(ax, 'ab'):
        a.set_xlabel(r'$\mu$ (passive)'); a.set_xlim(0.22, 1.03); a.grid(alpha=0.25); a.set_title(f'({l})', loc='left', fontsize=9)
    ax[0].set_ylim(0.3, 13)
    ax[0].legend(fontsize=6.2, frameon=False, loc='lower right', ncol=2)
    fig.savefig(os.path.join(DESTINO, 'fig1_core_indicators.png'), dpi=180); plt.close(fig)

    # 2. Densidad en el plano de fases junto a la órbita circular, en unidades de Fbar_max.
    R0, kap = 5.742, 0.07107
    casos = (('B2q', r'$M_{gas}/M=0.01$'), ('B015q', r'$M_{gas}/M=0.015$'), ('B02q', r'$M_{gas}/M=0.02$'),
             ('J20_0.04', r'$\sigma_J=0.2$, $M_{gas}/M=0.04$'))
    fig, ax = plt.subplots(1, 4, figsize=(7.6, 2.35), layout='constrained', sharey=True)
    for a, (n, tit) in zip(ax, casos):
        w = np.asarray(mapa(n, 'w'), float); sj = B.parametro(n, 'sr', 0.1)
        t, tau = np.asarray(mapa(n, 't')), B.banda(n)[3]
        marcos = np.where(t >= 20*tau)[0][::5]               # promedio temporal: una instantánea es ruido de partículas
        RR, PP = mapa(n, 'R'), mapa(n, 'P')
        b = np.linspace(-1.1, 1.1, 45)
        h = np.zeros((len(b) - 1, len(b) - 1))
        for m in marcos:
            r, p = np.asarray(RR[m], float), np.asarray(PP[m], float)
            h += np.histogram2d((r - R0)*np.sqrt(kap/(2*sj)), p/np.sqrt(2*sj*kap), bins=(b, b), weights=w)[0]
        h /= len(marcos)*w.sum()
        n0 = np.exp(-1.0)/(sj*np.sqrt(np.pi)/4)
        d = h/((b[1] - b[0])**2*2*sj)*2*np.pi/n0
        im = a.pcolormesh(b, b, np.where(d > 0, d, np.nan).T, cmap='Blues', norm=LogNorm(0.1, 30), rasterized=True)
        a.add_patch(Circle((0, 0), np.sqrt(0.05), fill=False, ec='#D55E00', lw=0.9))
        a.set_title(tit, fontsize=7.5); a.set_aspect('equal'); a.set_xlabel(r'$(r-R_0)\sqrt{\kappa_0/2\sigma_J}$', fontsize=7)
        a.tick_params(labelsize=6.5)
    ax[0].set_ylabel(r'$p_r/\sqrt{2\sigma_J\kappa_0}$', fontsize=7)
    c = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.01); c.set_label(r'DF / $\hat{F}_{max}$', fontsize=7); c.ax.tick_params(labelsize=6.5)
    fig.savefig(os.path.join(DESTINO, 'fig2_phase_space.png'), dpi=200); plt.close(fig)

    # 3. |h_1| dentro y fuera del núcleo, y sensibilidad de a_c a la selección de partículas.
    rs = [n for n in SERIE if n not in ('B015', 'B5')]
    mu = np.array([nu['mu'][i[n]] for n in rs]); hT, hN, hH = [], [], []
    for n in rs:
        d = serie(n); v = (d['t'] >= 10*float(d['tau'])) & (d['t'] <= 20*float(d['tau']))
        hT.append(np.abs(d['hn'][v] + d['hh'][v]).mean()); hN.append(np.abs(d['hn'][v]).mean()); hH.append(np.abs(d['hh'][v]).mean())
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.8), layout='constrained')
    ok = np.array([n not in ('B011', 'B012', 'B013', 'B014') for n in rs])    # ahí la elipse no está sobre el pico central
    for y, l, col, mk, q in ((hT, 'all the gas', 'k', 'o', np.ones_like(ok)), (hN, 'inside the ellipse', '#0072B2', 's', ok),
                             (hH, 'outside the ellipse', '#E69F00', '^', ok)):
        ax[0].semilogy(mu[q], np.array(y)[q], mk + '-', color=col, ms=3.5, lw=1, label=l)
    ax[0].set_ylabel(r'mean $|h_1|/h_0$, $10\leq t/\tau_1\leq 20$'); ax[0].legend(fontsize=6.5, frameon=False)
    rc = [n for n in ('B015q', 'B016', 'B017', 'B018', 'B019', 'B02q', 'B021', 'B022', 'B023', 'B024', 'B025h')]
    mc = [nu['mu'][i[n]] for n in rc]
    et = ('selected at the end', r'selected at $20\tau_1$', r'selected at $10\tau_1$', 'center of the ellipse')
    for k, (l, col, mk) in enumerate(zip(et, ('k', '#0072B2', '#009E73', '#D55E00'), 'osv^')):
        va, ti = [], []
        for n in rc:
            d = serie(n); x = d['rc'][:, k] if k < 3 else d['centro'][:, 0]
            _, a_, b_, _, _ = env_ac(x, d['t'], float(d['tau']))
            va.append(a_); ti.append(b_)
        ax[1].semilogy(mc, va, mk + '-', color=col, ms=3.5, lw=1, label=l)
        ax[2].plot(mc, ti, mk + '-', color=col, ms=3.5, lw=1)
    ax[1].set_ylabel('first minimum of $a_c$'); ax[2].set_ylabel(r'time of that minimum, $t/\tau_1$')
    ax[1].legend(fontsize=5.6, frameon=False, loc='upper center')
    ax[1].set_ylim(top=0.5)
    ax[2].set_ylim(0, 33); ax[2].axhline(30, color='0.5', lw=0.7, ls=':')
    for a, l in zip(ax, 'abc'):
        a.set_xlabel(r'$\mu$ (passive)'); a.grid(alpha=0.25); a.set_title(f'({l})', loc='left', fontsize=9); a.tick_params(labelsize=7)
    fig.savefig(os.path.join(DESTINO, 'fig3_h1_and_ac.png'), dpi=180); plt.close(fig)

    # 4. Lo que cambia con la masa en la retícula 400 x 100, de M_gas/M = 0.01 a 0.024.
    se = ['B2q', 'B011', 'B012', 'B013', 'B014', 'B015q', 'B016', 'B017', 'B018', 'B019', 'B02q', 'B021', 'B022', 'B023', 'B024']
    k = [i[n] for n in se]
    M = nu['M'][k]
    cn = np.load(os.path.join(SAL, 'centrado.npz'))
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.7), layout='constrained')
    ax[0].plot(M, nu['mc'][k, 1], 'o-', color='k', ms=3.5, lw=1, label='PIC')
    ax[0].plot(M, nu['nulo'][k, 1], '-', color='#0072B2', lw=1, label='free mixing')
    ax[0].fill_between(M, nu['nulo'][k, 1] - nu['nulo_de'][k, 1], nu['nulo_p99'][k, 1], color='#0072B2', alpha=0.15, lw=0)
    ax[0].set_ylabel('$m_c$'); ax[0].legend(fontsize=6.5, frameon=False)
    ax[1].semilogy(M, ce['datos'][k, 2], 'o-', color='k', ms=3.5, lw=1)
    ax[1].axhline(1, color='k', lw=0.7, ls=':')
    ax[1].set_ylabel(r'density in $J<0.025\,\sigma_J$ / $\hat{F}_{max}$')
    ax[2].plot(M, cn['datos'][k, 5], 'o-', color='k', ms=3.5, lw=1, label='ellipse at the densest point')
    ax[2].plot(M, ce['datos'][k, 5], 's--', color='#009E73', ms=3, lw=1, label=r'fixed region $J<0.05\,\sigma_J$')
    ax[2].set_ylim(-0.05, 1.3); ax[2].set_ylabel(r'fraction already there at $20\tau_1$'); ax[2].legend(fontsize=6, frameon=False, loc='upper left')
    ax[2].legend(fontsize=6, frameon=False, loc='center right')
    for a, l in zip(ax, 'abc'):
        a.set_xlabel(r'$M_{gas}/M$'); a.grid(alpha=0.25); a.set_title(f'({l})', loc='left', fontsize=9); a.tick_params(labelsize=7)
    fig.savefig(os.path.join(DESTINO, 'fig4_onset.png'), dpi=180); plt.close(fig)
    print('figuras en', DESTINO)
    print('Cociente máx/mín entre las cuatro selecciones (primer mínimo, rango tardío) y tiempos del mínimo:')
    for n in rc:
        d = serie(n)
        q = [env_ac(d['rc'][:, k] if k < 3 else d['centro'][:, 0], d['t'], float(d['tau']))[1:] for k in range(4)]
        q = np.array(q)
        print(f'  {n:6} mínimo {q[:, 0].min():.3f}-{q[:, 0].max():.3f} (x{q[:, 0].max()/q[:, 0].min():.1f}); t {q[:, 1].min():.1f}-{q[:, 1].max():.1f}; '
              f'tardío inferior {q[:, 2].min():.3f}-{q[:, 2].max():.3f}, superior {q[:, 3].min():.3f}-{q[:, 3].max():.3f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['nulo', 'perfil', 'oscilacion', 'halo', 'central', 'centrado', 'capa', 'figuras'])
    ap.add_argument('--solo', nargs='*')
    a = ap.parse_args()
    {'nulo': nulo, 'perfil': perfil, 'oscilacion': oscilacion, 'halo': halo, 'central': central, 'centrado': centrado, 'capa': capa,
     'figuras': figuras}[a.paso](a.solo)
