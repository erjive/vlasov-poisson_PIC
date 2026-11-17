"""¿Es lineal el frenado de la mezcla por el campo? Respuesta lineal alrededor del estado mezclado.

Hipótesis (2026-11-17). Pasadas unas tau_1, el gas de la nube está casi mezclado: su FD promediada en el ángulo,
F_0(J), es casi un estado estacionario del potencial Phi_0, y lo que queda, dF = F - F_0, es una perturbación cuyo
campo es pequeño (D ~ 1e-3). Si el frenado de la mezcla que se mide en las corridas (la corrida decae más despacio
que la mezcla sola, colas.py) es la respuesta lineal de F_0 al campo de dF, la ecuación linealizada

    d(dF)/dt + Omega_0(J) d(dF)/dQ = F_0'(J) d(dPhi)/dQ,       dPhi = potencial de dF,

con (Q, J) las variables de ángulo-acción de Phi_0, debe reproducir la corrida desde t0. Sin el término de la
derecha es la mezcla sola en Phi_0. Se comparan |h_n(t)| (con las variables de Phi_0 también en la corrida) y la
desviación D(t) del potencial del gas respecto de Phi_0, que no depende de ninguna variable de fase.

Pasos:
    python3 respuesta_mezcla.py extraer <bulto.npz> <salida.npz> --mgas M [--t0 5] [--tau T] [--procesos 4]
        Phi_0: media del potencial del gas en |t - t0| <= 0.5 tau_1. Variables de Phi_0 de las partículas en t0, y
        h_1..h_3 y D de la corrida en las instantáneas t >= t0 (una de cada dos). Corre donde esté el bulto.npz.
    python3 respuesta_mezcla.py lineal <extraido.npz> <salida.npz> [--nj 2400] [--nq 64] [--dt 0.5] [--t1 30]
        evoluciones lineal y libre desde t0 (unos minutos cada una).
    python3 respuesta_mezcla.py comparar <salida.npz> ...   tabla de la corrida frente a lineal y libre.
    python3 respuesta_mezcla.py resumen <salida.npz> ...    una línea por corrida, y figura (exe/demo_eta/respuesta/)

La discretización es la de lineal.py: rejilla fija de NJ x NQ nodos en (Q, J), transporte exacto en Q por FFT,
golpe del campo con separación de Strang y Poisson de capas esféricas, exacto en los nodos. dF(t0) se arma con los
armónicos n = 1..NQ/2 - 1 de las partículas en cada franja de J, y F_0' con un filtro de Savitzky y Golay sobre el
histograma de J.

Corrección de la rejilla. En la rejilla cada partícula gira con la frecuencia del centro de su franja de J; el error
de fase crece como |Omega'| (t - t0) dJ y deja un piso en |h_1| del orden de |Omega'| (t - t0) dJ/sqrt(12 N_ef)
(1e-3 en 20 tau_1 con 400 x 100 partículas y NJ = 2400). Como el problema es lineal, la respuesta al campo es la
diferencia entre las evoluciones lineal y libre en la rejilla, que parten del mismo estado, y la mezcla libre se
calcula exacta con las partículas: h_n = h_n^libre exacta + (h_n^lineal - h_n^libre) en la rejilla, y lo mismo para
el potencial, con los radios de las partículas en su órbita de Phi_0 interpolados de los nodos (Fourier en Q,
lineal en J). Con la corrida B2q (Gamma = 0.114) la evolución lineal así corregida sigue a la corrida a 2-10 % hasta
20 tau_1, y la mezcla sola no (factor 1.3 a 1.85).
"""
import os, sys, argparse, time
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from aa_numerico import MapaAA
from equilibrio import invertir
from landau_libre import omega_de_J

RVENT = (3.0, 15.0)


RMED = np.linspace(RVENT[0], RVENT[1], 121)


def potencial_capas(rr, mm, rq=RMED):
    """Potencial exacto de capas esféricas de radios rr y masas mm en los radios rq (G = 1)."""
    o = np.argsort(rr); a, mm = rr[o], mm[o]
    Mc = np.concatenate([[0.0], np.cumsum(mm)]); Pc = np.concatenate([np.cumsum((mm/a)[::-1])[::-1], [0.0]])
    j = np.searchsorted(a, rq)
    return -Mc[j]/rq - Pc[j]


def _mapear(tarea):
    rg, phi0, L, ks, Rs, Ps, w, mgas = tarea
    mm = mgas*w/w.sum()
    m = MapaAA(rg, phi0, L=L)
    sal = []
    for k, R, P in zip(ks, Rs, Ps):
        R = np.asarray(R, float)
        Q, J, _ = m(R, np.asarray(P, float))
        ok = np.isfinite(Q) & np.isfinite(J)
        sal.append((k, [np.sum(w[ok]*np.exp(-1j*n*Q[ok]))/w[ok].sum() for n in (1, 2, 3)], int((~ok).sum()),
                    potencial_capas(R, mm)))
    return sal


def extraer(npz, sal, mgas, T0=5.0, tau=327.15, ventana=0.5, cada=2, procesos=4, L=2.0):
    from multiprocessing import Pool
    d = np.load(npz)                            # un .npz no se mapea en memoria: cada d[...] lee el arreglo entero
    t, rg = np.asarray(d['t']), np.asarray(d['r'])
    ps = np.asarray(d['phi_self'], float); w = np.asarray(d['w'], float)
    Rall, Pall = d['R'], d['P']
    k0 = int(np.argmin(np.abs(t - T0*tau)))
    phi0 = ps[np.abs(t - t[k0]) <= ventana*tau].mean(0)
    m = MapaAA(rg, phi0, L=L)
    Q0, J0, _ = m(np.asarray(Rall[k0], float), np.asarray(Pall[k0], float))
    malos = int((~(np.isfinite(Q0) & np.isfinite(J0))).sum())
    ks = list(range(k0, len(t), cada))
    grupos = [(rg, phi0, L, ks[i:i + 10], Rall[ks[i:i + 10]], Pall[ks[i:i + 10]], w, mgas) for i in range(0, len(ks), 10)]
    t1 = time.time()
    with Pool(procesos) as pool:
        res = [x for b in pool.imap(_mapear, grupos) for x in b]
    res.sort(key=lambda x: x[0])
    hq = np.array([x[1] for x in res]); tq = t[[x[0] for x in res]]; Pq = np.array([x[3] for x in res])
    rv = (rg >= RVENT[0]) & (rg <= RVENT[1]); esc = np.abs(phi0[rv]).max()
    D = np.sqrt(np.mean((ps[:, rv] - phi0[rv])**2, axis=1))/esc
    np.savez(sal, T0=T0, tau=tau, t0=t[k0], rg=rg, phi0=phi0, Q0=Q0, J0=J0, m=mgas*w/w.sum(), mgas=mgas, L=L,
             tq=tq, hq=hq, t=t, D=D, esc=esc, malos=malos, R0=np.asarray(Rall[k0], float), P0=np.asarray(Pall[k0], float),
             Pq=Pq)
    print(f'extraído {sal}: t0 = {t[k0]:.0f} ({t[k0]/tau:.2f} tau_1), {len(w)} partículas ({malos} sin mapa en t0), '
          f'{len(ks)} instantáneas mapeadas en {time.time() - t1:.0f} s con {procesos} procesos', flush=True)


def tabla_E_J(m, jmax, n=4000):
    """E(J) en el potencial del mapa m, sobre órbitas que pasan por la órbita circular, hasta J >= jmax."""
    rc = float(m.rc); Ec = float(m.phi_ef(np.array([rc]))[0])
    dE = 1e-3
    while m(np.array([rc]), np.array([np.sqrt(2*dE)]))[1][0] < jmax:
        dE *= 1.5
    E = Ec + np.linspace(0.0, dE, n)
    J = m(np.full(n, rc), np.sqrt(2*(E - Ec)))[1]
    J[0] = 0.0
    return E, J


def evolucion(e, nj=2400, nq=64, dt=0.5, T1=30.0, cada=20.0, nsg=11, verboso=True, mascara=None):
    """Evoluciones lineal y libre desde t0 del estado extraído e. Devuelve un diccionario con t, h_1 y h_2 (pesos
    de la masa) y D de cada una, y los controles: autoconsistencia de Phi_0 y h_1(t0) de las partículas y de la
    rejilla."""
    from scipy.signal import savgol_filter
    L, mg, tau, t0 = float(e['L']), float(e['mgas']), float(e['tau']), float(e['t0'])
    m = MapaAA(e['rg'], e['phi0'], L=L)
    ok = np.isfinite(e['Q0']) & np.isfinite(e['J0'])
    Q0, J0, mi = e['Q0'][ok], e['J0'][ok], e['m'][ok]
    jmax = 1.02*J0.max(); dJ = jmax/nj
    Jn = (np.arange(nj) + 0.5)*dJ; Qn = (np.arange(nq) + 0.5)*2*np.pi/nq; dQ = 2*np.pi/nq
    E_t, J_t = tabla_E_J(m, jmax)
    om = omega_de_J(E_t, J_t)(Jn)
    b = np.clip((J0/dJ).astype(int), 0, nj - 1)
    c = 8*np.pi**2*L*2*np.pi*dJ                    # F_n(J_b) = sum_{i en b} m_i e^{-i n Q_i}/c
    F0 = np.bincount(b, weights=mi, minlength=nj)/c
    dF0 = savgol_filter(F0, nsg, 3, deriv=1, delta=dJ)
    if mascara is not None:                     # el campo actúa solo sobre el gas con jlo <= J < jhi
        dF0 = dF0*((Jn >= mascara[0]) & (Jn < mascara[1]))
    N = nq//2 - 1
    dF = np.zeros((nj, nq))
    for n in range(1, N + 1):
        Fn = (np.bincount(b, weights=mi*np.cos(n*Q0), minlength=nj) - 1j*np.bincount(b, weights=mi*np.sin(n*Q0), minlength=nj))/c
        dF += 2*np.real(Fn[:, None]*np.exp(1j*n*Qn)[None, :])
    JJ, QQ = np.meshgrid(Jn, Qn, indexing='ij')
    t1 = time.time(); r = np.empty(JJ.size)
    for s in range(0, JJ.size, 4000):
        r[s:s + 4000] = invertir(m, E_t, J_t, QQ.ravel()[s:s + 4000], JJ.ravel()[s:s + 4000])[0]
    if verboso:
        print(f'  radios de {JJ.size} nodos en {time.time() - t1:.0f} s', flush=True)
    orden = np.argsort(r); rs = r[orden]
    peso = 8*np.pi**2*L*dQ*dJ
    k = np.fft.fftfreq(nq, 1.0/nq)
    fase = np.exp(-1j*np.outer(om, k)*dt)
    rmed = np.linspace(RVENT[0], RVENT[1], 121)
    rv = (e['rg'] >= RVENT[0]) & (e['rg'] <= RVENT[1]); esc = float(e['esc'])

    def dPhi(f):
        mS = (f.ravel()*peso)[orden]
        phi = np.empty_like(r)
        phi[orden] = -(np.cumsum(mS) - 0.5*mS)/rs - (np.cumsum((mS/rs)[::-1])[::-1] - 0.5*mS/rs)
        return phi

    # Autoconsistencia: el potencial de F_0 (repartido igual en Q) frente a Phi_0, en la ventana radial.
    p0 = dPhi(np.repeat(F0[:, None], nq, axis=1))
    o = np.argsort(r)
    phiF0 = np.interp(rmed, r[o], p0[o]); phi0m = np.interp(rmed, e['rg'], e['phi0'])
    autoc = np.sqrt(np.mean((phiF0 - phi0m)**2))/esc
    eQ = np.exp(-1j*Qn)
    h1_part = np.sum(mi*np.exp(-1j*Q0))/mg
    h1_rej = np.sum(dF*eQ[None, :])*peso/mg

    def correr(libre):
        f = dF.copy()

        def golpe(f, h):
            if libre:
                return f
            P = dPhi(f).reshape(nj, nq)
            dP = np.real(np.fft.ifft(1j*k[None, :]*np.fft.fft(P, axis=1), axis=1))
            return f + h*dF0[:, None]*dP
        pasos = int(round((T1*tau - t0)/dt)); nsal = int(round(cada/dt))
        ts, h1, h2, D, PP = [], [], [], [], []
        t2 = time.time()
        for s in range(pasos + 1):
            if s % nsal == 0:
                ts.append(t0 + s*dt)
                h1.append(np.sum(f*eQ[None, :])*peso/mg); h2.append(np.sum(f*(eQ**2)[None, :])*peso/mg)
                P = np.interp(rmed, rs, dPhi(f)[orden])
                D.append(np.sqrt(np.mean(P**2))/esc); PP.append(P)
            if s == pasos:
                break
            f = golpe(f, 0.5*dt)
            f = np.real(np.fft.ifft(np.fft.fft(f, axis=1)*fase, axis=1))
            f = golpe(f, 0.5*dt)
        if verboso:
            print(f'  {"libre" if libre else "lineal"}: {pasos} pasos en {time.time() - t2:.0f} s', flush=True)
        return np.array(ts), np.array(h1), np.array(h2), np.array(D), np.array(PP)
    tl, h1l, h2l, Dl, Pl = correr(False)
    tf, h1f, h2f, Df, Pf = correr(True)
    # Mezcla libre exacta de las partículas en Phi_0: cada una con su ángulo, su acción y su frecuencia.
    Om = omega_de_J(E_t, J_t)(J0)
    h1x = np.array([np.sum(mi*np.exp(-1j*(Q0 + Om*(s - t0))))/mg for s in tl])
    h2x = np.array([np.sum(mi*np.exp(-2j*(Q0 + Om*(s - t0))))/mg for s in tl])
    # Radios de las partículas en su órbita: coeficientes de Fourier en Q de r(Q, J) de los nodos, interpolados
    # linealmente en s = sqrt(J), con un nodo más en la órbita circular (J = 0: r = r_c en todo Q). Cerca de J = 0 el
    # primer coeficiente crece como sqrt(J): en s es lineal.
    Rk = np.fft.fft(r.reshape(nj, nq), axis=1)/nq
    R0k = np.zeros(nq, complex); R0k[0] = float(m.rc)
    sn = np.concatenate([[0.0], np.sqrt(Jn)]); Rn = np.vstack([R0k[None, :], Rk])
    s = np.sqrt(np.maximum(J0, 0.0))
    b0 = np.clip(np.searchsorted(sn, s) - 1, 0, nj - 1); al = ((s - sn[b0])/(sn[b0 + 1] - sn[b0]))[:, None]
    Rp = (1 - al)*Rn[b0] + al*Rn[b0 + 1]

    def radio(fase):
        return np.real(np.sum(Rp*np.exp(1j*k[None, :]*(fase[:, None] - 0.5*dQ)), axis=1))
    ri0 = radio(Q0)
    if 'R0' in e.files:
        er = np.abs(ri0 - e['R0'][ok]); err_r = er.max()
        if verboso:
            print(f'  radios interpolados frente a los de la corrida en t0: mediana {np.median(er):.1e}, percentil 99 '
                  f'{np.percentile(er, 99):.1e}, máximo {err_r:.1e} (en J = {J0[np.argmax(er)]:.2e})', flush=True)
    else:
        err_r = np.nan

    potencial = potencial_capas
    # Potencial de F_0: cada partícula repartida en su órbita (NQ puntos), una vez.
    fases = (np.arange(nq) + 0.5)*dQ
    rF0 = np.concatenate([radio(np.full_like(Q0, f)) for f in fases]); mF0 = np.tile(mi, nq)/nq
    PF0 = potencial(rF0, mF0, rmed)
    Px = np.array([potencial(radio(Q0 + Om*(s - t0)), mi, rmed) - PF0 for s in tl])
    Dx = np.sqrt(np.mean(Px**2, axis=1))/esc
    Dc = np.sqrt(np.mean((Px + Pl - Pf)**2, axis=1))/esc
    # La corrida con el mismo estimador: capas exactas de sus partículas, respecto del potencial de F_0.
    Dq = np.sqrt(np.mean((e['Pq'] - PF0)**2, axis=1))/esc if 'Pq' in e.files else None
    return dict(t=tl, h1_lin=h1l, h2_lin=h2l, D_lin=Dl, h1_lib=h1f, h2_lib=h2f, D_lib=Df, autoc=autoc,
                h1_libx=h1x, h2_libx=h2x, h1_cor=h1x + h1l - h1f, h2_cor=h2x + h2l - h2f, D_libx=Dx, D_cor=Dc, err_r=err_r,
                **({} if Dq is None else {'D_runq': Dq}),
                h1_part=h1_part, h1_rej=h1_rej, Jn=Jn, F0=F0, dF0=dF0, om=om, nj=nj, nq=nq, dt=dt,
                mascara=np.array(mascara if mascara is not None else (0.0, np.inf)))


def lineal(ext, sal, **kw):
    e = np.load(ext)
    R = evolucion(e, **kw)
    np.savez(sal, **R, T0=e['T0'], tau=e['tau'], t0=e['t0'], tq=e['tq'], hq=e['hq'], t_run=e['t'], D_run=e['D'],
             mgas=e['mgas'], malos=e['malos'])
    print(f'escrito {sal}: autoconsistencia de Phi_0 {R["autoc"]:.1e} (D en t0: {e["D"][np.argmin(abs(e["t"] - e["t0"]))]:.1e}); '
          f'h_1(t0) partículas {abs(R["h1_part"]):.4e}, rejilla {abs(R["h1_rej"]):.4e}; error de los radios interpolados '
          f'{float(R["err_r"]):.1e}', flush=True)


def envolvente(x, t, ancho=50.0):
    n = max(1, int(round(ancho/(t[1] - t[0]))))
    return np.array([x[max(0, k - n):k + n + 1].max() for k in range(len(x))])


def pendiente(x, y, a, b):
    z = (x >= a) & (x <= b) & (y > 0)
    return np.polyfit(np.log(x[z]), np.log(y[z]), 1)[0]


def comparar(archivos, xs=(6, 8, 12, 16, 20, 25)):
    """Envolventes de |h_1|, |h_2| y D de la corrida, de la mezcla libre exacta y de la evolución lineal corregida."""
    for a in archivos:
        d = np.load(a); tau = float(d['tau'])
        tl = d['t']/tau; tq = d['tq']/tau; tr = d['t_run']/tau
        print(f'\n{os.path.basename(a)}: t0 = {float(d["T0"]):g} tau_1; autoconsistencia de Phi_0 {float(d["autoc"]):.1e}; '
              f'error de los radios {float(d["err_r"]):.1e}')
        Drun, xD = (d['D_runq'], tq) if 'D_runq' in d.files else (d['D_run'], tr)
        for nombre, run, lib, cor, xr in (('|h_1|/h_0', np.abs(d['hq'][:, 0]), np.abs(d['h1_libx']), np.abs(d['h1_cor']), tq),
                                         ('|h_2|/h_0', np.abs(d['hq'][:, 1]), np.abs(d['h2_libx']), np.abs(d['h2_cor']), tq),
                                         ('D', Drun, d['D_libx'], d['D_cor'], xD)):
            tt = d['tq'] if xr is tq else d['t_run']
            e = {'corrida': envolvente(run, tt), 'libre': envolvente(lib, d['t']), 'lineal': envolvente(cor, d['t'])}
            x = {'corrida': xr, 'libre': tl, 'lineal': tl}
            print(f'  {nombre} (variables de Phi_0), envolvente en t/tau_1 = ' + ', '.join(str(v) for v in xs) + ' | pendiente 10-20 tau_1')
            for q in ('corrida', 'libre', 'lineal'):
                print(f'    {q:8} ' + ' '.join(f'{e[q][np.argmin(abs(x[q] - v))]:8.1e}' for v in xs) + f' | {pendiente(x[q], e[q], 10, 20):6.2f}')
            for q in ('libre', 'lineal'):
                print(f'    corrida/{q:7}' + ' '.join(f'{e["corrida"][np.argmin(abs(xr - v))]/e[q][np.argmin(abs(tl - v))]:8.2f}' for v in xs))


def cocientes(d, xs=(8, 12, 16, 20, 25)):
    """Cocientes de las envolventes de |h_1| (corrida/lineal y corrida/libre) en t/tau_1 = xs, y de D hasta 12 tau_1."""
    tau = float(d['tau']); tl = d['t']/tau; tq = d['tq']/tau
    er = envolvente(np.abs(d['hq'][:, 0]), d['tq'])
    el, ef = envolvente(np.abs(d['h1_cor']), d['t']), envolvente(np.abs(d['h1_libx']), d['t'])
    v = lambda e, x, s: e[np.argmin(np.abs(x - s))]
    cl = [v(er, tq, s)/v(el, tl, s) for s in xs]; cf = [v(er, tq, s)/v(ef, tl, s) for s in xs]
    Dl = Df = None
    if 'D_runq' in d.files:
        eD = envolvente(d['D_runq'], d['tq']); eDl, eDf = envolvente(d['D_cor'], d['t']), envolvente(d['D_libx'], d['t'])
        Dl = [v(eD, tq, s)/v(eDl, tl, s) for s in (8, 12)]; Df = [v(eD, tq, s)/v(eDf, tl, s) for s in (8, 12)]
    return cl, cf, Dl, Df


def nombre_caso(a):
    b = os.path.basename(a)
    return os.path.basename(os.path.dirname(a)) if b.startswith('respuesta_') else b.split('_t5')[0]


def resumen(archivos, xs=(8, 12, 16, 20, 25)):
    print('Cocientes de las envolventes de |h_1| en t/tau_1 = ' + ', '.join(str(s) for s in xs) + ': corrida/lineal y corrida/libre;')
    print('de D en 8 y 12 tau_1 (después, piso de la reconstrucción). 1: la curva sigue a la corrida.')
    for a in archivos:
        d = np.load(a); cl, cf, Dl, Df = cocientes(d, xs)
        txt = f'{nombre_caso(a):13} lineal ' + ' '.join(f'{x:5.2f}' for x in cl) + '   libre ' + ' '.join(f'{x:5.2f}' for x in cf)
        if Dl is not None:
            txt += '   | D: lineal ' + ' '.join(f'{x:5.2f}' for x in Dl) + '  libre ' + ' '.join(f'{x:5.2f}' for x in Df)
        print(txt)
    figura(archivos)


def figura(archivos, sal=None):
    """|h_1| de cada corrida desde t0, frente a la evolución lineal corregida y la mezcla libre exacta."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    n = len(archivos); nc = 3; nf = (n + nc - 1)//nc
    fig, ax = plt.subplots(nf, nc, figsize=(11.5, 3.0*nf), squeeze=False)
    for i, a in enumerate(archivos):
        d = np.load(a); tau = float(d['tau']); c = ax[i//nc, i % nc]
        c.plot(d['tq']/tau, envolvente(np.abs(d['hq'][:, 0]), d['tq']), color=tinta, lw=1.6)
        c.plot(d['t']/tau, envolvente(np.abs(d['h1_cor']), d['t']), color='#2a78d6', lw=1.3, ls=(0, (4, 2)))
        c.plot(d['t']/tau, envolvente(np.abs(d['h1_libx']), d['t']), color='#eb6834', lw=1.3, ls=(0, (1, 1.2)))
        c.set_xscale('log'); c.set_yscale('log'); c.set_xlim(float(d['T0']), 30)
        c.set_title(nombre_caso(a), loc='left', fontsize=9.5, pad=4)
        c.grid(True, color=reja, lw=0.5); c.set_axisbelow(True)
        for lado in ('top', 'right'):
            c.spines[lado].set_visible(False)
        if i % nc == 0:
            c.set_ylabel(r'$|h_1|/h_0$')
        if i//nc == nf - 1:
            c.set_xlabel(r'$t/\tau_1$')
    for j in range(n, nf*nc):
        ax[j//nc, j % nc].axis('off')
    fig.legend(handles=[Line2D([], [], color=tinta, lw=1.6, label='run'),
                        Line2D([], [], color='#2a78d6', lw=1.3, ls=(0, (4, 2)), label='linear response about the mixed state'),
                        Line2D([], [], color='#eb6834', lw=1.3, ls=(0, (1, 1.2)), label='phase mixing alone')],
               loc='upper center', ncol=3, frameon=False, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.97 if nf > 1 else 0.92))
    sal = sal or os.path.join(os.path.dirname(AQUI), '..', 'exe', 'demo_eta', 'respuesta', 'respuesta.png')
    sal = os.path.abspath(sal)
    fig.savefig(sal, dpi=120); plt.close(fig)
    print('escrita', sal)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['extraer', 'lineal', 'comparar', 'resumen'])
    ap.add_argument('archivos', nargs='+')
    ap.add_argument('--mgas', type=float); ap.add_argument('--t0', type=float, default=5.0)
    ap.add_argument('--tau', type=float, default=327.15); ap.add_argument('--procesos', type=int, default=4)
    ap.add_argument('--nj', type=int, default=2400); ap.add_argument('--nq', type=int, default=64)
    ap.add_argument('--dt', type=float, default=0.5); ap.add_argument('--t1', type=float, default=30.0)
    ap.add_argument('--mascara', type=float, nargs=2, help='el campo actúa solo sobre el gas con J en [jlo, jhi)')
    a = ap.parse_args()
    if a.paso == 'extraer':
        extraer(a.archivos[0], a.archivos[1], a.mgas, a.t0, a.tau, procesos=a.procesos)
    elif a.paso == 'lineal':
        lineal(a.archivos[0], a.archivos[1], nj=a.nj, nq=a.nq, dt=a.dt, T1=a.t1, mascara=a.mascara)
    elif a.paso == 'resumen':
        resumen(a.archivos)
    else:
        comparar(a.archivos)
