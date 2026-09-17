"""h_1 descompuesto por filas de acción, con y sin autogravedad.

Con cuadratura, las partículas nacen en una rejilla de Nrc filas de J por Npc
columnas de Q. Cada fila tiene una sola acción inicial, así que bajo phase
mixing puro rota rígidamente a su frecuencia omega(J) y su contribución a h_1
nunca decae: el decaimiento de h_1 es la cancelación entre filas.

Escribiendo h_1(t) = sum_i g_i(t), con g_i la contribución de la fila i, se
ajusta cada fila como g_i(t) = S_i + A_i exp(-i w_i t) y se mide:

  1. el corrimiento de frecuencia de cada fila frente a omega_iso(J);
  2. la parte estática sum S_i, cuánto se cancelan las S_i entre filas y en qué
     J vive;
  3. la parte que gira del total frente al continuo C = sum A_i exp(-i w_i t) y
     frente a la suma de residuos del ajuste.

por_ventanas() repite el ajuste en ventanas sucesivas, para ver si las filas
siguen cambiando después de la mezcla inicial.

Uso:  python3 filas_J.py [directorio de exe/sg ...]
"""
import os, sys, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from df0 import rp_to_QJ, omega
from exact import ak_test

SG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'sg')
J1, SJ1 = 0.10, 0.10                     # función de prueba Phi_1
VENTANA = (1000.0, 2000.0)               # ya mezclado, antes del final
SQ1 = 0.40
# El código escribe h_k = 8 pi^2 L0 a_k sum f B e^{-ikQ}: h_1 lleva a_1 y h_0 lleva
# a_0. Aquí las sumas no llevan a_k, así que h_0 se divide por a_1/a_0 para que
# |S|/h_0 sea el mismo cociente que se lee en hk1.tl.
A1_A0 = ak_test(SQ1, 1)/ak_test(SQ1, 0)


def cargar(d):
    f = h5py.File(os.path.join(SG, d, 'vlasov_output.h5'), 'r')
    pasos = sorted([k for k in f if k.startswith('step_')],
                   key=lambda k: int(k.split('_')[1]))
    t = np.array([f[k].attrs['time'] for k in pasos])
    w = f[pasos[0]]['f'][:]
    Q = np.empty((len(pasos), len(w)))
    J = np.empty_like(Q)
    for n, k in enumerate(pasos):
        Q[n], J[n] = rp_to_QJ(f[k]['r_part'][:], f[k]['p_part'][:])
    return t, w, Q, J


def analizar(d):
    t, w, Q, J = cargar(d)
    # Filas: acción inicial, redondeada para agrupar la rejilla.
    J0 = J[0]
    claves, fila = np.unique(np.round(J0, 10), return_inverse=True)
    nf = len(claves)
    B = np.exp(-(J - J1)**2/SJ1**2)*J**2
    contrib = w[None, :]*B*np.exp(-1j*Q)          # (nt, Npart)
    g = np.zeros((len(t), nf), complex)
    for n in range(len(t)):
        g[n] = np.bincount(fila, weights=contrib[n].real, minlength=nf) \
             + 1j*np.bincount(fila, weights=contrib[n].imag, minlength=nf)
    h = g.sum(1)

    # Comprobación: la suma por filas reproduce el h_1 que escribe el código.
    a = np.loadtxt(os.path.join(SG, d, 'hk1_complex.tl'))
    href = np.interp(t, a[:, 0], a[:, 3]) + 1j*np.interp(t, a[:, 0], a[:, 4])
    m = np.abs(href) > 1e-3*np.abs(href).max()
    razon = href[m]/h[m]
    print(f'\n=== {d}: {nf} filas, {len(t)} instantáneas')
    print(f'  h_1 del código / suma por filas: {np.mean(razon.real):.6e}'
          f'  (dispersión relativa {np.std(np.abs(razon))/np.mean(np.abs(razon)):.1e},'
          f' fase {np.max(np.abs(np.angle(razon))):.1e})')
    h0 = (w*np.exp(-(J0 - J1)**2/SJ1**2)*J0**2).sum()/A1_A0

    # 1. frecuencia de cada fila
    v = (t >= 600) & (t <= 2000)
    peso = np.abs(g[0])
    principales = peso > 1e-3*peso.max()
    fase = np.unwrap(np.angle(g[v]), axis=0)
    w_med = -np.polyfit(t[v], fase, 1)[0]
    Jm = np.array([J[v][:, fila == i].mean() for i in range(nf)])
    dw = w_med - omega(claves)
    print('  1. corrimiento de frecuencia de las filas, omega_medida - omega_iso(J0):')
    for Jq in [0.05, 0.10, 0.15, 0.20, 0.30, 0.40]:
        i = np.argmin(np.abs(claves - Jq))
        print(f'     J0={claves[i]:.3f}: omega_iso={omega(claves[i]):.5f}  '
              f'dw={dw[i]:+.2e}  (dw/omega={dw[i]/omega(claves[i]):+.1e})'
              f'   <J_iso>-J0={Jm[i]-claves[i]:+.2e}')
    amp = np.abs(g[v]).mean(0)/np.maximum(np.abs(g[0]), 1e-300)
    print(f'     |g_i| tardío / inicial, filas con peso: mediana {np.median(amp[principales]):.4f}'
          f'  rango [{amp[principales].min():.4f}, {amp[principales].max():.4f}]')

    # 2. ajuste por fila: g_i(t) = S_i + A_i exp(-i w_i t) en t in [600,2000].
    #    La frecuencia sale de la pendiente de fase; S_i y A_i, de mínimos
    #    cuadrados lineales. Así la rotación propia de cada fila no contamina su
    #    parte estática, que es lo que falla al promediar en una ventana.
    tv = t[v]
    S_i = np.empty(nf, complex)
    A_i = np.empty(nf, complex)
    res = np.empty(nf)
    for i in range(nf):
        M = np.column_stack([np.ones_like(tv), np.exp(-1j*w_med[i]*tv)])
        coef, *_ = np.linalg.lstsq(M, g[v, i], rcond=None)
        S_i[i], A_i[i] = coef
        res[i] = np.linalg.norm(g[v, i] - M @ coef)/max(np.linalg.norm(g[v, i]), 1e-300)
    S = S_i.sum()
    print(f'  2. ajuste por fila: residuo relativo mediano {np.median(res[principales]):.1e},'
          f' máximo {res[principales].max():.1e}')
    print(f'     parte estática total |sum S_i|/h0 = {abs(S)/h0:.3e}, arg = {np.angle(S):+.3f}')
    print(f'     sum|S_i| / |sum S_i| = {np.abs(S_i).sum()/abs(S):.2f}  (1 = todas en fase)')
    orden = np.argsort(claves)
    acum = np.cumsum(S_i[orden].real)/S.real
    for q in [0.1, 0.5, 0.9]:
        print(f'     {int(q*100)}% de Re S acumulado en J0 <= '
              f'{claves[orden][min(np.searchsorted(acum, q), nf-1)]:.3f}')
    rel = np.abs(A_i)/np.maximum(np.abs(g[0]), 1e-300)
    print(f'     |A_i|/|g_i(0)| en filas con peso: mediana {np.median(rel[principales]):.4f},'
          f' rango [{rel[principales].min():.4f}, {rel[principales].max():.4f}]')
    for Jq in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40]:
        i = np.argmin(np.abs(claves - Jq))
        print(f'       J0={claves[i]:.3f}: |A|/|g(0)|={rel[i]:.4f}'
              f'  |S_i|/|A_i|={abs(S_i[i])/max(abs(A_i[i]),1e-300):.2e}')

    # 3. continuo extrapolado: C(t) = sum_i A_i exp(-i w_i t), hasta t = 20000.
    te = np.arange(2000.0, 20000.0 + 1, 2.0)
    C = np.array([np.sum(A_i*np.exp(-1j*w_med*tt)) for tt in te])
    for a_, b_ in [(2000, 15000), (15000, 20000)]:
        mm = (te >= a_) & (te <= b_)
        print(f'  3. continuo extrapolado, t in [{a_},{b_}]: rms|C|/h0 = '
              f'{np.sqrt(np.mean(np.abs(C[mm])**2))/h0:.2e}'
              f'  = {np.sqrt(np.mean(np.abs(C[mm])**2))/abs(S):.1%} de la parte estática')
    spc = np.abs(np.fft.fft(C))
    frc = -np.fft.fftfreq(len(te), 2.0)*2*np.pi
    print(f'     frecuencia dominante de C: {frc[np.argmax(spc)]:+.4f}')
    return dict(J0=claves, dw=dw, S_i=S_i, A_i=A_i, h0=h0, t=t, g=g)


if __name__ == '__main__':
    for d in (sys.argv[1:] or ['long_fino_nosg', 'long_fino']):
        analizar(d)


def por_ventanas(d, ancho=2000.0, t_ini=600.0):
    """Ajusta S_i + A_i exp(-i w_i t) por fila en ventanas sucesivas.

    Dice si las frecuencias y amplitudes de las filas siguen cambiando después
    de la mezcla inicial, y si la parte que gira del total la explica el
    continuo de filas (C) o los residuos, ventana por ventana.
    """
    t, w, Q, J = cargar(d)
    claves, fila = np.unique(np.round(J[0], 10), return_inverse=True)
    nf = len(claves)
    B = np.exp(-(J - J1)**2/SJ1**2)*J**2
    contrib = w[None, :]*B*np.exp(-1j*Q)
    g = np.zeros((len(t), nf), complex)
    for n in range(len(t)):
        g[n] = np.bincount(fila, weights=contrib[n].real, minlength=nf) \
             + 1j*np.bincount(fila, weights=contrib[n].imag, minlength=nf)
    h0 = (w*np.exp(-(J[0] - J1)**2/SJ1**2)*J[0]**2).sum()/A1_A0
    peso = np.abs(g[0]); P = peso > 1e-3*peso.max()
    rms = lambda x: np.sqrt(np.mean(np.abs(x)**2))
    print(f'\n=== {d}: ventanas de {ancho:g}')
    print(f"{'ventana':>16} {'|S|/h0':>10} {'arg S':>7} {'rms(h-S)/|S|':>13} {'rms C/|S|':>10}"
          f" {'rms res/|S|':>12} {'sum|S_i|/|S|':>13} {'<|A|/|g0|>':>11} {'dw medio':>10}")
    previo = None
    a = t_ini
    filas = []
    while a + ancho <= t[-1] + 1e-6:
        v = (t >= a) & (t <= a + ancho)
        tv = t[v]
        fase = np.unwrap(np.angle(g[v]), axis=0)
        w_i = -np.polyfit(tv, fase, 1)[0]
        S_i = np.empty(nf, complex); A_i = np.empty(nf, complex)
        for i in range(nf):
            M = np.column_stack([np.ones_like(tv), np.exp(-1j*w_i[i]*tv)])
            (S_i[i], A_i[i]), *_ = np.linalg.lstsq(M, g[v, i], rcond=None)
        S = S_i.sum()
        C = (A_i[None, :]*np.exp(-1j*np.outer(tv, w_i))).sum(1)
        res = g[v] - S_i[None, :] - A_i[None, :]*np.exp(-1j*np.outer(tv, w_i))
        h = g[v].sum(1)
        dw = np.average(w_i[P] - omega(claves[P]), weights=peso[P])
        print(f'{f"[{a:.0f},{a+ancho:.0f}]":>16} {abs(S)/h0:>10.4e} {np.angle(S):>+7.3f}'
              f' {rms(h - S)/abs(S):>13.3f} {rms(C)/abs(S):>10.3f} {rms(res.sum(1))/abs(S):>12.3f}'
              f' {np.abs(S_i).sum()/abs(S):>13.2f} {np.average(np.abs(A_i[P])/peso[P], weights=peso[P]):>11.4f}'
              f' {dw:>+10.2e}')
        filas.append((a, w_i, A_i, S_i))
        a += ancho
    # Cambio de las frecuencias y amplitudes de fila entre la primera y la última ventana.
    _, w1, A1, _ = filas[0]; _, w2, A2, _ = filas[-1]
    dW = (w2 - w1)[P]; dA = (np.abs(A2)/np.abs(A1) - 1)[P]
    print(f'  entre la primera y la última ventana, en filas con peso:')
    print(f'    cambio de w_i: medio {np.average(dW, weights=peso[P]):+.2e}, rms {np.sqrt(np.average(dW**2, weights=peso[P])):.2e}')
    print(f'    cambio relativo de |A_i|: medio {np.average(dA, weights=peso[P]):+.2e}, rms {np.sqrt(np.average(dA**2, weights=peso[P])):.2e}')
    # Rugosidad en J del cambio de frecuencia: lo que impide la cancelación suave.
    o = np.argsort(claves[P]); dWo = dW[o]
    liso = np.convolve(dWo, np.ones(9)/9, mode='same')
    print(f'    parte no suave del cambio de w_i (rms de w - media móvil de 9 filas): {np.sqrt(np.mean((dWo-liso)[4:-4]**2)):.2e}')
    return filas
