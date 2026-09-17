"""¿La componente que gira de h_1 es streaming libre después de la mezcla inicial?

En las variables ángulo-acción del potencial real (mapa numérico, potencial propio
promediado en t >= T0), cada partícula debería rotar libremente: J constante y
Q = Q_j(T0) + w_j (t - T0). Si es así, h_1 tardío es la suma de N rotadores libres,

    h_rec(t) = K1 sum_j f_j B(J_j) exp(-i (Q_j(T0) + w_j (t - T0))),

completamente determinada por el estado en T0. Se mide:

  1. rigidez de cada partícula: dispersión de J y residuo de fase de un ajuste lineal;
  2. si h_rec reproduce el h_1 medido en t > T0 (incluida la parte que gira);
  3. la estructura en J de G(J) = sum_{J_j en el bin} f_j B(J_j) e^{-i Q_j(T0)}, que
     es lo que la mezcla posterior convierte en h_1(t): una parte que gira persistente
     exige estructura fina en G.

Uso:  python3 rotadores.py [corrida]    (por omisión long20k_snap)
"""
import os, sys, time, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico import MapaAA, phi_iso

SG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'sg')
J1, SJ1 = 0.10, 0.10
T0 = 2000.0
B = lambda J: np.exp(-(J - J1)**2/SJ1**2)*J**2

corrida = sys.argv[1] if len(sys.argv) > 1 else 'long20k_snap'
f = h5py.File(os.path.join(SG, corrida, 'vlasov_output.h5'), 'r')
pasos = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
t = np.array([f[k].attrs['time'] for k in pasos])
rg = f['grid']['r'][:]
w = f[pasos[0]]['f'][:]
sel = np.where(t >= T0)[0]
ps = np.mean([f[pasos[n]]['potential'][:] - phi_iso(rg) for n in sel], axis=0)
mapa = MapaAA(rg, ps)

cache = os.path.join(SG, corrida, 'rotadores_QJ.npz')
if os.path.exists(cache):
    d = np.load(cache); Q, J = d['Q'], d['J']
else:
    Q = np.empty((len(sel), len(w))); J = np.empty_like(Q)
    t0 = time.time()
    for m, n in enumerate(sel):
        Q[m], J[m], _ = mapa(f[pasos[n]]['r_part'][:], f[pasos[n]]['p_part'][:])
        if m % 100 == 0:
            print(f'  {m}/{len(sel)} ({time.time()-t0:.0f} s)', flush=True)
    np.savez(cache, Q=Q, J=J)
ts = t[sel] - T0

d = np.load(os.path.join(SG, corrida, 'aa_meseta.npz'))
h_med = d['h1_promedio'][sel]
K1 = 1.688514e-02                         # prefactor del código (aa_meseta.py)
a = np.loadtxt(os.path.join(SG, corrida, 'hk1_complex.tl'))
h0_ref = a[0, 1]
c = K1*w*B(J[0])/h0_ref

# 1. rigidez
Jm = J.mean(0)
fase = np.unwrap(Q, axis=0)
coef = np.polyfit(ts, fase, 1)
w_j = coef[0]
res = fase - (np.outer(ts, coef[0]) + coef[1])
peso = np.abs(c); P = peso > 1e-3*peso.max()
print(f'{corrida}: {len(w)} partículas, {len(sel)} instantáneas en t >= {T0:g}')
print(f'1. rigidez (partículas con peso): dispersión relativa de J, mediana '
      f'{np.median((J.std(0)/Jm)[P]):.1e}, máx {np.max((J.std(0)/Jm)[P]):.1e};'
      f' residuo de fase rms, mediana {np.median(np.sqrt((res**2).mean(0))[P]):.1e},'
      f' máx {np.max(np.sqrt((res**2).mean(0))[P]):.1e} rad')

# 2. reconstrucción desde el estado en T0
rms = lambda x: np.sqrt(np.mean(np.abs(x)**2))
for nombre, wj in [('frecuencia ajustada', w_j)]:
    h_rec = np.array([np.sum(c*np.exp(-1j*(Q[0] + wj*tt))) for tt in ts])
    for lo, hi in [(2000, 6000), (6000, 12000), (12000, 20000)]:
        v = (t[sel] >= lo) & (t[sel] <= hi)
        Rm = h_med[v] - h_med[v].mean(); Rr = h_rec[v] - h_rec[v].mean()
        corr = abs(np.vdot(Rm, Rr))/np.linalg.norm(Rm)/np.linalg.norm(Rr)
        print(f'2. [{lo},{hi}] {nombre}: rms giro medido {rms(Rm):.3e}, reconstruido '
              f'{rms(Rr):.3e}, correlación {corr:.4f}, rms(medido-reconstruido) '
              f'{rms(h_med[v]-h_rec[v]):.2e}')

# 3. estructura en J de G(J) en T0
orden = np.argsort(Jm)
for nb in [100, 400, 1600]:
    bordes = np.linspace(Jm[P].min(), Jm[P].max(), nb + 1)
    k = np.clip(np.digitize(Jm, bordes) - 1, 0, nb - 1)
    G = np.bincount(k, weights=(c*np.exp(-1j*Q[0])).real, minlength=nb) \
      + 1j*np.bincount(k, weights=(c*np.exp(-1j*Q[0])).imag, minlength=nb)
    n_bin = np.bincount(k[P], minlength=nb)
    # Aspereza: G menos su versión suavizada (media móvil de 5 bins).
    liso = np.convolve(G, np.ones(5)/5, mode='same')
    print(f'3. G(J) en {nb} bins (ancho {bordes[1]-bordes[0]:.2e}, partículas por bin '
          f'mediana {np.median(n_bin[n_bin>0]):.0f}): sum|G| = {np.abs(G).sum():.3e},'
          f' |sum G| = {abs(G.sum()):.3e}, rms(G - suave)/rms(G) = '
          f'{rms((G-liso)[2:-2])/rms(G[2:-2]):.3f}')
np.savez(os.path.join(SG, corrida, 'rotadores.npz'), t=t[sel], h_med=h_med, h_rec=h_rec,
         Jm=Jm, w_j=w_j, c=c, Q0=Q[0])
