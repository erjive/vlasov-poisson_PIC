"""Prueba decisiva de la meseta: h_1 con las variables ángulo-acción del potencial real.

Si la parte estática de la meseta de |h_1| se debe a calcular (Q,J) con el mapa del
isócrono cuando el potencial es isócrono + autogravedad, al recalcular h_1 con el
mapa numérico del potencial real esa parte estática debe desaparecer.

Se calcula h_1 con tres mapas sobre las mismas partículas:
  iso        el analítico del isócrono (lo que hace el código);
  promedio   el numérico del potencial promediado en t >= T_PROM;
  instante   el numérico del potencial de cada instantánea.

El prefactor 8 pi^2 L0 a_1 se toma del propio código (cociente entre su h_1 y la
suma con el mapa isócrono), así que |h_1|/h_0 es directamente comparable con hk1.tl.

Uso:  python3 aa_meseta.py [corrida]      (por omisión long20k_snap)
Escribe exe/sg/<corrida>/aa_meseta.npz
"""
import os, sys, time, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from df0 import rp_to_QJ
from aa_numerico import MapaAA, phi_iso

SG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'sg')
J1, SJ1 = 0.10, 0.10
T_PROM = 2000.0

corrida = sys.argv[1] if len(sys.argv) > 1 else 'long20k_snap'
f = h5py.File(os.path.join(SG, corrida, 'vlasov_output.h5'), 'r')
pasos = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
t = np.array([f[k].attrs['time'] for k in pasos])
rg = f['grid']['r'][:]
w = f[pasos[0]]['f'][:]
PS = np.array([f[k]['potential'][:] - phi_iso(rg) for k in pasos])
mapa_prom = MapaAA(rg, PS[t >= T_PROM].mean(0))
print(f'{corrida}: {len(t)} instantáneas; potencial promediado en t >= {T_PROM:g} '
      f'({np.sum(t >= T_PROM)} instantáneas); r_c = {mapa_prom.rc:.4f}')

B = lambda J: np.exp(-(J - J1)**2/SJ1**2)*J**2
S = {k: np.empty(len(t), complex) for k in ('iso', 'promedio', 'instante')}
H0 = {k: np.empty(len(t)) for k in S}
idx = [4978, 7500, 9900]
Jpart = {k: np.empty((len(t), 3)) for k in S}
t0 = time.time()
for n, k in enumerate(pasos):
    r = f[k]['r_part'][:]; p = f[k]['p_part'][:]
    Qi, Ji = rp_to_QJ(r, p)
    Qp, Jp, _ = mapa_prom(r, p)
    Qn, Jn, _ = MapaAA(rg, PS[n])(r, p)
    for nombre, Q, J in (('iso', Qi, Ji), ('promedio', Qp, Jp), ('instante', Qn, Jn)):
        S[nombre][n] = np.sum(w*B(J)*np.exp(-1j*Q))
        H0[nombre][n] = np.sum(w*B(J))
        Jpart[nombre][n] = J[idx]
    if n % 50 == 0:
        print(f'  {n:4d}/{len(t)}  t={t[n]:7.0f}  ({time.time()-t0:.0f} s)', flush=True)

# Prefactores del código: h_1 = K1 * sum f B e^{-iQ}, h_0 = K0 * sum f B.
a = np.loadtxt(os.path.join(SG, corrida, 'hk1_complex.tl'))
h1c = np.interp(t, a[:, 0], a[:, 3]) + 1j*np.interp(t, a[:, 0], a[:, 4])
h0c = np.interp(t, a[:, 0], a[:, 1])
m = np.abs(h1c) > 1e-2*np.abs(h1c).max()
K1 = np.mean(h1c[m]/S['iso'][m]); K0 = np.mean(h0c/H0['iso'])
print(f'prefactores: K1 = {K1.real:.6e} (fase {np.angle(K1):.1e}), K0 = {K0:.6e};'
      f' comprobación h_0: {np.max(np.abs(K0*H0["iso"]/h0c - 1)):.1e}')
h0_ref = h0c[0]
h1 = {k: K1*S[k]/h0_ref for k in S}
h0 = {k: K0*H0[k]/h0_ref for k in S}
np.savez(os.path.join(SG, corrida, 'aa_meseta.npz'), t=t,
         **{f'h1_{k}': v for k, v in h1.items()}, **{f'h0_{k}': v for k, v in h0.items()},
         **{f'J_{k}': v for k, v in Jpart.items()})

rms = lambda x: np.sqrt(np.mean(np.abs(x)**2))
print(f"\n{'mapa':>9} {'ventana':>16} {'|<h1>|/h0':>11} {'arg':>7} {'rms(h1-<h1>)/h0':>16} {'mediana |h1|/h0':>16}")
for k in h1:
    for lo, hi in [(1300, 2000), (2000, 15000), (15000, 20000)]:
        v = (t >= lo) & (t <= hi)
        if not v.any():
            continue
        z = h1[k][v]; s_ = z.mean()
        print(f'{k:>9} {f"[{lo},{hi}]":>16} {abs(s_):>11.3e} {np.angle(s_):>+7.3f}'
              f' {rms(z - s_):>16.3e} {np.median(np.abs(z)):>16.3e}')
print(f"\n{'mapa':>9} {'h0(t=final)/h0(0)-1':>20}")
for k in h0:
    print(f'{k:>9} {h0[k][-1]/h0[k][0]-1:>+20.4e}')
print('\noscilación de J por partícula (t > 2000), pico a pico relativo:')
for k in Jpart:
    v = t > 2000
    print(f'{k:>9}: ' + '  '.join(f'{100*np.ptp(Jpart[k][v, i])/Jpart[k][v, i].mean():.3f}%'
                                  for i in range(3)))
