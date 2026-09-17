"""Prueba del paso de tiempo de la sección 9: la corrida de referencia con dt=0.05, 0.1, 0.2.

Imprime, para cada dt, el cambio de energía, el cambio de h_0, la parte estática de h_1
en t in [1600,2000], el coeficiente del ajuste gaussiano (ln|h_1| frente a t^2 en t<250)
y la diferencia máxima de h_1 con la corrida de dt=0.1 hasta t=800.

Uso:  python3 prueba_dt.py        (lee exe/sg/dt_c1, efix_quad, dt_c4)
"""
import os, numpy as np, h5py

SG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exe', 'sg')


def leer(d):
    a = np.loadtxt(os.path.join(SG, d, 'hk1_complex.tl'))
    return a[:, 0], a[:, 1], (a[:, 3] + 1j*a[:, 4])/a[0, 1]


def energia(d):
    f = h5py.File(os.path.join(SG, d, 'vlasov_output.h5'), 'r')
    st = sorted([s for s in f if s.startswith('step_')], key=lambda s: int(s.split('_')[1]))
    return f[st[-1]].attrs['total_energy']/f[st[0]].attrs['total_energy'] - 1


tr, _, zr = leer('efix_quad')
print(f"{'dt':>5} {'dE/E':>11} {'h0 final-1':>13} {'|S| [1600,2000]':>16} {'coef t^2':>12} {'max|dh1|/|h1| t<800':>20}")
for d, dt in [('dt_c1', 0.05), ('efix_quad', 0.1), ('dt_c4', 0.2)]:
    t, h0, z = leer(d)
    m = (t >= 1600) & (t <= 2000)
    mm = (t > 0) & (t < 250)
    c = np.polyfit(t[mm]**2, np.log(np.abs(z[mm])), 1)[0]
    zi = np.interp(t, tr, zr.real) + 1j*np.interp(t, tr, zr.imag)
    w = t < 800
    print(f'{dt:>5} {energia(d):>+11.3e} {h0[-1]/h0[0]-1:>+13.5e} {abs(z[m].mean()):>16.5e}'
          f' {c:>12.4e} {np.max(np.abs(z[w]-zi[w])/np.abs(zi[w])):>20.2e}')
