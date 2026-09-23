"""Parámetros de diseño de una corrida antes de lanzarla: banda de frecuencias,
corrimiento por autogravedad y el cociente eta entre ambos.

Con L fijo, cada órbita gira con su frecuencia radial Omega(J) = dE/dJ. El
soporte de F_eq(J) ocupa entonces una banda [Omega_min, Omega_max]:

  - Su ANCHO fija la mezcla de fases. Una perturbación con dependencia e^{ikQ}
    se desfasa en tau_k = 2 pi/(k dOmega), y eso ocurre con o sin autogravedad.
  - La autogravedad CORRE la banda entera hacia arriba (el potencial se hunde,
    las órbitas se aceleran) en una cantidad proporcional a la masa a0.

Para que un modo colectivo escape del continuo y deje de amortiguarse, el
corrimiento tiene que ser comparable al ancho:

  eta = (corrimiento de la frecuencia media) / (ancho de la banda).

eta es un criterio de DISEÑO, no un teorema: dice si la corrida está en el
rango donde la transición puede ocurrir. Quién decide es la simulación, que
mide la frecuencia de la cola y la compara con la banda.

Además reporta las dos escalas de tiempo que fijan cuánto hay que correr:

  tau_k  = 2 pi/(k dOmega)                 mezcla de fases; t_final ~ 10 tau_1
  T_rec  = 2 pi/(k |dOmega/dJ| J_max/Nrc)  recurrencia por aliasing de la malla
                                           en J: más allá, la meseta es un
                                           artefacto, y el límite baja con k

Uso:
    python3 eta.py --sigma 0.05 --a0 0.008 0.025 0.05 0.083 0.125
    python3 eta.py --sigma 0.03 --l0 1 --a0 0.01 0.04 --nrc 400
"""
import os, sys, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from equilibrio import Equilibrio, F_forma, SIGMA_J, J_MAX
from aa_numerico import L0 as L0_OMISION

UMBRAL = 0.01          # soporte: donde F_eq supera esta fracción de su máximo


def banda(a0, sigma, jmax, l0, nj=4000, verboso=False):
    """Banda de Omega sobre el soporte de F_eq, en el equilibrio autoconsistente."""
    eq = Equilibrio(a0, sigma=sigma, jmax=jmax, l0=l0)
    if a0 > 0:
        eq.iterar(tol=1e-10, maxit=25, verboso=verboso)
    else:                                   # a0 = 0: el isócrono desnudo
        m = eq.mapa()
        eq.E_t, eq.J_t = eq.tabla_J_de_E(m)
    J = np.linspace(1e-4, jmax, nj)
    E = np.interp(J, eq.J_t, eq.E_t)
    Om = np.gradient(E, J)                  # Omega = dE/dJ
    F = F_forma(J, sigma)
    sop = F >= UMBRAL*F.max()
    media = np.sum(F[sop]*Om[sop])/np.sum(F[sop])
    dOm = np.gradient(Om, J)
    return dict(J_lo=J[sop][0], J_hi=J[sop][-1], Om_min=Om[sop].min(),
                Om_max=Om[sop].max(), Om_media=media,
                dOmdJ=np.abs(np.interp(0.5*(J[sop][0] + J[sop][-1]), J, dOm)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a0', type=float, nargs='+', required=True)
    ap.add_argument('--sigma', type=float, default=SIGMA_J)
    ap.add_argument('--jmax', type=float, default=None,
                    help='por omisión, 6 sigma (el soporte de J^2 exp(-J^2/sigma^2))')
    ap.add_argument('--l0', type=float, default=L0_OMISION)
    ap.add_argument('--nrc', type=int, default=400, help='nodos en J, para la recurrencia')
    ap.add_argument('--kmax', type=int, default=4, help='modo más alto que se quiere fiable')
    ap.add_argument('--npart', type=int, default=None,
                    help='partículas, para estimar el costo (por omisión nrc*25)')
    arg = ap.parse_args()
    jmax = arg.jmax if arg.jmax is not None else 6*arg.sigma
    npart = arg.npart if arg.npart is not None else 25*arg.nrc

    ref = banda(0.0, arg.sigma, jmax, arg.l0)
    c = 0.5*(arg.l0 + np.sqrt(arg.l0**2 + 4))
    print(f'L0 = {arg.l0:g}   sigma_J = {arg.sigma:g}   J_max = {jmax:g}   '
          f'Nrc = {arg.nrc}   N = {npart}')
    print(f'soporte (F > {UMBRAL:g} F_max): J en [{ref["J_lo"]:.4f}, {ref["J_hi"]:.4f}];  '
          f'periodo radial 2 pi/Omega = {2*np.pi/ref["Om_media"]:.1f}')
    print()
    cab = (f'{"a0":>8} {"banda Omega":>21} {"ancho/media":>12} {"corrim.":>9} '
           f'{"eta":>7} {"tau_1":>8} {"t final":>9} {"T_rec(k=%d)" % arg.kmax:>11} '
           f'{"pasos":>9} {"costo 4 hilos":>14}')
    print(cab); print('-'*len(cab))
    for a0 in arg.a0:
        b = banda(a0, arg.sigma, jmax, arg.l0)
        ancho = b['Om_max'] - b['Om_min']
        rel = ancho/b['Om_media']
        corr = abs(b['Om_media'] - ref['Om_media'])/ref['Om_media']
        eta = corr/rel
        tau1 = 2*np.pi/ancho
        tfin = 10*tau1
        trec = 2*np.pi/(arg.kmax*b['dOmdJ']*jmax/arg.nrc)
        # dt = courant*dr/pmax con los valores de 11_landau (1, 0.1, 2).
        dt = 0.05
        pasos = tfin/dt
        # 1.67e-7 s por partícula y paso, medido con 4 hilos, yoshida4 y autogravedad.
        seg = 1.67e-7*pasos*npart
        print(f'{a0:8.4f} [{b["Om_min"]:.5f}, {b["Om_max"]:.5f}] {rel:12.3f} '
              f'{corr:9.4f} {eta:7.3f} {tau1:8.0f} {tfin:9.0f} {trec:11.0f} '
              f'{pasos:9.0f} {seg/60:11.1f} min')
    print()
    print(f'Avisos: t final = 10 tau_1; si T_rec(k={arg.kmax}) < t final, sube Nrc '
          f'(hace falta Nrc >~ {arg.kmax}*t_final*|dOmega/dJ|*J_max/2pi).')
    print('El eta de la tabla es de diseño; el de la corrida se mide de la cola '
          'de dPhi y se compara con la banda.')


if __name__ == '__main__':
    main()
