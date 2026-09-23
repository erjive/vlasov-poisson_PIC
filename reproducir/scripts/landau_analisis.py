"""Diagnósticos de una corrida que parte de un equilibrio autoconsistente.

Todo se mide en el marco FIJO del equilibrio (el potencial Phi_iso + Phi_self_eq de
equilibrio.py), que es el marco en que está formulada la teoría lineal:

  h_k(t)   = sum_j f_j B(J_j) exp(-i k Q_j) / sum_j f_j B(J_j) en t=0,
             con (Q,J) del mapa numérico del potencial de equilibrio;
  dPhi     = potencial del código - Phi_iso - Phi_self_eq, en la malla radial:
             no depende de ningún mapa, y es lo que el amortiguamiento de Landau
             hace decaer;
  deriva J = cambio de la acción de equilibrio de cada partícula.

Uso:  python3 landau_analisis.py <directorio de la corrida> <archivo _equilibrio.npz>
Escribe <directorio>/landau.npz
"""
import os, sys, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico import MapaAA, phi_iso

J1_OMISION, SJ1_OMISION = 0.10, 0.10   # los de las corridas de 11_landau


def pesos_prueba(corrida):
    """(j1, sj1) de la función de prueba Phi_1 con la que el código calcula h_k.

    No son el ancho de F_eq: son parámetros de entrada que la corrida deja
    escritos en params_usados.par. Se leen de ahí para que h_k medido aquí y
    h_k del código sean la misma cantidad aunque cambie la configuración."""
    p = os.path.join(corrida, 'params_usados.par')
    j1, sj1 = J1_OMISION, SJ1_OMISION
    if os.path.exists(p):
        for l in open(p):
            l = l.split('#')[0].split('!')[0]
            if '=' in l:
                k, v = l.split('=', 1)
                if k.strip().lower() == 'j1':
                    j1 = float(v)
                elif k.strip().lower() == 'sj1':
                    sj1 = float(v)
    return j1, sj1


def analizar(corrida, eqfile, ventana_r=(3.0, 15.0), verboso=True):
    eq = np.load(eqfile)
    L = float(eq['L0']) if 'L0' in eq.files else 2.0
    mapa = MapaAA(eq['r'], eq['phi_self'], L=L)
    j1, sj1 = pesos_prueba(corrida)
    B = lambda J: np.exp(-(J - j1)**2/sj1**2)*J**2
    if verboso:
        print(f'  L0 = {L:g}, función de prueba: j1 = {j1:g}, sj1 = {sj1:g}')
    f = h5py.File(os.path.join(corrida, 'vlasov_output.h5'), 'r')
    pasos = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    t = np.array([f[k].attrs['time'] for k in pasos])
    rg = f['grid']['r'][:]
    ps_eq = np.interp(rg, eq['r'], eq['phi_self'])
    zona = (rg >= ventana_r[0]) & (rg <= ventana_r[1])
    escala = np.max(np.abs(ps_eq[zona]))
    w = f[pasos[0]]['f'][:]
    hk = np.empty((len(t), 5), complex)
    dphi = np.empty((len(t), len(rg)))
    derivaJ = np.empty(len(t))
    J0 = None
    for n, k in enumerate(pasos):
        Q, J, _ = mapa(f[k]['r_part'][:], f[k]['p_part'][:])
        if J0 is None:
            J0 = J.copy()
        pesos = w*B(J)
        hk[n] = [np.sum(pesos*np.exp(-1j*m*Q)) for m in range(5)]
        dphi[n] = f[k]['potential'][:] - phi_iso(rg) - ps_eq
        derivaJ[n] = np.average(np.abs(J - J0), weights=w)
    hk = hk/hk[0, 0].real
    np.savez(os.path.join(corrida, 'landau.npz'), t=t, hk=hk, dphi=dphi, r=rg,
             derivaJ=derivaJ, escala=escala)
    if verboso:
        rms = np.sqrt(np.mean(dphi[:, zona]**2, axis=1))/escala
        print(f'{corrida}: {len(t)} instantáneas; |Phi_self_eq| máx en r in {ventana_r} = {escala:.3e}')
        print(f"{'t':>8} {'h0-1':>10} {'|h1|':>10} {'|h2|':>10} {'rms dPhi/|Phi_eq|':>18} {'<|J-J0|>':>10}")
        for n in list(range(0, len(t), max(1, len(t)//10))) + [len(t)-1]:
            print(f'{t[n]:>8.0f} {hk[n,0].real-1:>+10.2e} {abs(hk[n,1]):>10.3e} {abs(hk[n,2]):>10.3e}'
                  f' {rms[n]:>18.3e} {derivaJ[n]:>10.2e}')
    return t, hk, dphi, rg, derivaJ, escala


if __name__ == '__main__':
    analizar(sys.argv[1], sys.argv[2])
