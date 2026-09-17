"""Referencia de phase mixing puro en el potencial de equilibrio.

Las mismas partículas iniciales de la corrida, rotando libres en las variables
ángulo-acción del equilibrio: Q(t) = Q0 + omega_eq(J0) t, J constante. Es lo que
pasaría si la perturbación no generara campo propio. La diferencia entre la
simulación autoconsistente y esta referencia es el efecto colectivo.

omega_eq(J) = dE/dJ se obtiene derivando la tabla J(E) del equilibrio.

Uso:  python3 landau_libre.py <corrida> <archivo _equilibrio.npz> <archivo .dat de la IC>
Escribe <corrida>/libre.npz
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico import MapaAA
from equilibrio import Equilibrio

J1, SJ1 = 0.10, 0.10
B = lambda J: np.exp(-(J - J1)**2/SJ1**2)*J**2


def omega_de_J(E_t, J_t):
    """omega = dE/dJ por diferencias centradas sobre la tabla (J creciente)."""
    dEdJ = np.gradient(E_t, J_t)
    return lambda J: np.interp(J, J_t, dEdJ)


def validar():
    """Con potencial propio nulo, la tabla debe dar omega = 1/(J+c)^3."""
    eq = Equilibrio(1e-12)
    m = MapaAA()                           # isócrono puro
    E_t, J_t = eq.tabla_J_de_E(m)
    om = omega_de_J(E_t, J_t)
    c = 0.5*(2 + np.sqrt(8))
    J = np.linspace(0.02, 0.55, 200)
    err = np.max(np.abs(om(J)*(J + c)**3 - 1))
    print(f'validación omega(J) en el isócrono: max|omega_tabla/omega_exacta - 1| = {err:.1e}')
    return err


def libre(corrida, eqfile, icfile, t):
    eq = np.load(eqfile)
    m = MapaAA(eq['r'], eq['phi_self'])
    r, p, F = np.loadtxt(icfile, unpack=True)
    Q0, J0, _ = m(r, p)
    om = omega_de_J(eq['E_t'], eq['J_t'])(J0)
    c = F*B(J0)
    h = np.array([[np.sum(c*np.exp(-1j*k*(Q0 + om*tt))) for k in range(5)] for tt in t])
    h = h/h[0, 0].real
    np.savez(os.path.join(corrida, 'libre.npz'), t=t, hk=h, omega=om, J0=J0, Q0=Q0, F=F)
    return h


if __name__ == '__main__':
    validar()
    d = np.load(os.path.join(sys.argv[1], 'landau.npz'))
    libre(sys.argv[1], sys.argv[2], sys.argv[3], d['t'])
    print('escrito', os.path.join(sys.argv[1], 'libre.npz'))
