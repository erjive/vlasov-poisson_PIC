"""h_k(t) exacto para una F0(Q,J) arbitraria, separable o no.

Lo que analysish.f90 calcula es

    h_k = 8 pi^2 L0 * a_k * <b(J) e^{-ikQ}>_F ,

con a_k el coeficiente de Fourier en Q de la funcion de prueba y
b(J) su parte en J.  Bajo phase mixing F(Q,J,t) = F0(Q - omega(J) t, J), asi que

    h_k(t) = a0 * a_k * Int dJ b(J) c_k(J) e^{-i k omega(J) t} / Int dJ c_0(J),

donde c_k(J) es el coeficiente de Fourier en Q de F0 a cada J.  Para una F0 no
separable c_k depende de J y es complejo; se obtiene por FFT sobre Q, que para
una funcion periodica y suave converge espectralmente.
"""
import numpy as np
from df0 import df0, omega, Jrange

def ck_of_J(dftype, Js, kmax, nQ=256):
    """c_k(J) para k = 0..kmax, por FFT sobre Q."""
    Q = np.arange(nQ)*(2*np.pi/nQ)
    out = {k: np.empty(len(Js), complex) for k in range(kmax+1)}
    # Por trozos: la matriz (nJ, nQ) completa puede no caber en memoria.
    for a in range(0, len(Js), 20000):
        sl = slice(a, a+20000)
        C = np.fft.fft(df0(Q[None, :], Js[sl, None], dftype), axis=1)/nQ
        for k in range(kmax+1):
            out[k][sl] = C[:, k]
    return out

def ak_test(sq, k, nQ=200000):
    """Coeficiente de Fourier en Q de la funcion de prueba, como lo forma
    analysish.f90: (1/2pi) Int exp(-sin(Q/2)^2/sq^2) cos(kQ) dQ."""
    Q = np.arange(nQ)*(2*np.pi/nQ)
    return np.mean(np.exp(-np.sin(0.5*Q)**2/sq**2)*np.cos(k*Q))

def bJ_test(Js, j0, sj):
    """Parte en J de la funcion de prueba, como en analysish.f90."""
    return np.exp(-(Js-j0)**2/sj**2)*Js**2

def make_hk(dftype, a0, j0, sj, sq, modes=(0,1,2,3,4), nJ=20001, Jhi=None):
    Jlo0, Jhi0 = Jrange(dftype)
    Js = np.linspace(0.0, Jhi if Jhi is not None else max(Jhi0, 0.45), nJ)
    ck = ck_of_J(dftype, Js, max(modes))
    b  = bJ_test(Js, j0, sj)
    om = omega(Js)
    norm = np.trapezoid(ck[0].real, Js)
    ak = {k: ak_test(sq, k) for k in modes}
    def hk(k, t):
        t = np.atleast_1d(np.asarray(t, float))
        w = b*ck[k]
        out = np.array([np.trapezoid(w*np.exp(-1j*k*om*ti), Js) for ti in t])
        return a0*ak[k]*out/norm
    return hk
