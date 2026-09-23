"""Vlasov-Poisson linealizado alrededor del equilibrio, sin partículas.

En las variables ángulo-acción (Q,J) del equilibrio autoconsistente, con
F = F_eq(J) + dF(Q,J,t) y L fijo:

    d(dF)/dt + omega(J) d(dF)/dQ = F_eq'(J) d(dPhi)/dQ,

donde dPhi(r) es el potencial de la densidad de dF (solo la parte lineal: el campo
propio de la perturbación actúa sobre el equilibrio).

Discretización:
  * rejilla fija de NJ x NQ nodos en (Q,J); el radio r(Q,J) de cada nodo se calcula
    una vez con el mapa inverso del equilibrio (equilibrio.invertir);
  * transporte libre exacto en Q: en Fourier, f_k -> f_k exp(-i k omega dt);
  * golpe del campo: dF += dt F_eq'(J) dPhi/dQ, con la derivada en Q por FFT;
  * separación de Strang: medio golpe, transporte, medio golpe (segundo orden);
  * Poisson con la función de Green de capas esféricas, exacta para los nodos:
        dPhi(r_i) = -M(<r_i)/r_i - sum_{r_j > r_i} m_j/r_j,
    con m = 8 pi^2 L0 dF dQ dJ y la capa propia contada a medias. Los radios no
    cambian, así que el orden se calcula una vez.

h_1(t) = sum dF B(J) e^{-iQ} / sum F_eq B(J), con dF(t=0) = F_eq cos Q: es
directamente comparable con (h_1(eps) - h_1(0))/eps de las simulaciones.

Uso:  python3 lineal.py <archivo _equilibrio.npz> <salida.npz> [--nj 1600] [--nq 32]
      [--dt 0.5] [--tmax 3000] [--libre]
"""
import os, sys, argparse, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico import MapaAA, L0
from equilibrio import invertir, F_forma, J_MAX
from landau_libre import omega_de_J

J1, SJ1 = 0.10, 0.10
B = lambda J: np.exp(-(J - J1)**2/SJ1**2)*J**2


def params_eq(eqf):
    """(sigma_J, J_max, L0) del equilibrio; los omitidos son los de 11_landau."""
    d = np.load(eqf)
    sig = float(d['sigma_J']) if 'sigma_J' in d.files else 0.10
    jmx = float(d['J_max']) if 'J_max' in d.files else J_MAX
    l0 = float(d['L0']) if 'L0' in d.files else L0
    return sig, jmx, l0


def dF_forma(J, sigma):
    return np.exp(-J**2/sigma**2)*(2*J - 2*J**3/sigma**2)


def radios(eqf, nj, nq, cache):
    if os.path.exists(cache):
        d = np.load(cache)
        if d['r'].shape == (nj, nq):
            return d['r']
    eq = np.load(eqf)
    _, jmx, l0 = params_eq(eqf)
    m = MapaAA(eq['r'], eq['phi_self'], L=l0)
    Jn = (np.arange(nj) + 0.5)*jmx/nj
    Qn = (np.arange(nq) + 0.5)*2*np.pi/nq
    JJ, QQ = np.meshgrid(Jn, Qn, indexing='ij')
    r = np.empty(JJ.size); p = np.empty(JJ.size)
    t0 = time.time()
    for s in range(0, JJ.size, 4000):
        r[s:s+4000], p[s:s+4000] = invertir(m, eq['E_t'], eq['J_t'],
                                            QQ.ravel()[s:s+4000], JJ.ravel()[s:s+4000])
    print(f'  radios de {JJ.size} nodos en {time.time()-t0:.0f} s')
    r = r.reshape(nj, nq)
    np.savez(cache, r=r)
    return r


def resolver(eqf, nj=1600, nq=32, dt=0.5, tmax=3000.0, libre=False, cada=2.0, verboso=True):
    eq = np.load(eqf)
    A = float(eq['A'])
    sigma, jmx, l0 = params_eq(eqf)
    Jn = (np.arange(nj) + 0.5)*jmx/nj
    Qn = (np.arange(nq) + 0.5)*2*np.pi/nq
    dJ, dQ = jmx/nj, 2*np.pi/nq
    cache = os.path.splitext(eqf)[0] + f'_radios_{nj}x{nq}.npz'
    r = radios(eqf, nj, nq, cache).ravel()
    orden = np.argsort(r)
    rs = r[orden]
    om = omega_de_J(eq['E_t'], eq['J_t'])(Jn)
    Feq = A*F_forma(Jn, sigma)
    dFeq = A*dF_forma(Jn, sigma)
    k = np.fft.fftfreq(nq, 1.0/nq)                 # enteros 0..nq/2-1, -nq/2..-1
    fase = np.exp(-1j*np.outer(om, k)*dt)          # transporte exacto en Q
    peso = 8*np.pi**2*l0*dQ*dJ

    def dPhi(dF):
        mS = (dF.ravel()*peso)[orden]
        Mmenor = np.cumsum(mS) - 0.5*mS
        exterior = np.cumsum((mS/rs)[::-1])[::-1] - 0.5*mS/rs
        phi = np.empty_like(r)
        phi[orden] = -Mmenor/rs - exterior
        return phi.reshape(nj, nq)

    def golpe(dF, h):
        if libre:
            return dF
        dphidQ = np.real(np.fft.ifft(1j*k[None, :]*np.fft.fft(dPhi(dF), axis=1), axis=1))
        return dF + h*dFeq[:, None]*dphidQ

    dF = Feq[:, None]*np.cos(Qn)[None, :]
    norma = np.sum(Feq*B(Jn))*nq
    eQ = np.exp(-1j*Qn)
    pasos = int(round(tmax/dt)); nsal = int(round(cada/dt))
    ts, h1, h2 = [], [], []
    rmed = np.linspace(3, 15, 121)
    dphi_r = []
    t0 = time.time()
    for n in range(pasos + 1):
        if n % nsal == 0:
            ts.append(n*dt)
            h1.append(np.sum(dF*B(Jn)[:, None]*eQ[None, :])/norma)
            h2.append(np.sum(dF*B(Jn)[:, None]*(eQ**2)[None, :])/norma)
            ph = dPhi(dF).ravel()[orden]
            dphi_r.append(np.interp(rmed, rs, ph))
        if n == pasos:
            break
        dF = golpe(dF, 0.5*dt)
        dF = np.real(np.fft.ifft(np.fft.fft(dF, axis=1)*fase, axis=1))
        dF = golpe(dF, 0.5*dt)
    if verboso:
        print(f'  {pasos} pasos en {time.time()-t0:.0f} s (NJ={nj}, NQ={nq}, dt={dt}, libre={libre})')
    return np.array(ts), np.array(h1), np.array(h2), rmed, np.array(dphi_r)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('eq'); ap.add_argument('salida')
    ap.add_argument('--nj', type=int, default=1600)
    ap.add_argument('--nq', type=int, default=32)
    ap.add_argument('--dt', type=float, default=0.5)
    ap.add_argument('--tmax', type=float, default=3000.0)
    ap.add_argument('--libre', action='store_true')
    a = ap.parse_args()
    t, h1, h2, rmed, dphi = resolver(a.eq, a.nj, a.nq, a.dt, a.tmax, a.libre)
    np.savez(a.salida, t=t, h1=h1, h2=h2, r=rmed, dphi=dphi,
             nj=a.nj, nq=a.nq, dt=a.dt, libre=a.libre)
    print('escrito', a.salida)
