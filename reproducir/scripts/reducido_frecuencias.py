"""Pruebas estáticas del modelo reducido contra las corridas del PIC (tabla de frecuencias de las notas
docs/modelo_reducido). No corre nada: usa los mapas (Q, J) de bulto.py.

    python3 reducido_frecuencias.py B2 B02 B02q

Para cada corrida:
  sonda    frecuencia de cada partícula (pendiente del ángulo desenrollado entre 20 y 30 tau_1) y acción media;
  perfil   esa frecuencia frente a Omega(J) + (g/rho) dpsi/drho del estado mezclado (anillos con la acción de
           cada partícula; integral elíptica E), y la diferencia frente al corrimiento uniforme delta Omega_u;
  nucleo   frecuencia de las pequeñas oscilaciones en el potencial promediado en el tiempo, con todo el gas
           (órbitas casi circulares y partículas del núcleo) y solo con el halo suavizado con el ancho del núcleo
           (oscilación rígida del núcleo), frente a la frecuencia medida del radio medio del núcleo;
  hueco    masa por unidad de acción y frecuencia media cerca de la resonancia con la oscilación del núcleo.
Los resultados intermedios quedan en exe/demo_eta/reducido/sondas/.
"""
import os, sys, warnings
import numpy as np
from scipy.special import ellipe
from scipy.optimize import brentq
from scipy.ndimage import gaussian_filter1d
warnings.filterwarnings('ignore')
import reducido as R
import bulto as B

SONDAS = os.path.join(R.DATOS, 'sondas')
F0 = R.fondo()                                   # isócrono, L = 2
C, KAP, RC, G = F0['c'], F0['kappa'], F0['rc'], F0['g']
dU = lambda r: -R.LFIX**2/r**3 + r/(np.sqrt(1 + r*r)*(1 + np.sqrt(1 + r*r))**2)   # derivada del potencial efectivo
_h = 1e-3
D2U = (dU(RC + _h) - dU(RC - _h))/(2*_h); D3U = (dU(RC + _h) - 2*dU(RC) + dU(RC - _h))/_h**2
COEF_U = -(D3U/(RC**2*KAP**2) + 2/RC**3)/(2*KAP)  # corrimiento uniforme por unidad de masa encerrada


def sonda(n):
    sal = os.path.join(SONDAS, f'sonda_{n}.npz')
    if os.path.exists(sal):
        return np.load(sal)
    os.makedirs(SONDAS, exist_ok=True)
    d = np.load(B.ruta(n, 'bulto.npz'))
    t, w = d['t'], d['w'].astype(float); tau = B.banda(n)[3]
    Q, J = d['Q'].astype(float), d['J'].astype(float)
    ok = np.isfinite(J) & np.isfinite(Q)
    v = t >= 20*tau
    Qu = np.unwrap(np.where(ok, Q, 0.0)[v], axis=0); tt = t[v]
    Om = np.polyfit(tt - tt.mean(), Qu, 1)[0]
    np.savez(sal, t=t, tau=tau, w=w, Om=Om, Jf=np.where(ok, J, 0.0)[v].mean(0))
    return np.load(sal)


def perfil(n):
    s = sonda(n); a0 = B.info[n]['a0']; sj = B.parametro(n, 'sr', 0.1)
    m = a0*s['w']/s['w'].sum(); Om, Jf = s['Om'], s['Jf']
    ok = np.isfinite(Om) & np.isfinite(Jf) & (m > 0)
    m, Om, Jf = m[ok], Om[ok], Jf[ok]; rj = np.sqrt(2*Jf)
    rho = np.concatenate([np.linspace(1e-4, 0.2, 400), np.linspace(0.2, 1.2, 400)[1:]])
    psi = np.zeros_like(rho)
    for a, mj in zip(rj, m):
        psi += mj*(2/np.pi)*(rho + a)*ellipe(4*rho*a/(rho + a)**2)
    pred = (Jf + C)**-3 + np.interp(rj, rho, G*np.gradient(psi, rho)/rho)
    print(f'{n}: mu = {B.escalas(n)[2]:.2f}, M_gas = {a0}')
    print('   J/sigma_J        masa    medida    Omega + g psi\'/rho   diferencia')
    cortes = [0, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0]
    for a, b in zip(cortes[:-1], cortes[1:]):
        v = (Jf >= a*sj) & (Jf < b*sj)
        if m[v].sum() > 0:
            om, pr = np.average(Om[v], weights=m[v]), np.average(pred[v], weights=m[v])
            print(f'   {a:6.4f}-{b:6.4f} {m[v].sum()/a0:7.4f}   {om:.5f}   {pr:.5f}             {om - pr:+.2e}')


def nucleo(n, corte=0.0725):
    """corte: las partículas con frecuencia mayor forman el núcleo (el resultado cambia en 3e-5 entre 0.0722 y 0.0728)."""
    s = sonda(n); a0 = B.info[n]['a0']; tau = float(s['tau'])
    d = np.load(B.ruta(n, 'bulto.npz')); t = d['t']; m = a0*s['w']/s['w'].sum()
    v = t >= 20*tau; Rr = d['R'][v].astype(float)
    nu = B.nucleo(n); x = nu['rc'] - nu['rc'][t > 2*tau].mean(); dt = t[1] - t[0]
    F = np.fft.fft(x); fr = 2*np.pi*np.fft.fftfreq(len(x), dt)
    F[(np.abs(fr) < 0.05) | (np.abs(fr) > 0.10)] = 0; F[fr < 0] = 0
    ph = np.unwrap(np.angle(2*np.fft.ifft(F)))
    frec = lambda a, b: np.polyfit(t[(t >= a*tau) & (t <= b*tau)], ph[(t >= a*tau) & (t <= b*tau)], 1)[0]
    ancho = nu['dr'][v].mean()
    nuc = np.isfinite(s['Om']) & (s['Om'] >= corte)
    dr = 0.005; rg = np.arange(2.0, 14.0, dr) + 0.5*dr
    out = {}
    for nombre, sel, mself, sig in (('halo', ~nuc, m[nuc].sum(), ancho), ('total', np.ones(len(m), bool), 0.0, 0.0)):
        rho = np.histogram(Rr[:, sel].ravel(), bins=np.arange(2.0, 14.0 + dr/2, dr), weights=np.tile(m[sel], Rr.shape[0]))[0]/Rr.shape[0]/dr
        if sig > 0:
            rho = gaussian_filter1d(rho, sig/dr)
        Menc = np.cumsum(rho)*dr - 0.5*rho*dr
        Ff = -dU(rg) - (Menc + 0.5*mself)/rg**2
        i = np.where((Ff[:-1] > 0) & (Ff[1:] <= 0))[0]; i = i[np.argmin(np.abs(rg[i] - RC))]
        r0 = rg[i] + dr*Ff[i]/(Ff[i] - Ff[i + 1])
        k2 = D2U + D3U*(r0 - RC) + np.interp(r0, rg, rho)/r0**2 - 2*np.interp(r0, rg, Menc + 0.5*mself)/r0**3
        out[nombre] = (np.sqrt(k2), np.interp(r0, rg, Menc)/a0)
    print(f'{n}: frecuencia medida del radio medio del núcleo: {frec(2, 14):.5f} (2-14 tau_1), {frec(22, 30):.5f} (22-30); '
          f'ancho {ancho:.3f}; masa con frecuencia >= {corte}: {m[nuc].sum()/a0:.3f}')
    print(f'   campo medio total: {out["total"][0]:.5f} (masa encerrada {out["total"][1]:.3f} M_gas, corrimiento uniforme '
          f'{COEF_U*out["total"][1]*a0:.2e});  solo el halo: {out["halo"][0]:.5f}')


def hueco(n):
    s = sonda(n); sj = B.parametro(n, 'sr', 0.1)
    w, Om, Jf = s['w'], s['Om'], s['Jf']/sj
    ok = np.isfinite(Om) & np.isfinite(Jf); m = w/w[ok].sum()
    e = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.17, 0.20, 0.25, 0.30, 0.40, 0.55, 0.70, 1.0])
    print(f'{n}: J/sigma_J, masa por unidad de J/sigma_J, frecuencia media, partículas')
    for a, b in zip(e[:-1], e[1:]):
        v = ok & (Jf >= a) & (Jf < b)
        if v.sum():
            print(f'   {a:5.2f}-{b:5.2f}  {m[v].sum()/(b - a):7.3f}   {np.average(Om[v], weights=m[v]):.5f}  {v.sum():5d}')


if __name__ == '__main__':
    print(f'R_0 = {RC:.4f}, kappa_0 = {KAP:.5f}, g = {G:.5f}, A/M_gas = {np.pi*G*np.sqrt(0.2):.4f}, '
          f'corrimiento uniforme por unidad de masa encerrada = {COEF_U:.4f}')
    for n in sys.argv[1:]:
        perfil(n); nucleo(n); hueco(n)
