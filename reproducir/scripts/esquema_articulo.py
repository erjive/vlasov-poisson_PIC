"""Implementación independiente del cálculo del campo de la Sección 4 del artículo
(docs/articulo/sec4_numerics.tex, pasos (i)-(iv)), escrita a partir del texto y no del código
Fortran, para comprobar tres propiedades que el texto afirma:

  1. el campo ve exactamente la masa asignada, sea cual sea el radio de la capa;
  2. la fuerza de una capa sobre sí misma es -m/(2 r^2) salvo una corrección de orden dr/r;
  3. para una esfera homogénea la fuerza exterior es exacta y la interior converge con el
     número de capas.

    python3 esquema_articulo.py

Tarda unos segundos y no escribe archivos.
"""
import numpy as np

def campo(rp, mp, dr, Nr):
    R = (np.arange(1, Nr + 1) - 0.5)*dr
    # (i) asignación con pesos lineales (nube en celda)
    mu = np.zeros(Nr)
    u = rp/dr - 0.5; i0 = np.floor(u).astype(int); f = u - i0
    for idx, w in ((i0, 1 - f), (i0 + 1, f)):
        ok = (idx >= 0) & (idx < Nr)
        np.add.at(mu, idx[ok], (mp*w)[ok])
    # (ii) densidad lineal a trozos
    rho = mu/(4*np.pi*dr*(R**2 + dr**2/6))
    # (iii) masa encerrada y potencial, en forma cerrada en cada intervalo
    M = np.zeros(Nr); pot = np.zeros(Nr)
    M[0] = 4*np.pi/3*rho[0]*R[0]**3; pot[0] = 2*np.pi/3*rho[0]*R[0]**2
    for i in range(1, Nr):
        a, b = R[i-1], R[i]; s = (rho[i] - rho[i-1])/dr
        c3 = 4*np.pi*(rho[i-1] - s*a)/3; c4 = np.pi*s; c0 = M[i-1] - c3*a**3 - c4*a**4
        M[i] = c0 + c3*b**3 + c4*b**4
        pot[i] = pot[i-1] + (-c0/b + c3*b*b/2 + c4*b**3/3) - (-c0/a + c3*a*a/2 + c4*a**3/3)
    dpot = M/R**2
    pot -= pot[-1] + R[-1]*dpot[-1]
    return R, pot, dpot, M

def interp(R, g, r, dr):
    u = r/dr - 0.5; i0 = np.floor(u).astype(int); f = u - i0
    return g[i0]*(1 - f) + g[i0 + 1]*f

print('1) masa que ve el campo, una sola capa de masa 1 (dr = 0.1, malla hasta 20):')
for r0 in (0.55, 1.0, 2.03, 5.0, 12.34):
    R, pot, dpot, M = campo(np.array([r0]), np.array([1.0]), 0.1, 200)
    print(f'   r0 = {r0:6.2f}: M(r_max) = {M[-1]:.12f},  Phi(r_max) r_max = {pot[-1]*R[-1]:+.12f}')

print('2) fuerza de una capa sobre sí misma, frente a -m/(2 r^2), con r0 = 1:')
for dr in (0.2, 0.1, 0.05, 0.025, 0.0125):
    Nr = int(round(4/dr))
    R, pot, dpot, M = campo(np.array([1.0 + 0.3*dr]), np.array([1.0]), dr, Nr)
    r0 = 1.0 + 0.3*dr
    F = -interp(R, dpot, np.array([r0]), dr)[0]
    Fa = -1/(2*r0**2)
    print(f'   dr = {dr:7.4f}: F/F_auto = {F/Fa:.6f},  (1 - F/F_auto)/(dr/r) = {(1 - F/Fa)/(dr/r0):.3f}')

print('3) esfera homogénea de masa 1 y radio 5 con N capas (dr = 0.05): error relativo de la fuerza')
for N in (500, 2000, 8000):
    k = np.arange(1, N + 1)
    redge = 5*(k/N)**(1/3); rin = 5*((k - 1)/N)**(1/3)
    rp = (0.5*(redge**3 + rin**3))**(1/3)            # capas de igual masa
    R, pot, dpot, M = campo(rp, np.full(N, 1.0/N), 0.05, 200)
    dentro = (R > 0.5) & (R < 4.5); fuera = R > 5.2
    e1 = np.max(np.abs(dpot[dentro]/(R[dentro]/125) - 1)); e2 = np.max(np.abs(dpot[fuera]*R[fuera]**2 - 1))
    print(f'   N = {N:5d}: dentro (0.5 < r < 4.5) {e1:.2e},  fuera (r > 5.2) {e2:.2e}')
