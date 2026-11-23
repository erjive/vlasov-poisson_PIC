"""Las tres maneras de calcular mu que conviven en el artículo, para todas las corridas del PIC.

mu^2 = A/(Delta Omega Delta J) (Ec. 45), con la banda al 1 % del factor en J del dato. Las tres convenciones solo
difieren en A y en Delta Omega:

    corrida   A medida en la propia corrida con autogravedad: la mitad de la mayor variación de Phi_gas(t, r) a r
              fijo en la primera vuelta radial, 0 <= t <= 2 pi/Omega_min, maximizada en 3 <= r <= 15, con las
              instantáneas del PIC (cada 10 unidades de tiempo) y su retícula radial (bulto.escalas). Sección 5A,
              tablas I a III, secciones 6 y 8.
    pasiva    la misma A para la nube que se mueve sola en el potencial externo, calculada aquí con la solución
              exacta (Q = Q_0 + Omega(J) t) y el potencial exacto de las cáscaras. Es A = a M_gas, con a fijo por la
              forma del dato y el potencial. Sección 7E y tabla VIII (formas.amplitud_libre, solo isócrono y
              s_Q = s_J = 0.1; aquí para todas las familias).
    orden0    orden dominante en J (Ecs. 58 y 59): A = pi g M_gas (2 s_J)^(1/2) y Delta Omega = |Omega_0'| Delta J, es
              decir mu^2 = 0.69 Gamma. Solo está definida para el factor J^2 exp(-J^2/s_J^2).

No corre simulaciones: lee las corridas de exe/demo_eta/bulto (este repositorio) y el resumen.npz de las de
exe/rep/12_nube de VlasovPoisson_PIC_sp, y calcula la mezcla libre exacta del dato.

    python3 mu_convenciones.py convergencia   a = A/M_gas de la nube del artículo con varias retículas y pasos
    python3 mu_convenciones.py libre          A de B0 (sin autogravedad) medida como en la corrida frente a la exacta
    python3 mu_convenciones.py tabla          las tres convenciones y Gamma para todas las corridas (exe/demo_eta/mu)
"""
import os, sys
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import bulto as B
import reducido as R
import formas as F
from demo_eta import EXE

SP = os.environ.get('VP_PIC_SP', os.path.join(os.path.dirname(R.RAIZ), 'VlasovPoisson_PIC_sp'))
SAL = os.path.join(EXE, 'demo_eta', 'mu')
L0 = 2.0
RVENT = B.RVENT


def radio(Q, J, pot):
    """r(Q, J) con L = L0 en el isócrono (formas.radio_isocrono) o en Kepler (bulto.dato_kepler)."""
    if pot == 'isocrono':
        return F.radio_isocrono(Q, J, L0)
    Q = np.mod(Q, 2*np.pi)
    a = (J + L0)**2
    e = np.sqrt(np.clip(1 - L0**2/a, 0, None))
    eta = Q.copy()
    for _ in range(60):
        eta = eta - (eta - e*np.sin(eta) - Q)/(1 - e*np.cos(eta))
    return a*(1 - e*np.cos(eta))


def amplitud_pasiva(perfil=None, sj=0.1, sq=0.1, pot='isocrono', nj=400, nq=100, dt=0.25):
    """a = A/M_gas de la nube que se mueve sola en el potencial externo, con A = (G M_gas/2) Delta<1/r>: la mitad de la
    variación del potencial dentro de la nube, en 0 <= t <= 2 pi/Omega_min. Mezcla exacta en la retícula (nj, nq) recortada
    a 1e-7, para cualquier forma, anchura y potencial. En las nubes del isócrono coincide con el máximo en 3 <= r <= 15 (el
    menor pericentro es 3.4); en Kepler la ventana cortaba la meseta y daba mu 0.9 % menor."""
    Q0, J, w = R.reticula(nj, nq, sj=sj, sq=sq, perfil=perfil)
    w = w/w.sum()
    om = (J + R.c_de_L(L0, pot))**-3
    omin = R.banda_perfil(perfil, sj, pot)[0]
    v = np.array([np.sum(w/radio(Q0 + om*t, J, pot)) for t in np.arange(0.0, 2*np.pi/omin + dt, dt)])
    return 0.5*(v.max() - v.min())


def jmed(perfil=None, sj=0.1):
    """Mediana de la acción de la nube (el factor en Q no depende de J)."""
    J = np.linspace(0.0, 8*sj, 800001)
    f = R.perfil_J(J, sj, perfil)
    c = np.concatenate([[0.0], np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(J))])
    return np.interp(0.5*c[-1], c, J)


def gamma(m, perfil=None, sj=0.1, pot='isocrono'):
    """Gamma = g M_gas/(|Omega_0'| J_med^(3/2)), Ec. (57), con Omega_0' = -3 c^-4."""
    fo = R.fondo(pot, L0)
    return fo['g']*m/(3*fo['c']**-4*jmed(perfil, sj)**1.5)


def jmedia(perfil=None, sj=0.1):
    """Promedio de la acción sobre la masa de la nube; 2 s_J/sqrt(pi) para J^2 exp(-J^2/s_J^2)."""
    return R.J_promedio(sj, perfil)


def gamma_media(m, perfil=None, sj=0.1, pot='isocrono'):
    """Gamma = g M_gas/(|Omega_0'| <J>^(3/2)), la definición de la Sección V desde el 2026-11-23."""
    fo = R.fondo(pot, L0)
    return fo['g']*m/(3*fo['c']**-4*jmedia(perfil, sj)**1.5)


def mu_orden0(m, sj=0.1, pot='isocrono'):
    """Ecs. (58)-(59) con la banda al 1 % de J^2 exp(-J^2/s_J^2): A = pi g M (2 s_J)^(1/2), Delta Omega = |Omega_0'| Delta J."""
    fo = R.fondo(pot, L0)
    dj = R.banda_perfil(None, sj, pot)[2]
    return np.sqrt(np.pi*fo['g']*m*np.sqrt(2*sj)/(3*fo['c']**-4*dj**2))


_A = {}


def coef_pasivo(perfil=None, sj=0.1, sq=0.1, pot='isocrono'):
    """(a, c) con mu_pasiva = c sqrt(M_gas/M) y la banda del dato."""
    clave = (str(perfil), sj, sq, pot)
    if clave not in _A:
        a = amplitud_pasiva(perfil, sj, sq, pot)
        omin, omax, dj, _ = R.banda_perfil(perfil, sj, pot)
        _A[clave] = (a, np.sqrt(a/((omax - omin)*dj)))
    return _A[clave]


def convergencia():
    print('a = A/M_gas de la nube del artículo (s_Q = s_J = 0.1, isócrono); formas.coeficientes da 0.04967')
    for kw in (dict(), dict(nj=800), dict(nq=200), dict(dt=0.1), dict(dt=10.0)):
        print(f'  {str(kw) or "base":16} a = {amplitud_pasiva(**kw):.6f}')


def libre():
    """A/M_gas medida como en las corridas (bulto.escalas) frente a la de la nube pasiva exacta. B0 no sirve (sin
    autogravedad el PIC no guarda Phi_gas); con M_gas/M = 10^-3 (B1) la autogravedad es despreciable, y la diferencia
    con la exacta en la misma retícula (400 x 25) y con instantáneas cada 10 mide la discretización. El resto de la
    caída con la masa es el efecto de la autogravedad en la primera vuelta."""
    a = coef_pasivo()[0]
    a25t = amplitud_pasiva(nq=25, dt=10.0)
    print(f'exacta: 400 x 100 cada 0.25: {a:.5f};  400 x 25 cada 10 (como la corrida): {a25t:.5f}')
    for n in ('B1', 'B2', 'B2q', 'B015', 'B02', 'B02q', 'B025', 'B3', 'B04', 'B05', 'B07', 'B4', 'B4q', 'B5'):
        m = B.info[n]['a0']; A = B.escalas(n)[0]
        print(f'  {n:5} M_gas/M = {m:<6g} A/M = {A/m:.5f}  cociente con la exacta {A/m/a:.4f}  en mu {np.sqrt(A/m/a):.4f}')


def familias_bulto():
    """(nombre, M_gas, perfil, s_J, s_Q, potencial, sección) de las corridas con autogravedad de bulto.py."""
    out = []
    for n, a0, npc, sg, _ in B.CORRIDAS:
        if not sg:
            continue
        sj, sq = B.parametro(n, 'sr', 0.1), B.parametro(n, 'sp', 0.1)
        pot = 'kepler' if B.FONDO.get(n) == 'puntual' else 'isocrono'
        out.append((n, a0, None, sj, sq, pot, 'V'))
    return out


def corridas_sp():
    """Las corridas de L igual de VlasovPoisson_PIC_sp con su mu medido (resumen.npz)."""
    sys.path.insert(0, os.path.join(SP, 'reproducir', 'scripts'))
    base = B.BASE
    import nube_L as N                         # cambia B.BASE y añade sus corridas a B.CAMBIOS y B.FONDO
    out = []
    for n, c in N.CORRIDAS.items():
        if c['nl'] != 1 or not os.path.exists(N.ruta(n, 'resumen.npz')):
            continue
        sec = 'VI' if n.startswith('H') else 'VII.E' if n[0] in 'FG' else 'VIII' if '_s0' in n else 'B.4'
        out.append((n, float(c['m']), c['perfil'], c['sj'], 0.1, c['pot'], sec, float(N.resumen(n)['mu'])))
    B.BASE = base
    return out


def nombre_forma(perfil):
    if perfil is None:
        return 'J^2'
    t, x = perfil
    if t == 'p':
        return f'J^{x:g}'
    return f'(J-{x:g})^2' if t == 'hueco' else t


def tabla():
    os.makedirs(SAL, exist_ok=True)
    filas = []
    for n, m, pf, sj, sq, pot, sec in familias_bulto():
        filas.append((sec, n, m, pf, sj, sq, pot, B.escalas(n)[2]))
    filas += [(sec, n, m, pf, sj, sq, pot, mu) for n, m, pf, sj, sq, pot, sec, mu in corridas_sp()]
    cab = (f'{"sec.":5} {"corrida":14} {"forma":9} {"s_J":>4} {"s_Q":>4} {"pot.":8} {"M_gas/M":>8} {"Gamma":>6} | '
           f'{"corrida":>7} {"pasiva":>7} {"orden0":>7} | {"corr/pas":>8} {"mu_pas^2/Gamma":>14}')
    lineas = [cab, '-'*len(cab)]
    for sec, n, m, pf, sj, sq, pot, mu in filas:
        a, c = coef_pasivo(pf, sj, sq, pot)
        mp, g = c*np.sqrt(m), gamma(m, pf, sj, pot)
        m0 = f'{mu_orden0(m, sj, pot):7.3f}' if pf is None else f'{"-":>7}'
        lineas.append(f'{sec:5} {n:14} {nombre_forma(pf):9} {sj:4g} {sq:4g} {pot:8} {m:8.5g} {g:6.3f} | '
                      f'{mu:7.3f} {mp:7.3f} {m0} | {mu/mp:8.3f} {mp**2/g:14.3f}')
    lineas += ['', 'a = A/M_gas de la nube pasiva y c (mu_pasiva = c sqrt(M_gas/M)) por familia:']
    for (pf, sj, sq, pot), (a, c) in _A.items():
        lineas.append(f'  forma {pf:28} s_J {sj:4g} s_Q {sq:4g} {pot:8}  a = {a:.5f}  c = {c:.4f}')
    texto = '\n'.join(lineas) + '\n'
    open(os.path.join(SAL, 'mu_convenciones.txt'), 'w').write(texto)
    print(texto)


REF = ('B1', 'B2', 'B015', 'B02', 'B021', 'B022', 'B023', 'B024', 'B025', 'B3', 'B04', 'B05', 'B07', 'B4', 'B5')
FAMILIAS = {'s_J = 0.05': ('J05_0.005', 'J05_0.01', 'J05_0.02'), 's_J = 0.2': ('J20_0.02', 'J20_0.04', 'J20_0.08'),
            'Kepler': ('K_0.01', 'K_0.02', 'K_0.04', 'K_0.08'), 's_Q = 0.2': ('Q2_0.01', 'Q2_0.02', 'Q2_0.04', 'Q2_0.08'),
            's_Q = 0.4': ('Q4_0.01', 'Q4_0.02', 'Q4_0.04', 'Q4_0.08')}


def colapso():
    """Qué convención reúne mejor las familias de la sección 5: m_c de cada corrida con núcleo (m_c > 0.1) frente al de
    la familia de referencia interpolado (lineal) al mismo valor del parámetro. orden0 es proporcional a Gamma^(1/2),
    así que vale también para Gamma."""
    def x(n):
        m, sj, sq = B.info[n]['a0'], B.parametro(n, 'sr', 0.1), B.parametro(n, 'sp', 0.1)
        pot = 'kepler' if B.FONDO.get(n) == 'puntual' else 'isocrono'
        return dict(corrida=B.escalas(n)[2], pasiva=coef_pasivo(None, sj, sq, pot)[1]*np.sqrt(m), orden0=mu_orden0(m, sj, pot))
    mc = {n: float(B.nucleo(n)['mnuc']) for n in REF + sum(FAMILIAS.values(), ())}
    X = {n: x(n) for n in mc}
    conv = ('corrida', 'pasiva', 'orden0')
    print(f'{"familia":11} {"corrida":10} {"m_c":>6} | ' + ' '.join(f'{c:>15}' for c in conv))
    print(f'{"":11} {"":10} {"":>6} | ' + ' '.join(f'{"param.  desv.":>15}' for c in conv))
    for fam, ns in FAMILIAS.items():
        for n in ns:
            if mc[n] < 0.1:
                continue
            cols = []
            for c in conv:
                xr = np.array([X[r][c] for r in REF]); yr = np.array([mc[r] for r in REF]); o = np.argsort(xr)
                ref = np.interp(X[n][c], xr[o], yr[o])
                cols.append(f'{X[n][c]:6.3f} {100*(mc[n]/ref - 1):+7.1f}%')
            print(f'{fam:11} {n:10} {mc[n]:6.3f} | ' + ' '.join(f'{s:>15}' for s in cols))


if __name__ == '__main__':
    pasos = {'convergencia': convergencia, 'libre': libre, 'tabla': tabla, 'colapso': colapso}
    for p in sys.argv[1:] or ['tabla']:
        pasos[p]()
