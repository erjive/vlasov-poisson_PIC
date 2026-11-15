"""¿Aparece el núcleo de la nube en el mismo mu cuando cambia la forma del dato en J?

Todas las nubes del artículo tienen el factor J^2 exp(-J^2/s_J^2), que se anula como J^2 en la órbita circular.
Aquí se cambia esa forma (reducido.perfil_J), con s_Q = s_J = 0.1, L_0 = 2 y el isócrono:

    p0   exp(-J^2/s_J^2)               no se anula en la órbita circular (no es suave allí; solo con el modelo)
    p2   J^2 exp(-J^2/s_J^2)           la nube del artículo
    p4   J^4 exp(-J^2/s_J^2)           se anula como J^4
    h05  (J - 0.05)_+^2 exp(-J^2/s_J^2)  su soporte no llega a la órbita circular (J >= 0.05), como en Chaturvedi y Luk

Las formas se comparan a igual mu. mu^2 = A/(Delta Omega Delta J), con la banda al 1 % del factor en J
(reducido.banda_perfil) y la amplitud A del potencial de la nube en la primera vuelta radial (sección 6A del
artículo). Para masas pequeñas A es proporcional a la masa, A = a M_gas, y a se calcula aquí con la nube moviéndose
sola en el isócrono (amplitud_libre). Así mu = c sqrt(M_gas/M) con c = (a/(Delta Omega Delta J))^(1/2), y la masa de
cada forma con el mu de la nube del artículo de masa M_ref es M = M_ref (c_p2/c)^2.

    python3 formas.py masas     coeficientes y masas (las que llevan reducido.SERIES y nube_L.CORRIDAS)
    python3 formas.py tabla     resultados del modelo (y del PIC de VlasovPoisson_PIC_sp, si están)
"""
import os, sys
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import reducido as R

SJ, SQ, L0 = 0.1, 0.1, 2.0
PERFILES = {'p0': ('p', 0.0), 'p2': None, 'p4': ('p', 4.0), 'h05': ('hueco', 0.05)}
MASAS_REF = (0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15)   # modelo
MASAS_PIC = (0.02, 0.03, 0.05)                                              # PIC
JSTAR = 0.05*SJ                                  # radio (2 J_*)^(1/2) del disco del núcleo, como en el artículo


def radio_isocrono(Q, J, L=L0):
    """r(Q, J) en el isócrono con L fija: la inversión de bulto.dato_isocrono (Binney y Tremaine, 3.5.2)."""
    Q = np.mod(Q, 2*np.pi)
    E = -0.5/(J + 0.5*(L + np.sqrt(L*L + 4)))**2
    c = -0.5/E - 1.0
    e = np.sqrt(np.clip(1 - (L*L/c)*(1 + 1/c), 0, None))
    k = e*c/(c + 1)
    eta = Q.copy()
    for _ in range(60):
        eta = eta - (eta - k*np.sin(eta) - Q)/(1 - k*np.cos(eta))
    s = 2 + c*(1 - e*np.cos(eta))
    return np.sqrt(s*(s - 2))


def amplitud_libre(perfil, nj=400, nq=100, dt=0.25):
    """a = A/M_gas de la nube que se mueve sola en el isócrono: la mitad de la mayor variación de Phi_gas(t, r) a r
    fijo en 0 <= t <= 2 pi/Omega_min, maximizada en 3 <= r <= 15 (bulto.escalas)."""
    Q0, J, w = R.reticula(nj, nq, sj=SJ, sq=SQ, perfil=perfil)
    w = w/w.sum()
    om = (J + R.c_de_L(L0))**-3
    omin = R.banda_perfil(perfil, SJ)[0]
    ts = np.arange(0.0, 2*np.pi/omin + dt, dt)
    rg = np.linspace(3.0, 15.0, 241)
    P = np.empty((len(ts), len(rg)))
    for k, t in enumerate(ts):
        r = radio_isocrono(Q0 + om*t, J)
        o = np.argsort(r); rs, ws = r[o], w[o]
        mc = np.concatenate([[0.0], np.cumsum(ws)])
        pc = np.concatenate([np.cumsum((ws/rs)[::-1])[::-1], [0.0]])
        i = np.searchsorted(rs, rg)
        P[k] = -(mc[i]/rg + pc[i])
    return 0.5*(P.max(0) - P.min(0)).max()


def coeficientes():
    """Para cada forma: a, la banda, c (mu = c sqrt(M)), tau_1 y m_mix = J_* max f/int f dJ."""
    out = {}
    J = np.linspace(0.0, 8*SJ, 400001)
    for n, pf in PERFILES.items():
        omin, omax, dj, tau = R.banda_perfil(pf, SJ)
        a = amplitud_libre(pf)
        f = R.perfil_J(J, SJ, pf)
        out[n] = dict(a=a, omin=omin, omax=omax, dj=dj, tau=tau, c=np.sqrt(a/((omax - omin)*dj)),
                      mmix=JSTAR*f.max()/np.sum(0.5*(f[1:] + f[:-1])*np.diff(J)))   # trapecios (sin np.trapz)
    return out


def masa(n, mref, co=None):
    co = co or coeficientes()
    return mref*(co['p2']['c']/co[n]['c'])**2


def tfin(n, co=None):
    """t_fin = 10^4, o 30 tau_1 redondeado hacia arriba a un múltiplo de 200 si es mayor."""
    co = co or coeficientes()
    return max(1e4, 200.0*np.ceil(30*co[n]['tau']/200.0))


def nombre_modelo(n, mref):
    return f'rF_{n}_{int(round(mref*1000)):03d}'


def nombre_pic(n, mref):
    return f'F_{n}_{int(round(mref*1000)):03d}'


def masas():
    co = coeficientes()
    print('forma   a = A/M   banda de Omega     Delta J   tau_1   c (mu = c sqrt M)   m_mix')
    for n, d in co.items():
        print(f'{n:5}  {d["a"]:.5f}  [{d["omin"]:.4f}, {d["omax"]:.4f}]  {d["dj"]:.3f}  {d["tau"]:6.0f}   {d["c"]:.4f}'
              f'             {d["mmix"]:.4f}')
    print('\nmasas a igual mu (modelo; las del PIC son las de M_ref = ' + ', '.join(f'{m:g}' for m in MASAS_PIC) + '):')
    print('M_ref   mu    ' + '  '.join(f'{n:>7}' for n in PERFILES))
    for m in MASAS_REF:
        print(f'{m:5g}  {co["p2"]["c"]*np.sqrt(m):.3f}  ' + '  '.join(f'{masa(n, m, co):7.4f}' for n in PERFILES))
    print('t_fin: ' + ', '.join(f'{n} {tfin(n, co):.0f}' for n in PERFILES))
    # Comprueba que las masas de las series del modelo y de las corridas del PIC son estas.
    errores = 0
    for n in PERFILES:
        s = R.SERIES.get(f'formas_{n}', [])
        for nom, m, *_ in s:
            mref = int(nom.rsplit('_', 1)[1])/1000
            if abs(m - round(masa(n, mref, co), 5)) > 1e-12:
                print('  DISTINTA', nom, m, masa(n, mref, co)); errores += 1
    print('masas de reducido.SERIES:', 'bien' if errores == 0 else f'{errores} distintas')
    return co


def _pic(n):
    """(medidas.npz, resumen.npz) de la corrida n del PIC de VlasovPoisson_PIC_sp, o None si no está."""
    base = os.environ.get('VP_PIC_SP', os.path.join(os.path.dirname(R.RAIZ), 'VlasovPoisson_PIC_sp'))
    d = os.path.join(base, 'exe', 'rep', '12_nube', n)
    if not os.path.exists(os.path.join(d, 'resumen.npz')):
        return None
    return np.load(os.path.join(d, 'medidas.npz')), np.load(os.path.join(d, 'resumen.npz'))


def _numeros(M):
    """|h_1|/h_0 medio en 10 <= t/tau_1 <= 20, m_c (media en t >= 20 tau_1) y masa con J < 0.2 s_J en t >= 20 tau_1."""
    t, tau = M['t'], float(M['tau'])
    v = (t >= 10*tau) & (t <= 20*tau)
    return np.abs(M['h1'])[v].mean(), float(M['mnuc']), M['m02'][t >= 20*tau].mean()


def tabla():
    co = coeficientes()
    print('m_c es la masa en el disco del núcleo (radio (2 J_*)^(1/2), J_* = 0.05 s_J), m02 la masa con J < 0.2 s_J;')
    print('un gas mezclado a acción fija tiene m_c <= m_mix. mu del PIC: medido como en la sección 6A con la banda de la forma.')
    print(f'{"forma":5} {"M_ref":>5} {"mu":>6} {"M_gas":>8} | {"modelo |h1|":>11} {"m_c":>6} {"m02":>6} |'
          f' {"PIC mu":>6} {"|h1|":>8} {"m_c":>6} {"m02":>6}')
    for n in PERFILES:
        for mref in MASAS_REF:
            fila = f'{n:5} {mref:5g} {co["p2"]["c"]*np.sqrt(mref):6.3f} {masa(n, mref, co):8.5f} |'
            nm = R._ref(mref) if n == 'p2' else nombre_modelo(n, mref)
            if os.path.exists(os.path.join(R.DATOS, nm + '_diag.npz')) or os.path.exists(os.path.join(R.DATOS, nm + '.npz')):
                h1, mc, m02 = _numeros(R.diagnosticos(nm))
                fila += f' {h1:11.2e} {mc:6.3f} {m02:6.3f} |'
            else:
                fila += f' {"--":>11} {"--":>6} {"--":>6} |'
            p = _pic(nombre_pic(n, mref)) if mref in MASAS_PIC else None
            if p is not None:
                h1, mc, m02 = _numeros(p[0])
                fila += f' {float(p[1]["mu"]):6.3f} {h1:8.2e} {mc:6.3f} {m02:6.3f}'
            print(fila)
        print(f'      m_mix = {co[n]["mmix"]:.4f}, tau_1 = {co[n]["tau"]:.0f}')


if __name__ == '__main__':
    paso = sys.argv[1] if len(sys.argv) > 1 else 'masas'
    if paso == 'masas':
        masas()
    elif paso == 'tabla':
        tabla()
    else:
        sys.exit(__doc__)
