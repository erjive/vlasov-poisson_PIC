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
    python3 formas.py colapso   cruces de m_c y |h_1| como función de mu y de Gamma (con la mediana de J)
    python3 formas.py cruces    lo mismo para todas las formas (las dos tandas)
    python3 formas.py figura    figura del artículo (docs/articulo/figuras/formas_gamma.pdf)
    python3 formas.py articulo  tabla (docs/articulo/tab_formas.tex) y figura del artículo
    python3 formas.py familias  Gamma y m_c del modelo para las familias de la sección 6 (s_J, Kepler, s_Q)
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


_COEF = {}


def coeficientes():
    """Para cada forma: a, la banda, c (mu = c sqrt(M)), tau_1 y m_mix = J_* max f/int f dJ (se calcula una vez)."""
    if _COEF:
        return _COEF
    out = _COEF
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


OMEGA_P = -3.0*R.c_de_L(L0)**-4                 # dOmega/dJ en la órbita circular del isócrono con L_0 = 2 (-0.08831)


def J_mediana(perfil):
    """Mediana de la acción del dato (pesada con la masa) en la retícula de 400 x 100: la distancia típica de la nube
    a la órbita circular en el plano de zeta es (2 J_med)^(1/2)."""
    Q, J, w = R.reticula(400, 100, sj=SJ, sq=SQ, perfil=perfil)
    o = np.argsort(J); c = np.cumsum(w[o])/w.sum()
    return J[o][np.searchsorted(c, 0.5)]


def gamma(n, M):
    """Gamma = g M_gas/(|Omega_0'| <J>^(3/2)): el parámetro del modelo epicíclico con el promedio de la acción."""
    return R.G*M/(abs(OMEGA_P)*R.J_promedio(SJ, PERFILES[n])**1.5)


def _cruce(x, y, nivel):
    o = np.argsort(x); x, y = x[o], np.maximum.accumulate(y[o])
    for i in range(1, len(x)):
        if y[i - 1] < nivel <= y[i]:
            return np.exp(np.log(x[i - 1]) + (nivel - y[i - 1])/(y[i] - y[i - 1])*(np.log(x[i]) - np.log(x[i - 1])))
    return np.nan


def colapso():
    """Dónde cruzan m_c y |h_1|/h_0 ciertos niveles en las series del modelo, como función de mu y de Gamma."""
    co = coeficientes()
    datos = {}
    for n in PERFILES:
        fil = []
        for mref in MASAS_REF:
            nm = R._ref(mref) if n == 'p2' else nombre_modelo(n, mref)
            if os.path.exists(os.path.join(R.DATOS, nm + '_diag.npz')) or os.path.exists(os.path.join(R.DATOS, nm + '.npz')):
                h1, mc, _ = _numeros(R.diagnosticos(nm))
                M = masa(n, mref, co)
                fil.append((co[n]['c']*np.sqrt(M), gamma(n, M), mc, h1))
        datos[n] = np.array(fil)
    print('J_med: ' + ', '.join(f'{n} {J_mediana(pf):.4f}' for n, pf in PERFILES.items()))
    print('cruces (p0 ya está por encima en todo el barrido); cociente entre el mayor y el menor de p2, p4 y h05:')
    for col, nivel, txt in ((2, 0.1, 'm_c = 0.1'), (2, 0.2, 'm_c = 0.2'), (2, 0.3, 'm_c = 0.3'), (3, 0.1, '|h1| = 0.1')):
        for k, nom in ((0, 'mu'), (1, 'Gamma')):
            v = {n: _cruce(d[:, k], d[:, col], nivel) for n, d in datos.items()}
            r = np.array([v['p2'], v['p4'], v['h05']])
            print(f'  {txt:10} {nom:5}: ' + '  '.join(f'{n} {x:6.3f}' for n, x in v.items()) + f'   {r.max()/r.min():.2f}')
    print('corridas del PIC: Gamma, m_c, |h1|/h0')
    for n in ('p2', 'p4', 'h05'):
        for mref in MASAS_PIC:
            p = _pic(nombre_pic(n, mref))
            if p is not None:
                h1, mc, _ = _numeros(p[0])
                print(f'  {nombre_pic(n, mref):10} Gamma = {gamma(n, masa(n, mref, co)):.3f}  m_c = {mc:.3f}  |h1| = {h1:.2e}')


# Todas las formas (primera y segunda tanda): perfil, potencial, filas de la retícula, color y marcador (paleta
# categórica validada de la guía de gráficas, en este orden) y rótulo de la figura.
TODAS = {
    'p2':     (None, 'isocrono', 400, '#2a78d6', 'o', r'$J^2e^{-J^2/\sigma_J^2}$ (sections 6-8)'),
    'p4':     (('p', 4.0), 'isocrono', 400, '#eb6834', 's', r'$J^4e^{-J^2/\sigma_J^2}$'),
    'h05':    (('hueco', 0.05), 'isocrono', 400, '#1baf7a', '^', r'$(J-0.05)_+^2e^{-J^2/\sigma_J^2}$'),
    'anillo': (('anillo', (0.15, 0.03)), 'isocrono', 800, '#eda100', 'D', r'ring, $e^{-(J-0.15)^2/0.03^2}$'),
    'meseta': (('meseta', (0.05, 0.25)), 'isocrono', 400, '#e87ba4', 'v', r'water bag, $0.05\leq J\leq 0.25$'),
    'p0':     (('p', 0.0), 'isocrono', 400, '#008300', 'P', r'$e^{-J^2/\sigma_J^2}$ (not smooth at $J=0$)'),
    'h05K':   (('hueco', 0.05), 'kepler', 400, '#6250d6', 'X', r'$(J-0.05)_+^2e^{-J^2/\sigma_J^2}$, Kepler'),
}
_CACHE = {}


def escalas(n):
    """<J>, k (M = k Gamma) y kb (M = kb Gamma_b) de la forma n. Gamma_b = g M <J>^(1/2)/(|Omega_0'| Delta J^2) es Gamma con
    el ancho de la banda en lugar de la distancia <J> en el denominador."""
    if n not in _CACHE:
        pf, pot, nj = TODAS[n][:3]
        Q, J, w = R.reticula(nj, 100, sj=SJ, sq=SQ, perfil=pf)
        jm = R.J_promedio(SJ, pf)
        fo = R.fondo(pot); k = 3*fo['c']**-4*jm**1.5/fo['g']
        kb = 3*fo['c']**-4*R.banda_perfil(pf, SJ, pot)[2]**2/(fo['g']*jm**0.5)
        _CACHE[n] = (jm, k, kb)
    return _CACHE[n]


def corridas_modelo(n):
    """(nombre, M) de las corridas del modelo de la forma n."""
    if n == 'p2':
        s = R.SERIES['referencia'] + R.SERIES['formas_p2']
    else:
        s = R.SERIES.get(f'formas_{n}', []) + R.SERIES.get(f'formas2_{n}', [])
    return [(nom, m) for nom, m, *_ in s]


def corridas_pic(n):
    """(nombre, M) de las corridas del PIC (VlasovPoisson_PIC_sp, tandas formas y formas2) de la forma n."""
    out = []
    if n in ('p2', 'p4', 'h05'):
        out += [(nombre_pic(n, mr), masa(n, mr)) for mr in MASAS_PIC]
    for nom in ('G_anillo_150', 'G_anillo_300', 'G_anillo_450', 'G_meseta_150', 'G_meseta_450', 'G_p0_150'):
        f, gm = nom[2:].rsplit('_', 1)
        if f == n:
            out.append((nom, R.FORMAS2_MASAS[n][int(gm)/1000]))
    return out


def datos(n):
    """Arreglos (Gamma, Gamma_b, m_c, |h1|/h0) del modelo y del PIC de la forma n, ordenados en Gamma."""
    jm, k, kb = escalas(n)
    res = []
    for fuente in ('modelo', 'pic'):
        fil = []
        for nom, M in (corridas_modelo(n) if fuente == 'modelo' else corridas_pic(n)):
            if fuente == 'modelo':
                if not (os.path.exists(os.path.join(R.DATOS, nom + '_diag.npz')) or os.path.exists(os.path.join(R.DATOS, nom + '.npz'))):
                    continue
                h1, mc, _ = _numeros(R.diagnosticos(nom))
            else:
                p = _pic(nom)
                if p is None:
                    continue
                h1, mc, _ = _numeros(p[0])
            fil.append((M/k, M/kb, mc, h1))
        a = np.array(fil).reshape(-1, 4)
        res.append(a[np.argsort(a[:, 0])])
    return res


def figura():
    """docs/articulo/figuras/formas_gamma.pdf: m_c frente a Gamma_b y frente a Gamma, y |h_1|/h_0 frente a Gamma, para todas
    las formas; modelo en líneas, PIC en marcadores."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.6))
    paneles = ((1, 2, r'core mass $m_c$ against $\Gamma_b$', r'$\Gamma_b$', False),
               (0, 2, r'core mass $m_c$ against $\Gamma$', r'$\Gamma$', False),
               (0, 3, r'$|h_1|/h_0$, mean over $10\leq t/\tau_1\leq 20$', r'$\Gamma$', True))
    for a, (ix, iy, titulo, xl, log) in zip(ax, paneles):
        for n, (pf, pot, nj, col, mk, rot) in TODAS.items():
            mod, pic = datos(n)
            if len(mod):
                a.plot(mod[:, ix], mod[:, iy], color=col, lw=1.3, marker=mk, ms=3.2, mec=sup, mew=0.4, zorder=2)
            if len(pic):
                a.plot(pic[:, ix], pic[:, iy], ls='none', color=col, marker=mk, ms=7, mec=sup, mew=0.9, zorder=3)
        a.set_xscale('log')
        marcas = (0.02, 0.05, 0.1, 0.2, 0.5, 1) if ix == 1 else (0.05, 0.1, 0.2, 0.5, 1, 2)
        a.xaxis.set_major_locator(FixedLocator(marcas)); a.xaxis.set_minor_formatter(NullFormatter())
        a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
        if log:
            a.set_yscale('log')
        else:
            a.set_ylim(0, 0.8)
        a.grid(True, which='major', color=reja, lw=0.6); a.set_axisbelow(True)
        for lado in ('top', 'right'):
            a.spines[lado].set_visible(False)
        a.tick_params(length=3, width=0.6); a.tick_params(which='minor', length=2, width=0.4)
        a.set_title(titulo, loc='left', fontsize=9.5, color=tinta, pad=6)
        a.set_xlabel(xl)
    ax[0].set_xlim(0.008, 2.0); ax[1].set_xlim(0.04, 2.0); ax[2].set_xlim(0.04, 2.0)
    h = [Line2D([], [], color=c, lw=1.3, marker=mk, ms=5, mec=sup, mew=0.6, label=r) for (pf, pot, nj, c, mk, r) in TODAS.values()]
    h += [Line2D([], [], color=tinta2, lw=1.3, label='epicyclic model'),
          Line2D([], [], color=tinta2, ls='none', marker='o', ms=7, mec=sup, label='PIC')]
    fig.tight_layout(rect=(0, 0, 1, 0.80), w_pad=1.8)
    fig.legend(handles=h, loc='upper center', ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.0),
               handlelength=2.4, columnspacing=1.4)
    sal = os.path.join(R.RAIZ, 'docs', 'articulo', 'figuras', 'formas_gamma.pdf')
    fig.savefig(sal); fig.savefig(os.path.join(R.DATOS, 'formas_gamma.png'), dpi=110); plt.close(fig)
    print('escrita', sal)


def cruces():
    """Cruces de m_c = 0.1 y 0.3 y de |h1|/h0 = 0.1 en mu y en Gamma, por forma (modelo)."""
    for col, nivel, txt in ((2, 0.1, 'm_c = 0.1'), (2, 0.3, 'm_c = 0.3'), (3, 0.1, '|h1| = 0.1')):
        for k, nom in ((1, 'Gamma_b'), (0, 'Gamma')):
            v = {}
            for n in TODAS:
                mod = datos(n)[0]
                if len(mod) and np.isfinite(mod[:, k]).all():
                    v[n] = _cruce(mod[:, k], mod[:, col], nivel)
            print(f'  {txt:10} {nom:5}: ' + '  '.join(f'{n} {x:6.3f}' for n, x in v.items()))
    v = [_cruce(datos(n)[0][:, 0], datos(n)[0][:, 2], 0.1) for n in TODAS if n != 'p0']
    gc = float(np.mean(v))
    print(f'  m_c = 0.1 con las seis formas que se anulan en la órbita circular: Gamma = {min(v):.3f}-{max(v):.3f}, media {gc:.3f}')
    for pot in ('isocrono', 'kepler'):
        fo = R.fondo(pot)
        print(f'  masa umbral M_gas/M = Gamma_c |Omega_0\'|/g <J>^(3/2) = {gc*3*fo["c"]**-4/fo["g"]:.3f} <J>^(3/2) ({pot})')


def familias():
    """Gamma (con J_med = 1.088 s_J) y m_c del modelo para las familias de la sección 6 (series referencia y familias)."""
    jr = J_mediana(None)
    for nom, m, nj, nq, dt, tf, op, pic in R.SERIES['referencia'] + R.SERIES['familias']:
        if not os.path.exists(os.path.join(R.DATOS, nom + '_diag.npz')):
            continue
        sj, sq, pot = op.get('sj', 0.1), op.get('sq', 0.1), op.get('pot', 'isocrono')
        fo = R.fondo(pot)
        g = fo['g']*m/(3*fo['c']**-4*(jr*sj/0.1)**1.5)
        h1, mc, _ = _numeros(R.diagnosticos(nom))
        print(f'  {nom:9} s_J = {sj:4.2f}  s_Q = {sq:3.1f}  {pot:8}  M = {m:6.3f}  Gamma = {g:6.3f}  m_c = {mc:.3f}  |h1| = {h1:.1e}')

ROTULO_TABLA = {'p2': r'$J^2e^{-J^2/\sigma_J^2}$', 'p4': r'$J^4e^{-J^2/\sigma_J^2}$', 'h05': r'$(J-0.05)_+^2e^{-J^2/\sigma_J^2}$',
                'anillo': 'ring', 'meseta': 'water bag', 'p0': r'$e^{-J^2/\sigma_J^2}$'}


def _central(M, n):
    """F_c/Fhat_max de la forma n: densidad en J < 0.025 s_J (t >= 20 tau_1, una instantánea de cada cinco) en unidades del
    máximo de la masa por unidad de J de esa forma. rho_c viene en unidades del de J^2 exp(-J^2/s_J^2), cuyo m_mix es 0.2/(e sqrt(pi))."""
    t, tau = M['t'], float(M['tau'])
    return M['rho_c'][np.where(t >= 20*tau)[0][::5]].mean()*(0.2/(np.e*np.sqrt(np.pi)))/(JSTAR*_mmix(n))


def tabla_articulo():
    """docs/articulo/tab_formas.tex (macro TabFormas): las corridas del PIC de las dos tandas, con el modelo en la misma
    masa. Columnas: forma, M_gas/M, Gamma_b, Gamma, m_mix, F_c/Fhat_max (PIC, modelo), m_c (PIC, modelo), |h1|/h0 (PIC, modelo)."""
    filas = []
    for n in ('p2', 'p4', 'h05', 'anillo', 'meseta', 'p0'):
        jm, k, kb = escalas(n)
        mmix = JSTAR*_mmix(n)
        pic = corridas_pic(n)
        for i, (nom, M) in enumerate(pic):
            p = _pic(nom)
            if p is None:
                continue
            h1p, mcp, _ = _numeros(p[0])
            if nom.startswith('F_'):
                mref = int(nom.rsplit('_', 1)[1])/1000
                nm = R._ref(mref) if n == 'p2' else nombre_modelo(n, mref)
            else:
                nm = 'r' + nom
            h1m, mcm, _ = _numeros(R.diagnosticos(nm))
            Fp, Fm = _central(p[0], n), _central(R.diagnosticos(nm), n)
            filas.append(f'{ROTULO_TABLA[n] if i == 0 else ""} & ${M:.3g}$ & ${M/kb:.3f}$ & ${M/k:.2f}$ & '
                         + (f'${mmix:.3f}$' if i == 0 else '') +
                         f' & ${Fp:.2f}$ & ${Fm:.2f}$ & ${mcp:.3f}$ & ${mcm:.3f}$ & {R._sci(h1p)} & {R._sci(h1m)} \\\\')
        filas.append(r'\hline')
    filas = filas[:-1]
    R._macro('TabFormas', 'tab_formas.tex', filas, destino=os.path.join(R.RAIZ, 'docs', 'articulo'))
    print('escrita docs/articulo/tab_formas.tex con', sum(1 for f in filas if f != r'\hline'), 'filas')


def _mmix(n):
    """max f/int f dJ de la forma n (m_mix = J_* por esto)."""
    J = np.linspace(0.0, 8*SJ, 400001)
    f = R.perfil_J(J, SJ, TODAS[n][0])
    return f.max()/np.sum(0.5*(f[1:] + f[:-1])*np.diff(J))

if __name__ == '__main__':
    paso = sys.argv[1] if len(sys.argv) > 1 else 'masas'
    if paso == 'masas':
        masas()
    elif paso == 'tabla':
        tabla()
    elif paso == 'colapso':
        colapso()
    elif paso == 'figura':
        figura()
    elif paso == 'cruces':
        cruces()
    elif paso == 'articulo':
        tabla_articulo(); figura()
    elif paso == 'familias':
        familias()
    else:
        sys.exit(__doc__)
