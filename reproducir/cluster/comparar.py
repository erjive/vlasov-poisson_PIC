"""Compara corridas hechas en el cluster con las mismas corridas hechas en el portátil.

Los dos ejecutables salen de compiladores y procesadores distintos, así que los resultados
coinciden al redondeo solo al principio; después la diferencia crece. Este guion mide cuánto.

    python3 reproducir/cluster/comparar.py [--cluster exe/cluster/hadzic] [--local exe/hadzic]
                                           [--eps 0.03] <D> [<Z>]

Para cada corrida:
  - que el archivo de parámetros sea el mismo (suma del .meta) y el paso de tiempo del registro;
  - las series que escribe el código (hk1_complex.tl: h_n con el mapa del isócrono, cada
    instantánea, con todas las cifras): diferencia relativa por ventanas y el primer tiempo en
    que pasa de 1e-12, 1e-9, 1e-6 y 1e-3;
  - si las dos tienen serie.npz: lo mismo para h_1 en el mapa del equilibrio, y la energía.
Con dos corridas (la perturbada y su referencia) y --eps, además la respuesta
chi_1 = (h_1 de D - h_1 de Z)/eps de cada máquina: diferencia por ventanas y cociente complejo
de amplitudes (hadzic.ganancia) de la del cluster respecto de la del portátil.

Las corridas del cluster se traen con reproducir/cluster/traer.sh.
"""
import argparse, os, sys
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
sys.path.insert(0, os.path.join(RAIZ, 'reproducir', 'scripts'))
VENTANAS = ((0, 500), (500, 1000), (1000, 2000), (2000, 4000))
UMBRALES = (1e-12, 1e-9, 1e-6, 1e-3)


def meta(base, nombre):
    """El .meta de la corrida como diccionario (vacío si no existe)."""
    p = os.path.join(base, nombre + '.meta')
    if not os.path.exists(p):
        return {}
    return dict((l.split(None, 1) + [''])[:2] for l in open(p).read().splitlines() if l.strip())


def paso(base, nombre):
    p = os.path.join(base, nombre + '.log')
    if os.path.exists(p):
        for l in open(p, errors='replace'):
            if 'time step fixed at size' in l.lower():
                return float(l.split()[-1])
    return float('nan')


def por_ventanas(t, a, b):
    """max|a - b| / max|b| en cada ventana, y el primer tiempo en que |a - b|/max|b| pasa de
    cada umbral (max|b| de toda la serie)."""
    d = np.abs(a - b)
    fila = []
    for lo, hi in VENTANAS:
        v = (t >= lo) & (t <= hi)
        fila.append(d[v].max()/np.abs(b[v]).max() if v.any() else float('nan'))
    rel = d/np.abs(b).max()
    cruces = [t[np.argmax(rel > u)] if np.any(rel > u) else float('inf') for u in UMBRALES]
    return fila, cruces


def linea(etiqueta, t, a, b):
    fila, cruces = por_ventanas(t, a, b)
    return (f'{etiqueta:>22} ' + ' '.join(f'{x:11.1e}' for x in fila) + '   '
            + ' '.join(f'{("nunca" if np.isinf(c) else format(c, ".0f")):>7}' for c in cruces))


def comparar(nombre, bc, bl):
    print(f'\n=== {nombre}')
    mc, ml = meta(bc, nombre), meta(bl, nombre)
    for clave in ('maquina', 'VP_PIC', 'HEAD', 'src/'):
        print(f'  {clave:<8} cluster: {mc.get(clave, "-")[:70]}\n  {"":<8} portátil: {ml.get(clave, "-")[:70]}')
    pc, pl = mc.get('par', '').split()[:1], ml.get('par', '').split()[:1]
    mismo = 'sin .meta en alguna de las dos' if not (pc and pl) else 'el mismo' if pc == pl else 'DISTINTO'
    print(f'  archivo de parámetros: {mismo};  paso de tiempo: {paso(bc, nombre):g} y {paso(bl, nombre):g}')
    cabecera = (f'{"":>22} ' + ' '.join(f'{f"[{a},{b}]":>11}' for a, b in VENTANAS) + '   '
                + ' '.join(f'{f">{u:.0e}":>7}' for u in UMBRALES))
    print('  diferencia relativa por ventana, y primer tiempo en que pasa de cada umbral:')
    print(cabecera)
    series = {}
    tc = np.loadtxt(os.path.join(bc, nombre, 'hk1_complex.tl'))
    tl = np.loadtxt(os.path.join(bl, nombre, 'hk1_complex.tl'))
    n = min(len(tc), len(tl))
    if len(tc) != len(tl):
        print(f'  aviso: {len(tc)} y {len(tl)} tiempos; se comparan los primeros {n}')
    assert np.allclose(tc[:n, 0], tl[:n, 0]), 'los tiempos de salida no coinciden'
    t = tl[:n, 0]
    for k in (0, 1, 2):
        a, b = tc[:n, 1 + 2*k] + 1j*tc[:n, 2 + 2*k], tl[:n, 1 + 2*k] + 1j*tl[:n, 2 + 2*k]
        print(linea(f'h_{k} del código (.tl)', t, a, b))
    sc, sl = (os.path.join(b, nombre, 'serie.npz') for b in (bc, bl))
    if os.path.exists(sc) and os.path.exists(sl):
        c, l = np.load(sc), np.load(sl)
        m = min(len(c['t']), len(l['t']))
        print(linea('h_1 (serie.npz)', l['t'][:m], c['h1'][:m], l['h1'][:m]))
        print(f'  energía: max|E_cluster/E_portátil - 1| = {np.max(np.abs(c["E"][:m]/l["E"][:m] - 1)):.1e};  '
              f'max|E/E(0) - 1| = {np.max(np.abs(c["E"]/c["E"][0] - 1)):.1e} (cluster), '
              f'{np.max(np.abs(l["E"]/l["E"][0] - 1)):.1e} (portátil)')
        series = dict(c=c, l=l, m=m)
    else:
        print('  sin serie.npz en alguna de las dos: no se compara h_1 en el mapa del equilibrio')
    return series


def respuesta(d, z, eps):
    from hadzic import ganancia
    m = min(d['m'], z['m'])
    t = d['l']['t'][:m]
    xc = (d['c']['h1'][:m] - z['c']['h1'][:m])/eps
    xl = (d['l']['h1'][:m] - z['l']['h1'][:m])/eps
    print(f'\n=== respuesta chi_1 = (h_1 de D - h_1 de Z)/eps, eps = {eps:g}')
    print(f'{"":>22} ' + ' '.join(f'{f"[{a},{b}]":>11}' for a, b in VENTANAS) + '   '
          + ' '.join(f'{f">{u:.0e}":>7}' for u in UMBRALES))
    print(linea('cluster - portátil', t, xc, xl))
    T0 = (500, 1000, 1500, 2000, 2500, 3000, 3500)
    g = [ganancia(t, xc, xl, a) for a in T0]
    print('  cociente complejo de amplitudes del cluster respecto del portátil (ventana de Hann de ancho 200):')
    print('  ' + f'{"t0":>8} ' + ' '.join(f'{a:>8}' for a in T0))
    print('  ' + f'{"|kappa|":>8} ' + ' '.join(f'{abs(v):8.4f}' for v in g))
    print('  ' + f'{"arg":>8} ' + ' '.join(f'{np.angle(v):+8.4f}' for v in g))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('nombres', nargs='+', help='una corrida, o la perturbada y su referencia')
    ap.add_argument('--cluster', default=os.path.join(RAIZ, 'exe', 'cluster', 'hadzic'))
    ap.add_argument('--local', default=os.path.join(RAIZ, 'exe', 'hadzic'))
    ap.add_argument('--eps', type=float, help='amplitud de la corrida perturbada, para la respuesta')
    a = ap.parse_args()
    res = [comparar(n, a.cluster, a.local) for n in a.nombres]
    if len(res) == 2 and a.eps and all(res):
        respuesta(res[0], res[1], a.eps)
