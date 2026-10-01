"""Evolución con el código PIC de equilibrios construidos con la EDO
(equilibrio_edo.py), para comprobar que son estacionarios, y comparación con las
referencias Z de la demo eta, construidas con la iteración de punto fijo.

    python3 equilibrio_edo_pic.py preparar   datos iniciales (eps = 0) y .par
    python3 equilibrio_edo_pic.py correr     corridas PIC, una tras otra, 4 hilos
    python3 equilibrio_edo_pic.py analizar   deriva de J, dPhi y energía

Muestra: A1 (la masa menor), G1a (borde de King, g = 1), A4 (eta = 1) y M5
(a0 = 0.5, la mayor), con la resolución de la demo (400 x 25, dr = 0.1,
dt = 0.05); y M5 con cuatro veces más partículas (800 x 50) y con dr/2, para ver
de qué depende la deriva que quede. Salidas en exe/demo_eta/edo/.
"""
import os, sys, subprocess, time, argparse, numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from demo_eta import EXE, PARDIR, metadatos

BASE = os.path.join(EXE, 'demo_eta', 'edo')
IC_DEMO = os.path.join(EXE, 'demo_eta', 'ic')
JT, W0, DT = 0.138, 3.0, 0.05
SALIDA = 1000                                   # una instantánea cada 50 unidades de tiempo
# nombre, a0, g, referencia Z de la demo, t_fin, cambios en el .par, dato inicial
CORRIDAS = [
    ('E_A1', 0.0065, 2.0, 'Z_A1', 5550, {}, 'E_A1'),
    ('E_G1a', 0.0069, 1.0, 'Z_G1a', 5570, {}, 'E_G1a'),
    ('E_A4', 0.075, 2.0, 'Z_A4', 3950, {}, 'E_A4'),
    ('E_M5', 0.5, 2.0, 'Z_M5', 5000, {}, 'E_M5'),
    ('E_M5N', 0.5, 2.0, None, 5000, {'Nrc': '800', 'Npc': '50'}, 'E_M5N'),
    ('E_M5dr', 0.5, 2.0, None, 5000, {'dr': '0.05', 'courant': '2.0'}, 'E_M5'),
]


def ruta(*p):
    return os.path.join(BASE, *p)


def preparar():
    os.makedirs(ruta('ic'), exist_ok=True)
    for nombre, a0, g, ref, tfin, cambios, ic in CORRIDAS:
        dat = ruta('ic', f'{ic}.dat')
        nrc, npc = cambios.get('Nrc', '400'), cambios.get('Npc', '25')
        if not os.path.exists(dat):
            cmd = [sys.executable, os.path.join(AQUI, 'equilibrio.py'), '--metodo', 'edo',
                   '--forma', 'maxwell', '--jt', str(JT), '--k', str(g), '--w0', str(W0),
                   '--a0', str(a0), '--eps', '0', '--pert', 'suave3', '--nrc', nrc, '--npc', npc,
                   '--salida', dat]
            subprocess.run(cmd, check=True, stdout=open(ruta('ic', f'{ic}.log'), 'w'),
                           stderr=subprocess.STDOUT)
            print(f'  dato inicial {ic}', flush=True)
        valores = {'Nt': str(int(round(tfin/DT))), 'time_output': str(SALIDA),
                   'spatial_output': str(SALIDA), 'field_output': str(SALIDA),
                   'directory': f'demo_eta/edo/{nombre}', 'a0': str(a0),
                   'checkpointfile': f'demo_eta/edo/ic/{ic}.dat', **cambios}
        plantilla = open(os.path.join(PARDIR, 'demo__Z_A4.par')).read().split('\n')
        lineas = [f'# Equilibrio por EDO {nombre}: a0 = {a0}, g = {g:g}, eps = 0, t_fin = {tfin}.',
                  '# Generado por reproducir/scripts/equilibrio_edo_pic.py a partir de demo__Z_A4.par.']
        for l in plantilla:
            if l.startswith('#') or '=' not in l:
                continue
            clave = l.split('=')[0].strip()
            lineas.append(f'{clave:<16} = {valores[clave]}' if clave in valores else l)
        open(os.path.join(PARDIR, f'edo__{nombre}.par'), 'w').write('\n'.join(lineas) + '\n')
    print('preparado:', len(CORRIDAS), 'corridas')


def correr():
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    for nombre, *_ in CORRIDAS:
        if os.path.exists(ruta(f'{nombre}.ok')):
            continue
        par = os.path.join(PARDIR, f'edo__{nombre}.par')
        t0 = time.time()
        r = subprocess.run(['./VP_PIC', par], cwd=EXE, stdout=open(ruta(f'{nombre}.log'), 'w'),
                           stderr=subprocess.STDOUT, env=env)
        print(f'{"OK" if r.returncode == 0 else "FALLO":5} {nombre} ({time.time()-t0:.0f} s)', flush=True)
        if r.returncode == 0:
            open(ruta(f'{nombre}.ok'), 'w').write('')
            metadatos(nombre, par, ruta(f'{nombre}.meta'))


def medir(h5, npz, tmax, cada=50.0, mapa_pic=False):
    """Deriva de J (en el mapa del equilibrio inicial), dPhi respecto del potencial
    propio del equilibrio en r in [3, 15] y energía, en instantes múltiplos de cada.

    Con mapa_pic, J se calcula en el potencial que el propio código calcula en t = 0
    (el de su malla), no en el del equilibrio: así no entra la diferencia entre
    los dos, que hace oscilar la J medida sin que las partículas cambien de órbita."""
    import h5py
    from aa_numerico import MapaAA, phi_iso
    eq = np.load(npz)
    f = h5py.File(h5, 'r')
    pasos = sorted([k for k in f if k.startswith('step_')], key=lambda k: int(k.split('_')[1]))
    if mapa_pic:
        rg0 = f['grid']['r'][:]
        ps0 = f[pasos[0]]['potential'][:] - phi_iso(rg0)
        dentro = rg0 > 0
        mapa = MapaAA(rg0[dentro], ps0[dentro], L=float(eq['L0']))
    else:
        mapa = MapaAA(eq['r'], eq['phi_self'], L=float(eq['L0']))
    t = np.array([f[k].attrs['time'] for k in pasos])
    sel = [i for i, x in enumerate(t) if x <= tmax + 1e-6 and abs(x/cada - round(x/cada)) < 1e-6]
    rg = f['grid']['r'][:]
    zona = (rg >= 3) & (rg <= 15)
    phis = np.interp(rg, eq['r'], eq['phi_self'])
    escala = np.max(np.abs(phis[zona]))
    w = f[pasos[0]]['f'][:]
    E0 = f[pasos[0]].attrs['total_energy']
    J0 = None
    filas = []
    for i in sel:
        g = f[pasos[i]]
        _, J, _ = mapa(g['r_part'][:], g['p_part'][:])
        J0 = J if J0 is None else J0
        dJ = np.abs(J - J0)
        dphi = g['potential'][:] - phi_iso(rg) - phis
        filas.append((t[i], dJ.max(), np.average(dJ, weights=w),
                      np.sqrt(np.mean(dphi[zona]**2))/escala,
                      abs(g.attrs['total_energy'] - E0)/abs(E0)))
    f.close()
    return np.array(filas)


def analizar():
    lineas = []
    def w(x=''):
        lineas.append(x); print(x, flush=True)
    w('Estacionariedad: max|dJ| y <|dJ|> (pesado por la masa) al final, en unidades de J_t;')
    w('dPhi rms en 3 <= r <= 15 relativo a max|Phi_self| (al inicio y promedio en la segunda mitad);')
    w('error relativo máximo de la energía.')
    w(f'{"corrida":9} {"método":7} {"t_fin":>6}  {"max|dJ|/Jt":>10} {"<|dJ|>/Jt":>10}  '
      f'{"dPhi(0)":>8} {"dPhi 2a mitad":>13}  {"max|dE/E|":>9}')
    for nombre, a0, g, ref, tfin, cambios, ic in CORRIDAS:
        casos = [(nombre, 'EDO', ruta(nombre, 'vlasov_output.h5'), ruta('ic', f'{ic}_equilibrio.npz'))]
        if ref:
            casos.append((ref, 'Picard', os.path.join(EXE, 'demo_eta', ref, 'vlasov_output.h5'),
                          os.path.join(IC_DEMO, f'{ref}_equilibrio.npz')))
        for n, metodo, h5, npz in casos:
            if not os.path.exists(h5) or (metodo == 'EDO' and not os.path.exists(ruta(f'{n}.ok'))):
                w(f'{n:9} {metodo:7} (sin datos o en curso)')
                continue
            cache = ruta(f'medida_{n}_{tfin}.npy')        # se reutiliza si ya se midió
            if os.path.exists(cache):
                m = np.load(cache)
            else:
                m = medir(h5, npz, tfin)
                np.save(cache, m)
            seg = m[:, 0] >= 0.5*tfin
            w(f'{n:9} {metodo:7} {m[-1, 0]:6.0f}  {m[-1, 1]/JT:10.2e} {m[-1, 2]/JT:10.2e}  '
              f'{m[0, 3]:8.1e} {m[seg, 3].mean():13.1e}  {m[:, 4].max():9.1e}')
    open(ruta('resumen.txt'), 'w').write('\n'.join(lineas) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['preparar', 'correr', 'analizar'])
    a = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    globals()[a.paso]()
