"""Anima la evolucion de una distribucion en (r,p_r) y en (Q3,J3).

Un fichero por distribucion, con los dos espacios lado a lado y el mismo
instante en ambos.

Cada particula lleva dos codificaciones, porque con el esquema de cuadratura
hacen falta las dos: los nodos estan repartidos UNIFORMEMENTE en (Q,J) y la
distribucion vive en sus pesos f, no en la densidad de puntos.  Dibujar solo
posiciones mostraria la rejilla y no F0.  Asi que

  opacidad = peso de la particula   -> la forma de la distribucion
  tono     = angulo inicial Q3(0)   -> el enrollamiento por phase mixing,
                                       con un mapa ciclico porque Q es angular
"""
import sys, os, numpy as np, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from df0 import rp_to_QJ

DATA = '/home/erik/Documentos/Vlasov/vlasov-poisson_PIC/exe/dfvideo'
OUT  = '/home/erik/Documentos/Vlasov/vlasov-poisson_PIC/paper_runs/videos'

NPLOT = 30000     # particulas dibujadas, de las 50000
FPS   = 15
DPI   = 150       # 12.8 x 5.4 pulgadas a 150 ppp = 1920 x 810, dimensiones pares
FIGSZ = (12.8, 5.4)

# Calidad constante en vez de tasa de bits fija: una nube de puntos de un pixel
# es el peor caso para H.264, y a tasa fija el cuantizador la emborrona.  El
# preset lento y crf bajo cuestan tiempo de codificacion, no tamano excesivo.
X264 = ['-vcodec', 'libx264', '-crf', '15', '-preset', 'slow',
        '-pix_fmt', 'yuv420p', '-movflags', '+faststart']

TITULO = {
 'bimodal': 'bimodal  —  dos grupos de acciones, tres armonicos en $Q$',
 'spiral' : 'spiral  —  no separable, nace enrollada y se desenrolla en $t\\approx720$',
 'king'   : 'king  —  isoterma bajada con perturbacion de un armonico',
}

def cargar(df):
    f = h5py.File(os.path.join(DATA, df, 'vlasov_output.h5'), 'r')
    pasos = sorted([k for k in f if k.startswith('step_')],
                   key=lambda s: int(s.split('_')[1]))
    n0  = len(f[pasos[0]]['r_part'])
    sel = np.sort(np.random.default_rng(0).choice(n0, min(NPLOT, n0), replace=False))
    w = f[pasos[0]]['f'][:][sel]            # el peso de cada particula no cambia
    t, R, P, Q, J = [], [], [], [], []
    for k in pasos:
        g = f[k]
        r = g['r_part'][:][sel]; p = g['p_part'][:][sel]
        q, j = rp_to_QJ(r, p)
        t.append(float(g.attrs['time'])); R.append(r); P.append(p); Q.append(q); J.append(j)
    return np.array(t), np.array(R), np.array(P), np.array(Q), np.array(J), w

def animar(df):
    t, R, P, Q, J, w = cargar(df)

    # Color fijo por particula: tono del angulo inicial, opacidad proporcional
    # al peso.  La proporcion es directa, sin comprimir el rango: comprimirlo
    # hace visibles las alas y con ellas se pierde de vista donde esta de
    # verdad la masa, que es justo lo que hay que seguir.
    rgba = plt.get_cmap('twilight')(Q[0]/(2*np.pi))
    rgba[:, 3] = np.clip(w/w.max(), 0.0, 1.0)
    orden = np.argsort(w)                   # los pesos altos, encima

    rlo, rhi = np.nanmin(R), np.nanmax(R); dr = 0.04*(rhi-rlo)
    plo, phi = np.nanmin(P), np.nanmax(P); dp = 0.04*(phi-plo)
    jlo, jhi = np.nanmin(J), np.nanmax(J); dj = 0.04*(jhi-jlo)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=FIGSZ, dpi=DPI)
    fig.suptitle(TITULO[df], fontsize=13)
    kw = dict(s=2.0, c=rgba[orden], lw=0)
    s1 = a1.scatter(R[0][orden], P[0][orden], **kw)
    s2 = a2.scatter(Q[0][orden], J[0][orden], **kw)
    a1.set_xlim(rlo-dr, rhi+dr); a1.set_ylim(plo-dp, phi+dp)
    a1.set_xlabel('$r$', fontsize=12); a1.set_ylabel('$p_r$', fontsize=12)
    a1.set_title('espacio fisico', fontsize=11)
    a2.set_xlim(0, 2*np.pi); a2.set_ylim(jlo-dj, jhi+dj)
    a2.set_xticks([0, np.pi, 2*np.pi]); a2.set_xticklabels(['0', r'$\pi$', r'$2\pi$'])
    a2.set_xlabel('$Q_3$', fontsize=12); a2.set_ylabel('$J_3$', fontsize=12)
    a2.set_title('variables angulo-accion', fontsize=11)
    for a in (a1, a2):
        a.grid(alpha=.2); a.tick_params(labelsize=10)
    reloj = fig.text(0.5, 0.015, '', ha='center', fontsize=12)
    sm = plt.cm.ScalarMappable(cmap='twilight', norm=plt.Normalize(0, 2*np.pi))
    cb = fig.colorbar(sm, ax=a2, pad=0.02, ticks=[0, np.pi, 2*np.pi])
    cb.ax.set_yticklabels(['0', r'$\pi$', r'$2\pi$'])
    cb.set_label('$Q_3(0)$    (opacidad $\\propto$ peso)', fontsize=10)
    fig.tight_layout(rect=[0, 0.035, 1, 0.94])

    def frame(i):
        s1.set_offsets(np.column_stack([R[i][orden], P[i][orden]]))
        s2.set_offsets(np.column_stack([Q[i][orden], J[i][orden]]))
        reloj.set_text(f'$t = {t[i]:.0f}$        ({t[i]/100:.1f} periodos orbitales)')
        return s1, s2, reloj

    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, f'{df}_rp_QJ.mp4')
    FuncAnimation(fig, frame, frames=len(t), blit=False).save(
        dest, dpi=DPI,
        writer=FFMpegWriter(fps=FPS, extra_args=X264,
                            metadata=dict(title=f'VP_PIC {df}')))
    plt.close(fig)
    print(f'  {df:>8}: {len(t)} fotogramas -> {dest} ({os.path.getsize(dest)/1e6:.1f} MB)')

if __name__ == '__main__':
    for df in (sys.argv[1:] or ['bimodal', 'spiral', 'king']):
        animar(df)
