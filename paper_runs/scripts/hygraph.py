#!/usr/bin/env python3
"""
hygraph.py -- Visor/animador interactivo para la salida de VP_PIC
(equivalente a ygraph.py, pero para "vlasov_output.h5"/"vlasov_output.raw"
en vez de los .rl/.tl ASCII/xgraph que lee ygraph.py). Detecta el formato
por extensión del archivo (.h5 -> HDF5 vía h5py, .raw -> binario crudo vía
rawgraph_io.RawRun) y expone la misma interfaz de cuadros sin importar
cuál sea.

Cada corrida tiene una malla radial fija ("r") y un registro por
snapshot guardado (uno por cada spatial_output pasos), con datasets 1D
sobre la malla ("rho", "avg_rho", "curr", y "force"/"potential" si
autointeraction=.true.) y datasets por partícula ("r_part", "p_part",
"f"), más los atributos time/kinetic/potential/total_energy.

Dos modos:
  --mode line     (default) Curvas y(r) vs r -- una o más de
                  rho/avg_rho/curr/force/potential, superpuestas.
  --mode phase    Espacio de fase: scatter (r_part, p_part) coloreado
                  por f (como Fig. 1 del artículo).

Ninguno de los dos backends precarga todo el archivo en memoria (las
corridas grandes pesan cientos de MB a GB): cada cuadro se lee del
disco al vuelo, al mover el slider o durante la animación.

Uso:
    python3 hygraph.py corrida/vlasov_output.h5
    python3 hygraph.py corrida/vlasov_output.raw
    python3 hygraph.py corrida/vlasov_output.h5 --dataset rho avg_rho
    python3 hygraph.py corrida/vlasov_output.h5 --mode phase
    python3 hygraph.py corrida/vlasov_output.h5 --stride 10 -d 50
    python3 hygraph.py corrida/vlasov_output.h5 --static

Controles en la ventana: igual que ygraph.py -- slider de cuadro,
Play/Pause, ◀/▶, botón de velocidad (0.5x/1x/2x/4x/8x), zoom con la
barra de herramientas de matplotlib. Ejes fijos desde el primer cuadro
(recorridos una vez al arrancar) para que la escala no salte al animar.

Requiere: numpy, matplotlib, y h5py solo si vas a abrir un .h5.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


# ---------------------------------------------------------------------------
# Backends: cada uno expone n_steps, times[i], r (solo modo line), y
# get(i, name) -> np.ndarray para "rho"/"avg_rho"/.../"r_part"/"p_part"/"f".
# Así el visor de abajo no necesita saber si viene de HDF5 o de raw.
# ---------------------------------------------------------------------------

class HDF5Backend:
    def __init__(self, path, stride):
        import h5py
        self.h5 = h5py.File(path, "r")
        all_keys = sorted(
            (k for k in self.h5.keys() if k.startswith("step_")),
            key=lambda k: int(k.split("_")[1]),
        )
        if not all_keys:
            raise ValueError(f"'{path}' no tiene grupos /step_* -- ¿es un vlasov_output.h5?")
        self.keys = all_keys[::max(stride, 1)]
        if self.keys[-1] != all_keys[-1]:
            self.keys.append(all_keys[-1])
        self.times = [float(self.h5[k].attrs["time"]) for k in self.keys]
        self.r = self.h5["grid/r"][:]

    @property
    def n_steps(self):
        return len(self.keys)

    def get(self, i, name):
        return self.h5[self.keys[i]][name][:]

    def available(self, i):
        return list(self.h5[self.keys[i]].keys())

    def close(self):
        self.h5.close()


class RawBackend:
    def __init__(self, path, stride):
        from rawgraph_io import RawRun
        self.run = RawRun(path)
        if self.run.n_steps == 0:
            raise ValueError(f"'{path}' no tiene ningún registro -- ¿corrida vacía?")
        self.idx = list(range(0, self.run.n_steps, max(stride, 1)))
        if self.idx[-1] != self.run.n_steps - 1:
            self.idx.append(self.run.n_steps - 1)
        self.times = [self.run.times[i] for i in self.idx]
        self.r = self.run.r

    @property
    def n_steps(self):
        return len(self.idx)

    def get(self, i, name):
        return self.run.read(self.idx[i])[name]

    def available(self, i):
        return list(self.run.read(self.idx[i]).keys())

    def close(self):
        self.run.close()


def open_backend(path, stride):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".raw":
        return RawBackend(path, stride)
    return HDF5Backend(path, stride)  # default: .h5, or unrecognized -> try HDF5


class HYGraphViewer:
    def __init__(self, path, datasets=("rho",), mode="line",
                 stride=1, delay_ms=200, static=False):
        self.datasets = list(datasets)
        self.mode = mode
        self.delay_ms = delay_ms
        self.static = static

        self.backend = open_backend(path, stride)
        self.n_frames = self.backend.n_steps
        self.times = self.backend.times

        if self.mode == "line":
            avail = self.backend.available(0)
            for ds in self.datasets:
                if ds not in avail:
                    raise ValueError(f"dataset '{ds}' no existe -- disponibles: {avail}")

        self.fig, self.ax = plt.subplots(figsize=(9, 6))
        plt.subplots_adjust(bottom=0.22)

        if self.mode == "line":
            self._init_line()
        elif self.mode == "phase":
            self._init_phase()
        else:
            raise ValueError("--mode debe ser 'line' o 'phase'")

        self.title = self.ax.set_title("")

        self.frame_idx = 0
        self.playing = False
        self.speeds = [0.5, 1, 2, 4, 8]
        self.speed_idx = 1

        ax_slider = plt.axes([0.15, 0.1, 0.5, 0.03])
        self.slider = Slider(ax_slider, "Frame", 0, max(self.n_frames - 1, 0),
                              valinit=0, valstep=1)
        self.slider.on_changed(self._on_slider)

        ax_prev = plt.axes([0.67, 0.1, 0.045, 0.04])
        ax_play = plt.axes([0.72, 0.1, 0.07, 0.04])
        ax_next = plt.axes([0.795, 0.1, 0.045, 0.04])
        ax_speed = plt.axes([0.85, 0.1, 0.07, 0.04])
        self.btn_prev = Button(ax_prev, "◀")
        self.btn_play = Button(ax_play, "Play")
        self.btn_next = Button(ax_next, "▶")
        self.btn_speed = Button(ax_speed, f"{self.speeds[self.speed_idx]:g}x")
        self.btn_prev.on_clicked(lambda e: self.step(-1))
        self.btn_next.on_clicked(lambda e: self.step(1))
        self.btn_play.on_clicked(self._toggle_play)
        self.btn_speed.on_clicked(self._cycle_speed)

        self.timer = self.fig.canvas.new_timer(interval=self._current_interval())
        self.timer.add_callback(self._on_timer)

        self._draw_frame(0)

        if self.static or self.n_frames <= 1:
            for a in (ax_slider, ax_prev, ax_play, ax_next, ax_speed):
                a.set_visible(False)

    # ---- modo "line": rho/avg_rho/curr/force/potential vs r ------------

    def _init_line(self):
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(self.datasets), 1)))
        self.lines = []
        for i, ds in enumerate(self.datasets):
            (line,) = self.ax.plot([], [], "-", color=colors[i], lw=1.3, label=ds)
            self.lines.append(line)
        self.ax.legend(loc="best", fontsize=9)
        self.ax.set_xlabel("r")
        self.ax.set_ylabel("y")
        self._fix_axes_line()

    def _fix_axes_line(self):
        """Recorre todos los cuadros una vez (solo Nr valores, barato)
        para fijar xlim/ylim globales antes de animar."""
        gymin = gymax = None
        for ds in self.datasets:
            for i in range(self.n_frames):
                y = self.backend.get(i, ds)
                ymn, ymx = y.min(), y.max()
                gymin = ymn if gymin is None else min(gymin, ymn)
                gymax = ymx if gymax is None else max(gymax, ymx)
        pad = 0.05 * (gymax - gymin if gymax > gymin else 1.0)
        r = self.backend.r
        self.ax.set_xlim(r.min(), r.max())
        self.ax.set_ylim(gymin - pad, gymax + pad)

    def _draw_line(self, idx):
        r = self.backend.r
        for ds, line in zip(self.datasets, self.lines):
            line.set_data(r, self.backend.get(idx, ds))

    # ---- modo "phase": scatter (r_part, p_part) coloreado por f --------

    def _init_phase(self):
        # f se conserva a lo largo de cada trayectoria, así que el maximo
        # global es el del primer cuadro -- no hace falta recorrer todo.
        f0 = self.backend.get(0, "f")
        vmax = np.percentile(f0[f0 > 0], 99.5) if np.any(f0 > 0) else float(f0.max() or 1.0)

        self.scat = self.ax.scatter([], [], s=3, c=[], cmap="hot_r",
                                     vmin=0, vmax=vmax, linewidths=0)
        self.fig.colorbar(self.scat, ax=self.ax, label="f")
        self.ax.set_xlabel("r")
        self.ax.set_ylabel("$p_r$")
        self._fix_axes_phase()

    def _fix_axes_phase(self):
        gr0 = gr1 = gp0 = gp1 = None
        # Muestrear unos pocos cuadros (inicio/medio/fin) en vez de todos
        # para no cargar todas las partículas de cada snapshot al arrancar.
        sample_idx = sorted(set([0, self.n_frames // 2, self.n_frames - 1]))
        for idx in sample_idx:
            r = self.backend.get(idx, "r_part")
            p = self.backend.get(idx, "p_part")
            rmn, rmx = np.percentile(r, [0.1, 99.9])
            pmn, pmx = np.percentile(p, [0.1, 99.9])
            gr0 = rmn if gr0 is None else min(gr0, rmn)
            gr1 = rmx if gr1 is None else max(gr1, rmx)
            gp0 = pmn if gp0 is None else min(gp0, pmn)
            gp1 = pmx if gp1 is None else max(gp1, pmx)
        rpad = 0.05 * (gr1 - gr0)
        ppad = 0.05 * (gp1 - gp0)
        self.ax.set_xlim(gr0 - rpad, gr1 + rpad)
        self.ax.set_ylim(gp0 - ppad, gp1 + ppad)

    def _draw_phase(self, idx):
        r = self.backend.get(idx, "r_part")
        p = self.backend.get(idx, "p_part")
        f = self.backend.get(idx, "f")
        self.scat.set_offsets(np.column_stack([r, p]))
        self.scat.set_array(f)

    # ---- común -----------------------------------------------------------

    def _draw_frame(self, idx):
        idx = int(idx) % self.n_frames
        self.frame_idx = idx
        if self.mode == "line":
            self._draw_line(idx)
        else:
            self._draw_phase(idx)
        self.title.set_text(f"t = {self.times[idx]:.6g}")
        self.fig.canvas.draw_idle()

    def step(self, delta):
        self.slider.set_val((self.frame_idx + delta) % self.n_frames)

    def _on_slider(self, val):
        self._draw_frame(val)

    def _on_timer(self):
        self.step(1)

    def _toggle_play(self, event):
        self.playing = not self.playing
        if self.playing:
            self.btn_play.label.set_text("Pause")
            self.timer.start()
        else:
            self.btn_play.label.set_text("Play")
            self.timer.stop()

    def _current_interval(self):
        return max(int(self.delay_ms / self.speeds[self.speed_idx]), 1)

    def _cycle_speed(self, event):
        self.speed_idx = (self.speed_idx + 1) % len(self.speeds)
        self.btn_speed.label.set_text(f"{self.speeds[self.speed_idx]:g}x")
        was_running = self.playing
        if was_running:
            self.timer.stop()
        self.timer.interval = self._current_interval()
        if was_running:
            self.timer.start()

    def show(self):
        plt.show()
        self.backend.close()


def main():
    parser = argparse.ArgumentParser(
        description="Visor/animador interactivo para vlasov_output.h5/.raw (equivalente de ygraph.py)."
    )
    parser.add_argument("file", help="Ruta a vlasov_output.h5 o vlasov_output.raw")
    parser.add_argument("--dataset", "-y", nargs="+", default=["rho"],
                         help="Dataset(s) a graficar en modo 'line': rho avg_rho curr force potential (default: rho)")
    parser.add_argument("--mode", choices=["line", "phase"], default="line",
                         help="'line': y(r) vs r. 'phase': scatter (r_part,p_part) coloreado por f. (default: line)")
    parser.add_argument("--stride", type=int, default=1,
                         help="Usar 1 de cada N snapshots guardados (default: 1, todos)")
    parser.add_argument("-d", "--delay", type=int, default=200,
                         help="Delay entre cuadros de animación en ms (default: 200)")
    parser.add_argument("--static", action="store_true",
                         help="Mostrar solo el primer cuadro, sin controles de animación")
    args = parser.parse_args()

    try:
        viewer = HYGraphViewer(args.file, datasets=args.dataset, mode=args.mode,
                                stride=args.stride, delay_ms=args.delay, static=args.static)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    viewer.show()


if __name__ == "__main__":
    main()
