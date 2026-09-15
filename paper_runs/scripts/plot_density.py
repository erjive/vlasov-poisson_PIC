#!/usr/bin/env python3
"""Plot rho/avg_rho (r^2-weighted, as written by hdf5_io.f90) vs r,
overlaid at a few times, from a vlasov_output.h5 file.

Usage:
    python3 plot_density.py path/to/vlasov_output.h5 [out.png]
"""
import sys
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

path = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else 'density_check.png'

with h5py.File(path, 'r') as f:
    r = f['grid/r'][:]
    keys = sorted(
        (k for k in f.keys() if k.startswith('step_')),
        key=lambda k: int(k.split('_')[1])
    )

    # 4 snapshots spread across the run: start, 1/4, 1/2, end.
    # Swap for e.g. keys[::len(keys)//20] to get more curves.
    picks = [keys[0], keys[len(keys)//4], keys[len(keys)//2], keys[-1]]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for k in picks:
        g = f[k]
        t = g.attrs['time']
        axes[0].plot(r, g['rho'][:], label=f't={t:.0f}')
        axes[1].plot(r, g['avg_rho'][:], label=f't={t:.0f}')

    axes[0].set_title(r'rho  ($r^2\rho$, crudo)')
    axes[1].set_title(r'avg_rho  ($r^2\bar\rho$, promediado en celda)')
    for ax in axes:
        ax.set_xlabel('r')
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f'saved {out}')
