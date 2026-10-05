# Running on the cluster

Scripts to build VP_PIC on a Slurm cluster and to run and analyze there the runs that are
too large for the laptop (10⁵ to 10⁶ particles). The defaults are those of
`xook.lamod.unam.mx` and its partition `olin`: one node, `tochtli-24`, with 48 CPUs and
128 GB, no time limit, and 2 GB of memory per job unless more is requested. In the examples
the repository is in `$DEST`, for instance

```bash
DEST=xook.lamod.unam.mx:/storage/cactusolin/erik/newVlasov/vlasov-poisson_PIC
```

| file | purpose |
|---|---|
| `diagnostico.sh` | prints what the cluster offers: compilers, HDF5, Python, the partition, storage |
| `entorno.sh` | modules and compiler; the only file to edit for another machine |
| `compilar.sh` | builds `exe/VP_PIC` with that environment and runs `probar.sh` |
| `probar.sh` | a run of a few seconds, compared with the result of the laptop |
| `corrida.slurm` | one run: `sbatch -J name corrida.slurm file.par [name=value ...]` |
| `enviar.sh` | submits runs of `reproducir/scripts/hadzic.py`, one job each |
| `escala.slurm` | time per step against the number of threads, to choose `-c` |
| `serie.slurm` | computes `serie.npz` of finished runs, so that only that file is copied back |

All commands are given from the root of the repository.

## 1. First time

Copy the repository to the cluster, with `git clone` if the cluster reaches GitHub and with
`rsync` otherwise:

```bash
rsync -av --exclude exe --exclude objs --exclude .git ./ $DEST/
```

On the login node, collect the information and build:

```bash
bash reproducir/cluster/diagnostico.sh > diagnostico.txt 2>&1
# edit reproducir/cluster/entorno.sh if needed (the lines marked COMPLETAR)
bash reproducir/cluster/compilar.sh
```

The code needs HDF5 with its Fortran interface, built with the same compiler. On `xook`
that is gfortran 12.2 with the module `lamod/hdf5/1.10`, which provides `h5fc`; `entorno.sh`
loads both. The Intel module has no such HDF5. On another machine, load in `entorno.sh` a
module that provides `h5fc`, or set `VP_HDF5_INC` and `VP_HDF5_LIBS` there. `compilar.sh`
ends with `PASA la prueba del ejecutable` when the build works.

Then measure the scaling with threads, once:

```bash
sbatch reproducir/cluster/escala.slurm
```

The table in `log.vp_escala.<job>.out` gives the time per step with 1, 2, 4, ... threads.
Request for each run the number of cores beyond which the time stops falling. Several runs
with few cores each use a node better than one run with all of them.

## 2. A series of runs

The initial data are generated on the laptop, where Python is set up, and copied:

```bash
python3 reproducir/scripts/hadzic.py preparar                      # laptop
rsync -av exe/hadzic/ic/ $DEST/exe/hadzic/ic/
rsync -av reproducir/ $DEST/reproducir/
```

On the cluster:

```bash
reproducir/cluster/enviar.sh -n DP_k1.25_a1 ZP_k1.25_a1           # shows the commands
reproducir/cluster/enviar.sh -c 8 DP_k1.25_a1 ZP_k1.25_a1          # submits
squeue -u $USER
```

Each run writes to `exe/hadzic/<name>/`, its log to `exe/hadzic/<name>.log`, and
`exe/hadzic/<name>.ok` when it ends well, as `hadzic.py correr` does on the laptop.

## 3. Analysis

The snapshots of a run with 10⁶ particles take 17 GB. The analysis only needs the time series
of each run, 1.3 MB, which is computed on the cluster:

```bash
sbatch reproducir/cluster/serie.slurm DP_k1.25_a1 ZP_k1.25_a1
```

This needs Python 3 with `numpy` and `h5py`. Then, on the laptop:

```bash
for n in DP_k1.25_a1 ZP_k1.25_a1; do
  mkdir -p exe/hadzic/$n
  rsync -av $DEST/exe/hadzic/$n/serie.npz exe/hadzic/$n/
  rsync -av $DEST/exe/hadzic/$n.{ok,meta,log} exe/hadzic/
done
python3 reproducir/scripts/hadzic.py analizar
```

## Notes

- The node of `olin` has 24 cores with two threads each, and Slurm counts threads: `-c 16`
  gives 16 logical CPUs. The jobs start one OpenMP thread per physical core among the CPUs
  they were given (`vp_nucleos` in `entorno.sh`); `VP_HILOS` sets another number.
- That node also belongs to the partition `icn`. A job of `olin` waits while the CPUs it
  asks for are taken or held for other jobs. Ask for the time the run needs
  (`enviar.sh -t`), not for the maximum: a short limit lets Slurm fit the job in a gap.
- The code reads a parameter file given as an argument, `./VP_PIC file.par [name=value ...]`.
  The old form `./VP_PIC < input_parameters` no longer exists.
- Results obtained with another compiler agree with those of the laptop to rounding, not bit
  by bit. A run and its reference (D and Z) must be made with the same executable.
- `corrida.slurm` removes the stack limit: `analysish` keeps arrays of the size of the number
  of particles on the stack, 20 MB with 10⁶ particles.
- The acceleration criterion of `set_timestep` uses the cell size
  `drc = (rmaxc - rminc)/Nrc`, so with many rows in J it reduces the time step. In the runs
  of the Hadžić setting with mass one this happens above about 30 000 rows (the step was
  0.087 instead of 0.1 with `Nrc = 40000`). Keep `Nrc` below that, or check the line
  `Time step fixed at size` of the log.
