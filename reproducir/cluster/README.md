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
| `cola.sh` | why a job waits: partitions, the node, jobs ahead, estimated start, scheduler |
| `estado.sh` | what has run: recent jobs, the scaling table, and the state of each run |
| `serie.slurm` | computes `serie.npz` and `orbitas.npz` of finished runs, so that only those files are copied back |
| `traer.sh` | on the laptop: fetches the small files of finished runs |
| `comparar.py` | on the laptop: a run made on the cluster against the same run made on the laptop |

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

The runs listed in `EN_CLUSTER` of `hadzic.py` are meant for the cluster: `hadzic.py preparar`
writes their initial data and `.par` files, and `hadzic.py correr` skips them on the laptop.
At present they are the nine runs of series 6 (`DP`, `DPe1` and `ZP` for `k0.75_a1`, `k1_a1`
and `k2_a1`).

To see what has run, on the cluster:

```bash
bash reproducir/cluster/estado.sh > estado.txt 2>&1
```

It prints the recent jobs of the user, the scaling table and, for each run with recent
activity, whether it ended well, its number of particles, time step, final time, threads,
seconds, size and whether `serie.npz` exists. For the finished runs without `serie.npz` it
ends with the `sbatch` line of the next section.

## 3. Analysis

The snapshots of a run with 10⁶ particles take 17 GB. The analysis only needs two small files
of each run, which are computed on the cluster: `serie.npz`, the time series (1.3 MB), and
`orbitas.npz`, the angle and action of the particles near the edge every 100 time units
(10 MB for 10⁵ particles):

```bash
sbatch reproducir/cluster/serie.slurm DP_k1.25_a1 ZP_k1.25_a1
```

This needs Python 3 with `numpy` and `h5py`. Then, on the laptop:

```bash
export VP_CLUSTER=$DEST                       # with user@ in front if the user names differ
reproducir/cluster/traer.sh -d exe/hadzic DP_k1_a1 DPe1_k1_a1 ZP_k1_a1
python3 reproducir/scripts/hadzic.py analizar
```

`traer.sh` copies `serie.npz`, `orbitas.npz`, the `.tl` series, `params_usados.par`, and the `.ok`, `.meta`,
`.log` and Slurm output of each run, in one connection; `-g` adds `vlasov_output.h5`. With
`-d exe/hadzic` the runs land where `hadzic.py analizar` expects them, and the script refuses
to overwrite a run that already exists on the laptop. Without `-d` they go to
`exe/cluster/hadzic/`, which is the place for a run that was made on both machines:

```bash
reproducir/cluster/traer.sh DP_k1.25_a1 ZP_k1.25_a1
python3 reproducir/cluster/comparar.py --eps 0.03 DP_k1.25_a1 ZP_k1.25_a1
```

`comparar.py` gives the difference between the two machines for the series written by the
code and for `h_1`, the time at which it exceeds 10⁻¹², 10⁻⁹, 10⁻⁶ and 10⁻³, and the same for
the response `(D - Z)/eps`.

## Notes

- The node of `olin` has 24 cores with two threads each, and Slurm counts threads: `-c 16`
  gives 16 logical CPUs. The jobs start one OpenMP thread per physical core among the CPUs
  they were given (`vp_nucleos` in `entorno.sh`); `VP_HILOS` sets another number.
- That node also belongs to the partition `icn`, and the cluster schedules in order of
  arrival without filling gaps (`sched/builtin`, `priority/basic`). While an older job of
  `icn` waits for resources, every node of `icn`, this one included, is closed to newer
  jobs. A job of `olin` can then wait with free CPUs on its node, and neither a short time
  limit nor fewer CPUs change that. `cola.sh` shows the jobs ahead and the start time that
  Slurm estimates.
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
