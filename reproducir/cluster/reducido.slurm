#!/bin/bash
#SBATCH --job-name=reducido
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

# Series del modelo reducido de la nube (reproducir/scripts/reducido.py), una corrida detrás de
# otra, y sus medidas. Desde la raíz del repositorio:
#
#   sbatch reproducir/cluster/reducido.slurm <serie> [<serie> ...]      (Lb, control, ...)
#
# El costo es la suma directa sobre pares, que crece con el cuadrado del número de partículas y
# se reparte bien entre núcleos. Por omisión pide 16 CPUs, que en "olin" son 8 núcleos: la serie
# Lb, 4.7 h en los 4 núcleos del portátil, debería bajar a unas 2.5 h. Con "sbatch -c 48" pide el
# nodo entero, que es más rápido pero solo arranca cuando el nodo está vacío, y la cola no
# adelanta trabajos (README, Notes). Las corridas que ya existen en exe/demo_eta/reducido/ se saltan.
# Necesita Python 3 con numpy y scipy, y la biblioteca en C del modelo: la compila todo.sh de
# VlasovPoisson_PIC_sp en el nodo de entrada; si falta, se compila aquí, y entonces hace falta gcc.
# Fuera de Slurm corre igual, con los hilos de OMP_NUM_THREADS (4 por omisión).

RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$RAIZ/reproducir/cluster/entorno.sh"

# Un hilo por núcleo físico asignado (vp_nucleos, en entorno.sh). VP_HILOS fija otro número.
if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
  hilos=$(vp_nucleos)
  [ "$hilos" -gt "$SLURM_CPUS_PER_TASK" ] && hilos=$SLURM_CPUS_PER_TASK
else
  hilos="${OMP_NUM_THREADS:-4}"
fi
export OMP_NUM_THREADS="${VP_HILOS:-$hilos}" OMP_PLACES=cores OMP_PROC_BIND=close

[ $# -gt 0 ] || { echo "uso: reducido.slurm <serie> [<serie> ...]"; exit 1; }
cd "$RAIZ/reproducir/scripts" || exit 1
echo "máquina: $(hostname)   trabajo: ${SLURM_JOB_ID:-sin Slurm}   hilos: $OMP_NUM_THREADS"
echo "inicio: $(date)"
rc=0
for s in "$@"; do
  echo "== serie $s"
  python3 reducido.py serie "$s" || rc=1
done
# Las medidas de cada corrida (<nombre>_diag.npz y <nombre>_capas.npz, pequeños), con un solo hilo:
# así el resultado no depende del número de hilos del producto de matrices.
OMP_NUM_THREADS=1 python3 - "$@" <<'PY' || rc=1
import os, sys
import reducido as R
for s in sys.argv[1:]:
    for n, *_ in R.SERIES[s]:
        if os.path.exists(os.path.join(R.DATOS, n + '.npz')):
            R.diagnosticos(n); R.capas(n)
            print('medidas de', n, flush=True)
        else:
            print('FALTA', n, flush=True)
PY
echo "fin: $(date)   código de salida $rc"
exit $rc
