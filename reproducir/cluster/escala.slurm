#!/bin/bash
#SBATCH --job-name=vp_escala
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=4G

# Tiempo por paso frente al número de hilos, con autogravedad y 10^5 y 10^6 partículas.
# Sirve para elegir cuántos núcleos pedir por corrida (-c). Desde la raíz del repositorio:
#
#   sbatch reproducir/cluster/escala.slurm            (o con -c <núcleos del nodo>)
#
# Cada medida es la diferencia entre una corrida con pasos y otra con cero pasos, para quitar
# el arranque. No escribe instantáneas intermedias. Tarda unos minutos; el límite de media
# hora es corto a propósito, para que Slurm pueda colarlo entre trabajos de más prioridad.
# Fuera de Slurm usa OMP_NUM_THREADS como tope:
# "OMP_NUM_THREADS=4 bash reproducir/cluster/escala.slurm".

export LC_ALL=C
RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$RAIZ/reproducir/cluster/entorno.sh"
tope="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-1}}"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_STACKSIZE=256M
ulimit -s unlimited 2>/dev/null
cd "$RAIZ/exe" || exit 1
echo "máquina: $(hostname)   CPUs pedidas: $tope   $(lscpu | grep -E 'Model name|Nombre del modelo' | sed 's/.*: *//')"
echo "CPUs asignadas: $(taskset -cp $$ 2>/dev/null | sed 's/.*: *//')   núcleos físicos: $(vp_nucleos)"

mide () {       # mide <hilos> <Nrc> <Npc> <pasos>: segundos de la corrida
  local t0 t1
  t0=$(date +%s.%N)
  OMP_NUM_THREADS=$1 ./VP_PIC ../reproducir/corridas/base/base_selfgrav.par state=aa_quad \
      Nrc=$2 Npc=$3 Nt=$4 spatial_output=1000000 field_output=1000000 time_output=1000000 \
      output_format=raw directory=cluster_escala > cluster_escala.log 2>&1 || { echo "falló"; return; }
  t1=$(date +%s.%N)
  awk "BEGIN { print $t1 - $t0 }"
}

printf "%6s %14s %12s %14s %12s\n" hilos "ms/paso 1e5" aceleración "ms/paso 1e6" aceleración
base5=""; base6=""
for n in 1 2 4 8 12 16 24 32 48; do
  [ "$n" -le "$tope" ] || break
  a0=$(mide $n 4000 25 0);  a1=$(mide $n 4000 25 200)
  b0=$(mide $n 4000 250 0); b1=$(mide $n 4000 250 40)
  m5=$(awk "BEGIN { print ($a1 - $a0)*1000/200 }"); m6=$(awk "BEGIN { print ($b1 - $b0)*1000/40 }")
  [ -z "$base5" ] && { base5=$m5; base6=$m6; }
  awk "BEGIN { printf \"%6d %14.2f %12.2f %14.2f %12.2f\\n\", $n, $m5, $base5/$m5, $m6, $base6/$m6 }"
done
rm -rf cluster_escala cluster_escala.log
