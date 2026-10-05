#!/bin/bash
# Envía a Slurm corridas del escenario de Hadžić (reproducir/scripts/hadzic.py), un trabajo
# por corrida. Desde la raíz del repositorio:
#
#   reproducir/cluster/enviar.sh [-c hilos] [-t hh:mm:ss] [-p partición] [-n] <nombre> [<nombre> ...]
#
# <nombre> es el de la corrida (DP_k1.25_a1, ZP_k1.25_a1, ...). Hacen falta su .par en
# reproducir/corridas/13_hadzic/ y su dato inicial en exe/hadzic/ic/, que escribe
# "hadzic.py preparar". -n muestra los comandos sin enviar. Una corrida con .ok se salta.
# Las corridas son independientes, así que D y Z pueden ir a la vez en nodos distintos.

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
hilos=16; tiempo=180:00:00; part=olin; seco=0
while getopts "c:t:p:n" o; do
  case $o in
    c) hilos=$OPTARG ;; t) tiempo=$OPTARG ;; p) part=$OPTARG ;; n) seco=1 ;;
    *) sed -n '2,11p' "${BASH_SOURCE[0]}"; exit 1 ;;
  esac
done
shift $((OPTIND-1))
[ $# -gt 0 ] || { sed -n '2,11p' "${BASH_SOURCE[0]}"; exit 1; }

cd "$RAIZ" || exit 1
mkdir -p exe/hadzic
for n in "$@"; do
  par="reproducir/corridas/13_hadzic/hadzic__$n.par"
  [ -f "$par" ] || { echo "$n: falta $par"; continue; }
  ic=$(awk -F= 'tolower($1) ~ /^checkpointfile/ {gsub(/ /,"",$2); print $2}' "$par")
  [ -f "exe/$ic" ] || { echo "$n: falta el dato inicial exe/$ic"; continue; }
  [ -f "exe/hadzic/$n.ok" ] && { echo "$n: ya está hecha (exe/hadzic/$n.ok)"; continue; }
  cmd=(sbatch -J "$n" -c "$hilos" -t "$tiempo" -p "$part"
       -o "exe/hadzic/slurm_%x_%j.out" -e "exe/hadzic/slurm_%x_%j.err"
       --export="ALL,VP_OK=hadzic/$n,VP_RAIZ=$RAIZ" reproducir/cluster/corrida.slurm "$par")
  if [ $seco = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
done
