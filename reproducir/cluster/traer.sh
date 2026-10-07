#!/bin/bash
# Trae del cluster al portátil lo necesario para analizar corridas hechas allá. En el portátil,
# desde la raíz del repositorio:
#
#   reproducir/cluster/traer.sh [-o origen] [-d destino] [-g] [-n] <nombre> [<nombre> ...]
#
# origen   repositorio en el cluster, usuario@máquina:/ruta (o la variable VP_CLUSTER).
# destino  carpeta local; por omisión exe/cluster/hadzic, para no pisar una corrida del mismo
#          nombre hecha en el portátil. Con "-d exe/hadzic" quedan donde "hadzic.py analizar"
#          las espera: úsese solo con corridas que no existan en el portátil.
# -g       trae también vlasov_output.h5 (1.6 GB por corrida de 10^5 partículas).
# -n       muestra lo que traería, sin copiar.
#
# Sin -g trae los archivos pequeños: serie.npz, las series .tl, params_usados.par, y el .ok,
# el .meta, el .log y la salida de Slurm de cada corrida. Una sola conexión para todas.

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
origen="${VP_CLUSTER:-}"; destino="exe/cluster/hadzic"; grande=0; seco=""
while getopts "o:d:gn" o; do
  case $o in
    o) origen=$OPTARG ;; d) destino=$OPTARG ;; g) grande=1 ;; n) seco="--dry-run" ;;
    *) sed -n '2,15p' "${BASH_SOURCE[0]}"; exit 1 ;;
  esac
done
shift $((OPTIND-1))
[ $# -gt 0 ] || { sed -n '2,15p' "${BASH_SOURCE[0]}"; exit 1; }
[ -n "$origen" ] || { echo "falta el origen: -o usuario@máquina:/ruta/del/repositorio, o la variable VP_CLUSTER"; exit 1; }

cd "$RAIZ" || exit 1
case "$destino" in /*) ;; *) destino="$RAIZ/$destino" ;; esac
if [ "$destino" = "$RAIZ/exe/hadzic" ] && [ -z "$seco" ]; then
  for n in "$@"; do
    [ -e "exe/hadzic/$n.ok" ] && { echo "$n ya existe en exe/hadzic del portátil: no se pisa. Use otro destino."; exit 1; }
  done
fi
mkdir -p "$destino"

filtros=()
for n in "$@"; do
  filtros+=(--include="/$n.ok" --include="/$n.meta" --include="/$n.log" --include="/slurm_${n}_*"
            --include="/$n/" --include="/$n/serie.npz" --include="/$n/*.tl" --include="/$n/params_usados.par")
  [ $grande = 1 ] && filtros+=(--include="/$n/vlasov_output.h5")
done
rsync -av $seco --prune-empty-dirs "${filtros[@]}" --exclude='*' "$origen/exe/hadzic/" "$destino/"
rc=$?
echo "destino: $destino"
exit $rc
