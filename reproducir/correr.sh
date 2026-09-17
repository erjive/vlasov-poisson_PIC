#!/bin/bash
# Repite las simulaciones de las notas (docs/introduccion/vlasov_intro.tex).
#
#   reproducir/correr.sh <grupo> [patrón]      corre, una detrás de otra
#   reproducir/correr.sh <grupo> [patrón] -n   solo muestra los comandos
#
# <grupo>: 08_verificacion, 09_autogravedad u 11_landau (carpetas de corridas/).
# [patrón]: opcional, solo los .par cuyo nombre lo contiene (p. ej. "quad_500").
#
# Cada .par es el params_usados.par que escribió la corrida original: la
# configuración completa, con rutas relativas a exe/. Las corridas se ejecutan
# desde exe/ y escriben en el mismo directorio que la original (se sobrescribe).
#
# Las corridas de 11_landau leen su estado inicial de exe/landau/*.dat, que se
# genera aquí con equilibrio.py si no existe.
#
# Nunca corre dos simulaciones a la vez.

set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
grupo="${1:?uso: correr.sh <grupo> [patrón] [-n]}"
patron=""; seco=0
for a in "${@:2}"; do
  if [ "$a" = "-n" ]; then seco=1; else patron="$a"; fi
done
dir="$AQUI/corridas/$grupo"
[ -d "$dir" ] || { echo "no existe el grupo $grupo"; exit 1; }

cd "$RAIZ/exe" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OMP_PLACES=cores OMP_PROC_BIND=close

ejecuta () { if [ $seco = 1 ]; then echo "  $*"; else "$@"; fi; }

if [ "$grupo" = "11_landau" ]; then
  mkdir -p landau
  # Estados iniciales: equilibrio autoconsistente a0=1e-2 con perturbación eps cos Q.
  for spec in "eq_a1e-2_eps0.1 400 0.1" "eq_a1e-2_eps0.05 400 0.05" \
              "L_a1e-2_n400_e0 400 0" "L_a1e-2_n400_e0.1 400 0.1" \
              "L_a1e-2_n800_e0 800 0" "L_a1e-2_n800_e0.1 800 0.1"; do
    set -- $spec
    if [ ! -f "landau/$1.dat" ] || [ $seco = 1 ]; then
      ejecuta python3 "$AQUI/scripts/equilibrio.py" --a0 1e-2 --eps "$3" --nrc "$2" --npc 25 \
              --salida "landau/$1.dat"
    fi
  done
fi

for par in "$dir"/*"$patron"*.par; do
  [ -f "$par" ] || continue
  nombre=$(basename "$par" .par)
  destino=$(awk -F= '/^directory/ {gsub(/ /,"",$2); print $2}' "$par")
  if [ $seco = 1 ]; then
    echo "  ./VP_PIC $par    # -> exe/$destino"
    continue
  fi
  mkdir -p "$(dirname "$destino")"
  t0=$SECONDS
  if ./VP_PIC "$par" > "$destino.log" 2>&1; then
    echo "OK    $nombre ($((SECONDS-t0)) s)"
  else
    echo "FALLO $nombre"; tail -3 "$destino.log"
  fi
done
