#!/bin/bash
# Qué se ha hecho en el cluster: los trabajos recientes del usuario, la prueba de escala y el
# estado de las corridas del escenario de Hadžić. Desde la raíz del repositorio:
#
#   bash reproducir/cluster/estado.sh [días]        (por omisión, los últimos 4 días)
#
# Solo lee. La salida es corta y se puede pegar tal cual.

export LC_ALL=C
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
dias="${1:-4}"
cd "$RAIZ" || exit 1

titulo () { echo; echo "=== $*"; }

echo "=== $(date '+%Y-%m-%d %H:%M')  $(hostname)  $RAIZ"
echo "repositorio: $(git log -1 --format='%h %ad %s' --date=short 2>/dev/null)"
cambios=$(git status --porcelain 2>/dev/null | grep -v '^??' | head -n 5)
[ -n "$cambios" ] && { echo "archivos versionados con cambios:"; echo "$cambios"; }
[ -x exe/VP_PIC ] && echo "exe/VP_PIC: $(sha256sum exe/VP_PIC | cut -c1-16), $(date -r exe/VP_PIC '+%Y-%m-%d %H:%M')" \
                  || echo "no hay exe/VP_PIC"

# Trabajos recientes. Si el cluster no guarda la contabilidad (sacct), se listan las salidas de
# Slurm de los últimos días, con su fecha y su última línea (la de "fin:" si la tiene).
titulo "trabajos de $USER en los últimos $dias días"
historia=""
if command -v sacct > /dev/null 2>&1; then
  historia=$(sacct -X -u "$USER" -S "now-${dias}days" \
             --format=JobID,JobName%22,Partition%9,AllocCPUS%5,State%12,Elapsed,Start,End,ExitCode 2>&1 | head -n 80)
fi
if [ -n "$historia" ] && ! echo "$historia" | grep -qi 'disabled\|error'; then
  echo "$historia"
else
  echo "sin contabilidad de Slurm (${historia:-no hay sacct}); salidas de Slurm recientes:"
  salidas=$( { find . -maxdepth 1 -name 'log.*.out' -mtime "-$dias"; find exe/hadzic -maxdepth 1 -name 'slurm_*.out' -mtime "-$dias"; } 2>/dev/null | sort)
  [ -n "$salidas" ] || echo "  ninguna"
  for f in $salidas; do
    fin=$(grep '^fin:' "$f" | tail -n 1); [ -n "$fin" ] || fin=$(tail -n 1 "$f")
    printf "  %-44s %s  %s\n" "${f#./}" "$(date -r "$f" '+%m-%d %H:%M')" "$fin"
  done
fi

titulo "en cola o corriendo (squeue)"
if command -v squeue > /dev/null 2>&1; then
  squeue -u "$USER" -o "%.10i %.22j %.9P %.8T %.11M %.4C %R" 2>&1 | head -n 40
else
  echo "no hay squeue"
fi

titulo "prueba de escala (log.vp_escala.*.out)"
hay=0
for f in log.vp_escala.*.out; do
  [ -f "$f" ] || continue
  hay=1
  echo "--- $f ($(date -r "$f" '+%Y-%m-%d %H:%M'))"
  cat "$f"
  e="${f%.out}.err"
  [ -s "$e" ] && { echo "--- $e, últimas líneas"; tail -n 5 "$e"; }
done
[ $hay = 0 ] && echo "no hay"

# Corridas con actividad en los últimos días: nombre, estado, partículas, paso de tiempo,
# tiempo alcanzado, hilos y segundos (del .meta que escribe corrida.slurm), tamaño y si ya
# tiene serie.npz.
titulo "corridas en exe/hadzic con cambios en los últimos $dias días"
printf "%-18s %-10s %9s %7s %8s %6s %8s %7s %6s\n" corrida estado partic dt t hilos segundos tamaño serie
falladas=""; sin_serie=""; n_corridas=0
for log in $(find exe/hadzic -maxdepth 1 -name '*.log' -mtime "-$dias" 2>/dev/null | sort); do
  n=$(basename "$log" .log)
  d="exe/hadzic/$n"
  par="reproducir/corridas/13_hadzic/hadzic__$n.par"
  [ -f "$par" ] || continue                      # solo las corridas definidas en hadzic.py
  n_corridas=$((n_corridas + 1))
  [ -f "$d/params_usados.par" ] && par="$d/params_usados.par"
  val () { awk -F= -v k="$1" 'tolower($1) ~ "^[ \t]*" k "[ \t]*$" { v = $2; sub(/#.*/, "", v); gsub(/[ \t]/, "", v); print v; exit }' "$par" 2>/dev/null; }
  nrc=$(val nrc); npc=$(val npc)
  partic="-"; [ -n "$nrc" ] && [ -n "$npc" ] && partic=$((nrc * npc))
  dt=$(grep -i -m 1 'time step fixed at size' "$log" 2>/dev/null | awk '{ printf "%.4g", $NF }')
  t="-"; [ -f "$d/hk2.tl" ] && t=$(tail -n 1 "$d/hk2.tl" | awk '{printf "%.0f", $1}')
  hilos="-"; seg="-"
  if [ -f "exe/hadzic/$n.meta" ]; then
    m=$(grep '^maquina' "exe/hadzic/$n.meta")
    hilos=$(echo "$m" | grep -o '[0-9]* hilos' | grep -o '[0-9]*')
    seg=$(echo "$m" | grep -o '[0-9]* s$' | grep -o '[0-9]*')
  fi
  if [ -f "exe/hadzic/$n.ok" ]; then
    estado=ok
  elif squeue -h -u "$USER" -n "$n" 2>/dev/null | grep -q .; then
    estado=en_cola
  else
    estado=SIN_OK; falladas="$falladas $n"
  fi
  serie=no
  if [ -f "$d/serie.npz" ]; then serie=si; elif [ "$estado" = ok ]; then sin_serie="$sin_serie $n"; fi
  tam="-"; [ -d "$d" ] && tam=$(du -sh "$d" 2>/dev/null | cut -f1)
  printf "%-18s %-10s %9s %7s %8s %6s %8s %7s %6s\n" "$n" "$estado" "$partic" "${dt:--}" "$t" "${hilos:--}" "${seg:--}" "$tam" "$serie"
done
[ $n_corridas = 0 ] && echo "ninguna"

for n in $falladas; do
  titulo "$n no tiene .ok: final de su registro y de la salida de Slurm"
  tail -n 6 "exe/hadzic/$n.log" 2>/dev/null
  for f in $(ls -t exe/hadzic/slurm_"$n"_*.out exe/hadzic/slurm_"$n"_*.err 2>/dev/null | head -n 2); do
    [ -s "$f" ] && { echo "--- $f"; tail -n 6 "$f"; }
  done
done

titulo "trabajos de series (log.vp_serie.*.out)"
hay=0
for f in $(find . -maxdepth 1 -name 'log.vp_serie.*.out' -mtime "-$dias" 2>/dev/null | sort); do
  hay=1
  echo "--- $f"; tail -n 8 "$f"
  e="${f%.out}.err"
  [ -s "$e" ] && { echo "--- $e, últimas líneas"; tail -n 5 "$e"; }
done
[ $hay = 0 ] && echo "no hay"

titulo "datos iniciales y espacio"
if [ -d exe/hadzic/ic ]; then
  echo "exe/hadzic/ic: $(ls exe/hadzic/ic/*.dat 2>/dev/null | wc -l) archivos .dat"
  ls exe/hadzic/ic/*.dat 2>/dev/null | sed 's|.*/||; s|\.dat$||' | tr '\n' ' ' | fold -s -w 100; echo
else
  echo "exe/hadzic/ic no existe: faltan los datos iniciales (mkdir -p exe/hadzic/ic y rsync desde el portátil)"
fi
[ -d exe/hadzic ] && echo "exe/hadzic ocupa $(du -sh exe/hadzic 2>/dev/null | cut -f1)"

if [ -n "$sin_serie" ]; then
  titulo "siguiente paso"
  echo "corridas terminadas sin serie.npz:$sin_serie"
  echo "  sbatch reproducir/cluster/serie.slurm$sin_serie"
fi
