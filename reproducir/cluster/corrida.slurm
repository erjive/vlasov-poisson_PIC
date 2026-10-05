#!/bin/bash
#SBATCH --job-name=vp_pic
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=180:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=4G

# Una corrida de VP_PIC. Se envía desde la raíz del repositorio:
#
#   sbatch -J <nombre> [-c <hilos>] reproducir/cluster/corrida.slurm <archivo.par> [nombre=valor ...]
#
# Los argumentos son los de VP_PIC: el archivo de parámetros y, si se quiere, cambios
# nombre=valor (el código ya no lee "< input_parameters"). Las rutas que hay dentro del .par
# son relativas a exe/, que es desde donde corre. Las opciones de arriba (partición, tiempo,
# hilos, memoria) se pueden cambiar al enviar: las de la línea de sbatch mandan. La memoria
# se pide porque la partición da 2 GB por omisión; una corrida de 10^6 partículas usa menos
# de 1 GB.
#
# Con la variable VP_OK=<prefijo>, relativo a exe/, la salida de VP_PIC va a <prefijo>.log y,
# si la corrida termina bien, se escriben <prefijo>.ok y <prefijo>.meta. Es lo que hace
# "hadzic.py correr" en el portátil y lo que esperan sus pasos de análisis; enviar.sh la pone.
#
# Fuera de Slurm corre igual (bash corrida.slurm archivo.par), con un hilo o los de
# OMP_NUM_THREADS.

RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [ ! -x "$RAIZ/exe/VP_PIC" ]; then
  echo "no hay $RAIZ/exe/VP_PIC: envíe el trabajo desde la raíz del repositorio o fije VP_RAIZ"
  exit 1
fi
source "$RAIZ/reproducir/cluster/entorno.sh"

# Un hilo por núcleo pedido, cada uno fijo en su núcleo.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-1}}"
export OMP_PLACES=cores OMP_PROC_BIND=close

# analysish reserva en la pila arreglos del tamaño del número de partículas (20 MB con un
# millón). Con ifort eso pasa del límite usual de 8 MB y la corrida muere sin mensaje.
ulimit -s unlimited 2>/dev/null || echo "aviso: no se pudo quitar el límite de pila ($(ulimit -s) kB)"
export OMP_STACKSIZE=256M

par="${1:?uso: corrida.slurm <archivo.par> [nombre=valor ...]}"; shift
case "$par" in /*) ;; *) par="$RAIZ/$par" ;; esac
[ -f "$par" ] || { echo "no existe $par"; exit 1; }

cd "$RAIZ/exe" || exit 1
echo "máquina: $(hostname)   trabajo: ${SLURM_JOB_ID:-sin Slurm}   hilos: $OMP_NUM_THREADS"
echo "inicio:  $(date)"
echo "corrida: ./VP_PIC $par $*"
t0=$SECONDS
if [ -n "${VP_OK:-}" ]; then
  mkdir -p "$(dirname "$VP_OK")"
  ./VP_PIC "$par" "$@" > "$VP_OK.log" 2>&1
  rc=$?
else
  ./VP_PIC "$par" "$@"
  rc=$?
fi
echo "fin:     $(date)   ($((SECONDS-t0)) s, código de salida $rc)"

if [ $rc = 0 ] && [ -n "${VP_OK:-}" ]; then
  : > "$VP_OK.ok"
  g () { (cd "$RAIZ" && git "$@") 2>/dev/null; }      # el git del cluster no tiene -C
  src=$(g log -1 --format=%h -- src/)
  [ -n "$(g status --porcelain src/)" ] && src="$src+cambios"
  {
    echo "corrida    $(basename "$VP_OK")"
    echo "VP_PIC     $(sha256sum VP_PIC | cut -d' ' -f1)"
    echo "par        $(sha256sum "$par" | cut -d' ' -f1)  ${par#$RAIZ/}"
    echo "HEAD       $(g rev-parse --short HEAD)"
    echo "src/       $src"
    echo "maquina    $(hostname), $OMP_NUM_THREADS hilos, compilador $VP_FC, trabajo ${SLURM_JOB_ID:-sin Slurm}, $((SECONDS-t0)) s"
  } > "$VP_OK.meta"
fi
exit $rc
