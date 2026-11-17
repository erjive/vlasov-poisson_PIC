#!/bin/bash
#SBATCH --job-name=perfiles
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

# Prueba con F_0 lento (reproducir/scripts/respuesta_mezcla.py): F_0(J, t) de una corrida con las variables de
# Phi_0 (paso perfiles), la evolución lineal con la parte lenta de ese F_0 en el término F_0'(J) d(dPhi)/dQ y la
# cuasilineal, en que F_0 cambia por el flujo en J de la propia solución lineal (lineal --perfiles --cuasilineal,
# que repite además la lineal con F_0 fijo y la libre). Necesita el respuesta_t5_ext.npz que escribe
# respuesta.slurm en la carpeta de la corrida. Desde la raíz del repositorio:
#
#   sbatch -J perf_<corrida> reproducir/cluster/respuesta_perfiles.slurm <carpeta de la corrida>
#
# Escribe respuesta_t5_perf.npz y respuesta_t5_lin_F0t.npz (unos 10 MB cada uno) en la carpeta de la corrida.
# El mapeo de unas 500 instantáneas tarda unos minutos con 8 procesos; las cuatro evoluciones, de 15 a 35 min.

RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$RAIZ/reproducir/cluster/entorno.sh"
export OMP_NUM_THREADS=1
[ $# -eq 1 ] || { echo "uso: respuesta_perfiles.slurm <carpeta de la corrida>"; exit 1; }
cd "$RAIZ/reproducir/scripts" || exit 1
echo "máquina: $(hostname)   trabajo: ${SLURM_JOB_ID:-sin Slurm}   corrida: $1   inicio: $(date)"
python3 respuesta_mezcla.py perfiles "$1/bulto.npz" "$1/respuesta_t5_ext.npz" "$1/respuesta_t5_perf.npz" \
        --procesos "${SLURM_CPUS_PER_TASK:-4}" &&
python3 respuesta_mezcla.py lineal "$1/respuesta_t5_ext.npz" "$1/respuesta_t5_lin_F0t.npz" \
        --perfiles "$1/respuesta_t5_perf.npz" --cuasilineal
rc=$?
echo "fin: $(date)   código de salida $rc"
exit $rc
