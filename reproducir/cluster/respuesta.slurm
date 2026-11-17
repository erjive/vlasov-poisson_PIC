#!/bin/bash
#SBATCH --job-name=respuesta
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G

# Respuesta lineal alrededor del estado mezclado de una corrida (reproducir/scripts/respuesta_mezcla.py): extrae el
# estado en t0 = 5 tau_1 de su bulto.npz y resuelve las evoluciones lineal y libre. Desde la raíz del repositorio:
#
#   sbatch -J resp_<corrida> reproducir/cluster/respuesta.slurm <carpeta de la corrida> <M_gas> <tau_1>
#
# Escribe respuesta_t5_ext.npz y respuesta_t5_lin.npz en la carpeta de la corrida (unos MB). La extracción mapea
# unas 400 instantáneas con un proceso por CPU; la evolución lineal tarda de 6 a 15 min.

RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$RAIZ/reproducir/cluster/entorno.sh"
export OMP_NUM_THREADS=1
[ $# -eq 3 ] || { echo "uso: respuesta.slurm <carpeta de la corrida> <M_gas> <tau_1>"; exit 1; }
cd "$RAIZ/reproducir/scripts" || exit 1
echo "máquina: $(hostname)   trabajo: ${SLURM_JOB_ID:-sin Slurm}   corrida: $1   inicio: $(date)"
python3 respuesta_mezcla.py extraer "$1/bulto.npz" "$1/respuesta_t5_ext.npz" --mgas "$2" --tau "$3" --procesos "${SLURM_CPUS_PER_TASK:-4}" &&
python3 respuesta_mezcla.py lineal "$1/respuesta_t5_ext.npz" "$1/respuesta_t5_lin.npz"
rc=$?
echo "fin: $(date)   código de salida $rc"
exit $rc
