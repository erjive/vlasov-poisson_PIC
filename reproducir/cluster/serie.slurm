#!/bin/bash
#SBATCH --job-name=vp_serie
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

# Análisis en el cluster de corridas del escenario de Hadžić: calcula de cada una serie.npz
# (h_1 en el mapa del equilibrio, dPhi en la malla y la energía) y orbitas.npz (ángulo y acción
# de las partículas del borde cada 100 unidades de tiempo), que es lo único que hay que traer
# al portátil; el vlasov_output.h5 de una corrida de 10^6 partículas ocupa 17 GB.
# Desde la raíz del repositorio:
#
#   sbatch reproducir/cluster/serie.slurm <nombre> [<nombre> ...]
#
# Necesita Python 3 con numpy y h5py (ver entorno.sh). Reparte las instantáneas de cada
# corrida entre tantos procesos como núcleos pedidos.

RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$RAIZ/reproducir/cluster/entorno.sh"
cd "$RAIZ/reproducir/scripts" || exit 1
echo "máquina: $(hostname)   inicio: $(date)"
python3 hadzic.py serie -j "${SLURM_CPUS_PER_TASK:-4}" "$@"
rc=$?
echo "fin: $(date)   código de salida $rc"
exit $rc
