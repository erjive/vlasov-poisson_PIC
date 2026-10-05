# Entorno del cluster: módulos, compilador y rutas. Lo cargan compilar.sh y los trabajos de
# Slurm. Es el único archivo que hay que editar al cambiar de máquina o de compilador.
#
# Los módulos son los del trabajo anterior en la partición "olin". Las líneas marcadas con
# COMPLETAR dependen de lo que diga diagnostico.sh.

if type module > /dev/null 2>&1; then
  module purge
  module load lamod/intel/oneAPI
fi

# Compilador para "make FC=...". El Makefile conoce gfortran e ifort.
export VP_FC="${VP_FC:-ifort}"

# HDF5 con interfaz de Fortran para ese mismo compilador. COMPLETAR una de las dos formas:
#   1. un módulo que deje h5fc en el PATH (el Makefile compila a través de él):
#        module load <hdf5>
#   2. las rutas a mano, lo que imprime "h5fc -show" de esa instalación:
#        export VP_HDF5_INC="-I/ruta/include"
#        export VP_HDF5_LIBS="-L/ruta/lib -Wl,-rpath,/ruta/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5"

# Python con numpy y h5py, solo para el análisis en el cluster (serie.slurm). COMPLETAR si
# hace falta un módulo:
#   module load <python>
