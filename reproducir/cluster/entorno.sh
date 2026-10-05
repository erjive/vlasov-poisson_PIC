# Entorno del cluster: módulos, compilador y rutas. Lo cargan compilar.sh y los trabajos de
# Slurm. Es el único archivo que hay que editar al cambiar de máquina o de compilador.
#
# Valores para xook.lamod.unam.mx y la partición "olin" (diagnóstico del 2026-10-05):
# gfortran 12.2 con el HDF5 1.10 de los módulos, que trae h5fc y es con lo que el Makefile
# compila sin más. El módulo de Intel (lamod/intel/oneAPI: ifort 2021.7 e ifx 2022.2) no
# tiene un HDF5 con interfaz de Fortran para ese compilador; para usarlo habría que
# compilar HDF5 con ifort y dar aquí VP_HDF5_INC y VP_HDF5_LIBS.

if type module > /dev/null 2>&1; then
  module purge
  module load lamod/gcc/12.2 lamod/hdf5/1.10
  # Python 3 con numpy y h5py, solo para el análisis en el cluster (serie.slurm).
  module load lamod/python3/3.12
fi

# Compilador para "make FC=...". El Makefile conoce gfortran e ifort.
export VP_FC="${VP_FC:-gfortran}"

# Solo si no hay h5fc para ese compilador: las rutas de HDF5 a mano, lo que imprime
# "h5fc -show" de la instalación que se quiera usar.
#   export VP_HDF5_INC="-I/ruta/include"
#   export VP_HDF5_LIBS="-L/ruta/lib -Wl,-rpath,/ruta/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5"
