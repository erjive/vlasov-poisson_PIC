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

# Núcleos físicos entre las CPUs que el trabajo tiene asignadas. Los nodos de "olin" tienen
# dos hilos por núcleo y Slurm cuenta hilos: "-c 16" son 16 CPUs lógicas, que pueden ser
# 8 núcleos. Los bucles del código están limitados por la memoria, así que se pone un hilo
# de OpenMP por núcleo físico.
vp_nucleos () {
  local lista
  lista=$(taskset -cp $$ 2>/dev/null | sed 's/.*: *//')
  LC_ALL=C lscpu -p=CPU,CORE,SOCKET 2>/dev/null | awk -F, -v lista="$lista" '
    BEGIN { n = split(lista, a, ",")
            for (i = 1; i <= n; i++) { m = split(a[i], b, "-"); for (c = b[1]; c <= (m > 1 ? b[2] : b[1]); c++) ok[c] = 1 } }
    !/^#/ && ($1 in ok) { nuc[$3 "," $2] = 1 }
    END { k = 0; for (x in nuc) k++; print (k > 0 ? k : 1) }'
}
