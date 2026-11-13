# Entorno del cluster: módulos, compilador y rutas. Lo cargan compilar.sh y los trabajos de
# Slurm. Es el único archivo que hay que editar al cambiar de máquina o de compilador.
#
# Valores para el cluster de LAMOD con AlmaLinux 9 (diagnóstico del 2026-11-12). Se usan el
# gfortran 11.5 y el HDF5 1.12 del sistema (h5fc en /usr/bin), con los que el código compila y la
# prueba del ejecutable difiere del portátil en 1.4e-16, y el Python 3.9 del sistema, que trae
# numpy 1.26. scipy no viene: se instala una vez con "python3 -m pip install --user scipy".
# No se cargan módulos: no se ha comprobado que el hdf5/2.2 de los módulos tenga la interfaz de
# Fortran para su gcc/12.2, y los python/3.11 y python/3.14 no se han revisado. "module purge"
# solo quita los que el usuario tenga cargados (Slurm no es un módulo).
#
# Hasta octubre de 2026, con CentOS 7, se cargaban lamod/gcc/12.2, lamod/hdf5/1.10 (estático) y
# lamod/python3/3.12; el módulo de Intel no tenía HDF5 para Fortran.

if type module > /dev/null 2>&1; then
  module purge
fi

# Compilador para "make FC=...". El Makefile conoce gfortran e ifort.
export VP_FC="${VP_FC:-gfortran}"

# Partición de Slurm en la que enviar.sh y todo.sh ponen los trabajos (las líneas #SBATCH de los
# guiones dicen "olin"; la opción -p de sbatch manda sobre ellas).
export VP_PARTICION="${VP_PARTICION:-olin}"

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
