#!/bin/bash
# Reúne lo que hace falta saber del cluster para configurar la compilación y los trabajos:
# compiladores, HDF5, Python, la partición de Slurm y el almacenamiento. No cambia nada.
#
#   bash reproducir/cluster/diagnostico.sh [partición] > diagnostico.txt 2>&1
#
# Se corre en el nodo de entrada. La partición por omisión es "olin". Los módulos que se
# cargan para la prueba son los de entorno.sh.

export LC_ALL=C
PART="${1:-olin}"
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sec () { echo; echo "=== $* ==="; }
ver () { for x in "$@"; do if command -v "$x" > /dev/null 2>&1; then echo "$x: $(command -v "$x")"; else echo "$x: no está"; fi; done; }

sec "máquina"
date; hostname; whoami
head -2 /etc/os-release 2>/dev/null; uname -r
echo "núcleos (nproc): $(nproc)"
lscpu 2>/dev/null | grep -E "Model name|Nombre del modelo|^CPU\(s\)|Thread|Hilo|Core|Núcleo|Socket"
free -g 2>/dev/null | head -2

sec "módulos disponibles (compiladores, hdf5, python)"
if type module > /dev/null 2>&1; then
  module --version 2>&1 | head -2
  module avail 2>&1 | grep -i -E "intel|oneapi|gcc|gnu|hdf5|python|conda|lamod" | head -60
else
  echo "no hay comando module"
fi

sec "antes de cargar módulos"
ver gfortran ifort ifx h5fc python3 git make
gfortran --version 2>/dev/null | head -1
git --version

sec "con los módulos de entorno.sh"
source "$AQUI/entorno.sh"
type module > /dev/null 2>&1 && module list 2>&1 | head -30
ver gfortran ifort ifx h5fc h5pfc python3
ifort --version 2>/dev/null | head -1
ifx --version 2>/dev/null | head -1
gfortran --version 2>/dev/null | head -1
echo "VP_FC=$VP_FC"

sec "HDF5"
if command -v h5fc > /dev/null 2>&1; then echo "h5fc -show:"; h5fc -show; else echo "no hay h5fc"; fi
env | grep -i -E "^(HDF5|H5)[A-Z_]*=" | head
ldconfig -p 2>/dev/null | grep -i hdf5 | head
for d in /usr/lib64/gfortran/modules /usr/include /usr/include/hdf5/serial /usr/lib64/openmpi/include; do
  [ -f "$d/hdf5.mod" ] && echo "hdf5.mod en $d"
done
rpm -qa 2>/dev/null | grep -i hdf5 | head

sec "Python"
if command -v python3 > /dev/null 2>&1; then
  python3 --version
  for m in numpy scipy h5py matplotlib; do
    python3 -c "import $m; print('$m', $m.__version__)" 2>/dev/null || echo "$m: no está"
  done
fi
ver conda mamba

sec "Slurm, partición $PART"
if command -v sinfo > /dev/null 2>&1; then
  sinfo -p "$PART" -o "%P %a %l %D %c %m %f %N" 2>&1
  sinfo -p "$PART" -N -o "%N %c %m %O %T" 2>&1 | head -12
  scontrol show partition "$PART" 2>&1 | head -12
  for nodo in $(sinfo -h -p "$PART" -N -o "%N" 2>/dev/null | sort -u | head -4); do
    echo "--- nodo $nodo:"
    scontrol show node "$nodo" 2>&1 | grep -o -E "(CPUAlloc|CPUTot|CPULoad|Sockets|CoresPerSocket|ThreadsPerCore|RealMemory|AllocMem|FreeMem|State|Partitions|Reason)=[^ ]*" | tr '\n' ' '; echo
    echo "trabajos en $nodo, de cualquier partición:"
    squeue -w "$nodo" -o "%.9i %.10P %.14j %.9u %.2t %.11M %.4C %.8m %R" 2>&1 | head -25
  done
  echo "trabajos propios:"; squeue -u "$USER" -o "%.9i %.10P %.14j %.2t %.11M %.4C %.8m %R" 2>&1 | head -10
  sacctmgr -n show assoc user="$USER" format=account,partition,maxjobs,maxsubmit,grptres%40 2>/dev/null | head -5
else
  echo "no hay sinfo"
fi

RAIZ="$(cd "$AQUI/../.." && pwd)"
g () { (cd "$RAIZ" && git "$@"); }          # el git 1.8 de CentOS 7 no tiene -C

sec "almacenamiento"
df -h "$HOME" 2>/dev/null | tail -1
df -h "$RAIZ" 2>/dev/null | tail -1
quota -s 2>/dev/null | tail -3
ulimit -s | sed 's/^/pila (ulimit -s): /'

sec "copia del código en el cluster"
echo "$RAIZ"
g status -sb 2>&1 | head -3
g log -1 --format="%h %ad %s" --date=short 2>&1
g log -1 --format="src/: %h %ad" --date=short -- src/ 2>&1
ls -la "$RAIZ/exe/VP_PIC" 2>&1

sec "acceso a GitHub desde aquí"
url=$(g config --get remote.origin.url 2>/dev/null)
echo "remoto: ${url:-ninguno}"
[ -n "$url" ] && { timeout 15 git ls-remote "$url" HEAD 2>&1 | head -2 || echo "sin acceso"; }
