#!/bin/bash
# Compila VP_PIC en el cluster con el entorno de entorno.sh y lo prueba (probar.sh).
#
#   bash reproducir/cluster/compilar.sh
#
# El Makefile compila a través de h5fc si está en el PATH y envuelve al compilador elegido.
# Si no, usa las rutas VP_HDF5_INC y VP_HDF5_LIBS de entorno.sh.

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
source "$AQUI/entorno.sh"
set -e
cd "$RAIZ"

args=(FC="$VP_FC")
if [ -n "${VP_HDF5_INC:-}" ]; then
  args+=(HDF5_WRAPPER= "HDF5_INC=$VP_HDF5_INC" "HDF5_LIBS=$VP_HDF5_LIBS")
fi
echo "compilador: $(command -v "$VP_FC" || echo "$VP_FC no está en el PATH")"
echo "h5fc:       $(command -v h5fc || echo "no está; se usan VP_HDF5_INC y VP_HDF5_LIBS")"
echo "make ${args[*]}"
make clean
make "${args[@]}"
echo "bibliotecas de exe/VP_PIC que no se encuentran:"
ldd exe/VP_PIC | grep "not found" || echo "  ninguna"
bash "$AQUI/probar.sh"
