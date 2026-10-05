#!/bin/bash
# Prueba corta del ejecutable, de segundos: 5000 partículas con autogravedad, 400 pasos,
# 2 hilos. La corre compilar.sh; se puede repetir en un nodo de cálculo con
# "srun -p olin -c 2 bash reproducir/cluster/probar.sh".
#
#  1. salida ascii: |h_0| ... |h_4| del final frente a los del portátil (gfortran 13.3).
#     Otro compilador no da los mismos bits; se pide acuerdo relativo a 1e-10. En el cluster
#     de LAMOD (gfortran 12.2, CentOS 7) la diferencia fue de 4.4e-16.
#  2. salida hdf5: que la corrida termine y escriba el archivo (comprueba el enlace con HDF5).
#  3. salida raw: lo mismo.

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
source "$AQUI/entorno.sh"
cd "$RAIZ/exe" || exit 1
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited 2>/dev/null

# Última línea de hk1.tl (t y |h_0| ... |h_4|) de esta misma corrida en el portátil.
REF="4.0000000000000298E+01  2.3454497237236343E-05  1.8864224087602536E-05  1.0389749004067114E-05  4.1788631000868525E-06  1.2890760909246304E-06"

corre () {      # corre <formato>
  rm -rf "cluster_prueba_$1"
  ./VP_PIC ../reproducir/corridas/base/base_selfgrav.par state=aa_quad Nrc=200 Npc=25 a0=0.01 \
      Nt=400 spatial_output=100 field_output=400 time_output=1000000 output_format="$1" \
      directory="cluster_prueba_$1" > "cluster_prueba_$1.log" 2>&1
}

fallos=0
if corre ascii && [ "$REF" = "__REFERENCIA__" ]; then
  echo "  ascii: corrió; falta la referencia del portátil. Última línea de hk1.tl:"
  tail -1 cluster_prueba_ascii/hk1.tl
elif [ -s cluster_prueba_ascii/hk1.tl ]; then
  tail -1 cluster_prueba_ascii/hk1.tl | awk -v ref="$REF" '
    { n = split(ref, r, " "); peor = 0
      for (i = 1; i <= n; i++) { d = ($i - r[i])/(r[i] == 0 ? 1 : r[i]); if (d < 0) d = -d; if (d > peor) peor = d }
      printf "  ascii: diferencia relativa máxima con el portátil %.1e\n", peor
      exit (peor < 1e-10 ? 0 : 1) }' || fallos=$((fallos+1))
else
  echo "  ascii: la corrida falló"; tail -5 cluster_prueba_ascii.log; fallos=$((fallos+1))
fi
for f in hdf5 raw; do
  if corre $f && [ -s "cluster_prueba_$f/vlasov_output.${f/hdf5/h5}" ]; then
    echo "  $f: bien ($(du -h "cluster_prueba_$f/vlasov_output.${f/hdf5/h5}" | cut -f1))"
  else
    echo "  $f: la corrida falló"; tail -5 "cluster_prueba_$f.log"; fallos=$((fallos+1))
  fi
done
grep -h "Time step fixed" cluster_prueba_ascii.log
[ $fallos = 0 ] && echo "PASA la prueba del ejecutable" || echo "FALLA la prueba del ejecutable ($fallos)"
exit $fallos
