#!/bin/bash
cd /home/erik/Documentos/Vlasov/vlasov-poisson_PIC/exe || exit 1
export OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close
BASE=../reproducir/corridas/base/base_dftest.par
mkdir -p dfstudy/logs
fail=0
# cuadratura: rehacer, la ventana en J cambio
for df in bimodal spiral king; do
  for N in 500 5000 50000; do
    d=dfstudy/${df}_quad_${N}; rm -rf $d
    ./VP_PIC $BASE dftype=$df state=aa_quad Nrc=$((N/25)) Npc=25 directory=$d \
        > dfstudy/logs/${df}_quad_${N}.log 2>&1 \
      && echo "OK   quad ${df} N=${N}" || { echo "FALLO quad ${df} N=${N}"; fail=1; }
  done
done
d=dfstudy/bimodal_quad_50000_integratoranalytic; rm -rf $d
./VP_PIC $BASE dftype=bimodal state=aa_quad Nrc=2000 Npc=25 integrator=analytic \
    directory=$d > dfstudy/logs/control.log 2>&1 \
  && echo "OK   control analytic" || { echo "FALLO control"; fail=1; }
# Monte Carlo: cinco realizaciones por N
for df in bimodal spiral king; do
  for N in 500 5000 50000; do
    for sd in 1 2 3 4 5; do
      d=dfstudy/${df}_mcs${sd}_${N}
      [ -f $d/hk1_complex.tl ] && [ "$(awk '{print $1}' $d/hk1_complex.tl | sort -u | wc -l)" = 1001 ] && { echo "ya  ${df} N=${N} s${sd}"; continue; }
      rm -rf $d
      ./VP_PIC $BASE dftype=$df state=aa_random Nrc=$((N/25)) Npc=25 \
          seed=$((1000*sd+7)) directory=$d > dfstudy/logs/${df}_mcs${sd}_${N}.log 2>&1 \
        && echo "OK   mc ${df} N=${N} s${sd}" || { echo "FALLO mc ${df} N=${N} s${sd}"; fail=1; }
    done
  done
done
echo "TERMINADO fallos=$fail"
