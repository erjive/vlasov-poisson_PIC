#!/bin/bash
# Instantaneas de particulas cada 40 unidades para las capturas de la seccion 8.2
# de docs/introduccion. Las tres corridas, una detras de otra (~30 s cada una).
cd /home/erik/Documentos/Vlasov/vlasov-poisson_PIC/exe || exit 1
export OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close
mkdir -p dfsnap/logs
for df in bimodal spiral king; do
  d=dfsnap/$df; rm -rf $d; t0=$SECONDS
  if ./VP_PIC ../reproducir/corridas/base/base_dftest.par dftype=$df state=aa_quad Nrc=2000 Npc=25 \
       field_output=400 directory=$d > dfsnap/logs/$df.log 2>&1
  then echo "OK   $df ($((SECONDS-t0)) s, $(du -sh $d | cut -f1))"
  else echo "FALLO $df"; tail -3 dfsnap/logs/$df.log; fi
done
echo TERMINADO
