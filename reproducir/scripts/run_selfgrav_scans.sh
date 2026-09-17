#!/bin/bash
cd /home/erik/Documentos/Vlasov/vlasov-poisson_PIC/exe || exit 1
export OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close
B=../reproducir/corridas/base/base_selfgrav.par
mkdir -p sg/logs
go () { # $1=etiqueta  resto=overrides
  local tag=$1; shift
  local d=sg/$tag; rm -rf $d
  local t0=$SECONDS
  if ./VP_PIC $B state=aa_quad Nrc=400 Npc=25 directory=$d "$@" > sg/logs/$tag.log 2>&1
  then echo "OK   $tag  ($((SECONDS-t0)) s)"
  else echo "FALLO $tag"; tail -3 sg/logs/$tag.log; fi
}
echo "--- 1. barrido en dr, con dt=0.1 fijo (courant compensa) ---"
go scan_dr_0.200 dr=0.2   courant=1.0
go scan_dr_0.050 dr=0.05  courant=4.0
go scan_dr_0.025 dr=0.025 courant=8.0
echo "--- 2. orden del B-spline ---"
go scan_bspl_2 bsplineorder=2
go scan_bspl_3 bsplineorder=3
echo "--- 5. barrido en a0 ---"
go scan_a0_1e-4 a0=1.0e-4
go scan_a0_1e-2 a0=1.0e-2
echo "--- 3. escalado en N con cuadratura ---"
d=sg/scan_N_1000; rm -rf $d; t0=$SECONDS
./VP_PIC $B state=aa_quad Nrc=40 Npc=25 directory=$d > sg/logs/scan_N_1000.log 2>&1 \
  && echo "OK   scan_N_1000  ($((SECONDS-t0)) s)" || echo "FALLO scan_N_1000"
d=sg/scan_N_100000; rm -rf $d; t0=$SECONDS
./VP_PIC $B state=aa_quad Nrc=4000 Npc=25 directory=$d > sg/logs/scan_N_100000.log 2>&1 \
  && echo "OK   scan_N_100000  ($((SECONDS-t0)) s)" || echo "FALLO scan_N_100000"
echo TERMINADO
