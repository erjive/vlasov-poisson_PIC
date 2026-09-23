# Reproducir las notas

Todo lo necesario para regenerar los resultados, las figuras y los notebooks de
`docs/introduccion/vlasov_intro.tex` a partir del código de este repositorio.

```
reproducir/
├── correr.sh          corre las simulaciones de un grupo, una detrás de otra
├── corridas/
│   ├── base/          archivos base (base_dftest.par, base_selfgrav.par)
│   ├── 08_verificacion/   67 corridas: estudio agnóstico y capturas (sección 8)
│   ├── 09_autogravedad/   30 corridas: caso de estudio y mapa numérico (secciones 9-10)
│   └── 11_landau/         10 corridas: equilibrio autoconsistente (sección 11)
├── scripts/           análisis, estados iniciales, teoría lineal
├── notebooks/         df_agnostic.ipynb, selfgrav_mixing.ipynb
└── figuras/           generar_figuras.py (escribe los PDF en docs/introduccion/figuras)
```

Cada `.par` de `corridas/` es el `params_usados.par` que escribió la corrida original:
la configuración completa, con la línea de comando en el encabezado. Las rutas
(`directory`, `checkpointfile`) son relativas a `exe/`.

> **Nota sobre `dftype`.** Hasta el commit que corrigió `dump_parameters`, el código no
> escribía `dftype` en `params_usados.par`. En estos archivos se agregó a partir de la
> línea `Command`. Verificado: `king_quad_500` y `spiral_mcs3_500` se reproducen a
> 1e-22 con los archivos corregidos.

## Requisitos

- Compilar el código: `make` en la raíz (deja `exe/VP_PIC`).
- Python 3 con `numpy`, `h5py`, `matplotlib`; `nbconvert` para ejecutar los notebooks;
  TeX Live con `latexmk` para el documento.

## Orden

Todos los comandos desde la raíz del repositorio. Nada corre en paralelo.

### 1. Simulaciones

```bash
reproducir/correr.sh 08_verificacion
reproducir/correr.sh 09_autogravedad
reproducir/correr.sh 11_landau              # genera antes los estados iniciales
```

Un subconjunto: `reproducir/correr.sh 09_autogravedad scan_dr`. Ver sin ejecutar:
añadir `-n`. Tiempos medidos con 4 hilos: con autogravedad y N=1e4, unos 65 s hasta
t=2000 y 11 min hasta t=20000; con N=1e5 hasta t=2000, 11 min; las de Landau, 2 min
(N=1e4) y 4 min (N=2e4) hasta t=3000; sin autogravedad, segundos.

### 2. Análisis que alimentan figuras y tablas

```bash
cd reproducir/scripts
python3 aa_meseta.py long20k_snap          # mapa numérico sobre 501 instantáneas
python3 filas_J.py                         # análisis por filas de acción
python3 rotadores.py                       # streaming libre en variables verdaderas
python3 giro_resolucion.py                 # componente que gira frente a resolución, a0=1e-3
python3 giro_resolucion.py long_a0_1e-2 long20k_a0_1e-2_nrc800
python3 prueba_dt.py                       # tabla del paso de tiempo

cd ../../exe
S=../reproducir/scripts
for d in eq_a1e-2_eps0.1 eq_a1e-2_eps0.05; do
  python3 $S/landau_analisis.py landau/$d landau/${d}_equilibrio.npz
done
for d in L_a1e-2_n400_e0 L_a1e-2_n400_e0.1 L_a1e-2_n800_e0 L_a1e-2_n800_e0.1 \
         L_a1e-2_n400_e0_c1 L_a1e-2_n400_e0.1_c1 L_a1e-2_n400_e0_c4 L_a1e-2_n400_e0.1_c4; do
  eq=landau/$(echo $d | sed 's/_c[14]$//')_equilibrio.npz
  python3 $S/landau_analisis.py landau/$d $eq               # ~9 min cada una
done
for d in L_a1e-2_n400_e0.1 L_a1e-2_n800_e0.1 L_a1e-2_n400_e0.1_c1 L_a1e-2_n400_e0.1_c4; do
  eq=landau/$(echo $d | sed 's/_c[14]$//')_equilibrio.npz
  python3 $S/landau_libre.py landau/$d $eq landau/$(echo $d | sed 's/_c[14]$//').dat
done
EQ=landau/L_a1e-2_n400_e0.1_equilibrio.npz
mkdir -p landau/lineal
python3 $S/lineal.py $EQ landau/lineal/libre_400x25.npz --nj 400 --nq 25 --dt 0.5 --libre
python3 $S/lineal.py $EQ landau/lineal/lin_400x25.npz   --nj 400  --nq 25 --dt 0.5
python3 $S/lineal.py $EQ landau/lineal/lin_800x25.npz   --nj 800  --nq 25 --dt 0.5
python3 $S/lineal.py $EQ landau/lineal/lin_1600x32.npz  --nj 1600 --nq 32 --dt 0.5
python3 $S/lineal.py $EQ landau/lineal/lin_1600x32_dt025.npz --nj 1600 --nq 32 --dt 0.25
python3 $S/lineal.py $EQ landau/lineal/lin_1600x64.npz  --nj 1600 --nq 64 --dt 0.5
cd ../reproducir/scripts
python3 landau_cola.py L_a1e-2 400 800                    # omega y gamma de la cola
python3 lineal_compara.py L_a1e-2 800 ../../exe/landau/lineal/lin_*.npz
```

### 3. Notebooks

```bash
cd reproducir/scripts
python3 build_nb.py && python3 build_nb_sg.py
cd ../notebooks
python3 -m nbconvert --to notebook --execute --inplace df_agnostic.ipynb
python3 -m nbconvert --to notebook --execute --inplace selfgrav_mixing.ipynb
```

### 4. Figuras y documento

```bash
cd reproducir/figuras && python3 generar_figuras.py        # todas; o un nombre: "fase"
cd ../../docs/introduccion && latexmk -pdf vlasov_intro.tex
```

`generar_figuras.py` imprime, para cada figura, las cifras que cita el texto.

## Qué produce cada figura

| figura (PDF) | sección | datos | además |
|---|---|---|---|
| `liouville`, `potencial`, `frecuencias`, `orbita`, `enrollamiento`, `hk_exacto` | 2-5 | analíticas | — |
| `convergencia` | 8.2 | `dfstudy/*_quad_*`, `dfstudy/*_mcs*_*` | — |
| `capturas_bimodal`, `capturas_spiral`, `capturas_king` | 8.3 | `dfsnap/*` | — |
| `agnosticas` | 8.3 | `dfstudy/bimodal_quad_50000`, `spiral_quad_50000`, `king_quad_*` | — |
| `salud` | 9 | `sg/long_fino`, `sg/quad`, `sg/quad_nosg`, `sg/mc` | — |
| `envolvente` | 9 | `sg/quad`, `sg/quad_nosg`, `sg/mc` | — |
| `barridos` | 9 | `sg/quad`, `sg/quad_nosg`, `sg/scan_*` | — |
| `masa` | 9 | `sg/long_a0_*`, `sg/mc_1000`, `sg/mc`, `sg/mc_100000`, `sg/quad` | — |
| `fase` | 9 | `sg/long_a0_*` | — |
| `accion` | 9 | `sg/long_fino` | — |
| `mapa_numerico` | 9 | `sg/long20k_snap`, `sg/long_fino` | `aa_meseta.py` |
| `mapa_validacion` | 10 | `sg/long_fino_nosg`, `sg/long20k_snap` | — |
| `mapa_marco` | 10 | `sg/long20k_snap` | `aa_meseta.py` |
| `landau_respuesta` | 11 | `landau/L_a1e-2_n800_e0{,.1}` | `landau_analisis`, `landau_libre`, `lineal` (1600×64) |
| `landau_controles` | 11 | `landau/L_a1e-2_n400_*`, `landau/eq_a1e-2_eps{0.1,0.05}` | `landau_analisis`, `lineal` (todas) |

## Qué produce cada tabla o cifra que no está en una figura

| resultado | sección | corridas | script |
|---|---|---|---|
| energía antes y después de la corrección | 9 | `sg/quad`, `sg/efix_quad`, `sg/efix_nosg` | notebook `selfgrav_mixing` |
| prueba de Δt | 9 | `sg/efix_quad`, `sg/dt_c1`, `sg/dt_c4` | `prueba_dt.py` |
| análisis por filas | 9 | `sg/long_fino`, `sg/long_fino_nosg`, `sg/long20k_snap` | `filas_J.py` |
| parte que gira frente a resolución | 9 | `sg/long20k_{snap,nrc800,npc50,nosg}`, `sg/long_a0_1e-2`, `sg/long20k_a0_1e-2_nrc800` | `giro_resolucion.py` |
| marco del potencial (promedio frente a instantáneo) | 10 | `sg/long20k_snap` | `aa_meseta.py` |
| ω, γ de la cola; comparación con teoría lineal | 11 | `landau/*` | `landau_cola.py`, `lineal_compara.py` |
