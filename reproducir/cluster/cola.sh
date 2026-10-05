#!/bin/bash
# Por qué espera un trabajo: estado de la partición y de su nodo, las otras particiones que
# usan ese nodo, la estimación de inicio de los trabajos propios y la configuración del
# planificador. No cambia nada.
#
#   bash reproducir/cluster/cola.sh [partición] > cola.txt 2>&1
#
# La partición por omisión es "olin".

export LC_ALL=C
PART="${1:-olin}"
sec () { echo; echo "=== $* ==="; }

sec "usuario y grupos"
id

sec "particiones: CPUs asignadas/libres/otras/total"
sinfo -o "%.12P %.6a %.11l %.6D %.12T %.18C %N" 2>&1

sec "partición $PART"
scontrol show partition "$PART" 2>&1

for nodo in $(sinfo -h -p "$PART" -N -o "%N" 2>/dev/null | sort -u | head -4); do
  sec "nodo $nodo"
  scontrol show node "$nodo" 2>&1
  for p in $(scontrol show node "$nodo" 2>/dev/null | grep -o "Partitions=[^ ]*" | cut -d= -f2 | tr ',' ' '); do
    [ "$p" = "$PART" ] && continue
    sec "partición $p, que comparte $nodo"
    scontrol show partition "$p" 2>&1 | grep -o -E "(PartitionName|AllowGroups|AllowAccounts|PriorityTier|PriorityJobFactor|MaxTime|DefaultTime|TotalCPUs|TotalNodes|OverSubscribe|PreemptMode|State|Nodes)=[^ ]*" | tr '\n' ' '; echo
    sinfo -p "$p" -o "%.12P %.6D %.12T %.18C %N" 2>&1
  done
done

sec "trabajos que se ven en $PART"
squeue -p "$PART" -o "%.9i %.10P %.14j %.9u %.2t %.11M %.11l %.4C %.20S %R" 2>&1 | head -30

# Con sched/builtin (sin relleno de huecos) y priority/basic (orden de llegada), un trabajo
# en espera cierra todos los nodos de su partición a los trabajos más nuevos, también a los
# de otra partición que comparta nodos. Los que van delante son los de número menor.
sec "todos los trabajos en espera, del más antiguo al más nuevo"
squeue -a -t PD --sort=i -o "%.9i %.10P %.14j %.10u %.11l %.5C %.5D %.20S %R" 2>&1 | head -40

sec "permisos y nodos de todas las particiones"
for p in $(sinfo -h -o "%R" 2>/dev/null | sort -u); do
  scontrol show partition "$p" 2>&1 | grep -o -E "(PartitionName|AllowGroups|AllowAccounts|MaxTime|PriorityTier|Nodes)=[^ ]*" | tr '\n' ' '; echo
done

sec "trabajos propios, con la estimación de inicio"
squeue -u "$USER" -o "%.9i %.10P %.14j %.2t %.11M %.11l %.4C %.20S %R" 2>&1
squeue -u "$USER" --start 2>&1
for j in $(squeue -h -u "$USER" -t PD -o "%i" 2>/dev/null | head -5); do
  echo "--- trabajo $j:"
  scontrol show job "$j" 2>&1 | grep -o -E "(JobState|Reason|Priority|Partition|TimeLimit|SubmitTime|EligibleTime|StartTime|NumCPUs|CPUs/Task|MinMemoryNode|MinMemoryCPU|Features|ReqNodeList|ExcNodeList|SchedNodeList|Dependency|QOS)=[^ ]*" | tr '\n' ' '; echo
  sprio -j "$j" 2>&1 | head -3
done

sec "reservaciones"
scontrol show reservation 2>&1 | head -20

sec "planificador"
scontrol show config 2>&1 | grep -E "^(SchedulerType|SchedulerParameters|SelectType|SelectTypeParameters|PreemptType|PreemptMode|PriorityType|PriorityWeight[A-Za-z]*|PrivateData|EnforcePartLimits|SlurmctldParameters|SLURM_VERSION)"
