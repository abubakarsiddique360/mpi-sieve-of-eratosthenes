#!/bin/bash
# ============================================================================
#  run_repeat.sh -- repeat every measurement 3 times so that the minimum can
#                   be used in the report (the laboratory server is shared,
#                   the load average must be documented)
# ============================================================================
cd "$HOME/mpi_exp"
REP=${1:-3}

echo "=== machine : $(hostname) , $(nproc) cpus   date : $(date) ==="
echo "=== load before : $(uptime) ==="
echo

for r in $(seq 1 $REP); do
  echo "##################### REPETITION $r / $REP #####################"

  echo "===== TEST A : scaling , n = 1000000000 ====="
  for v in base optimize1 optimize2 optimize3 optimize4; do
    for np in 1 2 4 8; do
      echo "##A $v np=$np n=1000000000 rep=$r"
      mpirun -np $np ./$v 1000000000
    done
  done

  echo "===== TEST B : grading test , n = 4000000000 , np = 8 ====="
  for v in base optimize1 optimize2 optimize3 optimize4; do
    echo "##B $v np=8 n=4000000000 rep=$r"
    mpirun -np 8 ./$v 4000000000
  done

  echo "===== TEST C : scaling of the final version , n = 4000000000 ====="
  for np in 1 2 4 8 16; do
    echo "##C optimize4 np=$np n=4000000000 rep=$r"
    mpirun -np $np ./optimize4 4000000000
  done
done

echo "=== load after : $(uptime) ==="
echo "DONE $(date)"
