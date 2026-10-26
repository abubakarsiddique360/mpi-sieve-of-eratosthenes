#!/bin/bash
# ============================================================================
#  run_final.sh  --  final grading test:  n = 4 000 000 000 (all versions, np=8)
#                     plus the scaling of the final version with np = 1..16
#  (run on the laboratory server, results are stored in ./results)
# ============================================================================
cd "$HOME/mpi_exp"

echo "=== machine information ==="
hostname
nproc
echo "=== machine load before the tests ==="
uptime
free -h | head -2
echo

echo "=== final grading test : np = 8 , n = 4000000000 ==="
for v in base optimize1 optimize2 optimize3 optimize4; do
    echo "### VERSION=$v  np=8  n=4000000000  start=$(date +%T)"
    mpirun -np 8 ./$v 4000000000
    echo
done

echo "=== scaling of the final version : n = 4000000000 ==="
for np in 1 2 4 8 16; do
    echo "### VERSION=optimize4  np=$np  n=4000000000  start=$(date +%T)"
    mpirun -np $np ./optimize4 4000000000
    echo
done

echo "=== machine load after the tests ==="
uptime
echo "DONE $(date)"
