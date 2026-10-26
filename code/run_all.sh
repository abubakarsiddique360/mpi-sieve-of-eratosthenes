#!/bin/bash
# ============================================================================
#  run_all.sh -- build and test all five versions of the MPI sieve program
#  usage:  ./run_all.sh            (small + medium + final tests)
# ============================================================================
set -u
CC=mpic++
FLAGS="-O3"
mkdir -p results

echo "==============================================================="
echo " machine / environment"
echo "==============================================================="
hostname
nproc
uptime
free -h | head -2
$CC --version | head -2
mpirun --version | head -1
echo

# ---------------------------------------------------------------- compile ---
echo "==============================================================="
echo " compiling"
echo "==============================================================="
for v in base optimize1 optimize2 optimize3 optimize4; do
    echo "--- $CC $FLAGS -o $v $v.cpp"
    $CC $FLAGS -o $v $v.cpp || { echo "COMPILE ERROR in $v.cpp"; exit 1; }
done
echo "all five versions compiled successfully"
echo

# ------------------------------------------------- correctness (n = 1e8) ----
# reference: pi(10^8) = 5761455 , largest prime <= 1e8 = 99999989
echo "==============================================================="
echo " TEST 1 : correctness check,  n = 100000000 , np = 4"
echo "==============================================================="
for v in base optimize1 optimize2 optimize3 optimize4; do
    echo "### $v"
    mpirun -np 4 ./$v 100000000 | tee results/correct_${v}_np4.txt
    echo
done | tee results/correct_all.txt

# ------------------------------------------------- scaling   (n = 1e9) ------
echo "==============================================================="
echo " TEST 2 : scaling,  n = 1000000000 , np = 1,2,4,8"
echo "==============================================================="
for v in base optimize1 optimize2 optimize3 optimize4; do
  for np in 1 2 4 8; do
    echo "### $v   np=$np"
    mpirun -np $np ./$v 1000000000 | tee results/${v}_1e9_np${np}.txt
    echo
  done
done | tee results/scaling_1e9.txt

# ------------------------------------------------- big run   (n = 4e9) ------
echo "==============================================================="
echo " TEST 3 : final grading test,  n = 4000000000 , np = 8"
echo "==============================================================="
for v in base optimize1 optimize2 optimize3 optimize4; do
    echo "### $v   np=8   n=4000000000"
    mpirun -np 8 ./$v 4000000000 | tee results/${v}_4e9_np8.txt
    echo
done | tee results/final_4e9.txt

echo "==============================================================="
echo " all tests finished -- raw outputs are in ./results"
echo "==============================================================="
