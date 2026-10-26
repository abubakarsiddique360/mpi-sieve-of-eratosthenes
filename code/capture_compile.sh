#!/bin/bash
# ============================================================================
#  capture_compile.sh -- reproduce the build step and print a transcript of it
#                        (the transcript is used for the report screenshots)
# ============================================================================
cd "$HOME/mpi_exp"
PS1='a2023080146011@mpi25-1:~/mpi_exp$ '
set -x                      # print every command before executing it
rm -f base optimize1 optimize2 optimize3 optimize4
mpic++ --version | head -1
for v in base optimize1 optimize2 optimize3 optimize4; do
    mpic++ -O3 -o $v $v.cpp
done
ls -l *.cpp
ls -l base optimize1 optimize2 optimize3 optimize4
set +x
echo
echo "-> all five programs compiled without any warning or error"
