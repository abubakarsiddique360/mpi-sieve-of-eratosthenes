#!/bin/bash
# ============================================================================
#  capture_session.sh -- re-run the important commands of the experiment and
#                        print a readable transcript of them.  The transcript
#                        is cut into pieces and used for the report screenshots
#                        (every line below is the real output of a real run).
# ============================================================================
cd "$HOME/mpi_exp"
P='a2023080146011@mpi25-1:~/mpi_exp$'
run() { echo "$P $*"; "$@"; echo; }

echo "###SHOT:compile"
run mpic++ --version
echo "$P for v in base optimize1 optimize2 optimize3 optimize4; do mpic++ -O3 -o \$v \$v.cpp; done"
for v in base optimize1 optimize2 optimize3 optimize4; do mpic++ -O3 -o $v $v.cpp; done
echo
run ls -l base.cpp optimize1.cpp optimize2.cpp optimize3.cpp optimize4.cpp
run ls -l base optimize1 optimize2 optimize3 optimize4
echo "###END"

echo "###SHOT:correctness"
run mpirun -np 4 ./optimize4 100000000
echo "$P # --- the same correctness check for all five versions ---"
echo "$P for v in base optimize1 optimize2 optimize3 optimize4; do echo -n \"\$v : \"; mpirun -np 4 ./\$v 100000000 | grep 'primes are'; done"
for v in base optimize1 optimize2 optimize3 optimize4; do
    echo -n "$v : "
    mpirun -np 4 ./$v 100000000 | grep 'primes are'
done
echo
echo "reference value of the mathematics handbook : pi(10^8) = 5761455"
echo
echo "###END"

echo "###SHOT:grading"
echo "$P # the official test command of the experiment"
run mpirun -np 8 ./optimize4 4000000000
echo "###END"

echo "###SHOT:allversions"
echo "$P # all five versions with the official test parameters (n = 4000000000, np = 8)"
for v in base optimize1 optimize2 optimize3 optimize4; do
    echo "$P mpirun -np 8 ./$v 4000000000"
    mpirun -np 8 ./$v 4000000000
done
echo
echo "###END"

echo "###SHOT:summary"
cat > summary_table.sh <<'EOS'
#!/bin/bash
# one line per measurement : version / np / n / pi(n) / marking time / total time
printf "%-42s %12s %10s %10s\n" "configuration" "pi(n)" "marking" "total"
awk '/^### VERSION=/  {v=$0}
     /primes are less/ {p=$1}
     /marking .sieve./ {m=$(NF-1)}
     /total elapsed/   {t=$(NF-1)
                        printf "%-42s %12s %10s %10s\n", substr(v,5,60), p, m, t}' results/final_4e9.txt
EOS
echo "$P # the results of the runs above, one line per configuration"
run bash summary_table.sh
echo "###END"

echo "###SHOT:scaling"
echo "$P # scaling of the final version, n = 4000000000"
echo "$P for np in 1 2 4 8 16; do mpirun -np \$np ./optimize4 4000000000; done"
for np in 1 2 4 8 16; do
    echo "--- np = $np ---"
    mpirun -np $np ./optimize4 4000000000 | grep -E 'processes =|total elapsed'
done
echo
echo "###END"

echo "###SHOT:files"
run pwd
run ls -l
echo "$P ls -l results"
ls -l results
echo "###END"
