#!/bin/bash
# one line per measurement : version / np / n / pi(n) / marking time / total time
printf "%-42s %12s %10s %10s\n" "configuration" "pi(n)" "marking" "total"
awk '/^### VERSION=/  {v=$0}
     /primes are less/ {p=$1}
     /marking .sieve./ {m=$(NF-1)}
     /total elapsed/   {t=$(NF-1)
                        printf "%-42s %12s %10s %10s\n", substr(v,5,60), p, m, t}' results/final_4e9.txt
