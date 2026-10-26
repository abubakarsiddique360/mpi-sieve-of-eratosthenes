/* ============================================================================
 *  optimize1.cpp  --  MPI Sieve of Eratosthenes, OPTIMISATION 1
 *                     (optimisation: remove the even numbers)
 * ----------------------------------------------------------------------------
 *  Course  : Parallel Programming -- MPI Experiment 2
 *  Student : 2023080146011
 * ----------------------------------------------------------------------------
 *  Base of this version : base.cpp (block decomposition + MPI_Bcast + reduce)
 *
 *  OPTIMISATION 1 -- "delete the even integers"
 *    * 2 is the only even prime, therefore no even number > 2 can be a prime.
 *    * The list only stores the odd numbers 3, 5, 7, ..., n :
 *          local index j  <->  value  (low_value + 2*j)
 *    * => the sieve array (and the number of marking operations) is halved,
 *         twice as large n values fit into the same amount of memory.
 *    * The marks are still performed by stepping through the array with
 *      stride `prime` (in index space), because the multiples of the odd
 *      prime k that are odd are k*k, k*k+2k, k*k+4k, ...
 *    * The primes are still broadcast by process 0 (this is the next
 *      optimisation, see optimize2.cpp).
 *
 *  Compile :  mpic++ -O3 -o optimize1 optimize1.cpp
 *  Run     :  mpirun -np 8 ./optimize1 100000000
 * ============================================================================ */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))

/* ---------------------------------------------------------------------------
 * smallest odd multiple of the odd number k that is >= start
 * (only odd multiples of k are inside the odd-only list, they are 2k apart in
 *  value space ; the quotient must be odd, otherwise the product is even)
 * ------------------------------------------------------------------------ */
static inline long long first_odd_multiple (long long k, long long start)
{
    long long m = (start + k - 1) / k;      /* m = ceil(start / k)             */
    if ((m & 1LL) == 0) m++;                /* make the quotient odd           */
    return m * k;
}

int main (int argc, char *argv[])
{
    long long n;
    long long count_odd;         /* # odd numbers in the list (3..n)           */
    long long low_value, high_value, size;
    long long proc0_high;
    char     *marked;            /* 1 bit of information per ODD number        */
    long long i, index = 0, prime, first, q;
    long long count = 0, global_count = 0;    /* +1 : the prime 2 is not stored */
    long long last_prime = 0, global_last_prime = 0;

    int       id, p;
    double    t_mark, t_total, t_mark_max, t_total_max, t_start, t_mark_start;

    MPI_Init (&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);
    t_start = MPI_Wtime ();
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    if (argc != 2) {
        if (!id) printf ("Command line: %s <n>\n", argv[0]);
        MPI_Finalize ();
        exit (1);
    }
    n = atoll (argv[1]);
    if (n < 2) {
        if (!id) printf ("n must be >= 2\n");
        MPI_Finalize ();
        exit (1);
    }

    /* ---------------- list of the ODD numbers 3, 5, ..., n ----------------- */
    count_odd  = (n - 1) / 2;                       /* 3,5,...,n (or n-1)    */
    low_value  = 3 + 2 * BLOCK_LOW (id, p, count_odd);
    high_value = 3 + 2 * BLOCK_HIGH(id, p, count_odd);
    size       = BLOCK_SIZE (id, p, count_odd);

    /* process 0 must contain all odd sieving primes ( 3 .. sqrt(n) )        */
    if (p <= count_odd) {
        proc0_high = 3 + 2 * (count_odd / p - 1);
        if (proc0_high < (long long) sqrt ((double) n)) {
            if (!id) printf ("Too many processes\n");
            MPI_Finalize ();
            exit (1);
        }
    }

    marked = (char *) malloc (size > 0 ? (size_t) size : 1);
    if (marked == NULL) {
        printf ("Cannot allocate enough memory\n");
        MPI_Finalize ();
        exit (1);
    }
    for (i = 0; i < size; i++) marked[i] = 0;

    /* ======================= MARKING PHASE (parallel) ====================== */
    t_mark_start = MPI_Wtime ();
    if (!id) index = 0;            /* local index 0 <-> value 3 = first prime  */
    prime = 3;                     /* 2 is even -> not in the list at all      */
    do {
        /* first odd multiple of `prime` that is >= max(prime*prime, low_value) */
        first = first_odd_multiple (prime, (prime*prime > low_value)
                                             ? prime*prime : low_value);
        first = (first - low_value) / 2;          /* value -> local index      */
        /* odd multiples of an odd prime are 2*prime apart in value space,
           i.e. exactly `prime` apart in index space                            */
        for (i = first; i < size; i += prime) marked[i] = 1;

        /* process 0 searches the next (still unmarked) sieving prime          */
        if (!id) {
            do { ++index; } while (index < size && marked[index]);
            prime = 2 * index + 3;                /* local index -> value      */
        }
        MPI_Bcast (&prime, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    t_mark = MPI_Wtime () - t_mark_start;

    /* ---------------- count the primes of this block ----------------------- */
    if (!id) count = 1;                 /* the prime 2 is not in the odd list,
                                           it is counted by process 0 only     */
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;

    for (i = size - 1; i >= 0; i--)
        if (!marked[i]) { last_prime = low_value + 2 * i; break; }

    MPI_Reduce (&count,      &global_count,      1, MPI_LONG_LONG, MPI_SUM, 0,
                MPI_COMM_WORLD);
    MPI_Reduce (&last_prime, &global_last_prime, 1, MPI_LONG_LONG, MPI_MAX, 0,
                MPI_COMM_WORLD);

    MPI_Barrier (MPI_COMM_WORLD);
    t_total = MPI_Wtime () - t_start;
    MPI_Reduce (&t_mark,  &t_mark_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce (&t_total, &t_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!id) {
        printf ("=================================================================\n");
        printf (" Parallel Sieve of Eratosthenes -- OPTIMIZE 1\n");
        printf (" optimisation : EVEN NUMBERS REMOVED (only odd candidates kept)\n");
        printf (" ratio        : block decomposition (method 2)\n");
        printf (" commun.      : MPI_Bcast of the sieving prime + MPI_Reduce\n");
        printf (" student : 2023080146011      n = %lld      processes = %d\n", n, p);
        printf ("-----------------------------------------------------------------\n");
        printf (" %lld primes are less than or equal to %lld\n", global_count, n);
        printf (" largest prime found             : %lld\n", global_last_prime);
        printf (" marking (sieve) time    [max]   : %10.6f s\n", t_mark_max);
        printf (" total elapsed time      [max]   : %10.6f s\n", t_total_max);
        printf ("=================================================================\n");
        fflush (stdout);
    }

    free (marked);
    MPI_Finalize ();
    return 0;
}
