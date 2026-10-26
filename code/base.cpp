/* ============================================================================
 *  base.cpp  --  Parallel Sieve of Eratosthenes   (BASELINE, no optimisation)
 * ----------------------------------------------------------------------------
 *  Course  : Parallel Programming -- MPI Experiment 2
 *  Student : 2023080146011
 * ----------------------------------------------------------------------------
 *  This version is a literal implementation of the algorithm of Chapter 4
 *  "Message-Passing Programming - The Sieve of Eratosthenes" :
 *
 *    1. every process creates its own share of the list 2 .. n
 *       (block decomposition method #2, blocks differ by at most one element)
 *    2. every process marks the multiples of the current sieve prime k inside
 *       its own block, starting at k*k  (all smaller multiples are already
 *       marked by smaller primes)
 *    3. process 0 finds the next sieve prime (smallest unmarked number greater
 *       than the current one) and broadcasts it with MPI_Bcast
 *    4. when k*k > n the loop stops -- the unmarked numbers are the primes
 *    5. MPI_Reduce(MPI_SUM) collects the local prime counts on process 0
 *
 *  Compile :  mpic++ -O3 -o base base.cpp
 *  Run     :  mpirun -np 8 ./base 100000000
 * ============================================================================ */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ---- block decomposition (Chapter 4, decomposition method #2) --------------
 *   first element (index) owned by process id : BLOCK_LOW(id,p,n)
 *   last  element (index) owned by process id : BLOCK_HIGH(id,p,n)
 *   number of elements owned by process id    : BLOCK_SIZE(id,p,n)
 *   The first n%p blocks contain one element more than the remaining ones.
 * --------------------------------------------------------------------------*/
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))

int main (int argc, char *argv[])
{
    long long n;                 /* sieve limit: find all primes <= n          */
    long long low_value;         /* smallest value handled by this process     */
    long long high_value;        /* largest  value handled by this process     */
    long long size;              /* # candidates (elements) of this block      */
    long long proc0_size;        /* size of the block of process 0             */
    char     *marked;            /* marked[i]==1  => value (low_value+i) is composite */
    long long i, index = 0;
    long long prime;             /* current sieving prime                      */
    long long first;             /* local index of the first multiple to mark  */
    long long count = 0, global_count = 0;
    long long last_prime = 0, global_last_prime = 0;

    int       id, p;
    double    t_mark, t_total, t_mark_max, t_total_max, t_start, t_mark_start;

    /* ------------------------- start MPI ---------------------------------- */
    MPI_Init (&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);          /* make the timing meaningful      */
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

    /* --------------- decompose the list 2,3,...,n -------------------------- */
    /* n-1 elements are distributed over the p processes; local index i of
       process id represents the value  low_value + i                       */
    low_value  = 2 + BLOCK_LOW (id, p, n - 1);
    high_value = 2 + BLOCK_HIGH (id, p, n - 1);
    size       = BLOCK_SIZE (id, p, n - 1);

    /* the largest sieving prime is sqrt(n); process 0 has to own all of them */
    proc0_size = (n - 1) / p;
    if ((2 + proc0_size) < (long long) sqrt ((double) n)) {
        if (!id) printf ("Too many processes\n");
        MPI_Finalize ();
        exit (1);
    }

    /* --------------- allocate and initialise the local block --------------- */
    marked = (char *) malloc (size > 0 ? (size_t) size : 1);
    if (marked == NULL) {
        printf ("Cannot allocate enough memory\n");
        MPI_Finalize ();
        exit (1);
    }
    for (i = 0; i < size; i++) marked[i] = 0;

    /* ======================= MARKING PHASE (parallel) ====================== */
    t_mark_start = MPI_Wtime ();
    if (!id) index = 0;            /* process 0 searches the next prime       */
    prime = 2;                     /* first sieving prime                     */
    do {
        /* first local index to mark: multiples of `prime` start at prime*prime
           (everything below is already marked by smaller primes), and the
           search starts at the first multiple >= low_value                  */
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            long long r = low_value % prime;
            first = r ? prime - r : 0;
        }
        for (i = first; i < size; i += prime) marked[i] = 1;

        /* --- only process 0 searches the next prime ---------------------- */
        if (!id) {
            do { ++index; } while (index < size && marked[index]);
            prime = index + 2;     /* value of the found local index          */
        }
        /* --- everybody receives the next sieving prime ------------------- */
        MPI_Bcast (&prime, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    t_mark = MPI_Wtime () - t_mark_start;

    /* ---------------- count the primes of this block ----------------------- */
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;

    /* largest prime of this block (scan backwards, normally 0 or 1 steps) */
    for (i = size - 1; i >= 0; i--)
        if (!marked[i]) { last_prime = low_value + i; break; }

    MPI_Reduce (&count,      &global_count,      1, MPI_LONG_LONG, MPI_SUM, 0,
                MPI_COMM_WORLD);
    MPI_Reduce (&last_prime, &global_last_prime, 1, MPI_LONG_LONG, MPI_MAX, 0,
                MPI_COMM_WORLD);

    MPI_Barrier (MPI_COMM_WORLD);
    t_total = MPI_Wtime () - t_start;

    /* the maximum over all processes is the real parallel wall-clock time */
    MPI_Reduce (&t_mark,  &t_mark_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce (&t_total, &t_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!id) {
        printf ("=================================================================\n");
        printf (" Parallel Sieve of Eratosthenes -- BASE version (no optimisation)\n");
        printf (" ratio   : block decomposition (method 2)\n");
        printf (" commun. : MPI_Bcast of the sieving prime + MPI_Reduce\n");
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
