/* ============================================================================
 *  optimize2.cpp  --  MPI Sieve of Eratosthenes, OPTIMISATION 2
 *                     (even numbers removed + broadcast eliminated)
 * ----------------------------------------------------------------------------
 *  Course  : Parallel Programming -- MPI Experiment 2
 *  Student : 2023080146011
 * ----------------------------------------------------------------------------
 *  Base of this version : optimize1.cpp
 *
 *  OPTIMISATION 2 -- "eliminate the broadcast / each process finds its own
 *                     sieving primes"
 *    * In base.cpp / optimize1.cpp every iteration of the marking loop needs a
 *      MPI_Bcast of the next sieving prime.  The number of broadcasts is the
 *      number of primes <= sqrt(n) (for n = 4*10^9 that is 6'300 messages,
 *      each of them a synchronisation of all processes).
 *    * Instead of letting process 0 search the primes and send them around,
 *      EVERY process computes the list of all odd sieving primes 3..sqrt(n)
 *      itself, in a small private array.  sqrt(n) is tiny (6.3*10^4 for
 *      n = 4*10^9), so this replicated computation costs almost nothing.
 *    * Consequence: the main marking loop contains NO communication at all,
 *      only a single MPI_Reduce at the very end.  Every process also becomes
 *      independent of process 0 (better scalability, no synchronisation
 *      points, no "too many processes" restriction any more).
 *
 *  Compile :  mpic++ -O3 -o optimize2 optimize2.cpp
 *  Run     :  mpirun -np 8 ./optimize2 100000000
 * ============================================================================ */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))

/* smallest odd multiple of the odd number k that is >= start (see optimize1) */
static inline long long first_odd_multiple (long long k, long long start)
{
    long long m = (start + k - 1) / k;
    if ((m & 1LL) == 0) m++;
    return m * k;
}

int main (int argc, char *argv[])
{
    long long n;
    long long count_odd, low_value, high_value, size;
    char     *marked, *small;        /* `small` = local sieve for 3..sqrt(n)   */
    long long small_limit, nprimes, k, j, t;
    long long *prime_list;           /* all odd primes <= sqrt(n)              */
    long long i, first;
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
    count_odd  = (n - 1) / 2;
    low_value  = 3 + 2 * BLOCK_LOW (id, p, count_odd);
    high_value = 3 + 2 * BLOCK_HIGH(id, p, count_odd);
    size       = BLOCK_SIZE (id, p, count_odd);
    (void) high_value;

    /* ============ OPTIMISATION 2 : private sieving primes ================== */
    /* small limit: all sieving primes are <= sqrt(n)                          */
    small_limit = (long long) sqrt ((double) n) + 1;

    /* tiny sequential sieve (replicated on every process, no communication)  */
    small = (char *) calloc ((size_t) small_limit + 1, 1);
    if (small == NULL) {
        printf ("Cannot allocate enough memory\n"); MPI_Finalize (); exit (1);
    }
    for (k = 2; k * k <= small_limit; k++)
        if (!small[k])
            for (j = k * k; j <= small_limit; j += k) small[j] = 1;

    /* collect the odd sieving primes in a compact array                       */
    nprimes = 0;
    for (k = 3; k <= small_limit; k += 2) if (!small[k]) nprimes++;
    prime_list = (long long *) malloc ((size_t) (nprimes ? nprimes : 1) *
                                       sizeof (long long));
    if (prime_list == NULL) {
        printf ("Cannot allocate enough memory\n"); MPI_Finalize (); exit (1);
    }
    for (k = 3, t = 0; k <= small_limit; k += 2)
        if (!small[k]) prime_list[t++] = k;
    free (small);

    /* ---------------- the local block of candidates ------------------------ */
    marked = (char *) malloc (size > 0 ? (size_t) size : 1);
    if (marked == NULL) {
        printf ("Cannot allocate enough memory\n"); MPI_Finalize (); exit (1);
    }
    for (i = 0; i < size; i++) marked[i] = 0;

    /* ======================= MARKING PHASE (parallel) ====================== */
    t_mark_start = MPI_Wtime ();
    for (t = 0; t < nprimes; t++) {              /* no MPI call in this loop   */
        k = prime_list[t];
        first = first_odd_multiple (k, (k * k > low_value) ? k * k : low_value);
        first = (first - low_value) / 2;
        for (i = first; i < size; i += k) marked[i] = 1;
    }
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
        printf (" Parallel Sieve of Eratosthenes -- OPTIMIZE 2\n");
        printf (" optimisations: even numbers removed\n");
        printf ("                BROADCAST ELIMINATED (private sieving primes)\n");
        printf (" ratio        : block decomposition (method 2)\n");
        printf (" commun.      : only the final MPI_Reduce (no MPI_Bcast)\n");
        printf (" student : 2023080146011      n = %lld      processes = %d\n", n, p);
        printf ("-----------------------------------------------------------------\n");
        printf (" %lld primes are less than or equal to %lld\n", global_count, n);
        printf (" largest prime found             : %lld\n", global_last_prime);
        printf (" sieving primes per process      : %lld\n", nprimes);
        printf (" marking (sieve) time    [max]   : %10.6f s\n", t_mark_max);
        printf (" total elapsed time      [max]   : %10.6f s\n", t_total_max);
        printf ("=================================================================\n");
        fflush (stdout);
    }

    free (marked);
    free (prime_list);
    MPI_Finalize ();
    return 0;
}
