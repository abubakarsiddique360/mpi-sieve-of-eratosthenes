/* ============================================================================
 *  optimize3.cpp  --  MPI Sieve of Eratosthenes, OPTIMISATION 3
 *                     (even numbers removed + broadcast eliminated
 *                      + cache optimisation)
 * ----------------------------------------------------------------------------
 *  Course  : Parallel Programming -- MPI Experiment 2
 *  Student : 2023080146011
 * ----------------------------------------------------------------------------
 *  Base of this version : optimize2.cpp
 *
 *  OPTIMISATION 3 -- "reorganise the loops, increase the cache hit rate"
 *    * In optimize1/optimize2 the marking is done prime by prime over the WHOLE
 *      local block:
 *            for every prime k:  for i = ... ; i += k   marked[i] = 1
 *      For the small primes (3,5,7,...) the block (hundreds of MB) is streamed
 *      from / to the main memory once per prime, the same cache lines are
 *      loaded again and again -> the algorithm is limited by the memory
 *      bandwidth, not by the CPU.
 *    * Cache optimisation (loop interchange): the local block is split into
 *      tiles of TILE candidates that fit into the L1/L2 cache.  While one tile
 *      is in the cache, ALL sieving primes mark their multiples inside this
 *      tile, then the next tile is loaded:
 *            for every tile:  for every prime k with k*k <= tile_high:
 *                                 for i inside the tile; i += k   marked[i] = 1
 *      Every cache line is now loaded exactly once and is re-used by all
 *      primes  -> far less memory traffic, much better cache hit rate.
 *    * Correctness: inside a tile the primes are still processed in increasing
 *      order and each prime starts at max(k*k, tile_low); a number that is not
 *      marked by any smaller prime can therefore still be started at k*k, and
 *      primes with k*k > tile_high cannot mark anything unmarked in this tile
 *      (all their multiples are composed of smaller primes), so they are
 *      skipped -> this is exactly the sequential sieving order, only the
 *      traversal order of the memory is different.
 *
 *  Compile :  mpic++ -O3 -o optimize3 optimize3.cpp
 *  Run     :  mpirun -np 8 ./optimize3 100000000
 * ============================================================================ */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))

/* tile size in candidates: 2^18 odd numbers = 256 KB -> fits into the L2
   cache of the machine (Intel Xeon Gold 6342: 1.25 MB L2 per core)          */
#define TILE  (1LL << 18)

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
    char     *marked, *small;
    long long small_limit, nprimes, k, j, t;
    long long *prime_list;
    long long i, t0, t1, tile_low, tile_high, first;
    long long count = 0, global_count = 0;
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

    /* ============ private sieving primes (optimisation 2) ================== */
    small_limit = (long long) sqrt ((double) n) + 1;
    small = (char *) calloc ((size_t) small_limit + 1, 1);
    if (small == NULL) {
        printf ("Cannot allocate enough memory\n"); MPI_Finalize (); exit (1);
    }
    for (k = 2; k * k <= small_limit; k++)
        if (!small[k])
            for (j = k * k; j <= small_limit; j += k) small[j] = 1;

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

    marked = (char *) malloc (size > 0 ? (size_t) size : 1);
    if (marked == NULL) {
        printf ("Cannot allocate enough memory\n"); MPI_Finalize (); exit (1);
    }
    for (i = 0; i < size; i++) marked[i] = 0;

    /* =========== MARKING PHASE : tiled / cache-friendly loop =============== */
    t_mark_start = MPI_Wtime ();
    for (t0 = 0; t0 < size; t0 += TILE) {          /* one cache-sized tile   */
        t1 = t0 + TILE;
        if (t1 > size) t1 = size;
        tile_low  = low_value + 2 * t0;            /* smallest odd value     */
        tile_high = low_value + 2 * (t1 - 1);      /* largest  odd value     */

        for (t = 0; t < nprimes; t++) {            /* all primes in the tile */
            k = prime_list[t];
            if (k * k > tile_high) break;          /* k and all larger primes
                                                      cannot mark anything    */
            first = first_odd_multiple (k, (k * k > tile_low) ? k * k : tile_low);
            first = (first - low_value) / 2;
            for (i = first; i < t1; i += k) marked[i] = 1;
        }
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
        printf (" Parallel Sieve of Eratosthenes -- OPTIMIZE 3\n");
        printf (" optimisations: even numbers removed\n");
        printf ("                broadcast eliminated (private sieving primes)\n");
        printf ("                CACHE OPTIMISATION (%lld byte tiles, loop interchange)\n",
                TILE);
        printf (" ratio        : block decomposition (method 2)\n");
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
    free (prime_list);
    MPI_Finalize ();
    return 0;
}
