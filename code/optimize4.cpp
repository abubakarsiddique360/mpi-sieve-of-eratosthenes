/* ============================================================================
 *  optimize4.cpp  --  MPI Sieve of Eratosthenes, FINAL VERSION
 *                     (all optimisations + bit packing, the version that is
 *                      tested with:  mpirun -np 8 ./optimize4 4000000000 )
 * ----------------------------------------------------------------------------
 *  Course  : Parallel Programming -- MPI Experiment 2
 *  Student : 2023080146011
 * ----------------------------------------------------------------------------
 *  Base of this version : optimize3.cpp
 *
 *  This final version contains ALL optimisations of the four steps plus two
 *  additional (own) improvements:
 *
 *   (1) EVEN NUMBERS REMOVED
 *       only the odd numbers 3,5,7,...,n are stored, local index j <-> value
 *       (low_value + 2*j) ; the only even prime 2 is added to the result by
 *       process 0.  -> half of the candidates, half of the marking work.
 *
 *   (2) BROADCAST ELIMINATED
 *       every process computes the list of all odd sieving primes 3..sqrt(n)
 *       with a small private sequential sieve, so the marking loop needs no
 *       MPI_Bcast at all.  Only the final MPI_Reduce is a collective call.
 *
 *   (3) CACHE OPTIMISATION
 *       the local block is processed tile by tile (TILE = 2^18 candidates =
 *       32 KB of bits = L1 cache size) and inside a tile the loops are
 *       interchanged: all primes mark their multiples while the tile is in the
 *       cache, instead of streaming the whole block once per prime.
 *
 *   (4) BIT PACKING (extra optimisation)
 *       one BIT per odd number instead of one byte.  The array of 250 million
 *       odd candidates per process shrinks from 250 MB to 31 MB, therefore the
 *       memory traffic (which dominates the running time of the sieve) is
 *       reduced by a factor of 8 and much larger n values become possible
 *       ("frees storage for larger values of n").
 *
 *   (5) FAST COUNTING (extra optimisation)
 *       the number of unmarked numbers is computed with __builtin_popcount
 *       (one instruction for 8 candidates) instead of a byte-by-byte test;
 *       calloc() is used instead of malloc+memset for the initialisation.
 *
 *  Compile :  mpic++ -O3 -o optimize4 optimize4.cpp
 *  Run     :  mpirun -np 8 ./optimize4 4000000000
 * ============================================================================ */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))

/* tile size in candidates: 2^18 odd numbers = 32 KB of bits -> L1 cache      */
#define TILE  (1LL << 18)

/* bit i of the bitset represents the odd number (low_value + 2*i)            */
#define SETBIT(i)     (bits[(i) >> 3] |= (unsigned char) (1u << ((i) & 7)))
#define GETBIT(i)     ((bits[(i) >> 3] >> ((i) & 7)) & 1u)

/* smallest odd multiple of the odd number k that is >= start                 */
static inline long long first_odd_multiple (long long k, long long start)
{
    long long m = (start + k - 1) / k;      /* m = ceil(start / k)             */
    if ((m & 1LL) == 0) m++;                /* the quotient must be odd        */
    return m * k;
}

int main (int argc, char *argv[])
{
    long long n;
    long long count_odd, low_value, high_value, size;
    unsigned char *bits;             /* one bit per odd number -> 8x smaller   */
    char     *small;
    long long small_limit, nprimes, k, j, t, b, nbytes;
    long long *prime_list;
    long long i, t0, t1, tile_low, tile_high, first;
    long long ones = 0;
    long long count = 0, global_count = 0;  /* +1 : the prime 2 is not stored   */
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

    /* --------------- list of the ODD numbers 3, 5, ..., n ------------------ */
    count_odd  = (n - 1) / 2;
    low_value  = 3 + 2 * BLOCK_LOW (id, p, count_odd);
    high_value = 3 + 2 * BLOCK_HIGH(id, p, count_odd);
    size       = BLOCK_SIZE (id, p, count_odd);
    (void) high_value;

    /* ========== (2) private sieving primes 3..sqrt(n), no broadcast ======== */
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

    /* ========== (4) the bitset of this process (calloc -> zeroed) ========= */
    nbytes = (size + 7) / 8;
    bits = (unsigned char *) calloc ((size_t) (nbytes ? nbytes : 1), 1);
    if (bits == NULL) {
        printf ("Cannot allocate enough memory\n"); MPI_Finalize (); exit (1);
    }

    /* ==== (1)(3)(4) MARKING PHASE : odd only + bitpacked + cache tiling ==== */
    t_mark_start = MPI_Wtime ();
    for (t0 = 0; t0 < size; t0 += TILE) {        /* one L1-sized tile         */
        t1 = t0 + TILE;
        if (t1 > size) t1 = size;
        tile_low  = low_value + 2 * t0;          /* smallest odd value        */
        tile_high = low_value + 2 * (t1 - 1);    /* largest  odd value        */

        for (t = 0; t < nprimes; t++) {          /* all primes inside the tile */
            k = prime_list[t];
            if (k * k > tile_high) break;        /* nothing left to mark      */
            first = first_odd_multiple (k, (k * k > tile_low) ? k * k : tile_low);
            first = (first - low_value) / 2;     /* value -> local index      */
            for (i = first; i < t1; i += k) SETBIT (i);
        }
    }
    t_mark = MPI_Wtime () - t_mark_start;

    /* ========== (5) count the primes with popcount (8 candidates/op) ======= */
    if (!id) count = 1;                  /* the prime 2 is counted by process 0
                                            only (it is not in the odd list)    */
    for (b = 0; b < nbytes; b++)
        ones += __builtin_popcount ((unsigned) bits[b]);   /* marked bits       */
    /* `ones` counts the composites.  The bits behind `size` were never set
       (calloc + bounded loops), so the number of primes in this block is
       simply the number of remaining bits.                                  */
    count += size - ones;

    /* largest prime of this block: last unset bit                             */
    for (i = size - 1; i >= 0; i--)
        if (!GETBIT (i)) { last_prime = low_value + 2 * i; break; }

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
        printf (" Parallel Sieve of Eratosthenes -- OPTIMIZE 4 (FINAL)\n");
        printf (" optimisations: 1 even numbers removed (only odd candidates)\n");
        printf ("                2 broadcast eliminated (private sieving primes)\n");
        printf ("                3 cache optimisation (%lld byte tiles)\n", TILE / 8);
        printf ("                4 bit packing: 1 bit per candidate (8x less memory)\n");
        printf ("                5 popcount counting + calloc initialisation\n");
        printf (" student : 2023080146011      n = %lld      processes = %d\n", n, p);
        printf ("-----------------------------------------------------------------\n");
        printf (" %lld primes are less than or equal to %lld\n", global_count, n);
        printf (" largest prime found             : %lld\n", global_last_prime);
        printf (" sieving primes per process      : %lld\n", nprimes);
        printf (" bits (candidates) per process   : %lld  (%.1f MB)\n",
                size, (double) size / 8.0 / 1048576.0);
        printf (" marking (sieve) time    [max]   : %10.6f s\n", t_mark_max);
        printf (" total elapsed time      [max]   : %10.6f s\n", t_total_max);
        printf ("=================================================================\n");
        fflush (stdout);
    }

    free (bits);
    free (prime_list);
    MPI_Finalize ();
    return 0;
}
