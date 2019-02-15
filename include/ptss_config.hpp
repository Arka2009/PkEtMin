#ifndef __PTSS_CONFIG
#define __PTSS_CONFIG

// Global configuration
#define NPH     5
#define M       16
// #define LLIM    3
#define ULIM    M

/* Benchmarks (ID) */
#define BENCH_LACE_DFS                        0
#define BENCH_LACE_CILKSORT                   1
#define BENCH_LACE_FIB                        2
#define BENCH_LACE_PI                         3
#define BENCH_LACE_QUEENS                     4
#define BENCH_PRSC_BLACKSCHOLES               5
#define BENCH_PRSC_STREAMCLUSTER              6
#define BENCH_PRSC_CANNEAL                    7
#define BENCH_PRSC_FLUIDANIMATE               8
#define BENCH_PRSC_BODYTRACK                  9
#define BENCH_PRSC_DEDUP                      10
// #define BENCH_PRSC_SWAPTIONS               0x411

#endif