#ifndef __PTSS_CONFIG
#define __PTSS_CONFIG

// Global configuration
#define NPH     5
#define M       32
// #define D       216
#define LLIM    3
#define ULIM    M

/* Benchmarks (ID) */
#define BENCH_LACE_DFS                        400
#define BENCH_LACE_CILKSORT                   401
#define BENCH_LACE_FIB                        402
#define BENCH_LACE_PI                         403
#define BENCH_LACE_QUEENS                     404
#define BENCH_PRSC_BLACKSCHOLES               405
#define BENCH_PRSC_STREAMCLUSTER              406
#define BENCH_PRSC_CANNEAL                    407
#define BENCH_PRSC_FLUIDANIMATE               408
#define BENCH_PRSC_BODYTRACK                  409
#define BENCH_PRSC_DEDUP                      410
// #define BENCH_PRSC_SWAPTIONS                  0x411

#endif