#ifndef __PTSS_CONFIG
#define __PTSS_CONFIG

// Global configuration
#include "ptss_config_nph.hpp"

#define M       16
// #define         DEBUG1
#define ULIM    M
extern unsigned int LLIM[];
typedef unsigned long long int ptss_int_t;

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