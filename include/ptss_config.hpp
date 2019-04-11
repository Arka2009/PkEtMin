#ifndef __PTSS_CONFIG
#define __PTSS_CONFIG

// Global configuration
#include "ptss_config_nph.hpp"

#define M       8
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

/* Rule Graph Connections */
#define CONNECT_RULE1
#define CONNECT_RULE2
#define CONNECT_RULE3
// #define CONNECT_RULE4

/* Profile Timings */
#define TIMING1 /* Search Space Construction Timing */
#endif