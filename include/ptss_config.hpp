#ifndef __PTSS_CONFIG
#define __PTSS_CONFIG

// Global configuration
#include "ptss_config_nph.hpp"

#define M       3

// #define         DEBUG1
#define ULIM    M
#define TOT2     27
// #define TOT2     1048576    /* Size of the search space NPH = 5 */
// #define TOT2     16777216    /* Size of the search space NPH = 6 */
// #define TOT2     4294967296    /* Size of the search space NPH = 8 */

extern unsigned int LLIM[];
extern double AET[];
extern double BET[];
extern double AP[];
extern double BP[];

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
#define USEBGL /* USE BGL to construct a rule graph */
#define USEBITSET
#ifdef USEBGL
#define CONNECT_RULE1
#define CONNECT_RULE2
#define CONNECT_RULE3
#define CONNECT_RULE4
#endif

/* Profile Timings */
#define TIMING1 /* Search Space Construction Timing */
#define ENABLE_BRUTE_FORCE /* Enable Brute Force Search of Entire Space */
//#define DEBUGCHILDREN
#endif