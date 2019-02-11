#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <math.h>
#include <chrono>
#include <exception>
#include <queue>
#include <boost/math/distributions/normal.hpp>
#include <sys/time.h>

using namespace std;

// Global configuration
#define NPH 5
#define M   16
#define D   216

typedef vector<int> alloc_t;
typedef set<alloc_t> all_alloc_t;

class phase_t {
    public :
        unsigned int bench_id;
        unsigned int alloc;

        bool operator ==(const phase_t &b) const;
        bool operator !=(const phase_t &b) const;
        bool operator <(const phase_t &b) const;
        bool operator >(const phase_t &b) const;
        bool operator <=(const phase_t &b) const;
        bool operator >=(const phase_t &b) const;

        phase_t();
        phase_t(unsigned int bench_id,unsigned int alloc);
};

typedef vector<phase_t> alloc2_t;
typedef set<alloc2_t> all_alloc2_t;

bool is_eq(const alloc2_t &a, const alloc2_t &b);
bool is_gt(const alloc2_t &a, const alloc2_t &b);
bool is_lt(const alloc2_t &a, const alloc2_t &b);
bool is_incomparable(const alloc2_t &a, const alloc2_t &b);
bool lex_comp(const alloc2_t &a, const alloc2_t &b);


ostream& operator<<(ostream& os, const phase_t& pht);
ostream& operator<<(ostream& os, const alloc_t& vi);
ostream& operator<<(ostream& os, const all_alloc_t& vvi);
ostream& operator<<(ostream& os, const alloc2_t& vi);
ostream& operator<<(ostream& os, const all_alloc2_t& vvi);

double compute_risk(const alloc_t &x);
double compute_execution_time(const alloc_t &x);
double compute_execution_time(const alloc2_t &x);
double estimate_exec_time(const int x, const int bench);
int inv_estimate_exec_time(const double y, const int bench);
double estimate_power(const int x, const int bench);
int inv_estimate_power(const double y, const int bench);
double compute_estimated_util(const alloc_t &x);
double compute_pkpower(const alloc2_t &x);
unsigned int compute_bottleneck(const alloc2_t &x);
void balance_out(alloc2_t &x);
void construct_alloc(all_alloc_t &vvi, int ph);

class ptss_DSE {
    private :
        /*
         * The ptss_DSE partitions
         * the spaces into feasible
         * and discarded space
         */
        all_alloc_t search_space;
        all_alloc_t feasible_space;
        all_alloc_t discarded_space;

        double dmr;

        /* Optimal Point, just for testing */
        alloc_t opt_point;
        double opt_util;
        double opt_risk;

        // Lower and Upper "limits" points for DMR
        alloc_t lower;
        alloc_t upper;

        // Get the initial allocated point
        void init_point(double dmr);

        bool is_initialized;

    public :
        ptss_DSE();
        ptss_DSE(double dmr);

        // Evaluate all points in the Design Space
        void evaluate_all();

        // Comparison based elimination
        // void eliminate_points_rule12();
        void explore();

        // Display
        void display();
};

/* Benchmarks */
#define BENCH_LACE_DFS           0x400
#define BENCH_LACE_CILKSORT      0x401
#define BENCH_LACE_FIB           0x402
#define BENCH_LACE_PI            0x403
#define BENCH_LACE_QUEENS        0x404

class ptss_DSE_hrt {
    private :
        /*
         * The ptss_DSE partitions
         * the spaces into feasible
         * and discarded space
         */
        all_alloc2_t search_space;
        double deadline;

        /* Optimal Point, just for testing */
        alloc2_t opt_point;
        alloc2_t ext_point; /* Point of Highest Power Consumption */
        /* Corresponding Objective Function Value */
        double pkp_power;
        double exec_time;

    public :
        ptss_DSE_hrt();
        ptss_DSE_hrt(double);

        // Evaluate all points in the Design Space, and return the optimal value
        double evaluate_all();

        // Display
        void display();
};