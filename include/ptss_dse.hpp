#ifndef __PTSS_DSE
#define __PTSS_DSE

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
#include "ptss_config.hpp"

using namespace std;



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


// ostream& operator<<(ostream& os, const vector<double>& pht);
ostream& operator<<(ostream& os, const phase_t& pht);
ostream& operator<<(ostream& os, const alloc_t& vi);
ostream& operator<<(ostream& os, const all_alloc_t& vvi);
ostream& operator<<(ostream& os, const alloc2_t& vi);
ostream& operator<<(ostream& os, const all_alloc2_t& vvi);

double compute_risk(const alloc_t &x);
double compute_execution_time(const alloc_t &x);
double compute_execution_time(const alloc2_t &x);
double estimate_exec_time(const int x, const int bench);
// int inv_estimate_exec_time(const double y, const int bench);
double estimate_power(const int x, const int bench);
int inv_estimate_power(const double y, const int bench);
double compute_estimated_util(const alloc_t &x);
double compute_pkpower(const alloc2_t &x);
set<unsigned int> compute_bottleneck(const alloc2_t &x);
set<unsigned int> compute_maxgrad(const alloc2_t &x);
void balance_out(alloc2_t &x);
void construct_alloc_2();

/* Utility functions */
void construct_alloc(all_alloc2_t &vvi, \
                     const vector<int> &bench, \
                     double deadline, \
                     int ph,\
                     alloc2_t &opt_point,\
                     double &opt_pkp_power,\
                     double &opt_exec_time);
unsigned int gen_bench_id();

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

class ptss_DSE_hrt {
    private :
        /*
         * The ptss_DSE partitions
         * the spaces into feasible
         * and discarded space
         */
        all_alloc2_t search_space;
        double deadline;
        double pkp_cap;      /* Maximum peak power */

        /* Allocation Points */
        alloc2_t opt_point;  /* Oracle Computed pkp min point */
        alloc2_t ext_point;  /* Point of Highest Power Consumption, also serves as an init_point in some algorithms */
        alloc2_t ext_point2; /* Minimum allocation point*/
        alloc2_t cvx_point;  /* Optimal Point found using continuous relaxation and then dicretization */
        alloc2_t dggd_point; /* Point Found using discrete greedy gradient descent */
        alloc2_t opt_point2;  /* Oracle Computed et min point */
        alloc2_t dggd_point2; /* Dual of dggd_point */

        /* Corresponding Objective Function (and constraint) Value */
        // double opt_pkp_power;
        // double opt_exec_time;
        double cvx_pkp_min;   /* Continuous domain minimum obj function value */

        /* Parameter values for power and execution time characteristics */
        vector<int> bench; /* Benchmark Composition */
        vector<double> a_et;
        vector<double> b_et;
        vector<double> a_p;
        vector<double> b_p;

        double compute_cvx();
        double bench_create();

    public :
        ptss_DSE_hrt(double,double);
        ptss_DSE_hrt() : ptss_DSE_hrt(2000,15) {};
        void construct_alloc2();
        

        // Evaluate all points in the Design Space, and return the optimal value
        // double oracle();

        // Display
        void display();

        // Getters
        // alloc2_t& get_init_point();
        // double get_opt_pkp_power();
        // double get_opt_exec_time();

        // Membership testing
        bool contains_point(const alloc2_t&);

        /* The DGGD algorithm (Discrete Greedy Gradient Descent also referred to as
        pkmin in the paper) */
        double ptss_pkmin();
        double ptss_etmin(); /* Dual of previous problem */

};


/* Objective Function Constraints */
#endif