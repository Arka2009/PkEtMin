#ifndef __PTSS_DSE
#define __PTSS_DSE

#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <cmath>
#include <string>
#include <chrono>
#include <exception>
#include <queue>
#include <stack>
#include <bitset>
#include <boost/math/distributions/normal.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <sys/time.h>
#include <rapidcsv.h>
#include "ptss_config.hpp"

using namespace std;
using namespace boost;

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
enum  RulesEdgeType {RULE1,RULE2,RULE3,RULE4,INVALID};

#ifdef USEBGL
typedef adjacency_list<vecS,vecS,bidirectionalS,alloc2_t,RulesEdgeType> SRTDSEGraph;
typedef boost::graph_traits<SRTDSEGraph>::vertex_descriptor SRTDSEGraphVertex;
typedef boost::graph_traits<SRTDSEGraph>::edge_descriptor SRTDSEGraphEdge;
typedef boost::graph_traits<SRTDSEGraph>::out_edge_iterator SRTDSEGraphOutEdgeIterator;
typedef boost::graph_traits<SRTDSEGraph>::in_edge_iterator SRTDSEGraphInEdgeIterator;
typedef boost::graph_traits<SRTDSEGraph>::vertex_iterator SRTDSEGraphVertexIterator;
typedef boost::graph_traits<SRTDSEGraph>::degree_size_type SRTDSEGraphDegSizeType;
#endif
// typedef std::map<SRTDSEGraphVertex,int> SRTDSEComponent;

bool is_eq(const alloc2_t &a, const alloc2_t &b);
bool is_gt(const alloc2_t &a, const alloc2_t &b);
bool is_lt(const alloc2_t &a, const alloc2_t &b);
bool is_incomparable(const alloc2_t &a, const alloc2_t &b);
bool lex_comp(const alloc2_t &a, const alloc2_t &b);


// ostream& operator<<(ostream& os, const vector<double>& pht);
ostream& operator<<(ostream& os, const phase_t& pht);
ostream& operator<<(ostream& os, const alloc2_t& vi);
ostream& operator<<(ostream& os, const all_alloc2_t& vvi);
ostream& operator<<(ostream& os, const RulesEdgeType &rl);

double compute_risk(const alloc2_t &x, int deadline);
double compute_execution_time(const alloc2_t &x);
double compute_execution_time_std(const alloc2_t &x);
double estimate_exec_time(const int x, const int bench);
// int inv_estimate_exec_time(const double y, const int bench);
double estimate_power(const int x, const int bench);
int inv_estimate_power(const double y, const int bench);
double compute_pkpower(const alloc2_t &x);
set<unsigned int> compute_bottleneck(const alloc2_t &x);
set<unsigned int> compute_maxgrad(const alloc2_t &x);
void balance_out(alloc2_t &x);
ptss_int_t hashalloc_2(const alloc2_t &x);
alloc2_t invhashalloc_2(ptss_int_t id,const vector<int> &bench);
double eval_rel(const alloc2_t &pt1, const alloc2_t &pt2);
bool eval_manhtn(const alloc2_t &pt1, const alloc2_t &pt2);
unsigned int gen_bench_id();
std::set<ptss_int_t> union_sets(std::set<ptss_int_t> A, std::set<ptss_int_t> B);
std::pair<bool,double> checkFeasibilityAndPKP(const alloc2_t& v, double deadline, double dmr);
void does_allocation_exists(ptss_int_t allc, int rule, const vector<int> bench);

class ptss_DSE_srt {
    private :
        /* Benchmark Composition */
        vector<int> bench;

        #ifdef USEBGL
        /* use BGL */        
        SRTDSEGraph gSS;                      /* Directed Graph */
        // SRTDSEComponent gSSComponents;
        list<SRTDSEGraphVertex> validVertices;
        // list<SRTDSEGraphVertex> isolatedVertices; /* Vertices with outdegress zero */
        bool graphBuilt;
        ptss_int_t validOutDegrees(const SRTDSEGraphVertex&); /* Compute the number of valid Out Degrees */
        #endif

        /* Hold all the explored set of points */
        #ifdef USEBITSET
        bitset<TOT2> exploredSet; 
        #else
        vector<bool> exploredSet;
        #endif

        /* Constraints */
        double dmr;
        double deadline;

        /* Optimal Point and its corresponding Objective Function Value  */
        alloc2_t opt_point;
        double opt_pkp;
        double opt_risk;
        alloc2_t bktrack_point;

        /* Final point obtained by Random walker */
        alloc2_t rwalkerFinal;

        /* Get the initial allocated point */
        void construct_alloc2();

        /* Comparison operators */
        // bool operator<(const SRTDSEGraphVertex&, const SRTDSEGraphVertex&);
        // bool operator<=(const SRTDSEGraphVertex&, const SRTDSEGraphVertex&);

        /* Apply rules to get child points */
        set<ptss_int_t> applyRule1(const ptss_int_t pt1);
        set<ptss_int_t> applyRule2(const ptss_int_t pt1);
        set<ptss_int_t> applyRule3(const ptss_int_t pt1);
        set<ptss_int_t> applyRule4(const ptss_int_t pt1);

    public :
        ptss_DSE_srt();
        ptss_DSE_srt(double dmr);
        ptss_DSE_srt(double deadline, double dmr);
        ptss_DSE_srt(rapidcsv::Document &inDoc,unsigned int idx, bool &done); /* Read the workload from a CSV file */

        /* Brute Force Exploration */
        void exploreBruteForce();

        /* Obtain the DMR */
        double geDMR() {return this->dmr;}
        double getDeadline() {return this->deadline;}

        /* Display/Visualize the search space */
        void display();

        #ifdef USEBGL
        /* Randomly walk the graph */
        void random_walk();
        void random_walk2(); /*  Augmented with random walk and BnB */
        void random_walk3(ptss_int_t start_offset);
        ptss_int_t backtrackAll(const SRTDSEGraphVertex& start_node);
        #endif
        ptss_int_t backtrackAll(ptss_int_t start_node);
        void random_walk4(ptss_int_t start_offset);

        #ifdef USEBGL
        /* Visaulize the graph */
        void printRuleGraph();
        #endif
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
        double bench_create2(const vector<int>&);

    public :
        ptss_DSE_hrt(double,double);
        ptss_DSE_hrt() : ptss_DSE_hrt(2000,15) {};
        ptss_DSE_hrt(rapidcsv::Document &,ostream &,double,unsigned int,bool &);    /* Read the workload from a CSV file */
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

template <typename I>
I random_element(I begin, I end);
/* Objective Function Constraints */
#endif