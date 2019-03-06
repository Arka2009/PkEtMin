#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <math.h>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <exception>
#include <queue>
#include <cstdlib>
#include <boost/math/distributions/normal.hpp>
#include <sys/time.h>
#include <ecotools/roi_hooks.h>
#include <nlopt.hpp>
#include <algorithm>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
#include "ptss_nlopt.hpp"

using namespace std;

// AET 
static double AET[] = {\
0.0052912772458832795,\
0.004254273510905619,\
0.005266330309541947,\
0.003688517486111831,\
0.003959440926636103,\
0.009683975849832148,\
0.0024794976904501873,\
2.153404950919569e-05,\
0.0021638291897166048,\
0.003533669349872262,\
0.0002555000090127519};

static double BET[] = {\
0.002118298041182201,\
-0.0008304322052644017,\
0.0021923143697295602,\
0.0011143885733635485,\
0.001000559321893281,\
-0.002893833328927599,\
-0.0007496878381248601,\
0.0001677868522003634,\
-0.001202365402456029,\
-0.0027172901605883966,\
-7.198705194815577e-05};

static double AP[] = {\
2.9905877034358053,\
2.0279769439421345,\
3.3427689873417723,\
2.947167721518987,\
3.176307692307692,\
1.5027499999999998,\
2.299642857142857,\
0.6275357142857142,\
2.130357142857143,\
1.7288351648351648,\
0.8776666666666665};

static double BP[] = {\
1.9919258589511664,\
0.3829859855334483,\
2.140996835443037,\
1.758180379746836,\
1.7091208791208778,\
4.669250000000002,\
3.7145476190476217,\
1.810845238095241,\
0.834642857142855,\
3.2989230769230744,\
3.5970000000000013};

/* Lower Limit of core allocation for different phases */
unsigned int LLIM[] = {\
1,\
1,\
1,\
1,\
1,\
2,\
2,\
2,\
2,\
3,\
2};

bool phase_t::operator ==(const phase_t &b) const {
    return ((this->bench_id == b.bench_id) && (this->alloc == b.alloc));
}

bool phase_t::operator !=(const phase_t &b) const {
    return !((this->bench_id == b.bench_id) && (this->alloc == b.alloc));
}

bool phase_t::operator <(const phase_t &b) const {
    return (this->alloc < b.alloc);
}

bool phase_t::operator >(const phase_t &b) const {
    return (this->alloc > b.alloc);
}

bool phase_t::operator <=(const phase_t &b) const {
    return (this->alloc <= b.alloc);
}

bool phase_t::operator >=(const phase_t &b) const {
    return (this->alloc >= b.alloc);
}

ostream& operator<<(ostream& os, const phase_t& pht) {
    os << "{bench_id:" << pht.bench_id <<",alloc:" << pht.alloc << "}";
    return os;
}

phase_t::phase_t(unsigned int bench_id,unsigned int alloc) : bench_id{bench_id}, alloc{alloc} {}
phase_t::phase_t() : bench_id{0}, alloc{0} {}

double estimate_exec_time(const int x, const int bench) {
    if (bench > BENCH_PRSC_DEDUP || bench < BENCH_LACE_DFS) {
        throw invalid_argument( "Unrecognized Benchmark" );
    }
    double a = AET[bench];
    double b = BET[bench];
    double y;
    if (bench <= BENCH_LACE_QUEENS)
        y = 1/(a*x + b);
    else
        y = 1/(a*log(x) + b);
    return y;
}

double estimate_power(const int x, const int bench) {
    if (bench > BENCH_PRSC_DEDUP || bench < BENCH_LACE_DFS) {
        throw invalid_argument( "Unrecognized Benchmark" );
    }
    double a = AP[bench];
    double b = BP[bench];
    double y;
    y = a*x + b;
    return y;
}

int inv_estimate_power(const double y, const int bench) {
    if (bench > BENCH_PRSC_DEDUP || bench < BENCH_LACE_DFS) {
        throw invalid_argument( "Unrecognized Benchmark" );
    }
    double a = AP[bench];
    double b = BP[bench];
    double x;
    x = (1/a)*(y-b);

    /* The number of cores must be within [LLIM,ULIM] */
    if (x < LLIM[bench])
        x = LLIM[bench];
    if (x > ULIM)
        x = ULIM;
    return x;
}

/* 
 * Extrapolation of et can make the values often negative 
 * These need to be discarded
 */
double compute_execution_time(const alloc2_t &x) {
    double sum = 0.0;
    double et  = 0.0;
    for(auto jt = x.begin(); jt != x.end(); jt++) {
        et = estimate_exec_time(jt->alloc,jt->bench_id);
        if (et >= 0.0)
            sum += et;
        else
            return -HUGE_VAL;
        
    }
    return sum;
}


/* Compute the Peak Power */
double compute_pkpower(const alloc2_t &x) {
    double pkp = 0.0;
    for(auto jt = x.begin(); jt != x.end(); jt++)
        pkp = max(pkp,estimate_power(jt->alloc,jt->bench_id));
    return pkp;
}

/* Returns all the bottleneck phases (Those with highest power consumption) */
set<unsigned int> compute_bottleneck(const alloc2_t &x) {
    set<unsigned int> phs;
    double etp[NPH];

    /* Estimate the power of all phases */
    for (unsigned int i = 0; i < x.size();i++)
        etp[i] = estimate_power(x[i].alloc,x[i].bench_id);
    const int N = sizeof(etp) / sizeof(double);
    unsigned int least_idx = distance(etp, max_element(etp, etp + N));
    double maxp = etp[least_idx];

    /* Find all indices with maximum power consumption */
    for(unsigned int i = 0; i < x.size();i++) {
        // cout << "maxp : "<<maxp<<",power estimate : "<<estimate_power(x[i].alloc,x[i].bench_id)<<endl;
        if (maxp == etp[i]) {
            phs.insert(i);
        }
    }
    return phs;
}

/* Compute the nbd allocation with maximum decrease in et */
set<unsigned int> compute_maxgrad(const alloc2_t &x) {
    set<unsigned int> phs;
    double etp[NPH];

    /* Estimate the execution time of all phases */
    for (unsigned int i = 0; i < x.size();i++) {
        alloc2_t tmp = x;
        tmp[i].alloc += 1;
        etp[i] = compute_execution_time(x) - compute_execution_time(tmp);

        // if (etp[i]) {
        //     throw invalid_argument("Negative diff in et unexpected\n");
        // }
        // cout << "ph("<<i<<")-diff:"<<etp[i]<<endl;
    }
    const int N = sizeof(etp) / sizeof(double);
    unsigned int least_idx = distance(etp, max_element(etp, etp + N));
    double maxp = etp[least_idx];

    /* Find all indices with max decline in execution time */
    for(unsigned int i = 0; i < x.size();i++) {
        // cout << "maxp : "<<maxp<<",power estimate : "<<estimate_power(x[i].alloc,x[i].bench_id)<<endl;
        if (maxp == etp[i]) {
            phs.insert(i);
        }
    }
    // int i = 0;
    // for (set<unsigned int>::iterator it = phs.begin(); it != phs.end(); it++) {
    //     cout << "phs("<<i++<<")-diff:"<<*it<<endl;
    // }
    return phs;
}


/* Try to balance out all the power consumption in all the phases*/
void balance_out(alloc2_t &x) {
    double pkp = compute_pkpower(x);
    for(auto jt = x.begin(); jt != x.end(); jt++) {
        jt->alloc = inv_estimate_power(pkp,jt->bench_id);
    }
}


// /* Display Routines */
std::ostream& operator<<(std::ostream& os, const alloc2_t& vi) {
    unsigned int i;
    os << "<";
    for (i = 0; i < vi.size(); i++) {
        os << vi[i] ;
        os << ",";
    }
    os << ">";
    return os;
}

/* Display for vector */
ostream& operator<<(ostream& os, const std::vector<double> pht) {
    os << "Vector<" ;
    for (unsigned int i = 0; i < pht.size(); i++)
        os << pht[i] << "|";
    os << ">";
    return os;
}

std::ostream& operator<<(std::ostream& os, const all_alloc2_t& vvi) { 
  os << "{\n";
  for(auto it = vvi.begin();
      it != vvi.end();
      it++) {
      os << "  " << *it << "\n";
  }
  os << "}";
  return os;
}


/* An non-recursive version of the above */
static inline long long calculate_time_diff_spec(struct timespec t2, struct timespec t1) {
    long long elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000000000LL + t2.tv_nsec - t1.tv_nsec;
    return elapsedTime;                                                         
} 

void ptss_DSE_hrt::construct_alloc2() {
    alloc2_t vi2;
    double et, pkp;
    ptss_int_t r = M, j = 0, dec, cnt;
    ptss_int_t TOT = (ptss_int_t)pow(M,NPH);

    /* Temporary "Global Varirbles "*/
    double deadline = this->deadline;
    // unsigned int bench[] = {10,9,10,10,7};
    double opt_pkp_power = HUGE_VALF;
    double opt_exec_time = HUGE_VALF;
    //struct timespec t1, t2;

    for (cnt = 0; cnt < TOT; cnt++) {
        //clock_gettime(CLOCK_REALTIME,&t1);
        dec = cnt;

        /* Create the point */
        for (j = 0; j < NPH; j++) {
            phase_t phinfo;
            phinfo.alloc = dec%r + 1;
            phinfo.bench_id = this->bench[j];
            vi2.push_back(phinfo);
            dec = dec/r;
        }
        
        /* Find optimal pkp point */
        et = compute_execution_time(vi2);
        if (et <= deadline && et >= 0) {
            pkp = compute_pkpower(vi2);
            if (pkp <= opt_pkp_power) {
                this->opt_point     = vi2;
                opt_pkp_power       = pkp;
            }
        }
        // vi2.clear();

        /* Find optimal et point */
        pkp = compute_pkpower(vi2);
        if (pkp <= this->pkp_cap && pkp >= 0) {
            et = compute_execution_time(vi2);
            if (et >= 0 && et <= opt_exec_time) {
                this->opt_point2    = vi2;
                opt_exec_time       = et;
            }
        }
        vi2.clear();
    }

    // cout << "Optimal Point : "<<opt_point<<endl;
}

unsigned int gen_bench_id() {
    unsigned int a = BENCH_PRSC_BLACKSCHOLES + (rand() % 6);
    return a;
}


double ptss_DSE_hrt::compute_cvx() {
    nlopt::opt opt(nlopt::LD_MMA, NPH+1);
    // nlopt::opt opt(nlopt::LN_SBPLX, NPH+1);

    /* Set the Box Constraints */
    std::vector<double> lb(NPH+1);
    std::vector<double> ub(NPH+1);
    for (int i = 0; i < NPH; i++) {
        lb[i] = LLIM[this->bench[i]];
        ub[i] = ULIM;
    }
    lb[NPH] = 0.0;
    ub[NPH] = HUGE_VAL;
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    /* Set the Parameters for inequality constraints */
    // cout << "CVX Deadline " << this->deadline << endl;
    vector<ptss_constraint_param> param;
    param.push_back(ptss_constraint_param(this->a_et,this->b_et,this->deadline,-1));
    for (int i = 0; i < NPH; i++) {
        param.push_back(ptss_constraint_param(this->a_p,this->b_p,this->deadline,i));
    }
    opt.add_inequality_constraint(ptss_constraint_exectime, &param[0], 1e-8);
    for (int i = 1; i <= NPH; i++) {
        opt.add_inequality_constraint(ptss_constraint_power, &param[i], 1e-8);
    }

    /* Set the objective function */
    opt.set_min_objective(ptss_func_pkp, NULL);

    /* Stopping Criteria and Initial Point*/
    opt.set_ftol_rel(1e-4);
    vector<double> x(NPH+1);
    vector<phase_t> c(NPH);
    for (int i = 0; i < NPH; i++) {
        c[i].alloc    = M;
        c[i].bench_id = this->bench[i];
        x[i] = M;
    }
    x[NPH] = compute_pkpower(c);
    double minf = 0.0;
    
    try{
        //nlopt::result result = 
        opt.optimize(x, minf);
        // cout << "CVX-Opt f(" << x << ") = "<< std::setprecision(10) << minf << std::endl;
        // std::cout << "found minimum after " << count <<" evaluations\n";

        /* Update the cvx point */
        this->cvx_point.clear();
        for (int i = 0; i < NPH; i++) {
            phase_t tmp;
            tmp.alloc    = ceil(x[i]);
            tmp.bench_id = this->bench[i];
            this->cvx_point.push_back(tmp);
        }

        /* Display the continuous domain point */
        // cout <<"CVX-Opt continuous point x = <";
        // for (int i = 0; i < NPH+1; i++)
        //     cout <<x[i]<<",";
        // cout<<">"<<endl;
        // cout << "CVX-Opt minf = " << minf << endl;
        // cout << "CVX-Opt Relaxed Point = "<<this->cvx_point<<"\n";
        // cout << "CVX-Opt Power Consumption = "<<compute_pkpower(this->cvx_point)<<"\n";
        // cout << "CVX-Opt Execution Time = "<<compute_execution_time(this->cvx_point)<<"\n\n";
        
        return minf;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    return minf;
}

double ptss_DSE_hrt::bench_create() {
    /* Create a mix of benchmarks */
    int a;
    // unsigned int benchid[] = {10,9,10,10,7};
    for (int i = 0; i < NPH; i++) {
        a = gen_bench_id();
        // a = benchid[i];
        // cout << "Bench id Generated "<<a<<endl;
        this->bench.push_back(a);
        a_et.push_back(AET[a]);
        b_et.push_back(BET[a]);
        a_p.push_back(AP[a]);
        b_p.push_back(BP[a]);
    }
    
    double minf = this->compute_cvx();
    this->cvx_pkp_min = minf;
    return minf;
}

bool ptss_DSE_hrt::contains_point(const alloc2_t& a) {
    bool ret = this->search_space.find(a) != this->search_space.end();
    if (ret) {
        cout << "point:"<<a<<",pkp:"<<compute_pkpower(a)<<",et:"<<compute_execution_time(a)<<endl;
    } else {
        cout << "point not found : " << a << endl;
    }
    return ret;
}


ptss_DSE_hrt::ptss_DSE_hrt(double deadline,double pkp_cap) {
    this->deadline = deadline;
    this->pkp_cap  = pkp_cap;
    this->bench_create();
    /* Create an extreme point */
    for (int i = 0; i < NPH; i++) {
        phase_t p(this->bench[i],ULIM);
        this->ext_point.push_back(p);

        phase_t p2(this->bench[i],LLIM[this->bench[i]]);
        this->ext_point2.push_back(p2);
    }
#ifndef SHUTDOWN_ORACLE
    this->construct_alloc2();
#else
    /* Do not compute the Oracle*/
    this->opt_point     = ext_point;
#endif

    /* Use a DGGD Algorithm */
    this->ptss_pkmin();
    this->ptss_etmin();
}


// a == b ?
bool is_eq(const alloc2_t &a, const alloc2_t &b) {
    return (a == b);
}

// a > b ?
bool is_gt(const alloc2_t &a, const alloc2_t &b) {
    bool ret = !(is_eq(a,b));
    int idx = 0;
    for(auto it = a.begin();
        it != a.end();
        it++) {
            if (it->alloc >= b[idx++].alloc)
                ret = ret && true;
            else
                ret = ret && false;
    }
    return ret;
}

// a < b ?
bool is_lt(const alloc2_t &a, const alloc2_t &b) {
    bool ret = !(is_eq(a,b));
    int idx = 0;
    for(auto it = a.begin();
        it != a.end();
        it++) {
            if (it->alloc <= b[idx++].alloc)
                ret = ret && true;
            else
                ret = ret && false;
    }
    return ret;
}

// a and b are not comparable
bool is_incomparable(const alloc2_t &a, const alloc2_t &b) {
    bool ret = (!is_gt(a,b)) && \
               (!is_lt(a,b)) && \
               (!is_eq(a,b));
    return ret;
}

// lexicographic comparison of two allocations ( a < b ? )
bool lex_comp(const alloc2_t &a, const alloc2_t &b) {
    auto it2 = b.begin();
    for (auto it = a.begin(); it != a.end(); it++, it2++) {
        if (it2->alloc < it->alloc)
            return false;
    }
    return true;
}


void ptss_DSE_hrt::display() {
    if (this->cvx_point.size() < NPH) {
        throw invalid_argument("Did not create CVX Opt Point");
    }
    if (this->dggd_point.size() < NPH) {
        throw invalid_argument("Did not create DGGD Opt Point");
    }
    if (this->opt_point.size() < NPH) {
        cout << "NO FEASIBLE PKP POINT FOUND" << endl;
    }
    if (this->opt_point2.size() < NPH) {
        cout << "NO FEASIBLE ET POINT FOUND" << endl;
    }

#ifdef DEBUG1
    cout << "Optimal PKP Point:" << this->opt_point<< ",et:" \
    << compute_execution_time(this->opt_point) \
    <<",pkp:"<<compute_pkpower(this->opt_point)<<endl;
    
    cout << "Optimal ET Point:" << this->opt_point2<< ",et:" \
    << compute_execution_time(this->opt_point2) \
    <<",pkp:"<<compute_pkpower(this->opt_point2)<<endl;

    cout << "DGGD Point:" << this->dggd_point<<",et:" \
    << compute_execution_time(this->dggd_point) \
    <<",pkp:"<<compute_pkpower(this->dggd_point)<<endl;

    cout << "DGGD2 Point:" << this->dggd_point2<<",et:" \
    << compute_execution_time(this->dggd_point2) \
    <<",pkp:"<<compute_pkpower(this->dggd_point2)<<endl;
#endif
#ifndef SHUTDOWN_ORACLE
    // cout << "ufhew4r4{Oracle|CVX-Cont|CVX-Disc|DGGD},"<<endl;
    cout<<"ufhew4r4-pkpcap,"<<NPH<<","\
        <<this->pkp_cap<<","\
        <<this->deadline<<","\
        <<compute_pkpower(this->opt_point)<<","\
        <<compute_pkpower(this->ext_point)<<","\
        <<this->cvx_pkp_min<<","\
        <<compute_pkpower(this->cvx_point)<<","\
        <<compute_pkpower(this->dggd_point)<<","\
        <<compute_execution_time(this->opt_point2)<<","\
        <<compute_execution_time(this->dggd_point2)<<endl;
#else
    cout << "ufhew4r4{Worst|CVX-Cont|CVX-Disc|DGGD},";
#endif
    cout << "\n\n\n";
}
