#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <math.h>
#include <cmath>
#include <chrono>
#include <exception>
#include <queue>
#include <boost/math/distributions/normal.hpp>
#include <sys/time.h>
#include "ptss_dse.hpp"

using namespace std;

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
    double y;
    switch (bench) {
        case BENCH_LACE_CILKSORT : {
            y = 1/(0.004254273510905619*x - 0.0008304322052644017);
            break;
        }
        case BENCH_LACE_DFS : {
            y = 1/(0.0052912772458832795*x + 0.002118298041182201);
            break;
        }
        case BENCH_LACE_PI : {
            y = 1/(0.003688517486111831*x + 0.0011143885733635485);
            break;
        }
        case BENCH_LACE_QUEENS : {
            y = 1/(0.003959440926636103*x + 0.001000559321893281);
            break;
        }
        case BENCH_LACE_FIB : {
            y = 1/(0.005266330309541947*x + 0.0021923143697295602);
            break;
        }
        default : {
            cerr << "Unrecognized Benchmark";
            exit(EXIT_FAILURE);
        }
    }
    return y;
}

int inv_estimate_exec_time(const double y, const int bench) {
    double x;
    switch (bench) {
        case BENCH_LACE_CILKSORT : {
            x = ceil(((1/y) + 0.0008304322052644017)*(1/0.004254273510905619));
            break;
        }
        case BENCH_LACE_DFS : {
            x = ceil(((1/y) - 0.002118298041182201)*(1/0.0052912772458832795));
            break;
        }
        case BENCH_LACE_PI : {
            x = ceil(((1/y) - 0.0011143885733635485)*(1/0.003688517486111831));
            break;
        }
        case BENCH_LACE_QUEENS : {
            x = ceil(((1/y) - 0.001000559321893281)*(1/0.003959440926636103));
            break;
        }
        case BENCH_LACE_FIB : {
            x = ceil(((1/y) - 0.0021923143697295602)*(1/0.005266330309541947));
            break;
        }
        default : {
            cerr << "Unrecognized Benchmark";
            exit(EXIT_FAILURE);
        }
    }
    int x2 = (int)x;
    return x2;
}

double estimate_power(const int x, const int bench) {
    double y;
    switch (bench) {
        case BENCH_LACE_CILKSORT : {
            y = 2.0279769439421345*x + 0.3829859855334483;
            break;
        }
        case BENCH_LACE_DFS : {
            y = 2.9905877034358053*x + 1.9919258589511664;
            break;
        }
        case BENCH_LACE_PI : {
            y = 2.947167721518987*x + 1.758180379746836;
            break;
        }
        case BENCH_LACE_QUEENS : {
            y = 3.176307692307692*x + 1.7091208791208778;
            break;
        }
        case BENCH_LACE_FIB : {
            y = 3.3427689873417723*x + 2.140996835443037;
            break;
        }
        default : {
            cerr << "Unrecognized Benchmark";
            exit(EXIT_FAILURE);
        }
    }
    return y;
}

int inv_estimate_power(const double y, const int bench) {
    double x;
    switch (bench) {
        case BENCH_LACE_CILKSORT : {
            x = floor((y - 0.3829859855334483)*(1/2.0279769439421345));
            break;
        }
        case BENCH_LACE_DFS : {
            x = floor((y-1.9919258589511664)*(1/2.9905877034358053));
            break;
        }
        case BENCH_LACE_PI : {
            x = floor((y-1.758180379746836)*(1/2.947167721518987));
            break;
        }
        case BENCH_LACE_QUEENS : {
            x = floor((y-1.7091208791208778)*(1/3.176307692307692));
            break;
        }
        case BENCH_LACE_FIB : {
            x = floor((y-2.140996835443037)*(1/3.3427689873417723));
            break;
        }
        default : {
            cerr << "Unrecognized Benchmark";
            exit(EXIT_FAILURE);
        }
    }
    int x2  = (int)x;
    return x2;
}

/* Execution Time */
double compute_execution_time(const alloc2_t &x) {
    double sum = 0.0;
    for(auto jt = x.begin(); jt != x.end(); jt++) {
        sum += estimate_exec_time(jt->alloc,jt->bench_id);
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

/* Compute the bottleneck phase */
unsigned int compute_bottleneck(const alloc2_t &x) {
    unsigned int idx = 0, i;
    for(i = 0; i < x.size(); i++) {
        if (estimate_power(x[idx].alloc,x[idx].bench_id) <  estimate_power(x[i].alloc,x[i].bench_id))
            idx = i;
    }
    return idx;
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
//   os << "(|";
//   std::copy(vi.begin(), vi.end(), std::ostream_iterator<phase_t>(os, "|"));
//   os << ")";
    unsigned int i;
    os << "<";
    for (i = 0; i < vi.size(); i++) {
        os << vi[i] ;
        os << ",";
    }
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

void construct_alloc(all_alloc2_t &vvi, int ph) {
    alloc2_t vi2;
    static long unsigned int cnt = 0;
    long unsigned int r   = M;
    long unsigned int dec;
    long unsigned int j   = 0;

    if (ph == NPH) {
        //cout << "ph : " << ph << ", invoc : " << cnt << endl;
        dec = cnt++;

        /* Resolve cnt into M-radix number */
        for (j = 0; j < NPH; j++) {
            //cout << dec%r + 1 << ",";
            phase_t phinfo;
            phinfo.alloc = dec%r + 1;
            phinfo.bench_id = 0x400 + j;

            vi2.push_back(phinfo);
            dec = dec/r;
        }
        vvi.insert(vi2);
        //cout << "}\n";
        return;
    }
    for (int m = 1; m <= M; m++) {
        construct_alloc(vvi,ph+1);
    }
}

// Create All points
ptss_DSE_hrt::ptss_DSE_hrt() {
    construct_alloc(this->search_space,0);
    this->deadline = 200;
}

ptss_DSE_hrt::ptss_DSE_hrt(double deadline) {
    construct_alloc(this->search_space,0);
    this->deadline = deadline;
}

// Evaluate all points (Brute Force)
double ptss_DSE_hrt::evaluate_all() {
    this->pkp_power = 1e9;

    // cout << "Search space size " << search_space.size();
    for(auto it = search_space.begin();
        it != search_space.end();
        it++) {
        
        double et = compute_execution_time(*it);
        double pkp = compute_pkpower(*it);

        if (et <= this->deadline && pkp <= this->pkp_power) {
            this->opt_point = *it;
            this->pkp_power = pkp;
            this->exec_time = et;
        }
        // cout << *it << "," << risk << "," << util << "\n";
        // usleep(50);
    }
    cout << "Optimal Point : " << this->opt_point << "\n";
    cout << "Exec Time:" << this->exec_time << ",Power:" << this->pkp_power << "\n";
    return (this->pkp_power);
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
