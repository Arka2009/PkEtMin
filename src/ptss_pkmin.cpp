#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <math.h>
#include <chrono>
#include <exception>
#include <stdexcept>
#include <queue>
#include <boost/math/distributions/normal.hpp>
#include <sys/time.h>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
#include "ptss_pkmin.hpp"
#include "ptss_dse.hpp"

using namespace std;

/* Evaluate the relative difference of two objective functions */
double eval_rel(const alloc2_t &pt1, const alloc2_t &pt2) {
    double f1 = compute_pkpower(pt1);
    double f2 = compute_pkpower(pt2);

    if (f1 == 0 && f2 == 0)
        return 0.0;
    else if (f1 != 0) {
        double rel = (abs(f1-f2))/abs(f1);
        return rel;
    }
    else {
        double rel = (abs(f1-f2))/abs(f2);
        return rel;
    }
    return 0.0;
}

double ptss_DSE_hrt::ptss_pkmin() {
    alloc2_t tmp = this->cvx_point;
    alloc2_t opt_dggd;
    if (tmp.size() < NPH) {
        throw invalid_argument( "Invalid Initial Point for DGGD\n" );
    }
    unsigned int ph;
    bool stop = false; /* Stopping criteria */

    //  cout << "iter-" << i++ << ",point:" << tmp << ",power:" << compute_pkpower(tmp) << endl;
     do { /* Check after balancing */
        // cout << "iter-" << i++ << ",point:" << tmp << ",power:" << compute_pkpower(tmp) << endl;
        opt_dggd = tmp;
        ph = compute_bottleneck(tmp);

        /* Decrement the number of cores (if possible) */
        if (tmp[ph].alloc > LLIM[this->bench[ph]])
            tmp[ph].alloc--;
        
        stop = (compute_execution_time(tmp) > this->deadline) || (eval_rel(tmp,opt_dggd)<1e-3);
        if (!stop) continue;
        
        /* Balance the allocation (Will not change the power consumption) */
        balance_out(tmp);
        stop = (compute_execution_time(tmp) > this->deadline) || (eval_rel(tmp,opt_dggd)<1e-3);
    } while (!stop);

    cout << "DGGD Point : " << opt_dggd << "\n";
    cout << "Exec Time:" << compute_execution_time(opt_dggd) << ",Power:" << compute_pkpower(opt_dggd) << "\n";

    this->dggd_point = opt_dggd;
    return compute_pkpower(opt_dggd);
}