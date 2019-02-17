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

/* Nowhere else to go (lowest allocation possible) */
bool hit_the_wall(const alloc2_t &pt1) {
    for (unsigned int i = 0; i < pt1.size(); i++) {
        if (LLIM[pt1[i].bench_id] < pt1[i].alloc)
            return false;
    }
    return true;
}

double ptss_DSE_hrt::ptss_pkmin() {
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    alloc2_t tmp = this->cvx_point;
    alloc2_t opt_dggd;
    if (tmp.size() < NPH) {
        throw invalid_argument( "Invalid Initial Point for DGGD\n" );
    }
    set<unsigned int> phs;
    unsigned int ph;
    bool stop = false; /* Stopping criteria */
    unsigned int iter = 0;

    /* Test for negative execution time */
    // if(compute_execution_time(tmp) < 0) {
    //     cout <<"\n\n\n"<< "DGGD Init Value Illegal" <<tmp<<"et" << compute_execution_time(tmp)<<"\n\n\n";
    //     for (auto it = phs.begin(); it != phs.end(); it++) {
    //         cout << *it << endl;
    //     }
    //     throw invalid_argument("--------------------Illegal Value of tmp*****");
    // }

    //  cout << "iter-" << iter << ",point:" << tmp << ",power:" << compute_pkpower(tmp) << endl;
     do { /* Check after balancing */
        // cout << "iter-" << iter << ",point:" << tmp << ",power:" << compute_pkpower(tmp) << endl;
        iter++;
        opt_dggd = tmp;
        phs = compute_bottleneck(tmp);

         /* Test for negative execution time */
        // cout << "DGGD-Pre@"<<iter<<":"<<tmp<<",et:" << compute_execution_time(tmp)<<"\n";

        // cout << "DGGD:PHS SIZE : " <<phs.size()<<endl;
        for (auto it = phs.begin(); it != phs.end(); it++) {
            ph = *it;
            // cout << "DGGD"<<iter<<" : Peak Phase:"<<ph<<",power:"<<compute_pkpower(tmp)<<endl;
            /* Decrement the number of cores (if possible) (Guaranteed to reduce power if successful) */
            if (tmp[ph].alloc > LLIM[this->bench[ph]]) {
                // cout << "DGGD@"<<iter<<",ph:"<<ph<<",allocated:"<<tmp[ph].alloc<<",LLIM:"<<LLIM[this->bench[ph]]<< endl;
                tmp[ph].alloc--;
            } else if (tmp[ph].alloc < LLIM[this->bench[ph]]) {
                cout << "\n\n\nph:"<<ph<<" is already below the lower limit\n\n\n"<<endl;
                throw invalid_argument("already below the lower limit");
            }
            // cout << "DGGD"<<iter<<" : (New) Peak Phase:"<<ph<<",power:"<<compute_pkpower(tmp)<<endl;
        }
        
        /* Test for negative execution time */
        if(compute_execution_time(tmp) < 0) {
            cout << "DGGD@"<<iter<<":"<<tmp<<",et:" << compute_execution_time(tmp)<<"\n";
            throw invalid_argument("Negative Value of ET not allowed");
        }
        
        // cout << "iter-" << iter << ",point:" << tmp << ",power:" << compute_pkpower(tmp) << "\n\n";
        stop = (compute_execution_time(tmp) > this->deadline) || (eval_rel(tmp,opt_dggd) < 1e-3);
        if (!stop) continue;
        
        /* Balance the allocation (Will not change the power consumption) */
        balance_out(tmp);
        if(compute_execution_time(tmp) < 0) {
            cout << "DGGD-Balancing-Out@"<<iter<<":"<<tmp<<",et:" << compute_execution_time(tmp)<<"\n";
            throw invalid_argument("Negative Value of ET not allowed");
        }
        stop = (compute_execution_time(tmp) > this->deadline) || (eval_rel(tmp,opt_dggd) < 1e-3);
    } while (!stop);

    gettimeofday(&t2,NULL);
    long long diff = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
    cout << "dbcks{iter|elapsed time},"<<iter<<","<<diff<<endl;
    // cout << "DGGD Point : "<<opt_dggd<<"\n\n";
    this->dggd_point = opt_dggd;
    return compute_pkpower(opt_dggd);
}