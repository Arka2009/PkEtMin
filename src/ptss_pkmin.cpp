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
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
#include "ptss_pkmin.hpp"

using namespace std;

double ptss_pkmin(alloc2_t init, unsigned int deadline, alloc2_t *opt_dggd) {
    alloc2_t tmp = init;
    alloc2_t tmp2;
    unsigned int ph;

     do { /* Check after balancing */
        *opt_dggd = tmp;
        ph = compute_bottleneck(tmp);
        tmp[ph].alloc--;
        if (compute_execution_time(tmp) <= deadline) continue;
        /* Balance the allocation (Will not change the power consumption) */
        balance_out(tmp);
    } while (compute_execution_time(tmp) <= deadline);

    // /**/
    // ph = compute_bottleneck(tmp);
    // tmp[ph].alloc++;
    // balance_out(tmp);


    cout << "DGGD Point : " << *opt_dggd << "\n";
    cout << "Exec Time:" << compute_execution_time(*opt_dggd) << ",Power:" << compute_pkpower(*opt_dggd) << "\n";
    return compute_pkpower(*opt_dggd);
}