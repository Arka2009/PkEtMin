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

int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Enter The deadline" << endl;
        exit(EXIT_FAILURE);
    }
    int deadline = atoi(argv[1]);
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    ptss_DSE_hrt obj(deadline);
    // cout << obj;
    // obj.display();
    // obj.explore();
    double pkp_opt = obj.evaluate_all();
    alloc2_t init, dggd;
    int i;
    for (i = 0; i < NPH; i++) {
        phase_t p(0x400+i,M);
        init.push_back(p);
    }
    double pkp_opt2 = ptss_pkmin(init,deadline,&dggd);
    gettimeofday(&t2,NULL);
    double elapsed  = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);

    /* Percent Difference */
    double diff = (pkp_opt2-pkp_opt)*100/pkp_opt;
    cout << "deadline:" << deadline << ",diff:"<<diff<<"\n\n";
    // cout << elapsed << "us\n";
}