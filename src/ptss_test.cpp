#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <cmath>
#include <chrono>
#include <exception>
#include <queue>
#include <boost/math/distributions/normal.hpp>
#include <sys/time.h>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
#include "ptss_pkmin.hpp"

using namespace std;

void test_gen_bench_id() {
    int k = 10;
    while(k-- > 0) {
        cout << "Bench:"<<gen_bench_id()<<endl;
    }
}

void test_compute_execution_time() {
    /* Create an allocation */
    alloc2_t a;
    double s = 0.0;
    
    for (int i = 0; i < 6; i++) {
        a.push_back(phase_t(5+i,1));
        // a.push_back(phase_t(409,1));
        // a.push_back(phase_t(408,2));
        // a.push_back(phase_t(406,1));
        // a.push_back(phase_t(410,1));

        double d  = compute_execution_time(a);
        s += d;
        cout << "Execution Time "<< a <<", : " << d << endl;

        a.clear();
    }
    cout << "Total Time : " << s << endl;
}

void test_power() {
    ptss_int_t x = (ptss_int_t)pow(16.0,NPH);
    cout << x << endl;
    // construct_alloc_2();
}

int main() {
    //srand(time(NULL));
    // test_gen_bench_id();
    // test_compute_execution_time();
    test_power();
}
