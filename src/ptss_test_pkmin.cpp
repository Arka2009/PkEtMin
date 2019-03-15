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
#define NSAMPLES 100


int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Enter The deadline, Peak Power Cap" << endl;
        exit(EXIT_FAILURE);
    }
    double deadline   = atof(argv[1]);
    //vector<double> pkp_cap  = {8.5,9.22,9.94,10.66,11.38,12.11,12.83,13.55,14.27,15.0};
    vector<double> pkp_cap  = {9,10,11,12,13,14,15,16,17};

    // struct timeval t1, t2;
    // gettimeofday(&t1,NULL);
    int j = NSAMPLES,i = 0;
    for (i = 0; i < pkp_cap.size();i++) {
        j = NSAMPLES;
        while (j > 0) {
            cout << "("<<i<<","<<j<<") Iteration"<<endl;
            ptss_DSE_hrt obj(deadline,pkp_cap[i]);
            obj.display();
            j--;
        }
    }
    // gettimeofday(&t2,NULL);
    // double elapsed  = t2.tv_sec-t1.tv_sec;
}