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


int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Enter The deadline, Peak Power Cap" << endl;
        exit(EXIT_FAILURE);
    }
    double deadline = atof(argv[1]);
    double pkp_cap  = atof(argv[2]);
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    for (int i = 0; i < 100; i++) {
        ptss_DSE_hrt obj(deadline,pkp_cap);
        obj.display();
    }
    gettimeofday(&t2,NULL);
    double elapsed  = t2.tv_sec-t1.tv_sec;
}