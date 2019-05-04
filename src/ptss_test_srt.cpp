#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <cmath>
#include <chrono>
#include <exception>
#include <queue>
#include <rapidcsv.h>
#include <boost/math/distributions/normal.hpp>
#include <sys/time.h>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
#include "ptss_pkmin.hpp"

using namespace std;

int main() {

    char ifileName[BUFSIZ];
    sprintf(ifileName,"wkld_%d.csv",NPH);
    string idata(ifileName);
    rapidcsv::Document inDoc(idata,rapidcsv::LabelParams(-1,-1));
    int idx = 0;
    bool done = false;

    /* Read Each Row */
    // while (!done) {
        ptss_DSE_srt obj(inDoc,idx++,done);
        #ifdef USEBGL
        obj.printRuleGraph();
        obj.random_walk3(56771);
        #else
        obj.random_walk4(2766);
        #endif
        cout << "--------------------------------"<<endl;
    // }

    // ptss_int_t tmp = 8194;
    // vector<int> bench = {1,0,2,1,3};
    // alloc2_t vi2 = invhashalloc_2(tmp,bench);
    // bool t1;
    // double t2;
    // boost::tie(t1,t2) = checkFeasibilityAndPKP(vi2,400,0.11);
    // cout << t1 << "," << t2 << endl;
    // boost::tie(t1,t2) = checkFeasibilityAndPKP(vi2,400,0.11);
    // cout << t1 << "," << t2 << endl;
    // boost::tie(t1,t2) = checkFeasibilityAndPKP(vi2,400,0.11);
    // cout << t1 << "," << t2 << endl;
    return 0;
}
