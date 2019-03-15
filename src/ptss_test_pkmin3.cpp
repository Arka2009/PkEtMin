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
#include <rapidcsv.h>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
#include "ptss_pkmin.hpp"

using namespace std;
#define NSAMPLES 100


/**
 * Parse a CSV file containing the workloads
 * in the following format
 * Deadline,benchid1,benchid2,...
 */
int main(int argc, char **argv) {
    /* Setup the Input/Out File Names */
    char ifileName[BUFSIZ];
    char ofileName[BUFSIZ];
    sprintf(ifileName,"/home/amaity/Dropbox/workspace/islped_bonmip/workloads-exp1/wkld_%d.csv",NPH);
    sprintf(ofileName,"/home/amaity/Dropbox/workspace/islped_bonmip/workloads-exp1/wkld_%d_cpp.out.csv",NPH);

    string idata(ifileName);
    rapidcsv::Document inDoc(idata,rapidcsv::LabelParams(-1,-1));
    string odata(ofileName);
    ofstream ofs;
    ofs.open(ofileName,std::ofstream::out);

    bool done = false;
    unsigned int idx = 0;
    while (!done) {
        ptss_DSE_hrt obj(inDoc,ofs,9.94,idx++,done);
    }
    ofs.close();
}
