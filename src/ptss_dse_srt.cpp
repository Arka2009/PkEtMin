#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <exception>
#include <queue>
#include <cstdlib>
#include <rapidcsv.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <algorithm>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
extern "C" {
#include <ecotools/roi_hooks.h>
#include <sys/time.h>
}

// Profiled Data
static double mua  [] = {269.568837,	141.985532, 99.508883, 79.006262, 66.454092, 58.649068, 52.804705, 48.960141, 45.73406, 43.455542, 41.386889, 40.046129, 38.695398, 37.750183, 36.840822, 36.41925, 35.715777, 35.350802, 34.767695, 34.493228, 34.033093, 33.944396, 33.676262, 33.871417, 33.654013, 33.737325, 33.825452, 34.021759, 34.223725, 34.479067, 34.694037, 35.157394}; 
static double stda [] = {0.676473946, 0.335842225, 0.208458629, 0.319806191, 0.313754681, 0.340024999, 0.340440891, 0.441634464, 0.424892928, 0.45393722, 0.470011702, 0.479161768, 0.496199557, 0.493329504, 0.524060111, 0.523248507, 0.564837145, 0.590188106, 0.584813646, 0.602546264, 0.597015913, 0.600979201, 0.59552246, 0.595876665, 0.587283577, 0.645952011, 0.666749578, 0.732030737, 0.836999403, 0.942281274, 1.046115194, 1.103043517}; 
using namespace std;
using namespace boost;
using namespace boost::math;

/* Convert an allocation to equivalent big integer */
ptss_int_t hashalloc_2(const alloc2_t &x) {
    ptss_int_t sum = 0;
    unsigned int r = 1;
    for (int j = 0; j < NPH; j++) {
        sum += (x[j].alloc-1)*(r);
        r *= M;
    }
    return sum;
}

/* Inverse of invhashalloc_2 */
alloc2_t invhashalloc_2(ptss_int_t id, const vector<int> &bench) {
    ptss_int_t dec = id;
    ptss_int_t r = M;
    int j;
    /* Create the point */
    alloc2_t vi2;
    vi2.clear();
    for (j = 0; j < NPH; j++) {
        phase_t phinfo;
        phinfo.alloc = dec%r + 1;
        phinfo.bench_id = bench[j];
        vi2.push_back(phinfo);
        dec = dec/r;
    }
    return vi2;
}

ostream& operator<<(ostream& os, const RulesEdgeType &rl) {
    switch (rl) {
        case RULE1 : {os << "RULE1"; break;}
        case RULE2 : {os << "RULE2"; break;}
        case RULE3 : {os << "RULE3"; break;}
        case RULE4 : {os << "RULE4"; break;}
        default    : {os << "Invalid"; exit(127);}
    }
    return os;
}

/* Property Writer for Graphviz Vizualization */
 class SRTDSEGraphVertexPW {
     public:
       SRTDSEGraphVertexPW(SRTDSEGraph& g) : myGraph(g) {}
       template <class VertexOrEdge>
       void operator()(ostream& os, const VertexOrEdge& v) const {
            int i;
            os << "[label=\"{";
            for(i = 0; i < myGraph[v].size()-1; i++)
                os << myGraph[v][i].alloc << ",";
            os << myGraph[v][i].alloc << "}\"]";
       }
    private:
        SRTDSEGraph& myGraph;
 };

class SRTDSEGraphEdgePW {
    public:
       SRTDSEGraphEdgePW(SRTDSEGraph& g) : myGraph(g) {}
       template <class VertexOrEdge>
       void operator()(ostream& os, const VertexOrEdge& e) const {
            // string color;
            // switch(myGraph[e]) {
            //     case RULE1 : color = "red"; break;
            //     case RULE2 : color = "blue"; break;
            //     case RULE3 : color = "purple"; break;
            //     case RULE4 : color = "green"; break;
            // }
           os << "[label=\""<<myGraph[e]<<"\"]"<<endl;
       }
    private:
        SRTDSEGraph& myGraph;
};

double compute_risk(const alloc2_t &x, int D) {
    double mu = 0.0, var = 0.0;
    for(auto jt = x.begin();
        jt != x.end();
        jt++) {
                mu  += mua[jt->alloc-1];
                var += (stda[jt->alloc-1])*(stda[jt->alloc-1]);
    }
    /* Create a normal distribution and evaluate the risk */
    double sd = sqrt(var);
    normal dist(mu,sd);
    
    double risk = 1 - cdf(dist,D);
    return risk;
}

double compute_execution_time2(const alloc2_t &x) {
    double sum = 0.0;
    for(auto jt = x.begin(); jt != x.end(); jt++) {
        // cout << "jt->alloc : "<<jt->alloc<<endl;
        sum += mua[jt->alloc-1];
    }
    return sum;
}


void ptss_DSE_srt::construct_alloc2() {
    double risk, pkp;
    ptss_int_t r = M, j = 0, dec, cnt;
    ptss_int_t TOT = (ptss_int_t)pow(M,NPH);
#ifdef TIMING1
    __eco_roi_start_timer();
#endif
    /* Temporary "Global Variables "*/
    double deadline = this->deadline;
    double dmr      = this->dmr;

    /* Add dummy values for the benchmarks */
    if (this->bench.size() < NPH) {
        int i = NPH;
        while (i-- > 0)
            this->bench.push_back(0);
    }

    vector<SRTDSEGraphVertex> vdMap(TOT);
    for (cnt = 0; cnt < TOT; cnt++) {
        alloc2_t vi2 = invhashalloc_2(cnt,this->bench);
        vdMap[cnt]   = add_vertex(vi2,this->gSS);
        // vdMap[cnt]            = vd1;
    }
    
    /* Make edge connections */
    for (cnt = 0; cnt < TOT; cnt++) {
#ifdef CONNECT_RULE1
        for  (int idx = 0; idx < NPH; idx++) {
            alloc2_t rule1 = this->gSS[vdMap[cnt]];
            if (rule1[idx].alloc > 1) {
                rule1[idx].alloc--;
                ptss_int_t cnt2 = hashalloc_2(rule1);
                add_edge(vdMap[cnt],vdMap[cnt2],RULE1,this->gSS);
            }
        }
#endif
#ifdef CONNECT_RULE2
        for  (int idx = 0; idx < NPH; idx++) {
            alloc2_t rule2 = this->gSS[vdMap[cnt]];
            if (rule2[idx].alloc < M-1) {
                rule2[idx].alloc++;
                ptss_int_t cnt2 = hashalloc_2(rule2);
                add_edge(vdMap[cnt],vdMap[cnt2],RULE2,this->gSS);
            }
        }
#endif
#ifdef CONNECT_RULE3
        for (int i = 0; i < NPH; i++) {
            for (int j = 0; j < NPH; j++) {
                if (i != j) {
                    alloc2_t rule3 = this->gSS[vdMap[cnt]];

                    if (rule3[i].alloc > 1 && rule3[j].alloc < M-1) {
                        rule3[i].alloc--;
                        rule3[j].alloc++;

                        if (compute_execution_time2(rule3) > compute_execution_time2(this->gSS[vdMap[cnt]])) {
                            ptss_int_t cnt2 = hashalloc_2(rule3);
                            add_edge(vdMap[cnt],vdMap[cnt2],RULE3,this->gSS);
                        }
                    }
                }
            }
        }
#endif
    }
    
#ifdef TIMING1
    __eco_roi_stop_timer();
#endif
    /* Dump the graph */
    cout << "\n-- graphviz output START --" << endl;
    ofstream dotf;
    dotf.open("demo.dot",ios::out | ios::trunc);
    write_graphviz(dotf,this->gSS,SRTDSEGraphVertexPW(this->gSS),SRTDSEGraphEdgePW(this->gSS));
    dotf.close();
    cout << "\n-- graphviz output END --" << endl;

    /* Compute Degree properties */

 }

 ptss_DSE_srt::ptss_DSE_srt() {
    this->dmr = 0.9;
    this->deadline = 3000;
    this->construct_alloc2();
 }

 ptss_DSE_srt::ptss_DSE_srt(double dmr) {
    this->dmr = dmr;
    this->deadline = 3000;
    this->construct_alloc2();
 }

  ptss_DSE_srt::ptss_DSE_srt(double deadline, double dmr) {
    this->dmr = dmr;
    this->deadline = deadline;
    this->construct_alloc2();
 }


