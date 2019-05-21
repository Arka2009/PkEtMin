#include <iostream>
#include <vector>
#include <stack>
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
#include <boost/graph/strong_components.hpp>
#include <algorithm>
#include "ptss_dse.hpp"
#include "ptss_config.hpp"
extern "C" {
#include <ecotools/roi_hooks.h>
#include <sys/time.h>
}

// Profiled Data
// static double mua  [] = {269.568837,	141.985532, 99.508883, 79.006262, 66.454092, 58.649068, 52.804705, 48.960141, 45.73406, 43.455542, 41.386889, 40.046129, 38.695398, 37.750183, 36.840822, 36.41925, 35.715777, 35.350802, 34.767695, 34.493228, 34.033093, 33.944396, 33.676262, 33.871417, 33.654013, 33.737325, 33.825452, 34.021759, 34.223725, 34.479067, 34.694037, 35.157394}; 
// static double stda [] = {0.676473946, 0.335842225, 0.208458629, 0.319806191, 0.313754681, 0.340024999, 0.340440891, 0.441634464, 0.424892928, 0.45393722, 0.470011702, 0.479161768, 0.496199557, 0.493329504, 0.524060111, 0.523248507, 0.564837145, 0.590188106, 0.584813646, 0.602546264, 0.597015913, 0.600979201, 0.59552246, 0.595876665, 0.587283577, 0.645952011, 0.666749578, 0.732030737, 0.836999403, 0.942281274, 1.046115194, 1.103043517}; 
using namespace std;
using namespace boost;

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
        case RULE1   : {os << "RULE1"; break;}
        case RULE2   : {os << "RULE2"; break;}
        case RULE3   : {os << "RULE3"; break;}
        case RULE4   : {os << "RULE4"; break;}
        case INVALID : {os << "INVALID"; break;}
        default    : {os << "Invalid"; exit(127);}
    }
    return os;
}

#ifdef USEBGL
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
            if (myGraph[e] == INVALID)
                os << "[color=red]";
            else
                os << "[label=\""<<myGraph[e]<<"\"]"<<endl;
       }
    private:
        SRTDSEGraph& myGraph;
};
#endif

double compute_risk(const alloc2_t &x, int D) {
    double mu = compute_execution_time(x);
    double sd = compute_execution_time_std(x);
    /* Create a normal distribution and evaluate the risk */
    boost::math::normal dist(mu,sd);
    // cout <<"mu = " <<mu<<endl;
    // cout<<"sd = "<<sd<<endl;
    // cout<<"D = "<<D<<endl;
    // cout<<"dist value = "<<boost::math::cdf(dist,D)<<endl;
    double risk = 1 - boost::math::cdf(dist,D);
    return risk;
}

/* 
 * Compute the standard deviation in execution time 
 * Right now this is set to 0.01 of the mean, but
 * could be changed later, when profiling data becomes
 * available.
 */
double compute_execution_time_std(const alloc2_t &x) {
    double stdev = 0.0;
    double var   = 0.0;
    double et    = 0.0;
    for(auto jt = x.begin(); jt != x.end(); jt++) {
        et = estimate_exec_time(jt->alloc,jt->bench_id);
        if (et >= 0.0) {
            stdev  = 0.1*et;
            var   += stdev*stdev;
        }
        else {
            return -HUGE_VAL;
        }
        
    }
    return sqrt(var);
}

void ptss_DSE_srt::construct_alloc2() {
    // ptss_int_t TOT = M^NPH;
    // cout << "NPH = " << NPH << "M = " << M <<", TOT = " << TOT << endl;
    double risk, pkp;
    bool feasible;
    ptss_int_t r = M, j = 0, dec, cnt;
    #ifdef TIMING1
    __eco_roi_start_timer();
    #endif
    /* Temporary "Global Variables "*/
    double deadline = this->deadline;
    double dmr      = this->dmr;

    /* Add dummy values for the benchmarks (if not already created) */
    if (this->bench.size() < NPH) {
        int i = NPH;
        while (i-- > 0)
            this->bench.push_back(0);
    }

    #ifdef USEBGL
    /* Create the vertex set */
    vector<SRTDSEGraphVertex> vdMap(TOT2);
    #endif

    #ifdef ENABLE_BRUTE_FORCE
    /* Brute Force Search of the optimal point */
    double opt_pkp2  = HUGE_VAL;
    double opt_risk2 = 1;
    for (cnt = 0; cnt < TOT2; cnt++) {
        alloc2_t vi2 = invhashalloc_2(cnt,this->bench);
        // cout << "cntr : "<< cnt << ", Creating : " << vi2 << ", mean : "<<compute_execution_time(vi2) << ", std : " <<compute_execution_time_std(vi2)<<endl; 
        #ifdef USEBGL
        vdMap[cnt]   = add_vertex(vi2,this->gSS);
        validVertices.push_back(vdMap[cnt]);
        #endif
        

        /* Store the optimal point */
        boost::tie(feasible,pkp) = checkFeasibilityAndPKP(vi2,this->deadline,this->dmr);
        if (feasible) {
            if (pkp < opt_pkp2) {
                opt_pkp2        = pkp;
                this->opt_point = vi2;
            }
        }

        // cout << cnt << "|"<< cnt << "|" << vi2 << ",:" << opt_pkp2 << ",:"<<pkp <<",:" << feasible << endl;
        /* Unmark the exploredSet tag */
        #ifndef USEBITSET
        this->exploredSet.push_back(false);
        #else
        this->exploredSet[cnt] = false;
        #endif
    }
    cout << "Brute force search Optimal Point : " <<this->opt_point<<", pkp : "<<compute_pkpower(this->opt_point)<<endl;
    #else
    for (cnt = 0; cnt < TOT2; cnt++) {
        #ifndef USEBITSET
        this->exploredSet.push_back(false);
        #else
        this->exploredSet[cnt] = false;
        #endif
    }
    cout << "Brute force search disabled " <<endl;
    this->opt_point = invhashalloc_2(0,this->bench);
    #endif
    
    #ifdef USEBGL
    /* Make edge connections */
    for (cnt = 0; cnt < TOT2; cnt++) {
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
            if (rule2[idx].alloc < M) {
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


                        if (compute_execution_time(rule3) > compute_execution_time(this->gSS[vdMap[cnt]])) {
                            /* Add more edges */
                            ptss_int_t cnt2 = hashalloc_2(rule3);
                            add_edge(vdMap[cnt],vdMap[cnt2],RULE3,this->gSS);

                            ptss_int_t delta = std::min((this->gSS[vdMap[cnt]][i].alloc-1),(M-this->gSS[vdMap[cnt]][i].alloc));
                            for (int k = 1; k < delta; k++) {
                                alloc2_t tmp = rule3;
                                tmp[i].alloc--;
                                tmp[j].alloc++;
                                ptss_int_t cnt3 = hashalloc_2(tmp);
                                add_edge(vdMap[cnt],vdMap[cnt3],RULE3,this->gSS);
                            }
                        }
                    }
                }
            }
        }
        #endif
        #ifdef CONNECT_RULE4
        for (int i = 0; i < NPH; i++) {
            for (int j = 0; j < NPH; j++) {
                if (i != j) {
                    alloc2_t rule4 = this->gSS[vdMap[cnt]];
                    if (rule4[i].alloc > 1 && rule4[j].alloc < M-1) {
                        rule4[i].alloc--;
                        rule4[j].alloc++;
                    
                        if (compute_pkpower(rule4) > compute_pkpower(this->gSS[vdMap[cnt]])) {
                            /* Add more edges */
                            ptss_int_t cnt2 = hashalloc_2(rule4);
                            add_edge(vdMap[cnt],vdMap[cnt2],RULE4,this->gSS);

                            ptss_int_t delta = std::min((this->gSS[vdMap[cnt]][i].alloc-1),(M-this->gSS[vdMap[cnt]][i].alloc));
                            for (int k = 1; k < delta; k++) {
                                alloc2_t tmp2 = rule4;
                                tmp2[i].alloc--;
                                tmp2[j].alloc++;
                                ptss_int_t cnt6 = hashalloc_2(tmp2);
                                add_edge(vdMap[cnt],vdMap[cnt6],RULE4,this->gSS);
                            }
                        }
                    }
                }
            }
        }
        #endif
    }
    #endif

    #ifdef TIMING1
    __eco_roi_stop_timer();
    #endif
    #ifdef USEBGL
    this->graphBuilt = true;
    #endif
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

  ptss_DSE_srt::ptss_DSE_srt(double dmr,double deadline) {
    this->dmr = dmr;
    this->deadline = deadline;
    this->construct_alloc2();
 }



/* Read the workload from a File, compute and write the results to another CSV file */
ptss_DSE_srt::ptss_DSE_srt(rapidcsv::Document &inDoc,unsigned int idx, bool &done) {
    // vector<int> benchid;
    double deadline;

    /* Continuously Read the CSV Files */
    int i = 1, j = 0;
    try {
        std::vector<double> wkld = inDoc.GetRow<double>(idx);
        this->deadline = wkld[0];
        this->dmr      = wkld[1];
        for (i = 2; i < wkld.size(); i++) {
            this->bench.push_back((int)wkld[i]);
        }
        // struct timespec t1, t2;
        // clock_gettime(CLOCK_REALTIME,&t1);
    } catch(out_of_range) {
        done = true;
        std::cout << "Reached End of CSV File (numPhases = "<< NPH <<")" << std::endl;
        return;
    }
    this->construct_alloc2();
}


#ifdef USEBGL
ptss_int_t ptss_DSE_srt::validOutDegrees(const SRTDSEGraphVertex& v) {
    ptss_int_t k = 0;
    SRTDSEGraphOutEdgeIterator eiS, eiE;
    for (boost::tie(eiS, eiE) = out_edges(v,this->gSS); eiS != eiE;eiS++) {
        if (this->gSS[*eiS] != INVALID);
            k++;
    }
    return k;
}
#endif

template <class Q>
void clearQueue(Q & q) {
    q = Q();
}

/* Check if an allocation is valid, else return HUGE_VAL */
std::pair<bool,double> checkFeasibilityAndPKP(const alloc2_t& v, double deadline, double dmr) {
    if (compute_risk(v,deadline) > dmr) {
        return make_pair(false,compute_pkpower(v));
    } 
    return make_pair(true,compute_pkpower(v));
}

#ifdef USEBGL
/* Randomly walk the rule graph */
void ptss_DSE_srt::random_walk() {
    ptss_int_t TOT = num_vertices(this->gSS);
    
    SRTDSEGraphVertexIterator vertItBegin, verItEnd;
    tie(vertItBegin,verItEnd) = vertices(this->gSS);
    SRTDSEGraphOutEdgeIterator eiS, eiE;
    SRTDSEGraphInEdgeIterator eiIS, eiIE;

    /* Initialize with a randomly selected vertex */
    std::advance<SRTDSEGraphVertexIterator>(vertItBegin,2);
    SRTDSEGraphVertex rwalker     = *vertItBegin;
    SRTDSEGraphVertex prevrwalker = rwalker;
    priority_queue<pair<ptss_int_t,SRTDSEGraphVertex>> nbdQ;
    ptss_int_t deg;

    int j = 0;
    bool newset = false;
    while (j++ < 100000) {
        prevrwalker = rwalker;
        if (this->dmr <= compute_risk(this->gSS[rwalker],this->deadline)) {
            cout << "At("<<j<<") vertex (Infeasible)"<<this->gSS[rwalker] << endl;
            /* Invalidate all the incoming edges */
            for (boost::tie(eiIS, eiIE) = in_edges(rwalker,this->gSS); eiIS != eiIE;eiIS++) {
                this->gSS[*eiIS] = INVALID;
            }
            newset = false;
            for (boost::tie(eiS, eiE) = out_edges(rwalker,this->gSS); eiS != eiE;eiS++) {
                /* Invalidate Rule1 and Rule3 edges */      
                if(this->gSS[*eiS] == RULE1 || this->gSS[*eiS] == RULE3) {
                    this->gSS[*eiS] = INVALID;
                }
                /* Move along edge2/edge4  (which one) */
                else if ((this->gSS[*eiS] == RULE2 || this->gSS[*eiS] == RULE4)/* &&  !newset*/){
                    SRTDSEGraphVertex tmp = target(*eiS,this->gSS);
                    // newset = true;
                    nbdQ.push(pair<ptss_int_t,SRTDSEGraphVertex>(this->validOutDegrees(tmp),tmp));
                }
            }
        } else {
            cout << "At("<<j<<") vertex (Feasible)"<<this->gSS[rwalker] << endl;
            newset = false;
            for (boost::tie(eiS, eiE) = out_edges(rwalker,this->gSS); eiS != eiE;eiS++) {
                /* Invalidate Rule2 and Rule4 edges */
                if(this->gSS[*eiS] == RULE2 || this->gSS[*eiS] == RULE4) {
                    this->gSS[*eiS] = INVALID;
                }
                /* Move along edge1/edge3 (which one) */
                else if ((this->gSS[*eiS] == RULE1 || this->gSS[*eiS] == RULE3)/* &&  !newset*/){
                    SRTDSEGraphVertex tmp = target(*eiS,this->gSS);
                    // newset = true;
                    nbdQ.push(pair<ptss_int_t,SRTDSEGraphVertex>(this->validOutDegrees(tmp),tmp));
                }
            }
        }
        if (!nbdQ.empty()) {
            boost::tie(deg,rwalker) = nbdQ.top();
            nbdQ.pop();
            cout << "At("<<j<<") vertex (NBD)"<<this->gSS[rwalker] << " deg = "<<deg<< endl;
        } else {
            /* Ramove this vertex */
            cout << "Empty NBD" <<endl;
            this->validVertices.remove(rwalker);

            /* Random restart */
            ptss_int_t N = this->validVertices.size();
            if (N > 0) {
                auto it = this->validVertices.begin();
                std::advance(it,rand()%N);
                rwalker = *it;
                cout <<"At("<<j<<") vertex (Restart)"<<this->gSS[rwalker] << endl;
            }
        }
        // clearQueue(nbdQ);
        bool diffObj = eval_manhtn(this->gSS[prevrwalker],this->gSS[rwalker]);
        if (diffObj) {
            cout << "At("<<j<<") vertex (Break)"<<this->gSS[rwalker] << endl;
            break;
        }
    }
    this->rwalkerFinal = this->gSS[rwalker];
}


/* Randomly walk the rule graph (no restarts) */
void ptss_DSE_srt::random_walk2() {
    ptss_int_t TOT = num_vertices(this->gSS);
    
    SRTDSEGraphVertexIterator vertItBegin, verItEnd;
    tie(vertItBegin,verItEnd) = vertices(this->gSS);
    SRTDSEGraphOutEdgeIterator eiS, eiE;
    SRTDSEGraphInEdgeIterator eiIS, eiIE;

    /* Initialize with a randomly selected vertex */
    std::advance<SRTDSEGraphVertexIterator>(vertItBegin,2);
    SRTDSEGraphVertex rwalker     = *vertItBegin;
    SRTDSEGraphVertex prevrwalker = rwalker;
    priority_queue<pair<ptss_int_t,SRTDSEGraphVertex>> nbdQ;
    ptss_int_t deg;

    int j = 0;
    while (j++ < 100000) {
        prevrwalker = rwalker;
        if (this->dmr <= compute_risk(this->gSS[rwalker],this->deadline)) {
            cout << "At("<<j<<") vertex (Infeasible)"<<this->gSS[rwalker] << endl;
            /* Invalidate all the incoming edges */
            for (boost::tie(eiIS, eiIE) = in_edges(rwalker,this->gSS); eiIS != eiIE;eiIS++) {
                this->gSS[*eiIS] = INVALID;
            }
            for (boost::tie(eiS, eiE) = out_edges(rwalker,this->gSS); eiS != eiE;eiS++) {
                /* Invalidate Rule1 and Rule3 edges */      
                if(this->gSS[*eiS] == RULE1 || this->gSS[*eiS] == RULE3) {
                    this->gSS[*eiS] = INVALID;
                }
                /* Move along edge2/edge4  (which one) */
                else if ((this->gSS[*eiS] == RULE2 || this->gSS[*eiS] == RULE4)/* &&  !newset*/){
                    SRTDSEGraphVertex tmp = target(*eiS,this->gSS);
                    nbdQ.push(pair<ptss_int_t,SRTDSEGraphVertex>(this->validOutDegrees(tmp),tmp));
                }
            }
        } else {
            cout << "At("<<j<<") vertex (Feasible)"<<this->gSS[rwalker] << endl;
            for (boost::tie(eiS, eiE) = out_edges(rwalker,this->gSS); eiS != eiE;eiS++) {
                /* Invalidate Rule2 and Rule4 edges */
                if(this->gSS[*eiS] == RULE2 || this->gSS[*eiS] == RULE4) {
                    this->gSS[*eiS] = INVALID;
                }
                /* Move along edge1/edge3 (which one) */
                else if ((this->gSS[*eiS] == RULE1 || this->gSS[*eiS] == RULE3)/* &&  !newset*/){
                    SRTDSEGraphVertex tmp = target(*eiS,this->gSS);
                    nbdQ.push(pair<ptss_int_t,SRTDSEGraphVertex>(this->validOutDegrees(tmp),tmp));
                }
            }
        }
        if (!nbdQ.empty()) {
            boost::tie(deg,rwalker) = nbdQ.top();
            nbdQ.pop();
            // cout << "At("<<j<<") vertex (NBD)"<<this->gSS[rwalker] << " deg = "<<deg<< endl;
        }
        bool diffObj = eval_manhtn(this->gSS[prevrwalker],this->gSS[rwalker]);
        if (diffObj) {
            cout << "At("<<j<<") vertex (Break)"<<this->gSS[rwalker] << endl;
            break;
        }
    }
    this->rwalkerFinal = this->gSS[rwalker];
}

ptss_int_t ptss_DSE_srt::backtrackAll(const SRTDSEGraphVertex &start_node) {
    double pkp, feasible;
    ptss_int_t cntr;
    SRTDSEGraphOutEdgeIterator eiS, eiE;
    SRTDSEGraphInEdgeIterator eiIS, eiIE;
    SRTDSEGraphVertex node2     = start_node;
    double currentPkp           = HUGE_VAL;       /* Optimal Pkp found so far */
    stack<SRTDSEGraphVertex> stck;
    stck.push(node2);


    while (!stck.empty()) {
        SRTDSEGraphVertex node = stck.top();
        stck.pop();
        cntr++;
        boost::tie(feasible,pkp) = checkFeasibilityAndPKP(this->gSS[node],this->deadline,this->dmr);
        // cout << cntr << ",:" << this->gSS[node] << ",:" << currentPkp << ",:"<<pkp <<",:" << feasible << endl;

    
        /* Edge Invalidation */
        if (!feasible) {
            /* Incoming Edges */
            for (boost::tie(eiIS, eiIE) = in_edges(node,this->gSS); eiIS != eiIE;eiIS++) {
                this->gSS[*eiIS] = INVALID;
            }
            /* Rule 1 and rule 3 edges */
            for (boost::tie(eiS, eiE) = out_edges(node,this->gSS); eiS != eiE;eiS++) {
                /* Invalidate Rule1 and Rule3 edges */      
                if(this->gSS[*eiS] == RULE1 || this->gSS[*eiS] == RULE3) {
                    this->gSS[*eiS] = INVALID;
                }
            }
        } else {
            /* Incoming Edges */
            for (boost::tie(eiIS, eiIE) = in_edges(node,this->gSS); eiIS != eiIE;eiIS++) {
                this->gSS[*eiIS] = INVALID;
            }
            /* Rule 2 and rule 4 edges */
            for (boost::tie(eiS, eiE) = out_edges(node,this->gSS); eiS != eiE;eiS++) {
                /* Invalidate Rule1 and Rule3 edges */      
                if(this->gSS[*eiS] == RULE2 || this->gSS[*eiS] == RULE4) {
                    this->gSS[*eiS] = INVALID;
                }
            }
        }

        if (feasible) {
            if (pkp <= currentPkp) {
                this->bktrack_point = this->gSS[node];
                currentPkp = pkp;
            }
            /* Additional Edge Eliminations */
            for (boost::tie(eiS, eiE) = out_edges(node,this->gSS); eiS != eiE;eiS++) {
                if(this->gSS[*eiS]!=INVALID) {
                    SRTDSEGraphVertex tmp = target(*eiS,this->gSS);
                    if (compute_pkpower(this->gSS[tmp]) > currentPkp)
                        this->gSS[*eiS] = INVALID;
                }
            }
        }
    
        for (boost::tie(eiS, eiE) = out_edges(node,this->gSS); eiS != eiE;eiS++) {
            if(this->gSS[*eiS]!=INVALID) {
                SRTDSEGraphVertex tmp = target(*eiS,this->gSS);
                // this->backtrackAll(tmp, currentPkp,cntr);
                stck.push(tmp);
            }
        }
    }
    return cntr;
}

void ptss_DSE_srt::random_walk3(ptss_int_t start_offset) {
    /* Initialize with a randomly selected vertex */
    SRTDSEGraphVertexIterator vertItBegin, verItEnd;
    tie(vertItBegin,verItEnd) = vertices(this->gSS);
    std::advance<SRTDSEGraphVertexIterator>(vertItBegin,start_offset);
    SRTDSEGraphVertex rwalker     = *vertItBegin;
    double currentPkp             = HUGE_VAL;
    ptss_int_t cntr               = 0;
    stack<SRTDSEGraphVertex> stck;
    ptss_int_t numIter            = this->backtrackAll(rwalker);
    cout << "\nBacktrack Optimal Point="<<this->bktrack_point<<",pkp="<<compute_pkpower(this->bktrack_point) << ",numIter=" << numIter << endl;
}

void ptss_DSE_srt::printRuleGraph() {
    /* Dump the graph */
    if (TOT2 <= 200) {
        cout << "\n-- graphviz output START --" << endl;
        ofstream dotf;
        dotf.open("demo.dot",ios::out | ios::trunc);
        write_graphviz(dotf,this->gSS,SRTDSEGraphVertexPW(this->gSS),SRTDSEGraphEdgePW(this->gSS));
        dotf.close();
        cout << "\n-- graphviz output END --" << endl;
    }
    cout << "\nBrute Force Optimal Point="<<this->opt_point<<",pkp="<<compute_pkpower(this->opt_point) << endl;
    // cout << "\nRWalker Optimal Point="<<this->rwalkerFinal<<",pkp="<<compute_pkpower(this->rwalkerFinal) << endl;
}
#endif

std::set<ptss_int_t> union_sets(std::set<ptss_int_t> A, std::set<ptss_int_t> B) {
    A.insert(B.cbegin(), B.cend());
    return A;
}

void does_allocation_exists(ptss_int_t allc, int rule, const vector<int> bench) {
    if(allc >= TOT2) {
        cout << "Rule "<<rule<<" : Allocation does not exists : " << allc << "," << invhashalloc_2(allc,bench) << endl;
        exit(127);
    }
}

set<ptss_int_t> ptss_DSE_srt::applyRule1(const ptss_int_t pt1) {
    set<ptss_int_t> children;
    for (int idx = 0; idx < NPH; idx++) {
        alloc2_t rule1 = invhashalloc_2(pt1,this->bench);
        if (rule1[idx].alloc > 1) {
                rule1[idx].alloc--;
                ptss_int_t cnt2 = hashalloc_2(rule1);

                does_allocation_exists(cnt2, 1, this->bench);
                children.insert(cnt2);
        }
    }
    return children;
}

set<ptss_int_t> ptss_DSE_srt::applyRule2(const ptss_int_t pt1) {
    set<ptss_int_t> children;
    for (int idx = 0; idx < NPH; idx++) {
        alloc2_t rule1 = invhashalloc_2(pt1,this->bench);
        if (rule1[idx].alloc < M-1) { /* TODO : M-1 is a probable bug fix this */
                rule1[idx].alloc++;
                ptss_int_t cnt2 = hashalloc_2(rule1);

                // cout <<"pt1 : "<<pt1<<",allc : "<<invhashalloc_2(pt1,this->bench)<<endl;
                does_allocation_exists(cnt2, 2, this->bench);
                children.insert(cnt2);
        }
    }
    return children;
}

set<ptss_int_t> ptss_DSE_srt::applyRule3(const ptss_int_t pt1) {
    set<ptss_int_t> children;
    alloc2_t orig = invhashalloc_2(pt1,this->bench);
    for (int i = 0; i < NPH; i++) {
        for (int j = 0; j < NPH; j++) {
            if (i != j) {
                alloc2_t rule3 = orig;
                if (rule3[i].alloc > 1 && rule3[j].alloc < M-1) {
                    rule3[i].alloc--;
                    rule3[j].alloc++;
                    if (compute_execution_time(rule3) > compute_execution_time(orig)) {
                        /* Add more edges */
                        ptss_int_t cnt2 = hashalloc_2(rule3);
                        children.insert(cnt2);
                        does_allocation_exists(cnt2, 31, this->bench);
                        ptss_int_t delta = std::min((orig[i].alloc-1),(M-orig[i].alloc));
                        for (int k = 1; k < delta; k++) {
                            alloc2_t tmp = rule3;
                            tmp[i].alloc--;
                            tmp[j].alloc++;
                            ptss_int_t cnt3 = hashalloc_2(tmp);
                            does_allocation_exists(cnt3, 32, this->bench);
                            children.insert(cnt3);
                        }
                    }
                }
            }
        }
    }
    return children;
}

set<ptss_int_t> ptss_DSE_srt::applyRule4(const ptss_int_t pt1) {
    set<ptss_int_t> children;
    alloc2_t orig = invhashalloc_2(pt1,this->bench);
    for (int i = 0; i < NPH; i++) {
        for (int j = 0; j < NPH; j++) {
            if (i != j) {
                alloc2_t rule4 = orig;
                if (rule4[i].alloc > 1 && rule4[j].alloc < M-1) {
                    rule4[i].alloc--;
                    rule4[j].alloc++;
                
                    if (compute_pkpower(rule4) > compute_pkpower(orig)) {
                        /* Add more edges */
                        ptss_int_t cnt2 = hashalloc_2(rule4);
                        children.insert(cnt2);
                        does_allocation_exists(cnt2, 41, this->bench);
                        ptss_int_t delta = std::min((orig[i].alloc-1),(M-orig[i].alloc));
                        for (int k = 1; k < delta; k++) {
                            alloc2_t tmp2 = rule4;
                            tmp2[i].alloc--;
                            tmp2[j].alloc++;
                            ptss_int_t cnt6 = hashalloc_2(tmp2);
                            does_allocation_exists(cnt6, 42, this->bench);
                            children.insert(cnt6);
                        }
                    }
                }
            }
        }
    }
    return children;
}

ptss_int_t ptss_DSE_srt::backtrackAll(ptss_int_t start_node) {
    double pkp, feasible;
    double currentPkp           = HUGE_VAL;       /* Optimal Pkp found so far */
    ptss_int_t cntr = 0;
    stack<ptss_int_t> stck;
    stck.push(start_node);
    /* Temporary "Global Variables "*/
    const double deadline = this->deadline;
    const double dmr      = this->dmr;
    const vector<int> bench = this->bench;

    while (!stck.empty()) {
        ptss_int_t node = stck.top();
        alloc2_t node2  = invhashalloc_2(node,bench);
        stck.pop();
        /* Skip this node if it is already selected */
        if (this->exploredSet[node]) {
            continue;
        }
        cntr++;
        boost::tie(feasible,pkp) = checkFeasibilityAndPKP(node2,deadline,dmr);
        cout << cntr << "|"<< node << "|" << node2 << ",:" << currentPkp << ",:"<<pkp <<",risk:" << compute_risk(node2,deadline) <<",dmr:" << dmr<<",feasible:" << feasible << endl;
        cout << cntr << "|"<< node << "|" << ",dmr:" << this->dmr << endl; 

        /* Node Expansions */
        set<ptss_int_t> children1  = this->applyRule1(node);
        set<ptss_int_t> children2  = this->applyRule2(node);
        set<ptss_int_t> children3  = this->applyRule3(node);
        set<ptss_int_t> children4  = this->applyRule4(node);

        /* Invalidations */
        this->exploredSet[node] = true;
        if (!feasible) {
            for (auto it = children1.begin(); it != children1.end(); it++)
                this->exploredSet[*it] = true;
            for (auto it = children3.begin(); it != children3.end(); it++)
                this->exploredSet[*it] = true;
        } else {
            for (auto it = children2.begin(); it != children2.end(); it++)
                this->exploredSet[*it] = true;
            for (auto it = children4.begin(); it != children4.end(); it++)
                this->exploredSet[*it] = true;
        }

        /* Union all children */
        set<ptss_int_t> allchildren = union_sets(children1,children2);
        allchildren = union_sets(allchildren,children3);
        allchildren = union_sets(allchildren,children4);
        #ifdef DEBUGCHILDREN
        int k = 0;
        for (auto it = allchildren.begin(); it != allchildren.end(); it++) {
            cout << "children - "<< k++ << *it << ", alloc : " <<invhashalloc_2(*it,this->bench)<<endl;
        }
        #endif

        if (feasible) {
            if (pkp < currentPkp) {
                this->bktrack_point = node2;
                currentPkp = pkp;
            }
            /* Additional Edge Eliminations */
            for (auto it = allchildren.begin(); it != allchildren.end(); it++) {
                alloc2_t tmpNode = invhashalloc_2(*it,bench);
                if (compute_pkpower(tmpNode) > currentPkp) {
                    this->exploredSet[*it] = true;
                }
            }
        }

        /* push all unexplored children into the stack */
        for (auto it = allchildren.begin(); it != allchildren.end(); it++) {
            if (!this->exploredSet[*it])
               stck.push(*it); 
        }
    }
    return cntr;

}

void ptss_DSE_srt::random_walk4(ptss_int_t start_offset) {
    ptss_int_t numIter  = this->backtrackAll(start_offset);
    cout << "\nBacktrack Optimal Point="<<this->bktrack_point<<",pkp="<<compute_pkpower(this->bktrack_point) << ",numIter=" << numIter << endl;
}