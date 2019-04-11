/**
 * Implementation of Rules to eliminate
 * Design points. Check Risk-Base Allocation.pptx
 * for more details.
 */
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

/* HI - Test Editing */
using namespace std;
using boost::math::normal;


// a == b ?
bool is_eq(const alloc_t &a, const alloc_t &b) {
    return (a == b);
}

// a > b ?
bool is_gt(const alloc_t &a, const alloc_t &b) {
    bool ret = !(is_eq(a,b));
    int idx = 0;
    for(auto it = a.begin();
        it != a.end();
        it++) {
            if (*it >= b[idx++])
                ret = ret && true;
            else
                ret = ret && false;
    }
    return ret;
}

// a < b ?
bool is_lt(const alloc_t &a, const alloc_t &b) {
    bool ret = !(is_eq(a,b));
    int idx = 0;
    for(auto it = a.begin();
        it != a.end();
        it++) {
            if (*it <= b[idx++])
                ret = ret && true;
            else
                ret = ret && false;
    }
    return ret;
}

// a and b are not comparable
bool is_incomparable(const alloc_t &a, const alloc_t &b) {
    bool ret = (!is_gt(a,b)) && \
               (!is_lt(a,b)) && \
               (!is_eq(a,b));
    return ret;
}

// lexicographic comparison of two allocations ( a < b ? )
bool lex_comp(const alloc_t &a, const alloc_t &b) {
    auto it2 = b.begin();
    for (auto it = a.begin(); it != a.end(); it++, it2++) {
        if (*it2 < *it)
            return false;
    }
    return true;
}


// just for the sample -- print the data sets
std::ostream& operator<<(std::ostream& os, const alloc_t& vi) {
  os << "(|";
  std::copy(vi.begin(), vi.end(), std::ostream_iterator<int>(os, "|"));
  os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const all_alloc_t& vvi) { 
  os << "{\n";
  for(auto it = vvi.begin();
      it != vvi.end();
      it++) {
      os << "  " << *it << "\n";
  }
  os << "}";
  return os;
}

// Profiled Data
static double mua  [] = {269.568837,	141.985532, 99.508883, 79.006262, 66.454092, 58.649068, 52.804705, 48.960141, 45.73406, 43.455542, 41.386889, 40.046129, 38.695398, 37.750183, 36.840822, 36.41925, 35.715777, 35.350802, 34.767695, 34.493228, 34.033093, 33.944396, 33.676262, 33.871417, 33.654013, 33.737325, 33.825452, 34.021759, 34.223725, 34.479067, 34.694037, 35.157394}; 
static double stda [] = {0.676473946, 0.335842225, 0.208458629, 0.319806191, 0.313754681, 0.340024999, 0.340440891, 0.441634464, 0.424892928, 0.45393722, 0.470011702, 0.479161768, 0.496199557, 0.493329504, 0.524060111, 0.523248507, 0.564837145, 0.590188106, 0.584813646, 0.602546264, 0.597015913, 0.600979201, 0.59552246, 0.595876665, 0.587283577, 0.645952011, 0.666749578, 0.732030737, 0.836999403, 0.942281274, 1.046115194, 1.103043517}; 

double compute_risk(const alloc_t &x) {
    double mu = 0.0, var = 0.0;
    for(auto jt = x.begin();
        jt != x.end();
        jt++) {
                mu  += mua[*jt-1];
                var += (stda[*jt-1])*(stda[*jt-1]);
    }
    /* Create a normal distribution and evaluate the risk */
    double sd = sqrt(var);
    normal dist(mu,sd);
    
    double risk = 1 - cdf(dist,D);
    return risk;
}

double compute_execution_time(const alloc_t &x) {
    double sum = 0.0;
    for(auto jt = x.begin(); jt != x.end(); jt++) {
        sum += mua[*jt-1];
    }
    return sum;
}

double compute_estimated_util(const alloc_t &x) {
    double util = 0;
    for(auto jt = x.begin();
        jt != x.end();
        jt++) {
            util += (*jt) * mua[*jt-1];
    }
    return util;
}

/* recursively construct and add allocations */
void construct_alloc(all_alloc_t &vvi, int ph) {
    alloc_t vi2;
    static long unsigned int cnt = 0;
    long unsigned int r   = M;
    long unsigned int dec;
    long unsigned int j   = 0;

    if (ph == NPH) {
        //cout << "ph : " << ph << ", invoc : " << cnt << endl;
        dec = cnt++;

        /* Resolve cnt into M-radix number */
        for (j = 0; j < NPH; j++) {
            //cout << dec%r + 1 << ",";
            vi2.push_back(dec%r + 1);
            dec = dec/r;
        }
        vvi.insert(vi2);
        //cout << "}\n";
        return;
    }
    for (int m = 1; m <= M; m++) {
        construct_alloc(vvi,ph+1);
    }
}


void ptss_DSE::display() {
    cout << this->search_space << "\n";
    cout << "Lower : ";
    cout << this->lower << "\n";
    cout << "Upper : ";
    cout << this->upper << "\n";
}

// Create All points
ptss_DSE::ptss_DSE() {
    construct_alloc(this->search_space,0);
    this->is_initialized = true;
    this->dmr = 1.0;
}

ptss_DSE::ptss_DSE(double dmr) {
    construct_alloc(this->search_space,0);
    this->init_point(dmr);
    this->is_initialized = true;
    this->dmr = dmr;
}


// Evaluate all points
void ptss_DSE::evaluate_all() {
    this->opt_point = this->lower;
    this->opt_util  = 10e9;
    this->opt_risk  = 1;

    // cout << "Search space size " << search_space.size();
    for(auto it = search_space.begin();
        it != search_space.end();
        it++) {
        
        double risk = compute_risk(*it);
        double util = compute_estimated_util(*it);

        if (risk <= this->dmr && util <= this->opt_util) {
            this->opt_point = *it;
            this->opt_util  = util;
            this->opt_risk  = risk;
        }
        // cout << *it << "," << risk << "," << util << "\n";
        usleep(50);
    }
    cout << this->opt_point << "," << this->opt_risk << "," << this->opt_util << "\n";
}

void ptss_DSE::init_point(double dmr) {
    alloc_t tmp;
    double risk;
    int m;
    
    for (m = 1; m <= M; m++) {
        tmp.clear();

        // Create a uniform allocation
        for (int ph = 0; ph < NPH; ph++)
            tmp.push_back(m);
        
        // Compute the risk
        risk = compute_risk(tmp);

        if (risk <= dmr) {
            this->upper = tmp;
            break;
        }
    }

    tmp.clear();
    for (int ph = 0; ph < NPH; ph++)
        tmp.push_back(m-1);    
    this->lower = tmp;

    cout << "lower" << this->lower << ", dmr = " << compute_risk(this->lower) << endl;
    cout << "upper" << this->upper << ", dmr = " << compute_risk(this->upper) << endl;
}



/*
 * The action to generate a
 * set of all (child) points is to 
 * select a src and dst from 
 * startpoint and
 * transfer a single core
 * from one src to destitation.
 */
void epsilon_move2_rule1(all_alloc_t &children,\
                         const alloc_t startpoint,\
                         const all_alloc_t fbidden) {
    /* Expand all the elements of the frontier */
    for (int idx = 0; idx < NPH; idx++) {
        alloc_t tmp2  = startpoint;   
        if (tmp2[idx] > 1) {
            tmp2[idx]--;

            /* Inert only when it doesn't exist in forbidden set */
            // if (fbidden.find(tmp2) == fbidden.end())
                children.insert(tmp2);
        }
    }
}

void epsilon_move2_rule2(all_alloc_t &children,\
                         const alloc_t startpoint,\
                         const all_alloc_t fbidden) {
     /* Expand all the elements of the frontier */
    for (int idx = 0; idx < NPH; idx++) {
        alloc_t tmp2  = startpoint;   
        if (tmp2[idx] < M) {
            tmp2[idx]++;

            /* Inert only when it doesn't exist in forbidden set */
            // if (fbidden.find(tmp2) == fbidden.end())
                children.insert(tmp2);
        }
    }
}

void epsilon_move2_rule3(all_alloc_t &children,\
                         const alloc_t startpoint,\
                         const all_alloc_t fbidden) {
    int i, j;

    for (i = 0; i < NPH; i++) {
        for (j = 0; j < NPH; j++) {
            if (i != j) {
                if (startpoint[i] > 1) {
                    /* Transfer a core from i to j */
                    alloc_t new_point = startpoint;
                    new_point[i]--;
                    new_point[j]++;

                    /* Check Execution Time before insert */
                    if (compute_execution_time(new_point) > compute_execution_time(startpoint)) {
                        // if (fbidden.find(new_point) == fbidden.end())
                            children.insert(new_point);

                        /* Also insert other points derived from it */
                        alloc_t new_point2 = new_point;
                        while (--new_point2[i] > 0) {
                            if (++new_point2[j] <= NPH)
                                // if (fbidden.find(new_point2) == fbidden.end())
                                    children.insert(new_point2);
                        }
                    }
                }
            }
        }
    }
}



void epsilon_move2_test() {
    alloc_t startpoint = {5,5,6,4,4};
    all_alloc_t fbidden;
    all_alloc_t children = {};
    epsilon_move2_rule3(children, startpoint,fbidden);
    cout << "Start Point" << startpoint << endl;
    cout << "Rule-34 children" << children << endl;
    // cout << "Number of children " << children.size() << "\n\n";
    // cout << "Valid Actions " << valid_actions << "\n\n";
}

/* Comparison based elimination */
// void ptss_DSE::eliminate_points_rule12() {
//     // if (!this->is_initialized) {
//     //     throw domain_error("Object not initialized correctly");
//     // }
//     // /* Apply recursively */
//     // all_alloc_t children;
//     // set<alloc_t> fset1 = {this->upper};
//     // set<alloc_t> fset2 = {this->upper};
//     // set<alloc_t> aset;

//     // expand_rule1(children,fset1,aset,0);
//     // aset.clear(); 
//     // expand_rule2(children,fset2,aset,0);
//     // cout << children << endl;

//     // cout << "Eliminated " << children.size() << endl;
// }


/*
 * Whenever the DMR constraint is violated
 * by a point say "init", 
 * Rule1 and Rule3 will (recursively) be expanded
 * to create, sibling solutions which are guaranteed
 * to violate the DMR (Check the slides and the [PAPER])
 * 
 * fbidden : A set that holds all the discarded points.
 * aset    : Children of "init" point expanded according to rule 1 and rule 3
 * aset2   : An auxiliary aset. 
 * 
 */
void apply_action_rule13(all_alloc_t &fbidden,\
                         all_alloc_t &fset,\
                         all_alloc_t &aset,\
                         int actv_id) {
    bool no_child = true;
    int i, j;
    
    /* Expand the children and create the new frontier */
    for (set<alloc_t>::iterator it = fset.begin(); it != fset.end(); it++) {
        
        // cout << "Inserting to fbidden set (activ-"<<actv_id<<"): " << *it << endl;
        // fbidden.insert(*it); // Insert all the elements of frontier set into fbidden set.

        /* Expand all the elements of the frontier */
        // epsilon_move2_rule1(aset,*it,fbidden);
        // epsilon_move2_rule3(aset,*it,fbidden);
        
        alloc_t tmp = *it;
        fbidden.insert(tmp);

        /* Rule 1 Expansion */
        for (int idx = 0; idx < NPH; idx++) {
            alloc_t tmp2  = tmp;   
            if (tmp2[idx] > 1) {
                tmp2[idx]--;
                if (fbidden.count(tmp2) < 1)
                    aset.insert(tmp2);
            }
        }

        /* Rule 3 Expansion */
        for (i = 0; i < NPH; i++) {
            for (j = 0; j < NPH; j++) {
                // cout << "Expr : ("<<i<<","<<j<<") " << ((i != j) && (tmp[i] > 1) && (tmp[j] < M));
                if ((i != j) && (tmp[i] > 1) && (tmp[j] < M)) {
                    /* Transfer a core from i to j */
                    alloc_t new_point = tmp;
                    new_point[i]--;
                    new_point[j]++;
                    /* Check Execution Time before insert */
                    if (compute_execution_time(new_point) > compute_execution_time(tmp)) {
                        if (fbidden.count(new_point) < 1)
                            aset.insert(new_point);
                        /* Also insert other points derived from it */
                        alloc_t new_point2 = new_point;
                        while ((new_point2[i] > 1) && (new_point2[j] < M)) {
                                new_point2[i]--;
                                new_point2[j]++;
                                if (fbidden.count(new_point) < 1)
                                    aset.insert(new_point2);
                        }
                    }
                    // cout << "Blah : "<< "i = "<< i << ",j = " << j << ", tmp[i] = " << tmp[i] << ", tmp[j] = " << tmp[j] << endl; 
                } else {
                    // cout << "i = "<< i << ",j = " << j << ", tmp[i] = " << tmp[i] << ", tmp[j] = " << tmp[j] << endl; 
                }
            }
        }
    }
    fset.clear();
    no_child = no_child && (aset.empty()?true:false);
    
    if (!no_child) {
        apply_action_rule13(fbidden,aset,fset,++actv_id);
    }

}

void ptss_DSE::explore() {
    all_alloc_t tmp1 = {this->lower}, tmp2;
    apply_action_rule13(this->discarded_space,tmp1,tmp2,0);
    cout << "Search points discarded : "<<this->discarded_space.size()<<endl;
}
/********************************************************************/
int main () {

    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    ptss_DSE obj(0.25);
    // cout << obj;
    // obj.display();
    // obj.explore();
    obj.evaluate_all();
    gettimeofday(&t2,NULL);
    double elapsed  = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);

    cout << elapsed << "us\n";

    // alloc_t a = {1,2,3,4,5};
    // alloc_t b = {1,1,3,4,2};
    // cout << is_lt(b,a) << "\n";
}
