#ifndef __PTSS_NLOPT
#define __PTSS_NLOPT

#include <iomanip>
#include <iostream>
#include <vector>
#include <vector>
#include <nlopt.hpp>
#include "ptss_config.hpp"

using namespace std;

/**
 * The optimization space is x = <c,p>
 * where c is the vector of allocations (integer variables)
 * adjoined with the power consumption (p) (continuous variables)
 */
class ptss_constraint_param {
    public :
        vector<double> a;
        vector<double> b;
        double deadline;
        unsigned int sel_idx; /* Index to be selected for separable coefficients */
        ptss_constraint_param(vector<double> a,vector<double> b,double d,unsigned int idx) : a(a), b(b), deadline(d), sel_idx(idx) {}
} ;

double ptss_func_pkp(const std::vector<double> &x, \
                 std::vector<double> &grad, \
                 void *my_func_data);

double ptss_func_et(const std::vector<double> &x, \
                 std::vector<double> &grad, \
                 void *my_func_data);

/* Execution Time Constraints */
double ptss_constraint_exectime(const std::vector<double> &x, \
                                std::vector<double> &grad, \
                                void *param);



/* Minimax Power Objective converted into N additional constraints */
double ptss_constraint_power(const std::vector<double> &x, \
                             std::vector<double> &grad, \
                             void *param);


#endif