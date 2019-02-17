#include <iomanip>
#include <iostream>
#include <vector>
#include <vector>
#include <nlopt.hpp>
#include "ptss_config.hpp"
#include "ptss_nlopt.hpp"
// #include "ptss_dse.hpp"

using namespace std;

double ptss_func(const std::vector<double> &x, \
                 std::vector<double> &grad, \
                 void *my_func_data)
{
    if (!grad.empty()) {
        for (unsigned int i = 0; i < x.size()-1; i++) {
            grad[i] = 0.0;
        }
        grad[x.size()-1] = 1;
    }
    double f = x[x.size()-1];

    // static int count = 0;
    // cout <<"Obj_" << count++ << "(";
    // for (unsigned int i = 0; i < x.size(); i++)
    //     cout<<x[i]<<",";
    // cout<<") = "<< f << endl;
    return f;
}


/* Execution Time Constraints */
double ptss_constraint_exectime(const std::vector<double> &x, \
                                std::vector<double> &grad, \
                                void *param) {
    ptss_constraint_param *p = reinterpret_cast<ptss_constraint_param*>(param);
    vector<double> a = p->a, b = p->b;
    double d2 = p->deadline;
    if (!grad.empty()) {
        for (unsigned int i = 0; i < x.size()-1; i++) {
            grad[i] = (-a[i])/(x[i]*(a[i]*log(x[i])+b[i])*(a[i]*log(x[i])+b[i]));
        }
        grad[x.size()-1] = 0.0;
    }

    /* Compute the constraint function */
    double f = 0.0;
    double tmp = 0.0;
    for (unsigned int i = 0; i < x.size()-1; i++) {
        tmp = 1/(a[i]*log(x[i])+b[i]);
        // cout << "phase-et("<<i<<"):"<<tmp<<"|";
        f += tmp;
    }
    // cout << endl;
    f = f-d2;
    // static int count = 0;
    // cout << "Execution Time Constraint_"<<count++<<"("<<x<<") = "<<f<<endl;
    return f;
}

/* Minimax Power Objective converted into N additional constraints */
double ptss_constraint_power(const std::vector<double> &x, \
                             std::vector<double> &grad, \
                             void *param) {
    ptss_constraint_param *p = reinterpret_cast<ptss_constraint_param*>(param);
    vector<double> a = p->a, b = p->b;
    // static int count = 0;
    // double d2 = p->deadline;
    unsigned int idx = p->sel_idx;
    if (!grad.empty()) {
        for (unsigned int i = 0; i < x.size()-1; i++) {
            grad[i] = a[i];
        }
        grad[x.size()-1] = -1;
    }
    // double p2 = x[x.size()-1];
    double f = a[idx]*x[idx] + b[idx] - x[x.size()-1];
    // cout << "Power Constraint_"<<count++<<"("<<x<<") = "<<f<<endl;
    return f;
}