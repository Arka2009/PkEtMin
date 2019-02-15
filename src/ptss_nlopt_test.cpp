#include <iomanip>
#include <iostream>
#include <vector>
#include <vector>
#include <nlopt.hpp>
#include "ptss_config.hpp"
#include "ptss_nlopt.hpp"


int main() {
    nlopt::opt opt(nlopt::LD_MMA, 3);

    /* Set the Box Constraints */
    std::vector<double> lb(3);
    std::vector<double> ub(3);
    lb[0] = 3; lb[1] = 3, lb[2] = 0.0;
    ub[0] = ULIM; ub[1] = ULIM, ub[2] = HUGE_VAL;
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    
    /* Set the Parameters for inequality constraints */
    vector<double> a_et  = {0.009683975849832148,0.0024794976904501873};
    vector<double> b_et  = {-0.002893833328927599,-0.0007496878381248601};
    vector<double> a_p   = {1.5027499999999998,2.299642857142857};
    vector<double> b_p   = {4.669250000000002,3.7145476190476217};
    double deadline      = 700;
    vector<ptss_constraint_param> param;
    param.push_back(ptss_constraint_param(a_et,b_et,deadline,-1));
    param.push_back(ptss_constraint_param(a_p,b_p,deadline,0));
    param.push_back(ptss_constraint_param(a_p,b_p,deadline,1));
    opt.add_inequality_constraint(ptss_constraint_exectime, &param[0], 1e-8);
    opt.add_inequality_constraint(ptss_constraint_power, &param[1], 1e-8);
    opt.add_inequality_constraint(ptss_constraint_power, &param[2], 1e-8);

    /* Set the objective function */
    opt.set_min_objective(ptss_func, NULL);

    // opt.add_inequality_constraint(ptss_constraint_exectime, &data[1], 1e-8);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(3);
    x[0] = 3.234; x[1] = 5.678; x[2] = (a_p[0]*x[0]+b_p[0])+(a_p[1]*x[1]+b_p[1]);
    double minf;

    try{
        // nlopt::result result = 
        opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
            << std::setprecision(10) << minf << std::endl;
        // std::cout << "found minimum after " << count <<" evaluations\n";
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
}
