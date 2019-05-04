#include <cmath>
#include <vector>
#include <random>
#include "ptss_dse.hpp"
#include "ptss_StatisticalDistribution.hpp"
#include <boost/math/distributions/gamma.hpp>

PTSSETDistType::PTSSETDistType(const phase_t &ph) {
    double loc = 1/(0.02+(0.004*ph.alloc));     /* Compute the mean value of the benchmark (use regression to obtain an estimation) */
    this->loc  = loc; /* Mean Value of the distribution */
    #ifdef DIST_GEV
    double scale = 0.8*loc;
    this->dis = EtDistType(loc,scale);
    #endif
    #ifdef DIST_GAMMA
    double alpha = 3.2;
    double beta  = 2.8;
    this->dis = new EtDistType(alpha,beta);
    #endif
    
    /* Seed the URBG */
    this->u = new uniform_real_distribution<double>(0.0,1.0);
    this->gen.seed(this->rd());
}

PTSSETDistType::~PTSSETDistType() {
    std::cout << "Deleting all distribution objects " << std::endl;
    delete this->u;
    delete this->dis;
}

double PTSSETDistType::sample() {
    double p     = this->u->operator()(gen);
    double invTf = quantile(*this->dis,p);
    return this->loc + invTf;
}