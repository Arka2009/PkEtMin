#ifndef __STATISTICS_H
#define __STATISTICS_H

#include <cmath>
#include <vector>
#include <random>
#include <boost/math/distributions/gamma.hpp>
#include "ptss_dse.hpp"

using namespace std;

/* Distribution Selector */ 
#define DIST_GAMMA
// #define DIST_GEV
// #define DIST_FNORMAL
// #define DIST_WEIBULL
#define HISTBINSIZE 1000


#ifdef DIST_GEV
typedef extreme_value_distribution<double> EtDistType;
#endif
#ifdef DIST_GAMMA
typedef boost::math::gamma_distribution<double> EtDistType;
#endif
#ifdef DIST_WEIBULL
typedef weibull_distribution<double> EtDistType;
#endif

/* Generate a distrbution based on the resource allocation - n number of cores */
class PTSSETDistType {
    private:
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen;       //Standard mersenne_twister_engine seeded with rd()
        EtDistType *dis;
        std::uniform_real_distribution<double> *u;
        double loc;
    public:
        // PTSSETDistType() {this->gen.seed(this->rd());};
        PTSSETDistType(const phase_t&);
        ~PTSSETDistType();
        double sample();
};

/* Add two statistical distribution by convolving their PDFs * /
void add_update(PTSSETDistType&,const PTSSETDistType&);

/* Add two independant distributions */
// HistDistType operator+(const HistDistType &, const HistDistType &);
#endif