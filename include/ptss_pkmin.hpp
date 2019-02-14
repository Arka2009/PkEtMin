#ifndef __PTSS_PKMIN
#define __PTSS_PKMIN

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
#include "ptss_config.hpp"

double ptss_pkmin(alloc2_t init, unsigned int deadline, alloc2_t *opt_dggd);
#endif