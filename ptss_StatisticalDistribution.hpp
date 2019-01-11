#ifndef __STATISTICS_H
#define __STATISTICS_H

#include <cmath>
#include <vector>

class StatisticalDistribution {
 public:
  StatisticalDistribution();
  virtual ~StatisticalDistribution();

  // Distribution functions
  virtual double pdf(const double& x) const = 0;
  virtual double cdf(const double& x) const = 0;
  virtual double sf(const double& x) const = 0;

  // Inverse cumulative distribution functions (aka the quantile function)
  virtual double ppf(const double& quantile) const = 0;
  virtual double isf(const double& quantile) const = 0;
  
  // Descriptive stats
  virtual double mean() const = 0;
  virtual double var() const = 0;
  virtual double stdev() const = 0;

  // Obtain a sequence of random draws from this distribution
  virtual void random_draws(const std::vector<double>& uniform_draws,
                            std::vector<double>& dist_draws) = 0;
};

#endif