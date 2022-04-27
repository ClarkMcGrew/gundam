//
// Created by Shilin on 03/16/2022
//

#ifndef GUNDAM_MCMCENGINE_H
#define GUNDAM_MCMCENGINE_H

#include "FitterEngine.h"

class MCMCEngine : public FitterEngine {

public:
  // Keep the same, let the code do whatever is needed to intialize
  MCMCEngine();
  ~MCMCEngine() = default;

  // A function that performs the fit.
  // Defines a mcmc sampler and do the fit
  // MCMC sampler runs based on the json input information
  void fit();

protected:
  // MCMCEngine is estimating PDF [P(x)dx], and a transformation can change
  // the metric in a complicated way.  Transformations like renormalizing the
  // parameters, and eigen-value decomposition should not be done, or the
  // output results are very hard to interpret [You will get a point cloud,
  // but it will require a lot of internal information to interpret].
  bool allowFitSpaceTransformations() const override {return false;}
};
#endif
