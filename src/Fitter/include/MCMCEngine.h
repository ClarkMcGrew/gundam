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
  ~MCMCEngine();

  // A funtion that returns the dimension of fitting parameters.
 std::size_t GetDim() {
    return _nbFitParameters_;
  }
  
  // A function that performs the fit.
  // Defines a mcmc sampler and do the fit
  // MCMC sampler runs based on the json input information
  void fit(); 
};



#endif
