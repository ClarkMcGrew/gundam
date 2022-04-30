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

  // A function that performs the "fit".  This creates an MCMC sampler using
  // the likelihood and generates an MCMC chain based on settings in the JSON
  // configuration.  The chain will be saved into the output file.
  void fit();

 // Return true if the current fit parameters are at a valid point.
  bool ValidParameters();

protected:

  // Transform the accepted point into the fPoint vector so it can be saved.
  // The accepted point may be in a normalized fitting parameter value, or in
  // an eigenvalue decomposed state of the parameter space.  All of the
  // transformation applied to the accepted point MUST be linear, or the
  // metric of the likelihood has changed in a bad way and the chain isn't
  // valid.
  void FillPoints();

  // The point value that is associated with the last call to the likelihood.
  std::vector<double> fPoint;
};
#endif

// An MIT Style License

// Copyright (c) 2022 Shilin Liu

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/cmake/gundam-build.sh"
// End:
