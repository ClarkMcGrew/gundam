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
