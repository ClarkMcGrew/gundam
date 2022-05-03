//
// Created by Shilin on 03/18/2022
//

#ifndef GUNDAM_MCMCPROXYENGINE_H
#define GUNDAM_MCMCPROXYENGINE_H

#include "MCMCEngine.h"

struct MCMCProxyEngine {
    MCMCEngine *proxy;

    // An operator that is needed as a likelihood calculator for TSimpleMCMC.H
    // Will simply call FitterEngine::evalFit() method
    double operator() (const std::vector<double>& fCurrParams) {

        // FitterEngine::evalFit(const double* parArray_) does multiple
        // things to calculate the likelihood
        //
        // 1. It updates fit parameters _minimizerFitParameterPtr_ value with
        //    input parArray_
        // 2. It calculates likelihood by calling
        //    FitterEngine::updateChi2Cache() updateChi2Cache() does the
        //    "magic" to calculate _chi2Buffer_ with updated fit parameters.
        //
        // After the likelihood the parameter values have been updated by the
        // inversion of the eigenvalue decomposition, so the parameter values
        // change during the evaluation.  Note that this does not change the
        // values in parArray_.
        //
        // FitterEngine::evalFit(const double* parArray_) does a lot of things
        // to output timing info
        //
        // --- FitterEngine.cpp:653 _convergenceMonitor_, does it affect
        //     the likelihood calculation?
        //
        // --- FitterEngine.cpp:668 _chi2HistoryTree_, does it affect
        //     the likelihood calculation?
        double chiSquared = proxy->evalFit(fCurrParams.data());

        // If the parameters are invalid an infinite negative value, in
        // otherwords, a zero probability.  This will happen when one of the
        // parameters is out of a valid rage after eigenvalue decomposition.
        if (!proxy->ValidParameters()) {
            return -std::numeric_limits<double>::infinity();
        }

        return -0.5 * chiSquared;
    }
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
