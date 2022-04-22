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

	// FitterEngine::evalFit(const double* parArray_) does 2 things to calculate likelihood
	// 1. It updates fit parameters _minimizerFitParameterPtr_ value with input parArray_
	// 2. It calculates likelihood by calling FitterEngine::updateChi2Cache()
	//    updateChi2Cache() does the "magic" to calculate _chi2Buffer_ with updated fit 
	//    parameters _minimizerFitParameterPtr_.
	// FitterEngine::evalFit(const double* parArray_) does a lot of things to output timing info
	// --- FitterEngine.cpp:653 _convergenceMonitor_, does it affect likelihood calculation?
	// --- FitterEngine.cpp:668 _chi2HistoryTree_, does it affect likelihood calculation?
        double likelihood = proxy->evalFit(&fCurrParams[0]);
	return -0.5 * likelihood;
    }


};

#endif
