//
// Created by Shilin on Mar.22, 2022
//

#include "MCMCEngine.h"
#include "MCMCProxyEngine.h"

#include "GlobalVariables.h"
#include "Logger.h"
#include "JsonUtils.h"

// Tell TSimpleMCMC.H how much output to use and where to send it.
#define MCMC_DEBUG_LEVEL 1
#define MCMC_DEBUG(level) if (level <= (MCMC_DEBUG_LEVEL)) LogInfo
#define MCMC_ERROR (LogInfo << "ERROR: ")

#include "TSimpleMCMC.H"

LoggerInit([]{
  Logger::setUserHeaderStr("[MCMC]");
})

void MCMCEngine::FillPoints() {
    int count = 0;
    for (const FitParameterSet& parSet: _propagator_.getParameterSetsList()) {
        for (const FitParameter& iPar : parSet.getParameterList()) {
            fPoint[count++] = iPar.getParameterValue();
        }
    }
}

bool MCMCEngine::ValidParameters() {
    for (const FitParameterSet& parSet: _propagator_.getParameterSetsList()) {
        for (const FitParameter& iPar : parSet.getParameterList()) {
            if (iPar.isFixed()) continue;
            if (!iPar.isEnabled()) continue;
            if (std::isfinite(iPar.getMinValue())
                && iPar.getParameterValue() < iPar.getMinValue()) return false;
            if (std::isfinite(iPar.getMaxValue())
                && iPar.getParameterValue() > iPar.getMaxValue()) return false;
        }
    }
    return true;
}

void MCMCEngine::fit() {
    // Update "minimizer" type and algo. Although MCMC is a sampler, we still
    // use minimizer so we don't need to modify the base class FitterEngine
    // This will be printed in evalFit() function
    _minimizerType_ = "MCMC";
    _minimizerAlgo_ = "Metropolis";

    // Get output file name
    std::string outTreeName
        = JsonUtils::fetchValue(_minimizerConfig_, "mcmcOutputTree", "MCMC");

    // Create output tree in the existing
    LogInfo << "Adding MCMC tree " << outTreeName
            << " to file " << gFile->GetName()
            << std::endl;
    TTree *tree = new TTree(outTreeName.c_str(),"Tree of accepted points");

    // Storing parameter names
    TTree *parameterSetsTree = new TTree("parameterSets",
                                         "Tree of Parameter Set Information");
    std::vector<std::string> parameterSetNames;
    std::vector<int> parameterSetOffsets;
    std::vector<int> parameterSetCounts;
    std::vector<int> parameterIndex;
    std::vector<bool> parameterFixed;
    std::vector<bool> parameterEnabled;
    std::vector<std::string> parameterName;
    std::vector<double> parameterPrior;
    std::vector<double> parameterSigma;
    std::vector<double> parameterMin;
    std::vector<double> parameterMax;
    parameterSetsTree->Branch("parameterSetNames", &parameterSetNames);
    parameterSetsTree->Branch("parameterSetOffsets", &parameterSetOffsets);
    parameterSetsTree->Branch("parameterSetCounts", &parameterSetCounts);
    parameterSetsTree->Branch("parameterIndex", &parameterIndex);
    parameterSetsTree->Branch("parameterFixed", &parameterFixed);
    parameterSetsTree->Branch("parameterEnabled", &parameterEnabled);
    parameterSetsTree->Branch("parameterName", &parameterName);
    parameterSetsTree->Branch("parameterPrior", &parameterPrior);
    parameterSetsTree->Branch("parameterSigma", &parameterSigma);
    parameterSetsTree->Branch("parameterMin", &parameterMin);
    parameterSetsTree->Branch("parameterMax", &parameterMax);

    // Pull all of the parameters out of the parameter sets and save the
    // names, index, priors, and sigmas to the output.  The order in these
    // vectors define the parameters that are saved in the Points vector in
    // the output file.  Note that these parameters do not define the meaning
    // of the parameters in the call to the likelihood function.  Those
    // parameters are defined by the _minimizerFitParameterPtr_ vector.
    for (const FitParameterSet& parSet: _propagator_.getParameterSetsList()) {
      // Save name of parameter set
      parameterSetNames.push_back(parSet.getName());
      parameterSetOffsets.push_back(parameterIndex.size());

      int countParameters = 0;
      for (const FitParameter& iPar : parSet.getParameterList()) {
          ++countParameters;
          parameterIndex.push_back(iPar.getParameterIndex());
          parameterFixed.push_back(iPar.isFixed());
          parameterEnabled.push_back(iPar.isEnabled());
          parameterName.push_back(iPar.getTitle());
          parameterPrior.push_back(iPar.getPriorValue());
          parameterSigma.push_back(iPar.getStdDevValue());
          parameterMin.push_back(iPar.getMinValue());
          parameterMax.push_back(iPar.getMaxValue());
      }
      parameterSetCounts.push_back(countParameters);
    }
    parameterSetsTree->Fill();
    parameterSetsTree->Write();
    // End Storing parameter name informations

    // Create a vector to store the full output of the MCMCEngine chain.  This
    // will be filled based on the accepted states coming out of the MCMC
    // sampler.
    fPoint.resize(parameterName.size());
    tree->Branch("Points",&fPoint);

    TSimpleMCMC<MCMCProxyEngine> mcmc(tree);
    MCMCProxyEngine& like = mcmc.GetLogLikelihood();
    // Let the MCMCEngine proxy pointer in ProxyEngine point to 'this'
    // initialized object
    like.proxy = this;

    mcmc.GetProposeStep().SetDim(_minimizerFitParameterPtr_.size());

    // Create a fitting parameter vector and initialize it.  No need to worry
    // about resizing it or it moving, so be lazy and just use push_back.
    Vector p;
    for (const FitParameter* par : _minimizerFitParameterPtr_ ) {
        if (_useNormalizedFitSpace_) {
            // Changing the boundaries, change the value/step size?
            double val
                = FitParameterSet::toNormalizedParValue(
                    par->getParameterValue(), *par);
            double step
                = FitParameterSet::toNormalizedParRange(
                    par->getStepSize(), *par);
            mcmc.GetProposeStep().SetGaussian(p.size(),step);
            p.push_back(val);
            continue;
        }
        p.push_back(par->getPriorValue());
        switch (par->getPriorType()) {
        case PriorType::Flat: {
            // Gundam uses flat to mean "free", so this doesn't use a Uniform
            // step between bounds.
            double step = std::min(par->getStepSize(),par->getStdDevValue());
            if (step <= std::abs(1E-10*par->getPriorValue())) {
                step = std::max(par->getStepSize(),par->getStdDevValue());
            }
            step /= std::sqrt(_nbFitParameters_);
            mcmc.GetProposeStep().SetGaussian(p.size()-1,step);
            break;
        }
        default: {
            double step = par->getStdDevValue();
            mcmc.GetProposeStep().SetGaussian(p.size()-1,step);
            break;
        }
        }
    }

    LogInfo << "Chain is using " << p.size()
            << " effective dimensions for "
            << parameterName.size() << " parameters"
            << std::endl;

    // Set up the correlations in the priors.
    int count1 = 0;
    for (const FitParameter* par1 : _minimizerFitParameterPtr_ ) {
        ++count1;
        const FitParameterSet* set1 = par1->getParSetRef();
        if (!set1) {
            LogInfo << "Parameter set reference is not defined for"
                    << " " << par1->getName()
                    << std::endl;
            continue;
        }
        if (set1->isUseEigenDecompInFit()) {
            continue;
        }

        int count2 = 0;
        for (const FitParameter* par2 : _minimizerFitParameterPtr_ ) {
            ++count2;
            const FitParameterSet* set2 = par2->getParSetRef();
            if (!set2) {
                LogInfo << "Parameter set reference is not defined for"
                        << " " << par1->getName()
                        << std::endl;
                continue;
            }
            if (set2->isUseEigenDecompInFit()) {
                continue;
            }

            if (set1 != set2) continue;
            int in1 = par1->getParameterIndex();
            int in2 = par2->getParameterIndex();
            if (in2 <= in1) continue;
            const std::shared_ptr<TMatrixDSym>& corr
                = set1->getPriorCorrelationMatrix();
            if (!corr) continue;
            double correlation = (*corr)(in1,in2);
            // Don't impose very small correlations, let them be discovered.
            if (std::abs(correlation) < 0.01) continue;
            // Complain about large correlations.  When a correlation is this
            // large, then the user should (but probably won't) rethink the
            // parameter definitions!
            if (std::abs(correlation) > 0.98) {
                LogInfo << "VERY LARGE CORRELATION (" << correlation
                        << ") BETWEEN"
                        << " " << set1->getName() << "/" << par1->getName()
                        << " & " << set2->getName() << "/" << par2->getName()
                        << std::endl;
            }
            mcmc.GetProposeStep().SetCorrelation(count1-1,count2-1,
                                                 (*corr)(in1,in2));
        }
    }

    // Get MCMC burnin parameters.  Each burnin discards previous information
    // about the posterior and reset to the initial state (but starts from
    // the last accepted point.  The burnin will be skipped if the state
    // has been restored from a file.  The burnin can be skipped in favor of
    // discarding the initial parts of the MCMC chain.
    int burninCycle = JsonUtils::fetchValue(_minimizerConfig_,
                                            "mcmcBurninCycle", 0);
    int burninLength = JsonUtils::fetchValue(_minimizerConfig_,
                                             "mcmcBurninLength", 0);
    bool saveBurnin = JsonUtils::fetchValue(_minimizerConfig_,
                                            "mcmcSaveBurnin", true);

    // Get the MCMC chain parameters.  A run is broken into "mini-Chains"
    // called a "cycle" where the posterior covariance information is updated
    // after each mini-chain.  The cycle will have "mcmcRunLength" steps.
    int runCycle = JsonUtils::fetchValue(_minimizerConfig_,
                                         "mcmcRunCycle", 5);
    int runLength = JsonUtils::fetchValue(_minimizerConfig_,
                                          "mcmcRunLength", 100000);

    // Set the window to calculate the current acceptance value over.  If this
    // is set to short, the step size will fluctuate.  If this is set to long,
    // the step size won't be adjusted to match the target acceptance.
    int window = JsonUtils::fetchValue(_minimizerConfig_,
                                       "mcmcWindow", 5000);

    // Fill the initial point.
    FillPoints();

    // Initializing the mcmc sampler
    mcmc.Start(p, saveBurnin);
    mcmc.GetProposeStep().SetAcceptanceWindow(window);

    // Restore the chain if exist
    std::string restoreName = GlobalVariables::getRestoreName();
    if (!restoreName.empty()) {
        // Check for restore file
        LogInfo << "Restore from: " << restoreName << std::endl;
        std::unique_ptr<TFile> restoreFile
            (new TFile(restoreName.c_str(), "old"));
        if (!restoreFile) {
            LogInfo << "File to restore was not openned: "
                    << restoreName << std::endl;
            std::runtime_error("Old state file not open");
        }
        std::string treeName = "FitterEngine/fit/" + outTreeName;
        TTree* restoreTree = (TTree*) restoreFile->Get(treeName.c_str());
        if (!restoreTree) {
            LogInfo << "Tree to restore state is not found"
                    << treeName << std::endl;
            std::runtime_error("Old state tree not open");
        }
        mcmc.Restore(restoreTree);
        LogInfo << "State Restored" << std::endl;
    }
    else {
        // Burnin cycles
        for (int chain = 0; chain < burninCycle; ++chain){
            LogInfo << "Start Burnin chain " << chain << std::endl;
            // Reset the covariance to the initial state.  This forgets the
            // path of the previous cycle.
            mcmc.GetProposeStep().ResetProposal();
            // Burnin chain in each cycle
            for (int i = 0; i < burninLength; ++i) {
                // Run step
                if (mcmc.Step(false)) FillPoints();
                if (saveBurnin) mcmc.SaveStep(burninLength <= (i+1));

                if(burninLength > 100 && !(i%(burninLength/100))){
                    LogInfo << "Burn-in: " << chain
                            << " step: " << i << "/" << burninLength << " "
                            << i*100./burninLength << "%"
                            << " Trials: "
                            << mcmc.GetProposeStep().GetSuccesses()
                            << "/" << mcmc.GetProposeStep().GetTrials()
                            << " (" << mcmc.GetProposeStep().GetSigma() << ")"
                            << std::endl;
                }
            }
        }
        LogInfo << "Finished burnin chains" << std::endl;
    }

    // Run cycles
    for (int chain = 0; chain < runCycle; ++chain){
      LogInfo << "Start run chain " << chain << std::endl;
      // Update the covariance with the steps from the last cycle.  This
      // starts a new "reversible-chain".
      mcmc.GetProposeStep().UpdateProposal();
      // Run chain in each cycle
      for (int i = 0; i < runLength; ++i) {
          // Run step, but do not save the step.  The step isn't saved so the
          // accepted step can be copied into the points (which will have
          // any decomposition removed).
          if (mcmc.Step(false)) FillPoints();
          // Save the step.  Check to see if this is the last step of the run,
          // and if it is, then save the full state.
          mcmc.SaveStep(false);
          if(runLength > 100 && !(i%(runLength/100))){
              LogInfo << "Chain: " << chain
                      << " step: " << i << "/" << runLength << " "
                      << i*100./runLength << "%"
                      << " Trials: "
                      << mcmc.GetProposeStep().GetSuccesses()
                      << "/" << mcmc.GetProposeStep().GetTrials()
                      << " (" << mcmc.GetProposeStep().GetAcceptance()
                      << ":" << mcmc.GetProposeStep().GetSigma() << ")"
                      << std::endl;
          }
      }
      // Save the final state.  This step should be skipped when analyzing the
      // chain, the steps can be identified since the covariance is not empty.
      LogInfo << "Chain: " << chain << " complete"
              << " Run Length: " << runLength
              << " -- Saving state"
              << std::endl;
      mcmc.SaveStep(true);
    }
    LogInfo << "Finished Running chains" << std::endl;

    // Save the sampled points to the outputfile
    if (tree) tree->Write();

}

MCMCEngine::MCMCEngine() {}

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
