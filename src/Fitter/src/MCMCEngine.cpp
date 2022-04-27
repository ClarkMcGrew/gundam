//
// Created by Shilin on Mar.22, 2022
//

#include "MCMCEngine.h"
#include "MCMCProxyEngine.h"

#include "GlobalVariables.h"
#include "Logger.h"
#include "JsonUtils.h"

// Tell TSimpleMCMC.H how to send output.
#define MCMC_DEBUG_LEVEL 1
#define MCMC_DEBUG(level) if (level <= (MCMC_DEBUG_LEVEL)) LogInfo
#define MCMC_ERROR (LogInfo << "ERROR: ")

#include "TSimpleMCMC.H"

LoggerInit([]{
  Logger::setUserHeaderStr("[MCMCEngine]");
})

void MCMCEngine::fit() {

    // Update "minimizer" type and algo. Although MCMC is a sampler, we still
    // use minimizer so we don't need to modify the base class FitterEngine
    // This will be printed in evalFit() function
    _minimizerType_ = "MCMC";
    _minimizerAlgo_ = "Metropolis";

    LogThrowIf(_useNormalizedFitSpace_,"NormalizedFitSpace used in MCMC, it must be disabled.");

    // Get output file name
    std::string outFileName = JsonUtils::fetchValue(_minimizerConfig_, "mcmcOutputFile", "mcmc.root");
    std::string outTreeName = JsonUtils::fetchValue(_minimizerConfig_, "mcmcOutputTree", "accepted");

    // Check for restore file
    TFile *restoreFile{nullptr};
    TTree *restoreTree{nullptr};
    std::string restoreName = GlobalVariables::getRestoreName();
    LogThrowIf(restoreName == outFileName, "Writing the restore file! Need to change output file name in configBanffFit.yaml");
    if (not restoreName.empty()) {
      LogInfo << "Restore from: " << restoreName << std::endl;
      restoreFile = new TFile(restoreName.c_str(), "old");
      restoreTree = (TTree*) restoreFile->Get(outTreeName.c_str());
      LogThrowIf( not restoreTree, "No restore tree not found! Make sure outTreeName is the same as previous chain in configBanffFit.yaml");
    }

    // Create output file and a tree to save accepted points
    LogInfo << "Adding MCMC tree " << outTreeName
            << " to file " << gFile->GetName()
            << std::endl;
    TTree *tree = new TTree(outTreeName.c_str(),"Tree of accepted points");

    // Storing parameter names
    TTree *parName = new TTree("parameterSets", "Tree of parameterSets");
    std::vector<std::string> nameParameterSets;
    std::vector<int> nParameters;
    std::vector<int> parameter_index;
    std::vector<double> parameter_prior;
    std::vector<double> parameter_sigma;
    std::vector<std::string> parameter_name;
    parName->Branch("nameParameterSets", &nameParameterSets);
    parName->Branch("nParameters", &nParameters);
    parName->Branch("parameter_index", &parameter_index);
    parName->Branch("parameter_name", &parameter_name);
    parName->Branch("parameter_prior", &parameter_prior);
    parName->Branch("parameter_sigma", &parameter_sigma);

    for (const auto& parSet: _propagator_.getParameterSetsList()) {
      if (not parSet.isEnabled()) continue;
      LogThrowIf(parSet.isUseEigenDecompInFit(), "Eigen Decomp is used in MCMC! Don't run MCMC with Eigen Decomp");

      // Save name of parameter Set
      nameParameterSets.push_back(parSet.getName());
      int countParameters = 0;

      auto* parList = &parSet.getEffectiveParameterList();
      for (auto& iPar : *parList) {
	if (iPar.getParameterIndex()!= 9 && iPar.getParameterIndex()!= 10) {
            if (iPar.isFixed()) continue;
            if (!iPar.isEnabled()) continue;
        }
	countParameters ++;
	parameter_index.push_back(iPar.getParameterIndex());
	parameter_name.push_back(iPar.getTitle());
	parameter_prior.push_back(iPar.getPriorValue());
	parameter_sigma.push_back(iPar.getStdDevValue());
      }
      nParameters.push_back(countParameters);
    }
    parName->Fill();

    parName->Write();
    // End Storing parameter name informations

    // Create a mcmc sampler with this MCMCEngine serving as likelihood
    TSimpleMCMC<MCMCProxyEngine> mcmc(tree);
    MCMCProxyEngine& like = mcmc.GetLogLikelihood();
    // Let the MCMCEngine proxy pointer in ProxyEngine point to 'this'
    // initialized object
    like.proxy = this;

    // Set dimensiton for the MCMC sampler
    mcmc.GetProposeStep().SetDim(_nbFitParameters_);

    // Create a fitting parameter vector and initialize it
    Vector p(_nbFitParameters_);
    int count = 0;
    for(const FitParameter* par : _minimizerFitParameterPtr_ ){
        p[count] = par->getPriorValue();
        switch (par->getPriorType()) {
        case PriorType::Flat: {
            // Gundam uses flat to mean "free", so this doesn't use a Uniform
            // step between bounds.
            double step = std::min(par->getStepSize(),par->getStdDevValue());
            if (step <= std::abs(1E-10*par->getPriorValue())) {
                step = std::max(par->getStepSize(),par->getStdDevValue());
            }
            mcmc.GetProposeStep().SetGaussian(count,step);
            break;
        }
        default:
            mcmc.GetProposeStep().SetGaussian(count,par->getStdDevValue());
            break;
        }
        ++count;
    }

    // Get MCMC burnin and running info
    int burninCycle = JsonUtils::fetchValue(_minimizerConfig_, "mcmcBurninCycle", 0);
    int burninLength = JsonUtils::fetchValue(_minimizerConfig_, "mcmcBurninLength", 0);
    int runCycle = JsonUtils::fetchValue(_minimizerConfig_, "mcmcRunCycle", 1);
    int runLength = JsonUtils::fetchValue(_minimizerConfig_, "mcmcRunLength", 1000);
    bool saveBurnin = JsonUtils::fetchValue(_minimizerConfig_, "mcmcSaveBurnin", false);

    // Initializing the mcmc sampler
    mcmc.Start(p, saveBurnin);

    // Restore the chain if exist
    if (restoreTree) {
      mcmc.Restore(restoreTree);
      delete restoreFile;
      LogInfo << "State Restored" << std::endl;
    }

    // Burnin cycles
    for (int chain = 0; chain < burninCycle; ++chain){
      LogInfo << "Start Burnin chain " << chain << std::endl;
      // Reset the covariance to the initial state.  This forgets the path
      // of the previous cycle.
      mcmc.GetProposeStep().ResetProposal();
      // Burnin chain in each cycle
      for (int i = 0; i < burninLength; ++i) {
        // Run step
        mcmc.Step(saveBurnin);

        if(burninLength > 100 && !(i%(burninLength/100))){
          LogInfo << "Finished burn-in step: "
                  << i << "/" << burninLength << " "
                  << i*100./burninLength << "%" << std::endl;
        }
      }
    }
    LogInfo << "Finished burnin chains" << std::endl;

    // Run cycles
    for (int chain = 0; chain < runCycle; ++chain){
      LogInfo << "Start run chain " << chain << std::endl;
      // Update the covariance with the steps from the last cycle.  This
      // starts a new "reversible-chain".
      mcmc.GetProposeStep().UpdateProposal();
      // Run chain in each cycle
      for (int i = 0; i < runLength; ++i) {
        // Run step
        mcmc.Step(true);
        if(runLength > 100 && !(i%(runLength/100))){
          LogInfo << "Finished step: " << i << "/" << runLength << " "
                  << i*100./runLength << "%" << std::endl;
        }
      }
    }
    LogInfo << "Finished Running chains" << std::endl;

    // Save the final state. This is needed so that the current chain with the
    // current adaptive proposal can be restored and continued.
    mcmc.SaveStep();

    // Save the sampled points to the outputfile
    if (tree) tree->Write();

}

MCMCEngine::MCMCEngine() {}

// An MIT Style License

// Copyright (c) 2022 Shilin Liu (Preliminary: Checking with Shilin)

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
