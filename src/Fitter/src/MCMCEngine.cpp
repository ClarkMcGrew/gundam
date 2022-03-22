//
// Created by Shilin on Mar.22, 2022
//

#include "MCMCEngine.h"
#include "MCMCProxyEngine.h"
#include "TSimpleMCMC.H"

#include "JsonUtils.h"

void MCMCEngine::fit() {

    // Update "minimizer" type and algo. Although MCMC is a sampler, we still use
    // minimizer so we don't need to modify the base class FitterEngine
    // This will be printed in evalFit() function
    _minimizerType_ = "MCMC";
    _minimizerAlgo_ = "Metropolis";

    // Get output file name
    std::string outFileName = JsonUtils::fetchValue(_minimizerConfig_, "mcmcOutputFile", "mcmc.root");
    std::string outTreeName = JsonUtils::fetchValue(_minimizerConfig_, "mcmcOutputTree", "accepted");

    // Create output file and a tree to save accepted points
    TFile *outputFile = new TFile(outFileName.c_str(), "recreate");
    TTree *tree = new TTree(outTreeName.c_str(),"Tree of accepted points");

    // Create a mcmc sampler with this MCMCEngine serving as likelihood
    TSimpleMCMC<MCMCProxyEngine> mcmc(tree);
    MCMCProxyEngine& like = mcmc.GetLogLikelihood();
    // Let the MCMCEngine proxy pointer in ProxyEngine point to 'this' intialized object
    like.proxy = this;

    // Set dimensiton for the MCMC sampler
    mcmc.GetProposeStep().SetDim(like.proxy->_nbFitParameters_);

    // Create a fitting parameter vector and intialize it
    Vector p(like.proxy->_nbFitParameters_);
    for (int i = 0; i < _nbFitParameters_; i++) {
        p[i] = _minimizerFitParameterPtr_[i]->getParameterValue();
    }

    // Get MCMC burnin and running info
    int gBurninCycle = JsonUtils::fetchValue(_minimizerConfig_, "mcmcBurninCycle", 3);
    int gBurninLength = JsonUtils::fetchValue(_minimizerConfig_, "mcmcBurninLength", 100);
    int gRunCycle = JsonUtils::fetchValue(_minimizerConfig_, "mcmcRunCycle", 1);
    int gRunLength = JsonUtils::fetchValue(_minimizerConfig_, "mcmcRunLength", 1000);
    bool gSaveBurnin = JsonUtils::fetchValue(_minimizerConfig_, "mcmcSaveBurnin", false);

    // Initializing the mcmc sampler
    mcmc.Start(p, gSaveBurnin);

    // Burnin cycles
    for (int chain = 0; chain < gBurninCycle; ++chain){
      std::cout << "Start Burnin chain " << chain << std::endl;
      mcmc.GetProposeStep().UpdateProposal();
      // Burnin chain in each cycle
      for (int i = 0; i < gBurninLength; ++i) {
        // Run step
        mcmc.Step(gSaveBurnin);

        if(gBurninLength > 100 && !(i%(gBurninLength/100))){
          std::cout << "On step: " << i+1 << "/" << gBurninLength << " "
                    << i*100./gBurninLength << "%" << std::endl;
        }
      }
    }
    std::cout << "Finished burnin chains" << std::endl;

    // Run cycles
    for (int chain = 0; chain < gRunCycle; ++chain){
      std::cout << "Start run chain " << chain << std::endl;
      mcmc.GetProposeStep().UpdateProposal();
      // Run chain in each cycle
      for (int i = 0; i < gRunLength; ++i) {
        // Run step
        mcmc.Step(true);

        if(gRunLength > 100 && !(i%(gRunLength/100))){
          std::cout << "On step: " << i+1 << "/" << gRunLength << " "
                    << i*100./gRunLength << "%" << std::endl;
        }
      }
    }
    std::cout << "Finished Running chains" << std::endl;

    // Save the sampled points to the outputfile
    if (tree) tree->Write();
    if (outputFile) delete outputFile;

}
