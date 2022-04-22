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

      // Save name of parameter Set
      nameParameterSets.push_back(parSet.getName());
      int countParameters = 0;
      
      auto* parList = &parSet.getEffectiveParameterList();
      for (auto& iPar : *parList) {
	if (iPar.getParameterIndex()!= 9 && iPar.getParameterIndex()!= 10) {
	if (iPar.isFixed()) continue;
	if (!iPar.isEnabled()) continue;}
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
    // Let the MCMCEngine proxy pointer in ProxyEngine point to 'this' intialized object
    like.proxy = this;

    // Set dimensiton for the MCMC sampler
    mcmc.GetProposeStep().SetDim(_nbFitParameters_);

    // Create a fitting parameter vector and intialize it
    Vector p(_nbFitParameters_);
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
