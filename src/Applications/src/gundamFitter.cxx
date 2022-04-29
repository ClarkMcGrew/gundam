//
// Created by Nadrino on 01/06/2021.
//

#include "versionConfig.h"

#include "FitterEngine.h"
#include "MCMCEngine.h"
#include "MCMCProxyEngine.h"
#include "JsonUtils.h"
#include "GlobalVariables.h"
#include "GundamGreetings.h"
#ifdef GUNDAM_USING_CACHE_MANAGER
#include "CacheManager.h"
#endif
#include "CmdLineParser.h"
#include "Logger.h"
#include "GenericToolbox.h"
#include "GenericToolbox.Root.h"

#include <string>


LoggerInit([](){
  Logger::setUserHeaderStr("[gundamFitter.cxx]");
} )

int main(int argc, char** argv){

  GundamGreetings g;
  g.setAppName("GundamFitter");
  g.hello();

  // --------------------------
  // Read Command Line Args:
  // --------------------------
  CmdLineParser clParser;

  clParser.addTriggerOption("dry-run", {"--dry-run", "-d"},"Perform the full sequence of initialization, but don't do the actual fit.");
  clParser.addOption("cache", {"-C", "--cache-enabled"}, "Enable the event weight cache");
  clParser.addTriggerOption("generateOneSigmaPlots", {"--one-sigma"}, "Generate one sigma plots");

  clParser.addOption("configFile", {"-c", "--config-file"}, "Specify path to the fitter config file");
  clParser.addOption("nbThreads", {"-t", "--nb-threads"}, "Specify nb of parallel threads");
  clParser.addOption("outputFile", {"-o", "--out-file"}, "Specify the output file");
  clParser.addOption("scanParameters", {"--scan"}, "Enable parameter scan before and after the fit");
  clParser.addOption("toyFit", {"--toy"}, "Run a toy fit");
  clParser.addOption("restoreFile", {"-r", "--restore-file"}, "Specify the mcmc restore file");
  clParser.addOption("prngSeed", {"-s", "--seed"}, "Set the PRNG seed");
  clParser.addOption("fitterEngine", {"-E", "--engine"}, "Set fitter engine override");

  clParser.getOptionPtr("scanParameters")->setAllowEmptyValue(true); // --scan can be followed or not by the number of steps
  clParser.getOptionPtr("toyFit")->setAllowEmptyValue(true); // --toy can be followed or not by the number of steps

  LogInfo << "Usage: " << std::endl;
  LogInfo << clParser.getConfigSummary() << std::endl << std::endl;

  clParser.parseCmdLine(argc, argv);

  LogThrowIf(clParser.isNoOptionTriggered(), "No option was provided.");

  LogInfo << "Provided arguments: " << std::endl;
  LogInfo << clParser.getValueSummary() << std::endl << std::endl;
  LogInfo << clParser.dumpConfigAsJsonStr() << std::endl;

  std::string cacheDefault = "off";
#ifdef GUNDAM_USING_CACHE_MANAGER
  if (Cache::Manager::HasCUDA()) cacheDefault = "on";
#endif
  std::string cacheEnabled = clParser.getOptionVal("cache",cacheDefault);
  if (cacheEnabled != "on") {
      LogInfo << "Cache::Manager disabled" << std::endl;
      GlobalVariables::setEnableCacheManager(false);
  }
  else {
      LogInfo << "Enabling Cache::Manager" << std::endl;
      GlobalVariables::setEnableCacheManager(true);
  }

  auto configFilePath = clParser.getOptionVal("configFile", "");
  LogThrowIf(configFilePath.empty(), "Config file not provided.");

  GlobalVariables::setNbThreads(clParser.getOptionVal("nbThreads", 1));
  LogInfo << "Running the fitter with " << GlobalVariables::getNbThreads() << " parallel threads." << std::endl;

  // Make sure the seed is safely set.  It might be set on the command
  // line.  This uses the UUID of the process as the default seed.
  ULong_t seed = clParser.getOptionVal("prngSeed",0);
  LogInfo << "Pseudorandom sequence initialized with: " << seed << std::endl;
  GlobalVariables::setPrngSeed(seed);

  std::string engineOverride = clParser.getOptionVal("fitterEngine","");

  // --------------------------
  // Initialize the fitter:
  // --------------------------
  LogInfo << "Reading config file: " << configFilePath << std::endl;
  auto jsonConfig = JsonUtils::readConfigFile(configFilePath); // works with yaml

  bool isDryRun = clParser.isOptionTriggered("dry-run");
  bool enableParameterScan = clParser.isOptionTriggered("scanParameters") or JsonUtils::fetchValue(jsonConfig, "scanParameters", false);
  int nbScanSteps = clParser.getOptionVal("scanParameters", 100);

  bool isToyFit = clParser.isOptionTriggered("toyFit");
  int iToyFit = clParser.getOptionVal("toyFit", -1);

  std::string outFileName = configFilePath;
  if( isToyFit ){
    outFileName += "_toyFit";
    if( iToyFit != -1 ){ outFileName += "_" + std::to_string(iToyFit); }
  }
  outFileName = clParser.getOptionVal("outputFile", outFileName + ".root");
  LogWarning << "Creating output file: \"" << outFileName << "\"..." << std::endl;
  TFile* out = TFile::Open(outFileName.c_str(), "RECREATE");


  LogInfo << "Writing runtime parameters in output file..." << std::endl;

  auto restoreFileName = clParser.getOptionVal("restoreFile", "");
  if (not restoreFileName.empty()) {
    GlobalVariables::setRestoreName(restoreFileName.c_str());
    LogInfo << "MCMC continues to sample from: \"" << GlobalVariables::getRestoreName() << "\""<< std::endl;
  }

  // Gundam version?
  TNamed gundamVersionString("gundamVersion", getVersionStr().c_str());
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(out, "gundamFitter"), &gundamVersionString);

  // Command line?
  TNamed commandLineString("commandLine", clParser.getCommandLineString().c_str());
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(out, "gundamFitter"), &commandLineString);

  // Config unfolded ?
  auto unfoldedConfig = jsonConfig;
  JsonUtils::unfoldConfig(unfoldedConfig);
  std::stringstream ss;
  ss << unfoldedConfig << std::endl;
  TNamed unfoldedConfigString("unfoldedConfig", ss.str().c_str());
  GenericToolbox::writeInTFile(GenericToolbox::mkdirTFile(out, "gundamFitter"), &unfoldedConfigString);


  LogInfo << "FitterEngine setup..." << std::endl;

  // Fitter pointer
  FitterEngine *fitter{nullptr};

  std::string engine = JsonUtils::fetchValue(jsonConfig, "engine", "Fitter");
  if (not engineOverride.empty()) engine = engineOverride;

  if (engine == "Fitter") {
    fitter = new FitterEngine();
  }
  else if (engine == "MCMC") {
    fitter = new MCMCEngine();
  }
  else {
      throw std::runtime_error("Invalid fitting engine requested");
  }

  fitter->setConfig(JsonUtils::fetchSubEntry(jsonConfig, {"fitterEngineConfig"}));
  fitter->setSaveDir(GenericToolbox::mkdirTFile(out, "FitterEngine"));
  fitter->setNbScanSteps(nbScanSteps);
  fitter->setEnablePostFitScan(enableParameterScan);

  if( isToyFit ){ fitter->getPropagator().setThrowAsimovToyParameters(true); }

  fitter->initialize();

  fitter->updateChi2Cache();
  LogInfo << "Initial χ² = " << fitter->getChi2Buffer() << std::endl;
  LogInfo << "Initial χ²(stat) = " << fitter->getChi2StatBuffer() << std::endl;

  // --------------------------
  // Pre-fit:
  // --------------------------

  // LLH Visual Scan
  if( clParser.isOptionTriggered("generateOneSigmaPlots") or JsonUtils::fetchValue(jsonConfig, "generateOneSigmaPlots", false) ) fitter->generateOneSigmaPlots("preFit");
  if( clParser.isOptionTriggered("scanParameters") or JsonUtils::fetchValue(jsonConfig, "scanParameters", false) ) fitter->scanParameters(nbScanSteps, "preFit/scan");

  // Plot generators
  if( JsonUtils::fetchValue(jsonConfig, "generateSamplePlots", true) ) fitter->generateSamplePlots("preFit/samples");


  // --------------------------
  // Run the fitter:
  // --------------------------
  if( not isDryRun and JsonUtils::fetchValue(jsonConfig, "fit", true) ){
    fitter->fit();
  }

  LogWarning << "Closing output file \"" << out->GetName() << "\"..." << std::endl;
  out->Close();
  LogInfo << "Closed." << std::endl;

  // --------------------------
  // Goodbye:
  // --------------------------
  g.goodbye();

}
