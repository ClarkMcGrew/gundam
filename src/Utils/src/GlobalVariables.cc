//
// Created by Nadrino on 11/02/2021.
//

#include "GlobalVariables.h"

// INIT
bool GlobalVariables::_enableDevMode_{false};
int GlobalVariables::_nbThreads_ = 1;
std::string GlobalVariables::_restoreName_ = "";
std::mutex GlobalVariables::_threadMutex_;
std::map<std::string, bool> GlobalVariables::_boolMap_;
std::vector<TChain*> GlobalVariables::_chainList_;
GenericToolbox::ParallelWorker GlobalVariables::_threadPool_;
TRandom3 GlobalVariables::_prng_;
bool GlobalVariables::_enableCacheManager_{true};

void GlobalVariables::setNbThreads(int nbThreads_){
  _nbThreads_ = nbThreads_;
  _threadPool_.reset();
  _threadPool_.setCheckHardwareCurrency(false);
  _threadPool_.setNThreads(_nbThreads_);
  _threadPool_.setCpuTimeSaverIsEnabled(true);
  _threadPool_.initialize();
}
void GlobalVariables::setPrngSeed(ULong_t seed_){
  _prng_.SetSeed(seed_);
}

void GlobalVariables::setEnableCacheManager(bool enable) {_enableCacheManager_ = enable;}
bool GlobalVariables::getEnableCacheManager() {return _enableCacheManager_;}

void GlobalVariables::setRestoreName(std::string restoreName_) {
  _restoreName_ = restoreName_;
}

bool GlobalVariables::isEnableDevMode(){ return _enableDevMode_; }
const int& GlobalVariables::getNbThreads(){ return _nbThreads_; }
const std::string GlobalVariables::getRestoreName(){ return _restoreName_; }
std::mutex& GlobalVariables::getThreadMutex() { return _threadMutex_; }
std::map<std::string, bool>& GlobalVariables::getBoolMap() { return _boolMap_; }
std::vector<TChain*>& GlobalVariables::getChainList() { return _chainList_; }
GenericToolbox::ParallelWorker &GlobalVariables::getParallelWorker() {
  return _threadPool_;
}
TRandom3& GlobalVariables::getPrng(){
  return _prng_;
}
