//
// Created by Nadrino on 11/02/2021.
//

#ifndef GUNDAM_GLOBALVARIABLES_H
#define GUNDAM_GLOBALVARIABLES_H

#include "GenericToolbox.ParallelWorker.h"

#include <TTree.h>
#include <TChain.h>
#include <TRandom3.h>

#include <map>
#include <mutex>


class GlobalVariables{

public:

  // Setters
  static void setNbThreads(int nbThreads_);
  static void setPrngSeed(ULong_t seed_);
  static void setEnableCacheManager(bool enable = true);
  static void setRestoreName(std::string restoreName_);

  // Getters
  static bool isEnableDevMode();
  static const int& getNbThreads();
  static const std::string getRestoreName();
  static std::mutex& getThreadMutex();
  static std::map<std::string, bool>& getBoolMap();
  static std::vector<TChain*>& getChainList();
  static GenericToolbox::ParallelWorker &getParallelWorker();
  static TRandom* getPrng();
  static bool getEnableCacheManager();

private:

  static bool _enableDevMode_;
  static int _nbThreads_;
  static std::string _restoreName_;
  static std::mutex _threadMutex_;
  static std::map<std::string, bool> _boolMap_;
  static std::vector<TChain*> _chainList_;
  static GenericToolbox::ParallelWorker _threadPool_;
  static bool _enableCacheManager_;

};

#endif // XSLLHFITTER_GLOBALVARIABLES_H
