#ifndef __XsecFitter_hh__
#define __XsecFitter_hh__

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <TFile.h>
#include <TGraph.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TRandom3.h>
#include <TVectorT.h>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include "AnaFitParameters.hh"
#include "AnaSample.hh"
#include "ColorOutput.hh"
#include "OptParser.hh"

enum FitType
{
    kAsimovFit,
    kExternalFit,
    kDataFit,
    kToyFit
};

class XsecFitter
{
public:
    XsecFitter(TDirectory* dirout, const int seed);
    XsecFitter(TDirectory* dirout, const int seed, const int num_threads);
    ~XsecFitter();
    double CalcLikelihood(const double* par);
    void InitFitter(std::vector<AnaFitParameters*>& fitpara);
    void FixParameter(const std::string& par_name, const double& value);
    bool Fit(const std::vector<AnaSample*>& samples, int fit_type, bool stat_fluc);
    void ParameterScans(const std::vector<int>& param_list, unsigned int nsteps);

    void SetMinSettings(const MinSettings& ms);
    void SetSeed(int seed);
    void SetZeroSyst(bool flag) { m_zerosyst = flag; }
    void SetNumThreads(const unsigned int num) { m_threads = num; }
    void SetSaveFreq(int freq, bool flag = true)
    {
        m_freq = freq;
        m_save = flag;
    }
    void SetSaveEvents(bool flag = true) { m_save_events = flag; };

    // Declaration of leaf types
    int sample;
    int nutype;
    int topology;
    int reaction;
    int target;
    int sigtype;
    int sam_bin;
    float weight;
    float weightMC;
    std::vector<double> reco_var;
    std::vector<double> true_var;

    void InitOutputTree()
    {
        m_outtree->Branch("nutype", &nutype, "nutype/I");
        m_outtree->Branch("reaction", &reaction, "reaction/I");
        m_outtree->Branch("topology", &topology, "topology/I");
        m_outtree->Branch("target", &target, "target/I");
        m_outtree->Branch("sample", &sample, "sample/I");
        m_outtree->Branch("sigtype", &sigtype, "signal/I");
        m_outtree->Branch("sam_bin", &sam_bin, "sam_bin/I");
        m_outtree->Branch("reco_var", &reco_var);
        m_outtree->Branch("true_var", &true_var);
        m_outtree->Branch("weight", &weight, "weight_post/F");
        m_outtree->Branch("weightMC", &weightMC, "weight_nom/F");
    }

private:
    void GenerateToyData(int toy_type, bool stat_fluc = false);
    double FillSamples(std::vector<std::vector<double>>& new_pars, int datatype = 0);
    void SaveParams(const std::vector<std::vector<double>>& new_pars);
    void SaveEventHist(bool is_final = false);
    void SaveEventTree(std::vector<std::vector<double>>& par_results);
    void SaveChi2();
    void SaveResults(const std::vector<std::vector<double>>& parresults,
                     const std::vector<std::vector<double>>& parerrors);

    ROOT::Math::Minimizer* m_fitter;
    ROOT::Math::Functor* m_fcn;

    TTree* m_outtree;
    TRandom3* rng;
    TDirectory* m_dir;
    bool m_save;
    bool m_save_events;
    bool m_zerosyst;
    int m_threads;
    int m_npar, m_calls, m_freq;
    std::vector<std::string> par_names;
    std::vector<double> par_prefit;
    std::vector<double> vec_chi2_stat;
    std::vector<double> vec_chi2_sys;
    std::vector<double> vec_chi2_reg;
    std::vector<AnaFitParameters*> m_fitpara;
    std::vector<AnaSample*> m_samples;

    MinSettings min_settings;
    const std::string TAG = color::YELLOW_STR + "[XsecFitter]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + color::BOLD_STR + "[ERROR]: " + color::RESET_STR;
};
#endif
