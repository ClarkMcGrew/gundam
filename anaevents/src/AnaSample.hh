#ifndef __AnaSample_hh__
#define __AnaSample_hh__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>

#include "AnaEvent.hh"
#include "BinManager.hh"
#include "ColorOutput.hh"
#include "FitStructs.hh"
#include "Likelihoods.hh"
using xsllh::FitBin;

enum DataType
{
    kReset    = -1,
    kMC       = 0,
    kAsimov   = 1,
    kExternal = 2,
    kData     = 3
};

class AnaSample
{
public:
    AnaSample(int sample_id, const std::string& name, const std::string& detector,
              const std::string& binning, TTree* t_data);
    ~AnaSample();

    int GetN() const;
    AnaEvent* GetEvent(const unsigned int evnum);
    void ClearEvents();
    void AddEvent(const AnaEvent& event);
    void ResetWeights();
    void InitEventMap();

    void PrintStats() const;
    void SetNorm(const double val) { m_norm = val; }
    void MakeHistos();

    void SetBinning(const std::string& binning);
    int GetBinIndex(const double D1, const double D2) const;
    std::vector<FitBin> GetBinEdges() const { return m_bin_edges; }

    void SetLLHFunction(const std::string& func_name);
    double CalcLLH() const;
    double CalcChi2() const;

    void FillEventHist(bool reset_weights = false);
    void FillDataHist(int datatype, bool stat_fluc = false);

    void WriteEventHist(TDirectory* dirout, const std::string& bsname);
    void WriteDataHist(TDirectory* dirout, const std::string& bsname);

    double GetNorm() const { return m_norm; }
    int GetSampleID() const { return m_sample_id; }
    std::string GetName() const { return m_name; }
    std::string GetDetector() const { return m_detector; }
    std::string GetBinning() const { return m_binning; }

protected:
    int m_sample_id;
    int m_nbins;
    double m_norm;

    std::string m_name;
    std::string m_detector;
    std::string m_binning;
    std::vector<AnaEvent> m_events;
    std::vector<FitBin> m_bin_edges;

    BinManager bm;
    CalcLLHFunc* m_llh;

    TTree* m_data_tree;
    TH1D* m_hpred;
    TH1D* m_hdata;

    const std::string TAG = color::GREEN_STR + "[AnaSample]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
    const std::string WAR = color::RED_STR + "[WARNING]: " + color::RESET_STR;
};

#endif
