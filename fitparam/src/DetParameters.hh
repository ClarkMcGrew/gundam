#ifndef __DetParameters_hh__
#define __DetParameters_hh__

#include <iomanip>
#include <iostream>
#include <string>

#include "AnaFitParameters.hh"
#include "BinManager.hh"
#include "FitStructs.hh"

class DetParameters : public AnaFitParameters
{
public:
    DetParameters(const std::string& name);
    ~DetParameters();

    void InitParameters();
    void InitEventMap(std::vector<AnaSample*>& sample, int mode);
    void ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params);
    void AddDetector(const std::string& det, std::vector<AnaSample*>& v_sample);

private:
    std::map<int, BinManager> m_sample_bm;
    std::map<int, int> m_sample_offset;
    std::vector<int> v_samples;

    const std::string TAG = color::GREEN_STR + "[DetParameters]: " + color::RESET_STR;
};

#endif
