#ifndef EXTRAPARAMETERS_HH
#define EXTRAPARAMETERS_HH

#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

#include "AnaFitParameters.hh"
#include "FitStructs.hh"

class ExtraParameters : public AnaFitParameters
{
    public:
        ExtraParameters(const std::string& name);

        void InitParameters();
        void InitFixed();
        void InitEventMap(std::vector<AnaSample*>& sample, int mode);
        void ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params);
        void AddSample(std::vector<AnaSample*>& v_sample);
        std::vector<double> CalcConstraint(const std::vector<AnaSample*>& samples, std::vector<double>& v_pars);

    private:
        std::vector<int> v_samples;
        //std::vector<int> v_nevents;

        const std::string TAG = color::GREEN_STR + "[ExtraParameters]: " + color::RESET_STR;
};

#endif
