#include "FitParameters.hh"
using xsllh::FitBin;

FitParameters::FitParameters(const std::string& par_name)
{
    m_name = par_name;
    signal_id = 0;
}

FitParameters::~FitParameters()
{;}

void FitParameters::InitEventMap(std::vector<AnaSample*> &sample, int mode)
{
    for(const auto& s : sample)
    {
        if(!std::any_of(v_detectors.begin(), v_detectors.end(),
                        [&](std::string det){ return s->GetDetector() == det; }))
        {
            std::cerr << ERR << "In FitParameters::InitEventMap\n"
                      << ERR << "Detector " << s -> GetDetector() << " not part of fit parameters.\n"
                      << ERR << "Not building event map." << std::endl;
            return;
        }
    }

    InitParameters();
    m_evmap.clear();

    for(std::size_t s=0; s < sample.size(); s++)
    {
        std::vector<int> sample_map;
        for(int i=0; i < sample[s] -> GetN(); i++)
        {
            AnaEvent* ev = sample[s] -> GetEvent(i);

            // SIGNAL DEFINITION TIME
            // Warning, important hard coding up ahead:
            // This is where your signal is actually defined, i.e. what you want to extract an xsec for
            // N.B In Sara's original code THIS WAS THE OTHER WAY AROUND i.e. this if statement asked what was NOT your signal
            // Bare that in mind if you've been using older versions of the fitter.

            if(ev -> isSignalEvent())
            {
                const int sig = ev -> GetSignalType();
                const int bin = m_signal_bm[sig].GetBinIndex(ev->GetTrueVar());
#ifndef NDEBUG
                if(bin == BADBIN)
                {
                    std::cout << WAR << m_name << ", Event: " << i << std::endl
                              << WAR << "Falls outside bin ranges, and will be ignored in the analysis" << std::endl;

                    std::cout << WAR << "Event kinematics: " << std::endl;
                    for(const auto val : ev->GetTrueVar())
                        std::cout << "\t" << val << std::endl;
                }
#endif
                sample_map.push_back(bin);
            }
            else
            {
                sample_map.push_back(PASSEVENT);
                continue;
            }

        }
        m_evmap.push_back(sample_map);
    }
}

void FitParameters::ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double> &params)
{
#ifndef NDEBUG
    if(m_evmap.empty()) //need to build an event map first
    {
        std::cerr << ERR << "In FitParameters::ReWeight()\n"
                  << ERR << "Need to build event map index for " << m_name << std::endl;
        return;
    }
#endif

    const int bin = m_evmap[nsample][nevent];

    //skip event if not Signal
    if(bin == PASSEVENT || bin == BADBIN)
        return;

    else
    {
#ifndef NDEBUG
        if(bin > params.size())
        {
            std::cout << WAR << "In FitParameters::ReWeight()\n"
                      << WAR << "Number of bins in " << m_name << " does not match num of parameters.\n"
                      << WAR << "Setting event weight to zero." << std::endl;
            event -> AddEvWght(0.0);
        }
#endif

        event -> AddEvWght(params[bin + m_sig_offset.at(event->GetSignalType())]);

        /*
        std::cout << "-----------------" << std::endl;
        std::cout << "Ev D1: " << event -> GetTrueD1() << std::endl
                  << "Ev D2: " << event -> GetTrueD2() << std::endl;
        std::cout << "Ev ST: " << event -> GetSignalType() << std::endl;
        std::cout << "Ev TP: " << event -> GetTopology() << std::endl;
        std::cout << "Ev TG: " << event -> GetTarget() << std::endl;
        std::cout << "Bin  : " << bin << std::endl;
        std::cout << "Off  : " << m_sig_offset.at(event->GetSignalType()) << std::endl;
        */
    }
}

void FitParameters::InitParameters()
{
    unsigned int offset = 0;
    for(const auto& sig : v_signals)
    {
        m_sig_offset.emplace(std::make_pair(sig, offset));
        const int nbins = m_signal_bm.at(sig).GetNbins();
        for(int i = 0; i < nbins; ++i)
        {
            pars_name.push_back(Form("%s_sig%d_%d", m_name.c_str(), sig, i));
            pars_prior.push_back(1.0); //all weights are 1.0 a priori
            pars_step.push_back(0.05);
            pars_limlow.push_back(0.0);
            pars_limhigh.push_back(10.0);
            pars_fixed.push_back(false);
        }

        std::cout << TAG << "Total " << nbins << " parameters at "
                  << offset << " for signal ID " << sig << std::endl;
        offset += nbins;
    }

    Npar = pars_name.size();
    pars_original = pars_prior;
}

void FitParameters::AddDetector(const std::string& det, const std::vector<SignalDef>& v_input)
{
    std::cout << TAG << "Adding detector " << det << " for " << m_name << std::endl;

    for(const auto& sig : v_input)
    {
        if(sig.detector != det || sig.use_signal == false)
            continue;

        v_signals.emplace_back(signal_id);

        std::cout << TAG << "Adding signal " << sig.name << " with ID " << signal_id
                  << " to fit." << std::endl;

        m_signal_bm.emplace(std::make_pair(signal_id, BinManager(sig.binning)));
        v_detectors.emplace_back(det);
        signal_id++;
    }
}

double FitParameters::CalcRegularisation(const std::vector<double>& params) const
{
    return CalcRegularisation(params, m_regstrength, m_regmethod);
}

double FitParameters::CalcRegularisation(const std::vector<double>& params, double strength,
                                         RegMethod flag) const
{
    double L_reg = 0;
    unsigned int offset = 0;
    for(const auto& signal_id : v_signals)
    {
        const unsigned int nbins = m_signal_bm.at(signal_id).GetNbins();

        if(flag == kL1Reg)
        {
            for(int i = offset; i < offset+nbins-1; ++i)
                L_reg += std::fabs(params[i] - params[i+1]);
        }

        else if(flag == kL2Reg)
        {
            for(int i = offset; i < offset+nbins-1; ++i)
                L_reg += (params[i] - params[i+1]) * (params[i] - params[i+1]);
        }

        else
        {
            std::cout << WAR << "In CalcRegularisation(): "
                      << "Invalid regularisation method! Returning 0.\n";
            break;
        }

        offset += nbins;
    }

    return strength * L_reg;
}
