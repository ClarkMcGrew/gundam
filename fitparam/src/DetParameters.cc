#include "DetParameters.hh"

DetParameters::DetParameters(const std::string& name)
{
    m_name = name;
}

DetParameters::~DetParameters() { ; }

void DetParameters::InitEventMap(std::vector<AnaSample*>& sample, int mode)
{
    InitParameters();
    m_evmap.clear();

    if(mode == 2)
        std::cout << TAG << "Not using detector reweighting." << std::endl;

    for(std::size_t s = 0; s < sample.size(); ++s)
    {
        std::vector<int> sample_map;
        for(int i = 0; i < sample[s]->GetN(); ++i)
        {
            AnaEvent* ev = sample[s]->GetEvent(i);
            int bin = m_sample_bm[s].GetBinIndex(ev->GetRecoVar());
#ifndef NDEBUG
            if(bin < 0)
            {
                std::cout << WAR << m_name << ", Event: " << i << ", Sample: " << s <<  std::endl
                          << WAR << "Falls outside bin ranges, and will be ignored in the analysis" << std::endl;

                std::cout << WAR << "Event kinematics: " << std::endl;
                for(const auto val : ev->GetRecoVar())
                    std::cout << "\t" << val << std::endl;
            }
#endif
            // If event is signal let the c_i params handle the reweighting:
            if(mode == 1 && ev->isSignalEvent())
                bin = PASSEVENT;
            else if(mode == 2)
                bin = PASSEVENT;
            sample_map.push_back(bin);
        }
        m_evmap.push_back(sample_map);
    }
}

void DetParameters::ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params)
{
#ifndef NDEBUG
    if(m_evmap.empty()) // need to build an event map first
    {
        std::cerr << ERR << "In DetParameters::ReWeight()\n"
                  << ERR << "Need to build event map index for " << m_name << std::endl;
        return;
    }
#endif

    const int bin = m_evmap[nsample][nevent];

    if(bin == PASSEVENT || bin == BADBIN)
        return;
    else
    {
#ifndef NDEBUG
        if(bin > params.size())
        {
            std::cout << WAR << "In DetParameters::ReWeight()\n"
                      << WAR << "Number of bins in " << m_name << " does not match num of parameters.\n"
                      << WAR << "Setting event weight to zero." << std::endl;
            event->AddEvWght(0.0);
        }
#endif

        event->AddEvWght(params[bin + m_sample_offset.at(event->GetSampleType())]);
    }
}

void DetParameters::InitParameters()
{
    unsigned int offset = 0;
    for(const auto& sam : v_samples)
    {
        m_sample_offset.emplace(std::make_pair(sam, offset));
        const int nbins = m_sample_bm.at(sam).GetNbins();
        for(int i = 0; i < nbins; ++i)
        {
            pars_name.push_back(Form("%s_sam%d_%d", m_name.c_str(), sam, i));
            pars_prior.push_back(1.0);
            pars_step.push_back(0.05);
            pars_limlow.push_back(0.0);
            pars_limhigh.push_back(2.0);
            pars_fixed.push_back(false);
        }

        std::cout << TAG << "Total " << nbins << " parameters at "
                  << offset << " for sample ID " << sam << std::endl;
        offset += nbins;
    }

    Npar = pars_name.size();
    pars_original = pars_prior;

    if(m_decompose)
    {
        pars_prior = eigen_decomp -> GetDecompParameters(pars_prior);
        pars_limlow = std::vector<double>(Npar, -100);
        pars_limhigh = std::vector<double>(Npar, 100);

        const int idx = eigen_decomp -> GetInfoFraction(m_info_frac);
        for(int i = idx; i < Npar; ++i)
            pars_fixed[i] = true;

        std::cout << TAG << "Decomposed parameters.\n"
                  << TAG << "Keeping the " << idx << " largest eigen values.\n"
                  << TAG << "Corresponds to " << m_info_frac * 100.0
                  << "\% total variance.\n";
    }
}

void DetParameters::AddDetector(const std::string& det, std::vector<AnaSample*>& v_sample)
{
    std::cout << TAG << "Adding detector " << det << " for " << m_name << std::endl;

    for(const auto& sample : v_sample)
    {
        if(sample->GetDetector() != det)
            continue;

        const int sample_id = sample->GetSampleID();
        v_samples.emplace_back(sample_id);

        std::cout << TAG << "Adding sample " << sample->GetName()
                  << " with ID " << sample_id << " to fit." << std::endl;

        m_sample_bm.emplace(std::make_pair(sample_id, BinManager(sample->GetBinning())));
    }
}
