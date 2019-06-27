#include "ExtraParameters.hh"

ExtraParameters::ExtraParameters(const std::string& name)
{
    m_name = name;
}

void ExtraParameters::InitEventMap(std::vector<AnaSample*>& sample, int mode)
{
    InitParameters();
    m_evmap.clear();

    if(mode == 2)
        std::cout << TAG << "Not using reweighting for " << m_name << "." << std::endl;

    for(std::size_t s = 0; s < sample.size(); ++s)
    {
        std::vector<int> sample_map;
        for(int i = 0; i < sample[s]->GetN(); ++i)
        {
            AnaEvent* ev = sample[s]->GetEvent(i);

            int bin = 0;

            if(mode == 1 && ev->isSignalEvent())
                bin = PASSEVENT;
            else if(mode == 2)
                bin = PASSEVENT;
            sample_map.push_back(bin);
        }
        m_evmap.push_back(sample_map);
    }
}

void ExtraParameters::ReWeight(AnaEvent* event, const std::string& det, int nsample, int nevent, std::vector<double>& params)
{
    if(m_evmap.empty()) // need to build an event map first
    {
        std::cerr << ERR << "In ExtraParameters::ReWeight()\n"
                  << ERR << "Need to build event map index for " << m_name << std::endl;
        return;
    }

    const int bin = m_evmap[nsample][nevent];

    if(bin == PASSEVENT)
        return;
    if(bin == BADBIN)
        event->AddEvWght(0.0);
    else
    {
        if(bin > params.size())
        {
            std::cout << WAR << "In ExtraParameters::ReWeight()\n"
                      << WAR << "Number of bins in " << m_name << " does not match num of parameters.\n"
                      << WAR << "Setting event weight to zero." << std::endl;
            event->AddEvWght(0.0);
        }

        event->AddEvWght(params[bin + nsample]);
    }
}

void ExtraParameters::InitParameters()
{
    unsigned int offset = 0;
    for(const auto& sam : v_samples)
    {
        pars_name.push_back(Form("%s_sam%d_%d", m_name.c_str(), sam, offset++));
        pars_prior.push_back(1.0);
        pars_step.push_back(0.05);
        pars_limlow.push_back(0.0);
        pars_limhigh.push_back(4.0);
        pars_fixed.push_back(false);

        std::cout << TAG << "Added sample ID " << sam << " to " << m_name << std::endl;
    }

    Npar = pars_name.size();
    pars_original = pars_prior;
}

void ExtraParameters::AddSample(std::vector<AnaSample*>& v_sample)
{
    for(const auto& sample : v_sample)
    {
        const int sample_id = sample->GetSampleID();
        v_samples.emplace_back(sample_id);

        std::cout << TAG << "Adding sample " << sample->GetName()
                  << " with ID " << sample_id << " to " << m_name << std::endl;
    }
}
