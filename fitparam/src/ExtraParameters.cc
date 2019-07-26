#include "ExtraParameters.hh"

ExtraParameters::ExtraParameters(const std::string& name)
{
    m_name = name;
}

void ExtraParameters::InitEventMap(std::vector<AnaSample*>& sample, int mode)
{
    InitFixed();
    //InitParameters();
    m_evmap.clear();

    if(mode == 2)
        std::cout << TAG << "Not using reweighting for " << m_name << "." << std::endl;

    for(std::size_t s = 0; s < sample.size(); ++s)
    {
        std::vector<int> sample_map;
        for(int i = 0; i < sample[s]->GetN(); ++i)
        {
            AnaEvent* ev = sample[s]->GetEvent(i);

            int bin = s;

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
    std::vector<double> params_constrained = CalcConstraint(params);

    //std::cout << "---------------" << std::endl;
    //for(const auto& v : params_constrained)
    //    std::cout << "Par Extra: " << v << std::endl;

    if(bin == PASSEVENT)
        return;
    if(bin == BADBIN)
        event->AddEvWght(0.0);
    else
    {
        /*
        if(bin > params.size())
        {
            std::cout << WAR << "In ExtraParameters::ReWeight()\n"
                      << WAR << "Number of bins in " << m_name << " does not match num of parameters.\n"
                      << WAR << "Setting event weight to zero." << std::endl;
            event->AddEvWght(0.0);
        }
        */

        //event->AddEvWght(params[bin]);
        event->AddEvWght(params_constrained[bin]);
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

void ExtraParameters::InitFixed()
{
    const int num_samples = 8;
    for(int i = 0; i < num_samples; ++i)
    {
        pars_name.push_back(Form("%s_%d", m_name.c_str(), i));
        pars_prior.push_back(1.0);
        pars_step.push_back(0.05);
        pars_limlow.push_back(0.0);
        pars_limhigh.push_back(4.0);
        pars_fixed.push_back(false);

        std::cout << TAG << "Adding extra parameter " << i << std::endl;
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
        v_nevents.emplace_back(sample->GetN());

        std::cout << TAG << "Adding sample " << sample->GetName()
                  << " with ID " << sample_id << " to " << m_name << std::endl;
    }
}

std::vector<double> ExtraParameters::CalcConstraint(std::vector<double>& v_pars)
{
    //I hate this function so much.
    double Ntotal_nd280 = v_nevents[0]
                        + v_nevents[1]
                        + v_nevents[2]
                        + v_nevents[3]
                        + v_nevents[4]
                        + v_nevents[5]
                        + v_nevents[6]
                        + v_nevents[7];

    double constraint_nd280 = Ntotal_nd280
                            - v_pars[0] * v_nevents[0]
                            - v_pars[1] * v_nevents[1]
                            - v_pars[2] * v_nevents[2]
                            - v_pars[3] * v_nevents[3]
                            - v_pars[4] * v_nevents[4]
                            - v_pars[5] * v_nevents[5]
                            - v_pars[6] * v_nevents[6];
                            /* - v_pars[7] * v_nevents[7]; */
    constraint_nd280 = constraint_nd280 / v_nevents[7];

    double Ntotal_ingrid = v_nevents[8] + v_nevents[9];
    double constraint_ingrid = Ntotal_ingrid
                             - v_pars[7] * v_nevents[8];
                             /* - v_pars[9] * v_nevents[9]; */
    constraint_ingrid = constraint_ingrid / v_nevents[9];

    std::vector<double> new_pars = {v_pars[0],
                                    v_pars[1],
                                    v_pars[2],
                                    v_pars[3],
                                    v_pars[4],
                                    v_pars[5],
                                    v_pars[6],
                                    constraint_nd280,
                                    v_pars[7],
                                    constraint_ingrid};

    return new_pars;
}
