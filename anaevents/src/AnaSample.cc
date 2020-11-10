#include "AnaSample.hh"

AnaSample::AnaSample(int sample_id, const std::string& name, const std::string& detector,
                     const std::string& binning, TTree* t_data)
    : m_sample_id(sample_id)
    , m_name(name)
    , m_detector(detector)
    , m_binning(binning)
    , m_data_tree(t_data)
    , m_hpred(nullptr)
    , m_hdata(nullptr)
    , m_norm(1.0)
{
    TH1::SetDefaultSumw2(true);

    std::cout << TAG << "Sample: " << m_name << " (ID: " << m_sample_id << ")" << std::endl
              << TAG << "Detector: " << m_detector << std::endl;

    bm.SetBinning(m_binning);
    bm.Print();
    m_nbins = bm.GetNbins();

    m_llh = new PoissonLLH;

    MakeHistos(); // with default binning

    std::cout << TAG << "MakeHistos called." << std::endl;
}

AnaSample::~AnaSample()
{
    if(m_hpred != nullptr)
        delete m_hpred;

    if(m_hdata != nullptr)
        delete m_hdata;
}

AnaEvent* AnaSample::GetEvent(const unsigned int evnum)
{
#ifndef NDEBUG
    if(m_events.empty())
    {
        std::cerr << ERR << " In AnaSample::GetEvent()" << std::endl;
        std::cerr << ERR << " No events are found in " << m_name << " sample." << std::endl;
        return nullptr;
    }
    else if(evnum >= m_events.size())
    {
        std::cerr << ERR << " In AnaSample::GetEvent()" << std::endl;
        std::cerr << ERR << " Event number out of bounds in " << m_name << " sample." << std::endl;
        return nullptr;
    }
#endif

    return &m_events[evnum];
}

void AnaSample::ResetWeights()
{
    for(auto& event : m_events)
        event.SetEvWght(event.GetEvWghtMC());
}

void AnaSample::PrintStats() const
{
    const double mem_kb = sizeof(m_events) * m_events.size() / 1000.0;
    std::cout << TAG << "Sample " << m_name << " ID = " << m_sample_id << std::endl;
    std::cout << TAG << "Num of events = " << m_events.size() << std::endl;
    std::cout << TAG << "Memory used   = " << mem_kb << " kB." << std::endl;
}

void AnaSample::MakeHistos()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    m_hpred = new TH1D(Form("%s_pred_reco", m_name.c_str()),
                       Form("%s_pred_reco", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hpred->SetDirectory(0);

    if(m_hdata != nullptr)
        delete m_hdata;
    m_hdata = new TH1D(Form("%s_data", m_name.c_str()), Form("%s_data", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hdata->SetDirectory(0);
}

void AnaSample::InitEventMap()
{
    for(auto& e : m_events)
    {
        const int b = bm.GetBinIndex(e.GetRecoVar());
#ifndef NDEBUG
        if(b < 0)
        {
            std::cout << WAR << "In AnaSample::InitEventMap()\n"
                      << WAR << "No bin for current event." << std::endl;
            std::cout << WAR << "Event kinematics: " << std::endl;
            for(const auto val : e.GetRecoVar())
                std::cout << "\t" << val << std::endl;
        }
#endif
        e.SetSampleBin(b);
    }

    std::sort(m_events.begin(), m_events.end(), [](AnaEvent& a, AnaEvent& b){ return a.GetSampleBin() < b.GetSampleBin(); });
}

void AnaSample::FillEventHist(bool reset_weights)
{
#ifndef NDEBUG
    if(m_hpred == nullptr)
    {
        std::cerr << ERR << "In AnaSample::FillEventHist() h_pred is a nullptr!"
                         << "Returning from function." << std::endl;
        return;
    }
#endif
    m_hpred->Reset();

    for(const auto& e : m_events)
    {
        const double weight = reset_weights ? e.GetEvWghtMC() : e.GetEvWght();
        const int reco_bin  = e.GetSampleBin();
        m_hpred->Fill(reco_bin + 0.5, weight);
    }

    m_hpred->Scale(m_norm);
    return;
}

void AnaSample::FillDataHist(int datatype, bool stat_fluc)
{
    if(datatype == kAsimov)
    {
        m_hdata->Reset();
        if(stat_fluc)
            std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

        for(int j = 1; j <= m_hpred->GetNbinsX(); ++j)
        {
            double val = m_hpred->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val <= 0.0)
            {
                std::cout << WAR << "In AnaSample::FillEventHist()\n"
                          << WAR << "In Sample " <<  m_name << ", bin " << j
                          << " has 0 (or negative) entries. This may cause a problem with chi2 computations."
                          << std::endl;
            }
#endif
            m_hdata->SetBinContent(j, val);
        }
    }
    else if(datatype == kExternal || datatype == kData)
    {
        m_hdata->Reset();

        int sample = -99;
        float weight = 1.0;
        std::vector<double>* reco_var = 0;

        m_data_tree->SetBranchAddress("cut_branch", &sample);
        m_data_tree->SetBranchAddress("weight", &weight);
        m_data_tree->SetBranchAddress("reco_var", &reco_var);

        long int n_entries = m_data_tree->GetEntries();
        for(std::size_t i = 0; i < n_entries; ++i)
        {
            m_data_tree->GetEntry(i);
            if(sample != m_sample_id)
                continue;

            const int bin = bm.GetBinIndex(*reco_var);
            if(bin > -1)
            {
                m_hdata->Fill(bin + 0.5, weight);
            }
#ifndef NDEBUG
            else
            {
                std::cout << WAR << "In AnaSample::FillEventHist()\n"
                          << WAR << "No bin for current data event." << std::endl;
                std::cout << WAR << "Event kinematics: " << std::endl;
                for(const auto val : *reco_var)
                    std::cout << "\t" << val << std::endl;
            }
#endif
        }

        if(stat_fluc && datatype == kExternal)
        {
            std::cout << TAG << "Applying statistical fluctuations..." << std::endl;
            for(unsigned int i = 1; i <= m_hdata->GetNbinsX(); ++i)
            {
                const double val = gRandom->Poisson(m_hdata->GetBinContent(i));
                m_hdata->SetBinContent(i, val);
            }
        }
#ifndef NDEBUG
        std::cout << TAG << "Data histogram filled: " << std::endl;
        m_hdata->Print();
#endif
    }
    else
    {
        std::cout << WAR << "In AnaSample::FillDataHist()\n"
                  << WAR << "Invalid data type to fill histograms!\n";
    }
}

void AnaSample::SetLLHFunction(const std::string& func_name)
{
    if(m_llh != nullptr)
        delete m_llh;

    if(func_name.empty())
    {
        std::cout << TAG << "Likelihood function name empty. Setting to Poisson by default." << std::endl;
        m_llh = new PoissonLLH;
    }
    else if(func_name == "Poisson")
    {
        std::cout << TAG << "Setting likelihood function to Poisson." << std::endl;
        m_llh = new PoissonLLH;
    }
    else if(func_name == "Effective")
    {
        std::cout << TAG << "Setting likelihood function to Tianlu's effective likelihood." << std::endl;
        m_llh = new EffLLH;
    }
    else if(func_name == "Barlow")
    {
        std::cout << TAG << "Setting likelihood function to Barlow-Beeston." << std::endl;
        m_llh = new BarlowLLH;
    }
}

double AnaSample::CalcLLH() const
{
    const unsigned int nbins = m_hpred->GetNbinsX();
    double* exp_w  = m_hpred->GetArray();
    double* exp_w2 = m_hpred->GetSumw2()->GetArray();
    double* data   = m_hdata->GetArray();

    double chi2 = 0.0;
    for(unsigned int i = 1; i <= nbins; ++i)
        chi2 += (*m_llh)(exp_w[i], exp_w2[i], data[i]);

    return chi2;
}

double AnaSample::CalcChi2() const
{
    if(m_hdata == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                  << "[ERROR]: Need to define data histogram." << std::endl;
        return 0.0;
    }

    int nbins = m_hpred->GetNbinsX();
    if(nbins != m_hdata->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                  << "[ERROR]: Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    double chi2 = 0.0;
    for(int j = 1; j <= nbins; ++j)
    {
        double obs = m_hdata->GetBinContent(j);
        double exp = m_hpred->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }

    if(chi2 != chi2)
    {
        std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                  << "[WARNING]: Stat chi2 is nan, setting to 0." << std::endl;
        chi2 = 0.0;
    }

    return chi2;
}

void AnaSample::WriteEventHist(TDirectory* dirout, const std::string& bsname)
{
    dirout->cd();
    if(m_hpred != nullptr)
        m_hpred->Write(Form("evhist_sam%d_pred%s", m_sample_id, bsname.c_str()));
}

void AnaSample::WriteDataHist(TDirectory* dirout, const std::string& bsname)
{
    dirout->cd();
    if(m_hdata != nullptr)
        m_hdata->Write(Form("evhist_sam%d_data%s", m_sample_id, bsname.c_str()));
}
