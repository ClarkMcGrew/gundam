#include "AnaSample.hh"
using xsllh::FitBin;

// ctor
AnaSample::AnaSample(int sample_id, const std::string& name, const std::string& detector,
                     const std::string& binning, TTree* t_data)
    : m_sample_id(sample_id)
    , m_name(name)
    , m_detector(detector)
    , m_binning(binning)
    , m_data_tree(t_data)
    , m_norm(1.0)
{
    TH1::SetDefaultSumw2(true);
    SetBinning(m_binning);

    std::cout << TAG << m_name << ", ID " << m_sample_id << std::endl
              << TAG << "Detector: " << m_detector << std::endl
              << TAG << "Bin edges: " << std::endl;

    for(const auto& bin : m_bin_edges)
    {
        std::cout << bin.D2low << " " << bin.D2high << " " << bin.D1low << " " << bin.D1high
                  << std::endl;
    }

    m_hpred    = nullptr;
    m_hmc      = nullptr;
    m_hmc_true = nullptr;
    m_hsig     = nullptr;
    m_hdata    = nullptr;

    MakeHistos(); // with default binning

    std::cout << TAG << "MakeHistos called." << std::endl;
}

AnaSample::~AnaSample()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    if(m_hmc != nullptr)
        delete m_hmc;
    if(m_hmc_true != nullptr)
        delete m_hmc_true;
    if(m_hdata != nullptr)
        delete m_hdata;
    if(m_hsig != nullptr)
        delete m_hsig;
}

void AnaSample::SetBinning(const std::string& binning)
{
    m_binning = binning;
    m_nbins   = 0;

    std::ifstream fin(m_binning, std::ios::in);
    if(!fin.is_open())
    {
        std::cerr << ERR << "In AnaSample::SetBinning().\n"
                  << ERR << "Failed to open binning file: " << m_binning << std::endl;
    }
    else
    {
        std::string line;
        while(std::getline(fin, line))
        {
            std::stringstream ss(line);
            double D1_1, D1_2, D2_1, D2_2;
            if(!(ss >> D2_1 >> D2_2 >> D1_1 >> D1_2))
            {
                std::cerr << TAG << "Bad line format: " << line << std::endl;
                continue;
            }
            m_bin_edges.emplace_back(FitBin(D1_1, D1_2, D2_1, D2_2));
        }
        m_nbins = m_bin_edges.size();
    }
}

// ClearEvents -- clears all events from event vector
void AnaSample::ClearEvents() { m_events.clear(); }

// GetN -- get number of events stored
int AnaSample::GetN() const { return (int)m_events.size(); }

double AnaSample::GetIntegral() const { return m_hpred->Integral(); }

AnaEvent* AnaSample::GetEvent(int evnum)
{
    if(m_events.empty())
    {
        std::cerr << "[ERROR]: In AnaSample::GetEvent()" << std::endl;
        std::cerr << "[ERROR]: No events are found in " << m_name << " sample." << std::endl;
        return nullptr;
    }
    else if(evnum >= m_events.size())
    {
        std::cerr << "[ERROR]: In AnaSample::GetEvent()" << std::endl;
        std::cerr << "[ERROR]: Event number out of bounds in " << m_name << " sample." << std::endl;
        return nullptr;
    }

    return &m_events.at(evnum);
}

void AnaSample::AddEvent(const AnaEvent& event) { m_events.push_back(event); }

void AnaSample::ResetWeights()
{
    for(auto& event : m_events)
        event.SetEvWght(1.0);
}

void AnaSample::PrintStats() const
{
    double mem_kb = sizeof(m_events) * m_events.size() / 1000.0;
    std::cout << TAG << "Sample " << m_name << " ID = " << m_sample_id << std::endl;
    std::cout << TAG << "Num of events = " << m_events.size() << std::endl;
    std::cout << TAG << "Memory used   = " << mem_kb << " kB." << std::endl;
}

void AnaSample::MakeHistos()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    m_hpred = new TH1D(Form("%s_pred_recD1D2", m_name.c_str()),
                       Form("%s_pred_recD1D2", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hpred->SetDirectory(0);

    if(m_hmc != nullptr)
        delete m_hmc;
    m_hmc = new TH1D(Form("%s_mc_recD1D2", m_name.c_str()), Form("%s_mc_recD1D2", m_name.c_str()),
                     m_nbins, 0, m_nbins);
    m_hmc->SetDirectory(0);

    if(m_hmc_true != nullptr)
        delete m_hmc_true;
    m_hmc_true = new TH1D(Form("%s_mc_trueD1D2", m_name.c_str()),
                          Form("%s_mc_trueD1D2", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hmc_true->SetDirectory(0);

    if(m_hsig != nullptr)
        delete m_hsig;
    m_hsig = new TH1D(Form("%s_mc_trueSignal", m_name.c_str()),
                      Form("%s_mc_trueSignal", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hsig->SetDirectory(0);

    std::cout << TAG << m_nbins << " bins inside MakeHistos()." << std::endl;
}

void AnaSample::SetData(TObject* hdata)
{
    if(m_hdata != nullptr)
        delete m_hdata;
    m_hdata = (TH1D*)hdata->Clone(Form("%s_data", m_name.c_str()));
    m_hdata->SetDirectory(0);
}

int AnaSample::GetBinIndex(const double D1, const double D2) const
{
    for(int i = 0; i < m_bin_edges.size(); ++i)
    {
        if(D1 >= m_bin_edges[i].D1low && D1 < m_bin_edges[i].D1high && D2 >= m_bin_edges[i].D2low
           && D2 < m_bin_edges[i].D2high)
        {
            return i;
        }
    }
    return -1;
}

void AnaSample::FillEventHist(int datatype, bool stat_fluc)
{
    if(m_hpred != nullptr)
        m_hpred->Reset();
    if(m_hmc != nullptr)
        m_hmc->Reset();
    if(m_hmc_true != nullptr)
        m_hmc_true->Reset();
    if(m_hsig != nullptr)
        m_hsig->Reset();

    for(std::size_t i = 0; i < m_events.size(); ++i)
    {
        double D1_rec  = m_events[i].GetRecoD1();
        double D2_rec  = m_events[i].GetRecoD2();
        double D1_true = m_events[i].GetTrueD1();
        double D2_true = m_events[i].GetTrueD2();
        double wght    = datatype >= 0 ? m_events[i].GetEvWght() : m_events[i].GetEvWghtMC();

        int anybin_index_rec  = GetBinIndex(D1_rec, D2_rec);
        int anybin_index_true = GetBinIndex(D1_true, D2_true);

        m_hpred->Fill(anybin_index_rec + 0.5, wght);
        m_hmc->Fill(anybin_index_rec + 0.5, wght);
        m_hmc_true->Fill(anybin_index_true + 0.5, wght);

        if(m_events[i].isSignalEvent())
            m_hsig->Fill(anybin_index_true + 0.5, wght);
    }

    m_hpred->Scale(m_norm);
    m_hmc->Scale(m_norm);
    m_hsig->Scale(m_norm);

    if(datatype == 0 || datatype == -1)
        return;

    else if(datatype == 1)
    {
        SetData(m_hpred);
        m_hdata->Reset();

        if(stat_fluc)
            std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

        for(int j = 1; j <= m_hpred->GetNbinsX(); ++j)
        {
            double val = m_hpred->GetBinContent(j);
            if(stat_fluc)
                val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val == 0.0)
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] " << m_name << " bin " << j
                          << " has 0 entries. This may cause a problem with chi2 computations."
                          << std::endl;
                continue;
            }
#endif
            m_hdata->SetBinContent(j, val);
        }
    }

    else if(datatype == 2 || datatype == 3)
    {
        SetData(m_hpred);
        m_hdata->Reset();

        float D1_rec_tree, D2_rec_tree, wght;
        int cut_branch;

        m_data_tree->SetBranchAddress("cut_branch", &cut_branch);
        m_data_tree->SetBranchAddress("weight", &wght);
        m_data_tree->SetBranchAddress("D1Reco", &D1_rec_tree);
        m_data_tree->SetBranchAddress("D2Reco", &D2_rec_tree);

        long int n_entries = m_data_tree->GetEntries();
        for(std::size_t i = 0; i < n_entries; ++i)
        {
            m_data_tree->GetEntry(i);
            if(cut_branch != m_sample_id)
                continue;

            int anybin_index = GetBinIndex(D1_rec_tree, D2_rec_tree);
            if(anybin_index != -1)
            {
                m_hdata->Fill(anybin_index + 0.5, wght);
            }
#ifndef NDEBUG
            else
            {
                std::cout << "[WARNING] In AnaSample::FillEventHist()\n"
                          << "[WARNING] No bin for current data event.\n"
                          << "[WARNING] D1 Reco: " << D1_rec_tree << std::endl
                          << "[WARNING] D2 Reco: " << D2_rec_tree << std::endl;
            }
#endif
        }

        if(stat_fluc && datatype == 2)
        {
            if(stat_fluc)
                std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

            for(unsigned int i = 1; i <= m_hdata->GetNbinsX(); ++i)
            {
                double val = gRandom->Poisson(m_hdata->GetBinContent(i));
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
        std::cout << "[WARNING]: In AnaSample::FillEventHist()\n"
                  << "[WARNING]: Invalid data type to fill histograms!\n";
    }
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

double AnaSample::CalcEffLLH() const
{
    const unsigned int nbins = m_hpred->GetNbinsX();
    double* exp_w  = m_hpred->GetArray();
    double* exp_w2 = m_hpred->GetSumw2()->GetArray();
    double* data   = m_hdata->GetArray();

    //std::cout << m_name << std::endl;

    double llh_eff = 0.0;
    for(unsigned int i = 1; i <= nbins; ++i)
    {
        if(exp_w[i] <= 0.0)
            continue;

        const double b = exp_w[i] / exp_w2[i];
        const double a = (exp_w[i] * b) + 1.0;
        const double k = data[i];

        //std::cout << "--------------" << std::endl;
        //std::cout << "i  : " << i << std::endl;
        //std::cout << "w  : " << exp_w[i] << std::endl;
        //std::cout << "w2 : " << exp_w2[i] << std::endl;
        //std::cout << "a  : " << a << std::endl;
        //std::cout << "b  : " << b << std::endl;
        //std::cout << "k  : " << data[i] << std::endl;

        llh_eff += a * std::log(b) + std::lgamma(k+a) - std::lgamma(k+1) - ((k+a) * std::log1p(b)) - std::lgamma(a);
    }

    return -2 * llh_eff;
}

void AnaSample::Write(TDirectory* dirout, const std::string& bsname, int fititer)
{
    dirout->cd();
    m_hpred->Write(Form("%s_pred", bsname.c_str()));
    m_hmc_true->Write(Form("%s_true", bsname.c_str()));
    if(fititer == 0)
    {
        m_hmc->Write(Form("%s_mc", bsname.c_str()));
        if(m_hdata != nullptr)
            m_hdata->Write(Form("%s_data", bsname.c_str()));
    }
}

void AnaSample::GetSampleBreakdown(TDirectory* dirout, const std::string& tag,
                                   const std::vector<std::string>& topology, bool save)
{
    const int ntopology = topology.size();
    int compos[ntopology];
    std::vector<TH1D> hAnybin_rec;
    std::vector<TH1D> hAnybin_true;

    for(int i = 0; i < ntopology; ++i)
    {
        compos[i] = 0;
        hAnybin_rec.emplace_back(
            TH1D(Form("%s_Anybins_rec_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 Form("%s_Anybins_rec_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 m_nbins, 0, m_nbins));
        hAnybin_rec[i].SetDirectory(0);
        hAnybin_rec[i].GetXaxis()->SetTitle("Bin Index");

        hAnybin_true.emplace_back(
            TH1D(Form("%s_Anybins_true_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 Form("%s_Anybins_true_%s_%s", m_name.c_str(), topology[i].c_str(), tag.c_str()),
                 m_nbins, 0, m_nbins));
        hAnybin_true[i].SetDirectory(0);
        hAnybin_true[i].GetXaxis()->SetTitle("Bin Index");
    }

    int Ntot = GetN();
    for(std::size_t i = 0; i < m_events.size(); ++i)
    {
        double D1_rec    = m_events[i].GetRecoD1();
        double D2_rec    = m_events[i].GetRecoD2();
        double D1_true   = m_events[i].GetTrueD1();
        double D2_true   = m_events[i].GetTrueD2();
        double wght      = m_events[i].GetEvWght();
        int evt_topology = m_events[i].GetTopology();

        compos[evt_topology]++;
        int anybin_index_rec  = GetBinIndex(D1_rec, D2_rec);
        int anybin_index_true = GetBinIndex(D1_true, D2_true);
        hAnybin_rec[evt_topology].Fill(anybin_index_rec + 0.5, wght);
        hAnybin_true[evt_topology].Fill(anybin_index_true + 0.5, wght);
    }

    dirout->cd();
    for(int i = 0; i < ntopology; ++i)
    {
        hAnybin_true[i].Scale(m_norm);
        hAnybin_rec[i].Scale(m_norm);

        if(save == true)
        {
            hAnybin_true[i].Write();
            hAnybin_rec[i].Write();
        }
    }

    std::cout << TAG << "GetSampleBreakdown()\n"
              << "============ Sample " << m_name << " ===========" << std::endl;

    for(int j = 0; j < ntopology; ++j)
    {
        std::cout << std::setw(10) << topology[j] << std::setw(5) << j << std::setw(5) << compos[j]
                  << std::setw(10) << ((1.0 * compos[j]) / Ntot) * 100.0 << "%" << std::endl;
    }

    std::cout << std::setw(10) << "Total" << std::setw(5) << " " << std::setw(5) << Ntot
              << std::setw(10) << "100.00%" << std::endl;
}
