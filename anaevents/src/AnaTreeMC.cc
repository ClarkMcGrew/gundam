// This is the code that actually reads in the MC tree and fills the event info.
// The tree should be produced by feeding a HL2 microtree into the treeconvert macro.

#include "AnaTreeMC.hh"

AnaTreeMC::AnaTreeMC(const std::string& file_name, const std::string& tree_name, bool extra_var)
    : read_extra_var(extra_var)
{
    fChain = new TChain(tree_name.c_str());
    fChain->Add(file_name.c_str());
    SetBranches();
}

AnaTreeMC::~AnaTreeMC()
{
    if(fChain != nullptr)
        delete fChain->GetCurrentFile();
}

long int AnaTreeMC::GetEntry(long int entry) const
{
    // Read contents of entry.
    if(fChain == nullptr)
        return -1;
    else
        return fChain->GetEntry(entry);
}

void AnaTreeMC::SetBranches()
{
    // Set branch addresses and branch pointers
    fChain->SetBranchAddress("nutype", &nutype);
    fChain->SetBranchAddress("beammode", &beammode);
    fChain->SetBranchAddress("cut_branch", &sample);
    fChain->SetBranchAddress("topology", &topology);
    fChain->SetBranchAddress("reaction", &reaction);
    fChain->SetBranchAddress("target", &target);
    fChain->SetBranchAddress("D1True", &D1True);
    fChain->SetBranchAddress("D1Reco", &D1Reco);
    fChain->SetBranchAddress("D2True", &D2True);
    fChain->SetBranchAddress("D2Reco", &D2Reco);
    fChain->SetBranchAddress("D3True", &D3True);
    fChain->SetBranchAddress("D3Reco", &D3Reco);
    fChain->SetBranchAddress("D4True", &D4True);
    fChain->SetBranchAddress("D4Reco", &D4Reco);
    fChain->SetBranchAddress("enu_true", &enu_true);
    fChain->SetBranchAddress("weight", &weight);

    if(read_extra_var)
    {
        //Put extra variables here.
    }
}

void AnaTreeMC::GetEvents(std::vector<AnaSample*>& ana_samples,
                          const std::vector<SignalDef>& v_signal, const bool evt_type)
{
    if(fChain == nullptr || ana_samples.empty())
        return;

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Reading Events "));

    long int nentries = fChain->GetEntries();
    long int nbytes   = 0;

    std::cout << TAG << "Reading events...\n";

    // Loop over all events:
    for(long int jentry = 0; jentry < nentries; jentry++)
    {
        nbytes += fChain->GetEntry(jentry);
        AnaEvent ev(jentry);
        ev.SetTrueEvent(evt_type);
        ev.SetFlavor(nutype);
        ev.SetBeamMode(beammode);
        ev.SetSampleType(sample);
        ev.SetTopology(topology); // mectopology (i.e. CC0Pi,CC1Pi etc)
        ev.SetReaction(reaction); // reaction (i.e. CCQE,CCRES etc)
        ev.SetTarget(target);
        ev.SetTrueEnu(enu_true);
        ev.SetTrueD1(D1True);
        ev.SetRecoD1(D1Reco);
        ev.SetTrueD2(D2True);
        ev.SetRecoD2(D2Reco);
        ev.SetTrueD3(D3True);
        ev.SetRecoD3(D3Reco);
        ev.SetTrueD4(D4True);
        ev.SetRecoD4(D4Reco);
        ev.SetEvWght(weight);
        ev.SetEvWghtMC(weight);

        if(read_extra_var)
        {
            //Put extra variables here.
        }

        // Check for all events if they are signal events:

        // Number of current signal type:
        int signal_type = 0;
        // Loop over all sets of temlate parameters as defined in the .json config file in the "template_par" entry for this detector:
        for(const auto& sd : v_signal)
        {
            bool sig_passed = true;

            // Loop over all the different signal definitions for this set of template parameters (e.g. topology, target, nutype, etc.):
            for(const auto& kv : sd.definition)
            {
                bool var_passed = false;
                
                // Loop over all the values for the current signal definition (e.g. all the different topology integers):
                for(const auto& val : kv.second)
                {
                    if(ev.GetEventVar(kv.first) == val)
                        var_passed = true;
                }
                sig_passed = sig_passed && var_passed;
            }
            if(sig_passed)
            {
                ev.SetSignalType(signal_type);
                ev.SetSignalEvent();
                break;
            }
            signal_type++;
        }

        // Loop over all samples:
        for(auto& s : ana_samples)
        {
            // Check if current event falls into current sample:
            if(s->GetSampleID() == sample)
            {
                //If yes, add event to current AnaSample:
                s->AddEvent(ev);
            }
        }

        // Print current progress:
        if(jentry % 2000 == 0 || jentry == (nentries - 1))
            pbar.Print(jentry, nentries - 1);
    }

    // Print out number of events for each sample and how much memory it uses:
    for(auto& sample : ana_samples)
        sample->PrintStats();
}
