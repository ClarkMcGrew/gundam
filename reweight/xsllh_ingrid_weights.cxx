#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "TAxis.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"

#include "T2KReWeight.h"
#include "T2KSyst.h"

#include "T2KGenieReWeight.h"
#include "T2KGenieUtils.h"

#include "T2KNeutReWeight.h"
#include "T2KNeutUtils.h"

#include "T2KNIWGReWeight.h"
#include "T2KNIWGUtils.h"

#include "T2KJNuBeamReWeight.h"

#include "ND__GRooTrackerVtx.h"
#include "ND__NRooTrackerVtx.h"

#include "T2KWeightsStorer.h"

#include "SK__h1.h"

#include "neutrootTreeSingleton.h"

const std::string RESET("\033[0m");
const std::string RED("\033[31;1m");
const std::string GREEN("\033[92m");
const std::string COMMENT_CHAR("#");

const std::string TAG = GREEN + "[xsReWeight]: " + RESET;
const std::string ERR = RED + "[ERROR]: " + RESET;

int GetDialIndex(const std::string& dial_name, int step);
int GetIngridNutype(bool anti, bool nue);
int GetIngridReaction(int code);
int GetIngridTopology(int code);

int main(int argc, char** argv)
{
    std::cout << "--------------------------------------------------------\n"
              << TAG << "Welcome to the Super-xsLLhReWeight Interface.\n"
              << TAG << "Initializing the reweight machinery..." << std::endl;

    bool is_SK_tree = false;
    bool limits = false;
    unsigned int dial_steps = 7;
    double nominal = 0;
    double error_neg = 0;
    double error_pos = 0;
    std::string fname_input;
    std::string fname_output;
    std::string fname_dials;
    std::string dial_name;
    std::string dial_name_enum;

    char option;
    while((option = getopt(argc, argv, "i:o:r:d:n:SLh")) != -1)
    {
        switch(option)
        {
            case 'i':
                fname_input = optarg;
                break;
            case 'o':
                fname_output = optarg;
                break;
            case 'r':
                dial_name = optarg;
                break;
            case 'd':
                fname_dials = optarg;
                break;
            case 'n':
                break;
            case 'S':
                is_SK_tree = true;
                break;
            case 'L':
                limits = true;
                break;
            case 'h':
                std::cout << "USAGE: " << argv[0] << "\nOPTIONS\n"
                          << "-i : Input list of Highland trees (.txt)\n"
                          << "-o : Output ROOT filename\n"
                          << "-r : Systematic parameter to reweight\n"
                          << "-d : Dial value file (.txt)\n"
                          << "-S : Read as SK file\n"
                          << "-L : Read errors as upper and lower limits\n"
                          << "-h : Display this help message\n";
            default:
                return 0;
        }
    }

    if(fname_input.empty() || fname_output.empty() || dial_name.empty())
    {
        std::cout << ERR << "Missing necessary command line arguments.\n" << std::endl;
        return 1;
    }

    std::cout << TAG << "Input INGRID file: " << fname_input << std::endl
              << TAG << "Output weights file: " << fname_output << std::endl
              << TAG << "Input dial value file: " << fname_dials << std::endl
              << TAG << "Generating weights for " << dial_name << std::endl
              << TAG << "Number of dial steps: " << dial_steps << std::endl;

    std::ifstream fdial(fname_dials, std::ios::in);
    if(!fdial.is_open())
    {
        std::cerr << ERR << "Failed to open " << fname_dials << std::endl;
        return 1;
    }
    else
    {
        std::string line;
        while(std::getline(fdial, line))
        {
            std::stringstream ss(line);
            std::string dial;
            std::string dial_enum;
            double nom{0}, err_n{0}, err_p{0};

            ss >> dial >> nom >> err_n >> err_p >> dial_enum;
            if(ss.str().front() == COMMENT_CHAR)
                continue;

            if(dial == dial_name)
            {
                nominal = nom;
                error_neg = err_n;
                error_pos = err_p;
                dial_name_enum = dial_enum;
                break;
            }
        }
        fdial.close();
    }

    std::cout << "\n" << TAG << "Dial Name: " << dial_name
              << "\n" << TAG << "Dial enum: " << dial_name_enum
              << "\n" << TAG << "Nominal: " << nominal
              << "\n" << TAG << (limits ? "Limit -: " : "Error -: ") << error_neg
              << "\n" << TAG << (limits ? "Limit +: " : "Error +: ") << error_pos << std::endl;

    int rejected_weights = 0;
    int total_weights = 0;
    const float mu_mass = 105.658374;
    const float deg_to_rad = TMath::Pi() / 180;
    float enu_true, enu_reco;
    float pmu_true, pmu_reco;
    float cosmu_true, cosmu_reco;
    float weight_nom;
    float weight_syst[dial_steps];
    int topology, reaction, sample, target;
    std::string weight_str = "weight_syst[" + std::to_string(dial_steps) + "]/F";

    std::vector<TH1D> v_hists;
    std::vector<double> v_dial_val;
    for(unsigned int d = 0; d < dial_steps; ++d)
    {
        std::stringstream ss;
        ss << dial_name << "_step" << d;
        v_hists.emplace_back(TH1D(ss.str().c_str(), ss.str().c_str(), 25, 0, 5000));

        int step = -1.0 * std::floor(dial_steps / 2.0) + d;
        if(limits)
        {
            double inc = (error_pos - error_neg) / (dial_steps - 1);
            v_dial_val.emplace_back(step * inc / nominal);
        }
        else
        {
            if(step < 0)
                v_dial_val.emplace_back(step * error_neg / nominal);
            else
                v_dial_val.emplace_back(step * error_pos / nominal);
        }
    }

    TFile* file_output = TFile::Open(fname_output.c_str(), "RECREATE");
    TTree* tree_output = new TTree("selectedEvents", "selectedEvents");

    tree_output->Branch("Enutrue", &enu_true, "enu_true/F");
    tree_output->Branch("Enureco", &enu_reco, "enu_reco/F");
    tree_output->Branch("sample", &sample, "sample/I");
    tree_output->Branch("target", &target, "target/I");
    tree_output->Branch("reaction", &reaction, "reaction/I");
    tree_output->Branch("topology", &topology, "topology/I");
    tree_output->Branch("Pmutrue", &pmu_true, "pmu_true/F");
    tree_output->Branch("Pmureco", &pmu_reco, "pmu_reco/F");
    tree_output->Branch("CTHmutrue", &cosmu_true, "cosmu_true/F");
    tree_output->Branch("CTHmureco", &cosmu_reco, "cosmu_reco/F");
    tree_output->Branch("weight", &weight_nom, "weight_nom/F");
    tree_output->Branch("weight_syst", weight_syst, weight_str.c_str());

    std::cout << TAG << "Initializing T2KReWeight and reweight engines." << std::endl;
    t2krew::T2KReWeight rw;
    t2krew::T2KSyst_t t2ksyst = t2krew::T2KSyst::FromString(dial_name_enum);
    rw.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());

    rw.Systematics().Include(t2ksyst);
    rw.Systematics().SetAbsTwk(t2ksyst);
    rw.Systematics().PrintSummary();

    TFile* file_input = TFile::Open(fname_input.c_str(), "READ");
    TTree* tree_input = (TTree*)file_input -> Get("wtree");
    std::cout << TAG << "Reading file: " << fname_input << std::endl;

    //Set up INGRID variables
    int interaction_type, fsi_int;
    int muon_track, selected_sample;
    int track_sample;
    float enu, event_weight;
    float mom_true, angle_true;
    float mom_reco, angle_reco;
    float reweight_array[147];
    bool is_anti, is_nue;
    bool new_event;

    tree_input -> SetBranchAddress("InteractionType", &interaction_type);
    tree_input -> SetBranchAddress("FSIInt", &fsi_int);
    tree_input -> SetBranchAddress("SelectedSample", &selected_sample);
    tree_input -> SetBranchAddress("MuonCandidateTrack", &muon_track);
    tree_input -> SetBranchAddress("Enu", &enu);
    tree_input -> SetBranchAddress("IsAnti", &is_anti);
    tree_input -> SetBranchAddress("IsNuE", &is_nue);
    tree_input -> SetBranchAddress("NewEvent", &new_event);
    tree_input -> SetBranchAddress("weight", &event_weight);
    tree_input -> SetBranchAddress("TrueMomentumMuon", &mom_true);
    tree_input -> SetBranchAddress("TrueAngleMuon", &angle_true);
    tree_input -> SetBranchAddress("ReWeight[147]", reweight_array);

    const float iron_carbon_ratio = 7.640777;
    int sample_track0, sample_track1;
    float iron_track0, plastic_track0, angle_track0;
    float iron_track1, plastic_track1, angle_track1;

    tree_input -> SetBranchAddress("Sample_track0", &sample_track0);
    tree_input -> SetBranchAddress("TrackAngle_track0", &angle_track0);
    tree_input -> SetBranchAddress("IronDistance_track0", &iron_track0);
    tree_input -> SetBranchAddress("PlasticDistance_track0", &plastic_track0);
    tree_input -> SetBranchAddress("Sample_track1", &sample_track1);
    tree_input -> SetBranchAddress("TrackAngle_track1", &angle_track1);
    tree_input -> SetBranchAddress("IronDistance_track1", &iron_track1);
    tree_input -> SetBranchAddress("PlasticDistance_track1", &plastic_track1);

    const unsigned int nevents = tree_input -> GetEntries();
    for(unsigned int i = 0; i < nevents; ++i)
    {
        tree_input -> GetEntry(i);
        if(i % 5000 == 0)
            std::cout << TAG << "Processing event " << i << std::endl;

        if(selected_sample != 0)
            continue;

        //if(mom_true == 0 || angle_true == 0)
        //    continue;

        //Perform INGRID selection
        switch(muon_track)
        {
            case 0:
                pmu_reco = iron_track0 + (plastic_track0 / iron_carbon_ratio);
                angle_reco = angle_track0;
                cosmu_reco = TMath::Cos(angle_reco * deg_to_rad);
                track_sample = sample_track0;
                break;
            case 1:
                pmu_reco = iron_track1 + (plastic_track1 / iron_carbon_ratio);
                angle_reco = angle_track1;
                cosmu_reco = TMath::Cos(angle_reco * deg_to_rad);
                track_sample = sample_track1;
                break;
            default:
                pmu_reco = -999;
                angle_reco = -999;
                cosmu_reco = -999;
                break;
        }

        if(track_sample >= 6)
            continue;

        //Scale to correct variables and translate INGRID codes to fit codes
        enu_true = enu * 1000;
        enu_reco = enu * 1000;
        sample = 10;
        target = 6;
        reaction = GetIngridReaction(interaction_type);
        topology = GetIngridTopology(fsi_int);
        pmu_true = mom_true * 1000.0;
        cosmu_true = TMath::Cos(angle_true * deg_to_rad);
        weight_nom = event_weight;

        /*
        if(pmu_true == 0 && angle_true == 0)
        {
            pmu_true = 0;
            cosmu_true = 0;
        }
        */

        double emu_true = std::sqrt(pmu_true * pmu_true + mu_mass * mu_mass);
        double q2_true = 2.0 * enu_true * (emu_true - pmu_true * TMath::Cos(angle_true * deg_to_rad))
                         - mu_mass * mu_mass;
        double q0 = enu_true - emu_true;
        double q3 = std::sqrt(q2_true + q0*q0);

        /*
        std::cout << "Pmu : " << pmu_true << std::endl
                  << "Cos : " << TMath::Cos(angle_true * deg_to_rad) << std::endl
                  << "Emu : " << emu_true << std::endl
                  << "Enu : " << enu_true << std::endl;

        TLorentzVector nu(0,0,enu_true,enu_true);
        TLorentzVector mu(0, pmu_true * TMath::Sin(angle_true * deg_to_rad), pmu_true * TMath::Cos(angle_true * deg_to_rad), emu_true);
        std::cout << "Q2: " << (nu - mu)*(nu - mu) << std::endl;
        std::cout << "q0: " << (nu - mu).E() << std::endl;
        std::cout << "q3: " << (nu - mu).Vect().Mag() << std::endl;
        */

        //Create by hand our NIWG event for NIWGRW dials (e.g. 2p2h shape)
        niwg::rew::NIWGEvent nwg_event;
        nwg_event.detid = 0;
        nwg_event.targetA = 12;
        nwg_event.neutmode = interaction_type;
        nwg_event.nuid = GetIngridNutype(is_anti, is_nue);
        nwg_event.truenu = enu;
        nwg_event.Q2 = q2_true / 1.0E6;
        nwg_event.q0 = q0 / 1.0E3;
        nwg_event.q3 = q3 / 1.0E3;

        //NIWGRW expects a valid particle stack in the event, but the INGRID file has no such
        //information. So we create an empty "fake" stack and NIWGRW is happy.
        std::vector<niwg::rew::NIWGPartStack> fake_particle_stack = {niwg::rew::NIWGPartStack()};
        nwg_event.part_stack = fake_particle_stack;
        //nwg_event.Print();

        for(unsigned int step = 0; step < dial_steps; ++step)
        {
            const int idx = GetDialIndex(dial_name, step);
            if(idx >= 0)
            {
                //Undo the 2015 NIWG tuning by dividing by reweight_array[145]
                //This is mainly for removing the SF->RFG weighting
                if(reweight_array[idx] == 0 || reweight_array[146] == 0)
                    weight_syst[step] = 1;
                else
                    weight_syst[step] = reweight_array[idx] / reweight_array[145];
                //std::cout << i << ": " << idx << " " << reweight_array[idx] << std::endl;

                /*
                std::cout << "----------" << std::endl
                          << "Evt: " << i << std::endl
                          << "IDX: " << idx << std::endl
                          << "R  : " << reaction << std::endl
                          << "T  : " << topology << std::endl
                          //<< "W  : " << reweight_array[idx] << std::endl
                          << "W  : " << reweight_array[idx] / reweight_array[146] << std::endl
                          << "146: " << reweight_array[146] << std::endl;
                */
            }
            else
            {
                rw.Systematics().SetTwkDial(t2ksyst, v_dial_val.at(step));
                rw.Reconfigure();

                t2krew::T2KNIWGReWeight* niwg_rw = (t2krew::T2KNIWGReWeight*)rw.WghtEngine("niwg_rw");
                weight_syst[step] = niwg_rw -> CalcWeight(&nwg_event);
                //weight_syst[step] = 1.0;
            }

            if(weight_syst[step] > 5.0)
            {
                weight_syst[step] = 5.0;
                rejected_weights++;
            }

            //weight_syst[step] = 1.0;
            total_weights++;
        }

        tree_output -> Fill();
    }
    file_input -> Close();

    std::cout << TAG << "Rejected weights: " << rejected_weights << std::endl;
    std::cout << TAG << "Total weights: " << total_weights << std::endl;
    std::cout << TAG << "Saving results." << std::endl;

    file_output -> cd();
    tree_output -> Write();
    file_output -> Close();

    return 0;
}

int GetDialIndex(const std::string& dial_name, int step)
{
    //Programmer's note: This hurt to write as much as it does to look at.
    int offset = -100;
    if(dial_name == "MACCQE")
        offset = 42;
    else if(dial_name == "MEC_C")
        offset = 28;
    else if(dial_name == "PF_C")
        offset = 0;
    else if(dial_name == "CA5")
        offset = 49;
    else if(dial_name == "MARES")
        offset = 56;
    else if(dial_name == "I12RES")
        offset = 63;
    else if(dial_name == "CCCOH")
        offset = 77;
    else if(dial_name == "NCCOH")
        offset = 84;
    else if(dial_name == "NCOTH")
        offset = 91;
    else if(dial_name == "DIS")
        offset = 70;
    else if(dial_name == "FSI_PI_ABS")
        offset = 98;
    else if(dial_name == "FSI_PI_PROD")
        offset = 133;
    else if(dial_name == "FSI_CEX_LO")
        offset = 105;
    else if(dial_name == "FSI_CEX_HI")
        offset = 112;
    else if(dial_name == "FSI_INEL_LO")
        offset = 119;
    else if(dial_name == "FSI_INEL_HI")
        offset = 126;

    return offset + step;
}

int GetIngridNutype(bool anti, bool nue)
{
    if(nue)
        return anti ? -12 : 12;
    else
        return anti ? -14 : 14;
}

int GetIngridReaction(int code)
{
    int reaction = -1;

    if(code == 1)
        reaction = 0;
    else if(code == 2)
        reaction = 9;
    else if(code == 11 || code == 12 || code == 13)
        reaction = 1;
    else if(code == 21 || code == 26)
        reaction = 2;
    else if(code == 16)
        reaction = 3;
    else if(code > 30 && code < 53)
        reaction = 4;
    else if(code < 0)
        reaction = 5;

    return reaction;
}

int GetIngridTopology(int code)
{
    int topology;
    switch(code)
    {
        case 0:
            topology = 0;
            break;
        case 1:
            topology = 1;
            break;
        case 2:
            topology = 2;
            break;
        case 3:
            topology = 3;
            break;
        case 4:
            topology = 5;
            break;
        case 5:
            topology = 4;
            break;
        default:
            topology = 5;
            break;
    }

    return topology;
}

