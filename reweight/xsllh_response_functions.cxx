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

int RenameHLReaction(int code);
int RenameHLTopology(int code);
int RenameHLReaction_anti(int code);
int RenameHLTopology_anti(int code);

int main(int argc, char** argv)
{
    std::cout << "--------------------------------------------------------\n"
              << TAG << "Welcome to the Super-xsLLhReWeight Interface.\n"
              << TAG << "Initializing the reweight machinery..." << std::endl;

    bool limits = false;
    bool use_truth_tree = false;
    bool do_dummy_splines = false;
    unsigned int dial_steps = 15;
    double nominal = 0;
    double error_neg = 0;
    double error_pos = 0;
    std::string fname_input_numu;
    std::string fname_input_antinumu;
    std::string fname_output;
    std::string fname_dials;
    std::string dial_name;
    std::string dial_name_enum;

    const int num_samples = 12;
    const int cut_level[num_samples] = {7, 9, 8, 10, 9, 6, 4, 7, 7, 8, 10, 9};
    const int num_samples_anti = 12;
    const int cut_level_anti[num_samples_anti] = {8, 10, 7, 8, 10, 7, 8, 10, 7, 8, 10, 7};

    char option;
    while((option = getopt(argc, argv, "i:a:o:r:n:d:LTDh")) != -1)
    {
        switch(option)
        {
            case 'i':
                fname_input_numu = optarg;
                break;
            case 'a':
                fname_input_antinumu = optarg;
                break;
            case 'o':
                fname_output = optarg;
                break;
            case 'r':
                dial_name = optarg;
                break;
            case 'n':
                dial_steps = std::stoi(optarg);
                break;
            case 'd':
                fname_dials = optarg;
                break;
            case 'L':
                limits = true;
                break;
            case 'T':
                use_truth_tree = true;
                break;
            case 'D':
                do_dummy_splines = true;
                break;
            case 'h':
                std::cout << "USAGE: " << argv[0] << "\nOPTIONS\n"
                          << "-i : Input list of Highland trees (.txt)\n"
                          << "-a : Input list of antinumu Highland trees (.txt)\n"
                          << "-o : Output ROOT filename\n"
                          << "-r : Systematic parameter to reweight\n"
                          << "-n : Number of dial steps\n"
                          << "-d : Dial value file (.txt)\n"
                          << "-L : Read errors as upper and lower limits\n"
                          << "-T : Use truth tree for reweighting\n"
                          << "-D : Make inputs for dummy splines\n"
                          << "-h : Display this help message\n";
            default:
                return 0;
        }
    }

    if(fname_input_numu.empty() || fname_input_antinumu.empty() || fname_output.empty() || dial_name.empty())
    {
        std::cout << ERR << "Missing necessary command line arguments.\n" << std::endl;
        return 1;
    }

    std::cout << TAG << "Input numu Highland list: " << fname_input_numu << std::endl
              << TAG << "Input antinumu Highland list: " << fname_input_antinumu << std::endl
              << TAG << "Output weights file: " << fname_output << std::endl
              << TAG << "Reading tree: " << (use_truth_tree ? "truth" : "default") << std::endl
              << TAG << "Input dial value file: " << fname_dials << std::endl
              << TAG << "Generating weights for " << dial_name << std::endl
              << TAG << "Number of dial steps: " << dial_steps << std::endl;

    if(do_dummy_splines)
        std::cout << TAG << "Producing dummy spline inputs." << std::endl;

    #ifdef OA_FAST
    std::cout << TAG << "Using RooVtxIndex for vertex finding." << std::endl;
    #else
    std::cout << TAG << "Using TruthVtxID for vertex finding." << std::endl;
    #endif

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

    std::vector<std::string> hl_filenames;
    std::ifstream fin(fname_input_numu, std::ios::in);
    std::cout << TAG << "Reading list of Highland files." << std::endl;
    if(!fin.is_open())
    {
        std::cerr << ERR << "Failed to open " << fname_input_numu << std::endl;
        return 1;
    }
    else
    {
        std::string name;
        while(std::getline(fin, name))
        {
            hl_filenames.emplace_back(name);
        }
        for(const auto& fname : hl_filenames)
            std::cout << TAG << "Input: " << fname << std::endl;
    }

    std::vector<std::string> hl_filenames_anti;
    std::ifstream fin_anti(fname_input_antinumu, std::ios::in);
    std::cout << TAG << "Reading list of antinumu Highland files." << std::endl;
    if(!fin_anti.is_open())
    {
        std::cerr << ERR << "Failed to open " << fname_input_antinumu << std::endl;
        return 1;
    }
    else
    {
        std::string name;
        while(std::getline(fin_anti, name))
        {
            hl_filenames_anti.emplace_back(name);
        }
        for(const auto& fname : hl_filenames_anti)
            std::cout << TAG << "Input: " << fname << std::endl;
    }

    float enu_true, enu_reco;
    float pmu_true, pmu_reco;
    float cosmu_true, cosmu_reco;
    float weight_nom;
    float weight_syst[dial_steps];
    std::string weight_str = "weight_syst[" + std::to_string(dial_steps) + "]/F";

    int topology, topology_mod, reaction, reaction_mod, sample, target;
    int accum_level[num_samples][num_samples];
    int input_RooVtxIndex, input_RooVtxEntry;
    int TruthVertexID;
    float input_weight;
    float pmu_range;
    float s_mu_true_mom;
    float s_mu_true_costheta;
    float s_selmu_true_mom;
    float s_selmu_true_dir[3];

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
            //double inc = (error_pos - error_neg) / (dial_steps - 1);
            double lim = step < 0 ? error_neg : error_pos;
            double inc = 2 * std::abs(lim - nominal) / (dial_steps - 1);
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

    std::cout << TAG << "Dial Values: " << std::endl;
    for(const auto& v : v_dial_val)
        std::cout << nominal * (1 + v) << std::endl;

    TFile* file_output = TFile::Open(fname_output.c_str(), "RECREATE");
    TTree* tree_output = new TTree("selectedEvents", "selectedEvents");

    tree_output->Branch("Enutrue", &enu_true, "enu_true/F");
    tree_output->Branch("Enureco", &enu_reco, "enu_reco/F");
    tree_output->Branch("sample", &sample, "sample/I");
    tree_output->Branch("target", &target, "target/I");
    tree_output->Branch("reaction", &reaction_mod, "reaction/I");
    tree_output->Branch("topology", &topology_mod, "topology/I");
    tree_output->Branch("Pmutrue", &pmu_true, "pmu_true/F");
    tree_output->Branch("Pmureco", &pmu_reco, "pmu_reco/F");
    tree_output->Branch("CTHmutrue", &cosmu_true, "cosmu_true/F");
    tree_output->Branch("CTHmureco", &cosmu_reco, "cosmu_reco/F");
    tree_output->Branch("weight", &weight_nom, "weight_nom/F");
    tree_output->Branch("weight_syst", weight_syst, weight_str.c_str());

    std::cout << TAG << "Initializing T2KReWeight and reweight engines." << std::endl;
    t2krew::T2KNeutUtils::SetCardFile("/afs/cern.ch/work/c/cschloes/T2K_final_for_real/t2kreweight/T2KReWeight/neut.d.card");
    t2krew::T2KReWeight rw;
    t2krew::T2KSyst_t t2ksyst = t2krew::T2KSyst::FromString(dial_name_enum);
    rw.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
    rw.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());

    if(dial_name == "MACCQE")
    {
        //rw.Systematics().Include(t2krew::kNXSec_VecFFCCQE);
        //rw.Systematics().SetTwkDial(t2krew::kNXSec_VecFFCCQE, 2);
    }

    rw.Systematics().Include(t2ksyst);
    rw.Systematics().SetAbsTwk(t2ksyst);
    rw.Systematics().PrintSummary();

    for(unsigned int f = 0; f < hl_filenames.size(); ++f)
    {
        TFile* file_input = TFile::Open(hl_filenames.at(f).c_str(), "READ");

        TTree* tree_RooVtx = (TTree*)file_input -> Get("NRooTrackerVtx");
        TClonesArray* nRooVtxs = new TClonesArray("ND::NRooTrackerVtx");
        int NRooVtx = 0;

        tree_RooVtx -> SetBranchAddress("Vtx", &nRooVtxs);
        tree_RooVtx -> SetBranchAddress("NVtx", &NRooVtx);

        std::string tree_name = use_truth_tree ? "truth" : "default";
        TTree* tree_input = (TTree*)file_input -> Get(tree_name.c_str());
        std::cout << TAG << "Reading file: " << hl_filenames.at(f) << std::endl;

        if(use_truth_tree)
        {
            tree_input->SetBranchAddress("topology", &topology);
            tree_input->SetBranchAddress("reaction", &reaction);
            tree_input->SetBranchAddress("target", &target);
            tree_input->SetBranchAddress("truelepton_mom", &s_mu_true_mom);
            tree_input->SetBranchAddress("truelepton_costheta", &s_mu_true_costheta);
            tree_input->SetBranchAddress("accum_level", accum_level);
            tree_input->SetBranchAddress("nu_trueE", &enu_true);
            tree_input->SetBranchAddress("weight", &input_weight);
            tree_input->SetBranchAddress("RooVtxIndex", &input_RooVtxIndex);
            tree_input->SetBranchAddress("RooVtxEntry", &input_RooVtxEntry);
            tree_input->SetBranchAddress("TruthVertexID", &TruthVertexID);
        }

        else
        {
            tree_input->SetBranchAddress("topology", &topology);
            tree_input->SetBranchAddress("reaction", &reaction);
            tree_input->SetBranchAddress("target", &target);
            tree_input->SetBranchAddress("truelepton_mom", &s_mu_true_mom);
            tree_input->SetBranchAddress("truelepton_costheta", &s_mu_true_costheta);
            tree_input->SetBranchAddress("selmu_truemom", &s_selmu_true_mom);
            tree_input->SetBranchAddress("selmu_truedir", s_selmu_true_dir);
            tree_input->SetBranchAddress("selmu_mom", &pmu_reco);
            tree_input->SetBranchAddress("selmu_costheta", &cosmu_reco);
            tree_input->SetBranchAddress("selmu_mom_range_oarecon", &pmu_range);
            tree_input->SetBranchAddress("accum_level", accum_level);
            tree_input->SetBranchAddress("nu_e", &enu_reco);
            tree_input->SetBranchAddress("nu_trueE", &enu_true);
            tree_input->SetBranchAddress("weight_corr_total", &input_weight);
            tree_input->SetBranchAddress("RooVtxIndex", &input_RooVtxIndex);
            tree_input->SetBranchAddress("RooVtxEntry", &input_RooVtxEntry);
            tree_input->SetBranchAddress("TruthVertexID", &TruthVertexID);
        }

        unsigned int roovtx_events = tree_RooVtx->GetEntries();
        unsigned int select_events = tree_input->GetEntries();
        std::cout << TAG << "Default Tree Events  : " << select_events << std::endl
                  << TAG << "NRooTrackerVtx Events: " << roovtx_events << std::endl;

        for(unsigned int i = 0; i < select_events; ++i)
        {
            bool event_pass = false;
            tree_input -> GetEntry(i);
            tree_RooVtx -> LoadTree(input_RooVtxEntry);
            tree_RooVtx -> GetEntry(input_RooVtxEntry);

            #ifdef OA_FAST
            ND::NRooTrackerVtx* vtx = (ND::NRooTrackerVtx*)nRooVtxs->At(input_RooVtxIndex);
            #else
            ND::NRooTrackerVtx* vtx = nullptr;
            for(int v = 0; v < NRooVtx; ++v)
            {
                vtx = (ND::NRooTrackerVtx*)nRooVtxs->At(v);
                if(vtx->TruthVertexID == TruthVertexID)
                {
                    #ifdef DEBUG_MSG
                    std::cout << "NVtx = " << NRooVtx << std::endl;
                    std::cout << "Vtx->" << vtx->TruthVertexID << " vs. " << TruthVertexID << std::endl;
                    std::cout << "RooVtxIdx: " << input_RooVtxIndex << " vs. " << v << std::endl;
                    #endif
                    break;
                }
            }
            #endif

            if(i % 1000 == 0)
                std::cout << TAG << "Processing event " << i << std::endl;

            for(unsigned int s = 0; s < num_samples; ++s)
            {
                if(accum_level[0][s] > cut_level[s])
                {
                    event_pass = true;
                    sample = s;
                    break;
                }
            }

            event_pass = event_pass || use_truth_tree;

            if(pmu_reco < 0 || pmu_reco > 30000 || event_pass == false)
                continue;

            weight_nom = input_weight;

            if(s_mu_true_mom > 0 && s_mu_true_costheta > -2)
            {
                pmu_true = s_mu_true_mom;
                cosmu_true = s_mu_true_costheta;
            }
            else
            {
                pmu_true = s_selmu_true_mom;
                cosmu_true = s_selmu_true_dir[2]
                             / (std::sqrt(s_selmu_true_dir[0] * s_selmu_true_dir[0]
                                          + s_selmu_true_dir[1] * s_selmu_true_dir[1]
                                          + s_selmu_true_dir[2] * s_selmu_true_dir[2]));
            }

            for(unsigned int step = 0; step < dial_steps; ++step)
            {
                rw.Systematics().SetTwkDial(t2ksyst, v_dial_val.at(step));
                rw.Reconfigure();

                if(!vtx)
                {
                    if(step == 0)
                    {
                        std::cout << TAG << "No vertex for event " << i << std::endl;
                        std::cout << TAG << "NVtx = " << NRooVtx << std::endl;
                        std::cout << TAG << "RooVtxIndex = " << input_RooVtxIndex << std::endl;
                    }

                    continue;
                }
                else
                {
                    weight_syst[step] = do_dummy_splines ? 1.0 : rw.CalcWeight(vtx);
                    v_hists.at(step).Fill(pmu_reco, weight_nom * weight_syst[step]);
                }
            }

            reaction_mod = RenameHLReaction(reaction);
            topology_mod = RenameHLTopology(topology);

            tree_output -> Fill();
        }

        delete tree_input;
        delete tree_RooVtx;
        delete nRooVtxs;
        file_input -> Close();
    }




for(unsigned int f = 0; f < hl_filenames_anti.size(); ++f)
    {
        TFile* file_input = TFile::Open(hl_filenames_anti.at(f).c_str(), "READ");

        TTree* tree_RooVtx = (TTree*)file_input -> Get("NRooTrackerVtx");
        TClonesArray* nRooVtxs = new TClonesArray("ND::NRooTrackerVtx");
        int NRooVtx = 0;

        tree_RooVtx -> SetBranchAddress("Vtx", &nRooVtxs);
        tree_RooVtx -> SetBranchAddress("NVtx", &NRooVtx);

        std::string tree_name = use_truth_tree ? "truth" : "default";
        TTree* tree_input = (TTree*)file_input -> Get(tree_name.c_str());
        std::cout << TAG << "Reading file: " << hl_filenames_anti.at(f) << std::endl;

        if(use_truth_tree)
        {
            tree_input->SetBranchAddress("topology", &topology);
            tree_input->SetBranchAddress("reaction", &reaction);
            tree_input->SetBranchAddress("target", &target);
            tree_input->SetBranchAddress("truelepton_mom", &s_mu_true_mom);
            tree_input->SetBranchAddress("truelepton_costheta", &s_mu_true_costheta);
            tree_input->SetBranchAddress("accum_level", accum_level);
            tree_input->SetBranchAddress("nu_trueE", &enu_true);
            tree_input->SetBranchAddress("weight", &input_weight);
            tree_input->SetBranchAddress("RooVtxIndex", &input_RooVtxIndex);
            tree_input->SetBranchAddress("RooVtxEntry", &input_RooVtxEntry);
            tree_input->SetBranchAddress("TruthVertexID", &TruthVertexID);
        }

        else
        {
            tree_input->SetBranchAddress("topology", &topology);
            tree_input->SetBranchAddress("reaction", &reaction);
            tree_input->SetBranchAddress("target", &target);
            tree_input->SetBranchAddress("Truemu_mom", &s_mu_true_mom);
            tree_input->SetBranchAddress("Truemu_costheta", &s_mu_true_costheta);
            tree_input->SetBranchAddress("selmu_truemom", &s_selmu_true_mom);
            tree_input->SetBranchAddress("selmu_truedir", s_selmu_true_dir);
            tree_input->SetBranchAddress("selmu_mom", &pmu_reco);
            tree_input->SetBranchAddress("selmu_costheta", &cosmu_reco);
            tree_input->SetBranchAddress("selmu_MomRangeMuon", &pmu_range);
            tree_input->SetBranchAddress("accum_level", accum_level);
            tree_input->SetBranchAddress("true_Q2", &enu_reco);
            tree_input->SetBranchAddress("nu_trueE", &enu_true);
            tree_input->SetBranchAddress("weight_corr_total", &input_weight);
            tree_input->SetBranchAddress("RooVtxIndex", &input_RooVtxIndex);
            tree_input->SetBranchAddress("RooVtxEntry", &input_RooVtxEntry);
            tree_input->SetBranchAddress("TruthVertexID", &TruthVertexID);
        }

        unsigned int roovtx_events = tree_RooVtx->GetEntries();
        unsigned int select_events = tree_input->GetEntries();
        std::cout << TAG << "Default Tree Events  : " << select_events << std::endl
                  << TAG << "NRooTrackerVtx Events: " << roovtx_events << std::endl;

        for(unsigned int i = 0; i < select_events; ++i)
        {
            bool event_pass = false;
            tree_input -> GetEntry(i);
            tree_RooVtx -> LoadTree(input_RooVtxEntry);
            tree_RooVtx -> GetEntry(input_RooVtxEntry);

            #ifdef OA_FAST
            ND::NRooTrackerVtx* vtx = (ND::NRooTrackerVtx*)nRooVtxs->At(input_RooVtxIndex);
            #else
            ND::NRooTrackerVtx* vtx = nullptr;
            for(int v = 0; v < NRooVtx; ++v)
            {
                vtx = (ND::NRooTrackerVtx*)nRooVtxs->At(v);
                if(vtx->TruthVertexID == TruthVertexID)
                {
                    #ifdef DEBUG_MSG
                    std::cout << "NVtx = " << NRooVtx << std::endl;
                    std::cout << "Vtx->" << vtx->TruthVertexID << " vs. " << TruthVertexID << std::endl;
                    std::cout << "RooVtxIdx: " << input_RooVtxIndex << " vs. " << v << std::endl;
                    #endif
                    break;
                }
            }
            #endif

            if(i % 1000 == 0)
                std::cout << TAG << "Processing event " << i << std::endl;

            for(unsigned int s = 0; s < num_samples_anti; ++s)
            {
                if(accum_level[0][s] > cut_level_anti[s])
                {
                    event_pass = true;
                    sample = s;
                    break;
                }
            }

            event_pass = event_pass || use_truth_tree;

            if(pmu_reco < 0 || pmu_reco > 30000 || event_pass == false)
                continue;

            weight_nom = input_weight;

            if(s_mu_true_mom > 0 && s_mu_true_costheta > -2)
            {
                pmu_true = s_mu_true_mom;
                cosmu_true = s_mu_true_costheta;
            }
            else
            {
                pmu_true = s_selmu_true_mom;
                cosmu_true = s_selmu_true_dir[2]
                             / (std::sqrt(s_selmu_true_dir[0] * s_selmu_true_dir[0]
                                          + s_selmu_true_dir[1] * s_selmu_true_dir[1]
                                          + s_selmu_true_dir[2] * s_selmu_true_dir[2]));
            }

            for(unsigned int step = 0; step < dial_steps; ++step)
            {
                rw.Systematics().SetTwkDial(t2ksyst, v_dial_val.at(step));
                rw.Reconfigure();

                if(!vtx)
                {
                    if(step == 0)
                    {
                        std::cout << TAG << "No vertex for event " << i << std::endl;
                        std::cout << TAG << "NVtx = " << NRooVtx << std::endl;
                        std::cout << TAG << "RooVtxIndex = " << input_RooVtxIndex << std::endl;
                    }

                    continue;
                }
                else
                {
                    weight_syst[step] = do_dummy_splines ? 1.0 : rw.CalcWeight(vtx);
                    v_hists.at(step).Fill(pmu_reco, weight_nom * weight_syst[step]);
                }
            }

            reaction_mod = RenameHLReaction_anti(reaction);
            topology_mod = RenameHLTopology_anti(topology);

            tree_output -> Fill();
        }

        delete tree_input;
        delete tree_RooVtx;
        delete nRooVtxs;
        file_input -> Close();
    }





    file_output -> cd();
    tree_output -> Write();

    for(const auto& hist : v_hists)
        hist.Write();
    file_output -> Close();

    std::cout << TAG << "Finished." << std::endl;

    return 0;
}

int RenameHLReaction(int code)
{
    int reaction;
    switch(code)
    {
        case 0: // CCQE
            reaction = 0;
            break;
        case 1: // RES
            reaction = 1;
            break;
        case 2: // DIS
            reaction = 2;
            break;
        case 3: // COH
            reaction = 3;
            break;
        case 4: // NC
            reaction = 4;
            break;
        case 5: // CC Numubar
            reaction = 5;
            break;
        case 6: // CC Nuebar
            reaction = 6;
            break;
        case 7: // CC Nue
            reaction = 7;
            break;
        case 777: // Sand muon
            reaction = 8;
            break;
        case 9: // 2p2h
            reaction = 9;
            break;
        case 999: // other
            reaction = 10;
            break;
        default: // other
            reaction = 10;
            break;
    }
    return reaction;
}

int RenameHLReaction_anti(int code)
{
    int reaction;
    switch(code)
    {
        case 0: // CCQE
            reaction = 11;
            break;
        case 1: // RES
            reaction = 12;
            break;
        case 2: // DIS
            reaction = 13;
            break;
        case 3: // COH
            reaction = 14;
            break;
        case 4: // NC
            reaction = 15;
            break;
        case 5: // CC Numubar
            reaction = 16;
            break;
        case 6: // CC Nuebar
            reaction = 17;
            break;
        case 7: // CC Nue
            reaction = 18;
            break;
        case 777: // Sand muon
            reaction = 19;
            break;
        case 9: // 2p2h
            reaction = 20;
            break;
        case 999: // other
            reaction = 21;
            break;
        default: // other
            reaction = 21;
            break;
    }
    return reaction;
}

int RenameHLTopology(int code)
{
    int topology;
    switch(code)
    {
        case 0: // CC0pi
            topology = 0;
            break;
        case 1: // CC1pi
            topology = 1;
            break;
        case 2: // CCOther
            topology = 2;
            break;
        case 3: // BKG
            topology = 3;
            break;
        case 7: // OOFV
            topology = 4;
            break;
        case 777: // Sand muon
            topology = 5;
            break;
        default: // BKG
            topology = 3;
            break;
    }
    return topology;
}

int RenameHLTopology_anti(int code)
{
    int topology;
    switch(code)
    {
        case 0: // CC0pi
            topology = 6;
            break;
        case 1: // CC1pi
            topology = 7;
            break;
        case 2: // CCOther
            topology = 8;
            break;
        case 3: // BKG
            topology = 9;
            break;
        case 7: // OOFV
            topology = 10;
            break;
        case 777: // Sand muon
            topology = 11;
            break;
        default: // BKG
            topology = 9;
            break;
    }
    return topology;
}
