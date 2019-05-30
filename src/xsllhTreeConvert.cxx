#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>

#include "json.hpp"
using json = nlohmann::json;

#include "ColorOutput.hh"
#include "ProgressBar.hh"

struct HL2TreeVar
{
    std::string reaction;
    std::string topology;
    std::string target;
    std::string nutype;
    std::string enu_true;
    std::string enu_reco;
    std::string weight;
    std::string D1True;
    std::string D1Reco;
    std::string D2True;
    std::string D2Reco;
};

struct HL2FileOpt
{
    std::string fname_input;
    std::string sel_tree;
    std::string tru_tree;
    unsigned int file_id;
    unsigned int num_branches;
    std::vector<int> cuts;
    std::map<int, std::vector<int>> samples;

    HL2TreeVar sel_var;
    HL2TreeVar tru_var;
};

struct INGFileOpt
{
    std::string fname_input;
    std::string event_tree;
    unsigned int file_id;
    unsigned int sample;
    unsigned int branch;
};

template <typename T>
HL2TreeVar ParseHL2Var(T json_obj, bool flag);

int GetIngridNutype(bool anti, bool nue);
int GetIngridReaction(int code);
int GetIngridTopology(int code);

double GetTestEnuWeight(double enu);

int main(int argc, char** argv)
{
    const std::string TAG = color::GREEN_STR + "[xsTreeConvert]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + color::BOLD_STR
                            + "[ERROR]: " + color::RESET_STR;

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Reading Events "));

    std::cout << "-------------------------------------------------\n"
              << TAG << "Welcome to the Super-xsLLh Tree Converter.\n"
              << TAG << "Initializing the tree machinery..." << std::endl;

    const double POT_norm = 0.0669;
    bool do_apply_weights = false;
    std::string json_file;

    char option;
    while((option = getopt(argc, argv, "j:Wh")) != -1)
    {
        switch(option)
        {
            case 'j':
                json_file = optarg;
                break;
            case 'W':
                do_apply_weights = true;
                break;
            case 'h':
                std::cout << "USAGE: "
                          << argv[0] << "\nOPTIONS:\n"
                          << "-j : JSON input\n"
                          << "-W : Enable extra/arbitrary event weights.\n";
            default:
                return 0;
        }
    }

    std::fstream f;
    f.open(json_file, std::ios::in);
    std::cout << TAG << "Opening " << json_file << std::endl;
    if(!f.is_open())
    {
        std::cout << ERR << "Unable to open JSON configure file." << std::endl;
        return 1;
    }

    json j;
    f >> j;

    std::cout << TAG << "Reading configuration options..." << std::endl;
    std::string out_fname = j["output"]["fname"];
    std::string out_seltree_name = j["output"]["sel_tree"];
    std::string out_trutree_name = j["output"]["tru_tree"];

    std::cout << TAG << "Out File: " << out_fname << std::endl
              << TAG << "Out Selection Tree: " << out_seltree_name << std::endl
              << TAG << "Out Truth Tree    : " << out_trutree_name << std::endl;

    TFile* out_file = TFile::Open(out_fname.c_str(), "RECREATE");
    TTree* out_seltree = new TTree(out_seltree_name.c_str(), out_seltree_name.c_str());
    TTree* out_trutree = new TTree(out_trutree_name.c_str(), out_trutree_name.c_str());

    const float mu_mass = 105.658374;
    int nutype, nutype_true;
    int reaction, reaction_true;
    int topology, topology_true;
    int target, target_true;
    int cut_branch;
    float enu_true, enu_reco;
    float q2_true, q2_reco;
    float D1True, D1Reco;
    float D2True, D2Reco;
    float weight, weight_true;

    float selmu_mom_range;

    out_seltree -> Branch("nutype", &nutype, "nutype/I");
    out_seltree -> Branch("reaction", &reaction, "reaction/I");
    out_seltree -> Branch("topology", &topology, "topology/I");
    out_seltree -> Branch("target", &target, "target/I");
    out_seltree -> Branch("cut_branch", &cut_branch, "cut_branch/I");
    out_seltree -> Branch("enu_true", &enu_true, "enu_true/F");
    out_seltree -> Branch("enu_reco", &enu_reco, "enu_reco/F");
    out_seltree -> Branch("q2_true", &q2_true, "q2_true/F");
    out_seltree -> Branch("q2_reco", &q2_reco, "q2_reco/F");
    out_seltree -> Branch("D1True", &D1True, "D1True/F");
    out_seltree -> Branch("D1Reco", &D1Reco, "D1Reco/F");
    out_seltree -> Branch("D2True", &D2True, "D2True/F");
    out_seltree -> Branch("D2Reco", &D2Reco, "D2Reco/F");
    out_seltree -> Branch("weight", &weight, "weight/F");

    out_trutree -> Branch("nutype", &nutype_true, "nutype/I");
    out_trutree -> Branch("reaction", &reaction_true, "reaction/I");
    out_trutree -> Branch("topology", &topology_true, "topology/I");
    out_trutree -> Branch("target", &target_true, "target/I");
    out_trutree -> Branch("cut_branch", &cut_branch, "cut_branch/I");
    out_trutree -> Branch("enu_true", &enu_true, "enu_true/F");
    out_trutree -> Branch("q2_true", &q2_true, "q2_true/F");
    out_trutree -> Branch("D1True", &D1True, "D1True/F");
    out_trutree -> Branch("D2True", &D2True, "D2True/F");
    out_trutree -> Branch("weight", &weight_true, "weight/F");

    std::vector<HL2FileOpt> v_files;
    for(const auto& file : j["highland_files"])
    {
        if(file["use"])
        {
            HL2FileOpt f;
            f.fname_input = file["fname"];
            f.sel_tree = file["sel_tree"];
            f.tru_tree = file["tru_tree"];
            f.file_id = file["file_id"];
            f.num_branches = file["num_branches"];
            f.cuts = file["cut_level"].get<std::vector<int>>();

            std::map<std::string, std::vector<int>> temp_json = file["samples"];
            for(const auto& kv : temp_json)
                f.samples.emplace(std::make_pair(std::stoi(kv.first), kv.second));

            f.sel_var = ParseHL2Var(file["sel_var"], true);
            f.tru_var = ParseHL2Var(file["tru_var"], false);

            v_files.emplace_back(f);
        }
    }

    for(const auto& file : v_files)
    {
        std::cout << TAG << "Reading file: " << file.fname_input << std::endl
                  << TAG << "File ID: " << file.file_id << std::endl
                  << TAG << "Selected tree: " << file.sel_tree << std::endl
                  << TAG << "Truth tree: " << file.tru_tree << std::endl
                  << TAG << "Num. Branches: " << file.num_branches << std::endl;

        std::cout << TAG << "Branch to Sample mapping:" << std::endl;
        for(const auto& kv : file.samples)
        {
            std::cout << TAG << "Sample " << kv.first << ": ";
            for(const auto& b : kv.second)
                std::cout << b << " ";
            std::cout << std::endl;
        }

        TFile* hl2_file = TFile::Open(file.fname_input.c_str(), "READ");
        TTree* hl2_seltree = (TTree*)hl2_file -> Get(file.sel_tree.c_str());
        TTree* hl2_trutree = (TTree*)hl2_file -> Get(file.tru_tree.c_str());

        int accum_level[1][file.num_branches];

        hl2_seltree -> SetBranchAddress("accum_level", &accum_level);
        hl2_seltree -> SetBranchAddress(file.sel_var.nutype.c_str(), &nutype);
        hl2_seltree -> SetBranchAddress(file.sel_var.reaction.c_str(), &reaction);
        hl2_seltree -> SetBranchAddress(file.sel_var.topology.c_str(), &topology);
        hl2_seltree -> SetBranchAddress(file.sel_var.target.c_str(), &target);
        hl2_seltree -> SetBranchAddress(file.sel_var.D1Reco.c_str(), &D1Reco);
        hl2_seltree -> SetBranchAddress(file.sel_var.D2Reco.c_str(), &D2Reco);
        hl2_seltree -> SetBranchAddress(file.sel_var.D1True.c_str(), &D1True);
        hl2_seltree -> SetBranchAddress(file.sel_var.D2True.c_str(), &D2True);
        hl2_seltree -> SetBranchAddress(file.sel_var.enu_true.c_str(), &enu_true);
        hl2_seltree -> SetBranchAddress(file.sel_var.enu_reco.c_str(), &enu_reco);
        hl2_seltree -> SetBranchAddress(file.sel_var.weight.c_str(), &weight);

        hl2_seltree -> SetBranchAddress("selmu_mom_range_oarecon", &selmu_mom_range);

        long int npassed = 0;
        long int nevents = hl2_seltree -> GetEntries();
        std::cout << TAG << "Reading selected events tree." << std::endl
                  << TAG << "Num. events: " << nevents << std::endl;

        for(int i = 0; i < nevents; ++i)
        {
            hl2_seltree -> GetEntry(i);

            bool event_passed = false;
            for(const auto& kv : file.samples)
            {
                for(const auto& branch : kv.second)
                {
                    if(accum_level[0][branch] > file.cuts[branch])
                    {
                        cut_branch = kv.first;
                        event_passed = true;
                        npassed++;
                        break;
                    }
                }
            }

            if(reaction < 0 || topology < 0)
            {
                npassed--;
                continue;
            }

            if(cut_branch == 3 || cut_branch == 4)
                D1Reco = selmu_mom_range;

            float selmu_mom = D1Reco;
            float selmu_cos = D2Reco;
            float selmu_mom_true = D1True;
            float selmu_cos_true = D2True;

            double emu_true = std::sqrt(selmu_mom_true * selmu_mom_true + mu_mass * mu_mass);
            q2_true = 2.0 * enu_true * (emu_true - selmu_mom_true * selmu_cos_true)
                - mu_mass * mu_mass;

            double emu_reco = std::sqrt(selmu_mom * selmu_mom + mu_mass * mu_mass);
            q2_reco = 2.0 * enu_reco * (emu_reco - selmu_mom * selmu_cos)
                - mu_mass * mu_mass;

            if(do_apply_weights)
                weight = weight * GetTestEnuWeight(enu_true) * POT_norm;

            if(event_passed)
                out_seltree -> Fill();

            if(i % 2000 == 0 || i == (nevents-1))
                pbar.Print(i, nevents-1);
        }
        std::cout << TAG << "Selected events passing cuts: " << npassed << std::endl;

        hl2_trutree -> SetBranchAddress(file.tru_var.nutype.c_str(), &nutype_true);
        hl2_trutree -> SetBranchAddress(file.tru_var.reaction.c_str(), &reaction_true);
        hl2_trutree -> SetBranchAddress(file.tru_var.topology.c_str(), &topology_true);
        hl2_trutree -> SetBranchAddress(file.tru_var.target.c_str(), &target_true);
        hl2_trutree -> SetBranchAddress(file.tru_var.D1True.c_str(), &D1True);
        hl2_trutree -> SetBranchAddress(file.tru_var.D2True.c_str(), &D2True);
        hl2_trutree -> SetBranchAddress(file.tru_var.enu_true.c_str(), &enu_true);
        hl2_trutree -> SetBranchAddress(file.tru_var.weight.c_str(), &weight_true);

        nevents = hl2_trutree -> GetEntries();
        std::cout << TAG << "Reading truth events tree." << std::endl
                  << TAG << "Num. events: " << nevents << std::endl;
        for(int i = 0; i < nevents; ++i)
        {
            hl2_trutree -> GetEntry(i);
            cut_branch = -1;

            if(reaction_true < 0 || topology_true < 0)
                continue;

            float selmu_mom_true = D1True;
            float selmu_cos_true = D2True;

            double emu_true = std::sqrt(selmu_mom_true * selmu_mom_true + mu_mass * mu_mass);
            q2_true = 2.0 * enu_true * (emu_true - selmu_mom_true * selmu_cos_true)
                - mu_mass * mu_mass;

            if(do_apply_weights)
                weight_true = weight_true * GetTestEnuWeight(enu_true) * POT_norm;

            out_trutree -> Fill();

            if(i % 2000 == 0 || i == (nevents-1))
                pbar.Print(i, nevents-1);
        }

        hl2_file -> Close();
    }

    std::cout << TAG << "Reading INGRID files..." << std::endl;
    std::vector<INGFileOpt> v_files_ing;
    for(const auto& file : j["ingrid_files"])
    {
        if(file["use"])
        {
            INGFileOpt f;
            f.fname_input = file["fname"];
            f.event_tree = file["event_tree"];
            f.file_id = file["file_id"];
            f.sample = file["sample"];
            f.branch = file["branch"];

            v_files_ing.emplace_back(f);
        }
    }

    for(const auto& file : v_files_ing)
    {
        std::cout << TAG << "Reading file: " << file.fname_input << std::endl
                  << TAG << "File ID: " << file.file_id << std::endl
                  << TAG << "Events tree: " << file.event_tree << std::endl
                  << TAG << "Sample: " << file.sample << std::endl
                  << TAG << "Branch: " << file.branch << std::endl;

        TFile* ing_file = TFile::Open(file.fname_input.c_str(), "READ");
        TTree* ing_tree = (TTree*)ing_file -> Get(file.event_tree.c_str());

        int interaction_type, fsi_int;
        int muon_track, selected_sample;
        int track_sample;
        float enu, event_weight;
        float pmu_true, angle_true;
        float pmu_reco, angle_reco;
        bool is_anti, is_nue, is_fv;
        bool new_event;

        ing_tree -> SetBranchAddress("InteractionType", &interaction_type);
        ing_tree -> SetBranchAddress("FSIInt", &fsi_int);
        ing_tree -> SetBranchAddress("SelectedSample", &selected_sample);
        ing_tree -> SetBranchAddress("MuonCandidateTrack", &muon_track);
        ing_tree -> SetBranchAddress("Enu", &enu);
        ing_tree -> SetBranchAddress("IsAnti", &is_anti);
        ing_tree -> SetBranchAddress("IsNuE", &is_nue);
        ing_tree -> SetBranchAddress("IsFV", &is_fv);
        ing_tree -> SetBranchAddress("NewEvent", &new_event);
        ing_tree -> SetBranchAddress("weight", &event_weight);
        ing_tree -> SetBranchAddress("TrueMomentumMuon", &pmu_true);
        ing_tree -> SetBranchAddress("TrueAngleMuon", &angle_true);

        const float deg_to_rad = TMath::Pi() / 180;
        const float iron_carbon_ratio = 7.640777;
        int sample_track0, sample_track1, sample_track2;
        float iron_track0, plastic_track0, angle_track0;
        float iron_track1, plastic_track1, angle_track1;
        float iron_track2, plastic_track2, angle_track2;

        ing_tree -> SetBranchAddress("Sample_track0", &sample_track0);
        ing_tree -> SetBranchAddress("TrackAngle_track0", &angle_track0);
        ing_tree -> SetBranchAddress("IronDistance_track0", &iron_track0);
        ing_tree -> SetBranchAddress("PlasticDistance_track0", &plastic_track0);
        ing_tree -> SetBranchAddress("Sample_track1", &sample_track1);
        ing_tree -> SetBranchAddress("TrackAngle_track1", &angle_track1);
        ing_tree -> SetBranchAddress("IronDistance_track1", &iron_track1);
        ing_tree -> SetBranchAddress("PlasticDistance_track1", &plastic_track1);
        ing_tree -> SetBranchAddress("Sample_track2", &sample_track2);
        ing_tree -> SetBranchAddress("TrackAngle_track2", &angle_track2);
        ing_tree -> SetBranchAddress("IronDistance_track2", &iron_track2);
        ing_tree -> SetBranchAddress("PlasticDistance_track2", &plastic_track2);

        const unsigned int nevents = ing_tree -> GetEntries();
        std::cout << TAG << "Reading events tree." << std::endl
                  << TAG << "Num. events: " << nevents << std::endl;

        double sel_events = 0;
        double gen_events = 0;
        for(unsigned int i = 0; i < nevents; ++i)
        {
            ing_tree -> GetEntry(i);

            switch(muon_track)
            {
                case 0:
                    pmu_reco = iron_track0 + (plastic_track0 / iron_carbon_ratio);
                    angle_reco = angle_track0;
                    track_sample = sample_track0;
                    break;
                case 1:
                    pmu_reco = iron_track1 + (plastic_track1 / iron_carbon_ratio);
                    angle_reco = angle_track1;
                    track_sample = sample_track1;
                    break;
                case 2:
                    pmu_reco = iron_track2 + (plastic_track2 / iron_carbon_ratio);
                    angle_reco = angle_track2;
                    track_sample = sample_track2;
                    break;
                default:
                    pmu_reco = -999;
                    angle_reco = -999;
                    track_sample = 999;
                    break;
            }

            //D1Reco = (pmu_reco * 0.0114127 + 0.230608) * 1000;
            D1Reco = pmu_reco;
            D1True = pmu_true * 1000;
            D2Reco = angle_reco;
            //D2Reco = TMath::Cos(angle_reco * deg_to_rad);
            D2True = TMath::Cos(angle_true * deg_to_rad);

            nutype = GetIngridNutype(is_anti, is_nue);
            topology = GetIngridTopology(fsi_int);
            reaction = GetIngridReaction(interaction_type);
            target = 6;

            enu_true = enu * 1000.0;
            enu_reco = enu * 1000.0;

            double emu_true = std::sqrt(pmu_true * pmu_true + mu_mass * mu_mass);
            q2_true = 2.0 * enu_true * (emu_true - pmu_true * TMath::Cos(angle_true * deg_to_rad))
                - mu_mass * mu_mass;

            double emu_reco = std::sqrt(pmu_reco * pmu_reco + mu_mass * mu_mass);
            q2_reco = 2.0 * enu_reco * (emu_reco - pmu_reco * TMath::Cos(angle_reco * deg_to_rad))
                - mu_mass * mu_mass;

            cut_branch = file.branch;
            weight = event_weight;

            if(do_apply_weights)
                weight = weight * GetTestEnuWeight(enu_true);

            if(selected_sample == 1)
                weight *= 1.16;

            if(selected_sample == file.sample && track_sample < 6)
            {
                out_seltree -> Fill();
                sel_events += weight;
            }

            nutype_true = nutype;
            reaction_true = reaction;
            topology_true = topology;
            target_true = target;
            cut_branch = file.branch * -1;
            weight_true = weight;

            if(is_fv && !is_anti && !is_nue && fsi_int < 3 && new_event == 1)
            {
                out_trutree -> Fill();
                gen_events += weight;
            }

            if(i % 2000 == 0 || i == (nevents-1))
                pbar.Print(i, nevents-1);
        }

        std::cout << TAG << "INGRID Selected Ev : " << sel_events << std::endl;
        std::cout << TAG << "INGRID Generated Ev: " << gen_events << std::endl;
    }

    out_file -> cd();
    out_file -> Write();
    out_file -> Close();
    std::cout << TAG << "Finished." << std::endl;

    return 0;
}

template <typename T>
HL2TreeVar ParseHL2Var(T j, bool reco_info)
{
    HL2TreeVar v;
    v.reaction = j["reaction"];
    v.topology = j["topology"];
    v.target   = j["target"];
    v.nutype   = j["nutype"];
    v.enu_true = j["enu_true"];
    v.weight = j["weight"];

    v.D1True = j["D1True"];
    v.D2True = j["D2True"];

    if(reco_info)
    {
        v.enu_reco = j["enu_reco"];
        v.D1Reco = j["D1Reco"];
        v.D2Reco = j["D2Reco"];
    }

    return v;
}

int GetIngridNutype(bool anti, bool nue)
{
    if(nue)
        return anti ? 3 : 2;
    else
        return anti ? 1 : 0;
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
    else if(code > 30 || code < 53)
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

double GetTestEnuWeight(double enu)
{
    double weight = 1.0;

    if(enu < 500)
        weight = 1 + (0.5 * enu) / 500.0;
    else if(enu > 500 && enu < 1500)
        weight = 2 - enu / 1000.0;
    else if(enu > 1500 && enu < 2000)
        weight = (0.5 * enu) / 500.0 - 1;
    else
        weight = 1.0;

    return weight;
}
