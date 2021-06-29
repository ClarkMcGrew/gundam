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
#include <TH1.h>
#include <TH2.h>
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
    int beammode;
    unsigned int num_branches;
    double pot_norm;
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

class TemplateWeights
{
    public:
        TemplateWeights() {};
        TemplateWeights(const std::string& filename, int num_templates)
        {
            ReadFile(filename, num_templates);
        };

        double operator()(double q0, double q3, int reaction)
        {
            if(q0 < 0 || q0 > 5000 || q3 < 0 || q3 > 5000)
                return 1.0;
            else if(reaction == 6 || reaction == 7 || reaction == 8)
                return 1.0;
            else
                return v_templates.at(reaction).Interpolate(q0, q3);
        };

        bool ReadFile(const std::string& filename, int num_templates)
        {
            std::cout << "Reading " << filename << " for templates." << std::endl;
            TFile* temp_file = TFile::Open(filename.c_str(), "READ");
            if(temp_file == nullptr)
                return false;

            for(unsigned int i = 0; i < num_templates; ++i)
            {
                std::string name = "ratio_template_" + std::to_string(i);
                TH2D* temp_hist = (TH2D*)temp_file->Get(name.c_str());
                v_templates.emplace_back(*temp_hist);
            }
            temp_file->Close();

            return true;
        }

    private:
        std::vector<TH2D> v_templates;
};

template <typename T>
HL2TreeVar ParseHL2Var(T json_obj, bool flag);

int GetIngridNutype(bool anti, bool nue);
int GetIngridReaction(int code);
int GetIngridTopology(int code);

double GetTestEnuWeight(double enu);
double GetBeRPAWeight(double q2);
double GetCCResQ2Weight(double q2);
double GetCCZeroPiQ2Weight(double q2);
double GetFluxWeightND5(double enu);
double GetFluxWeightND2(double enu);

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

    bool do_apply_weights = false;
    bool do_apply_templates = false;
    std::string json_file;

    char option;
    while((option = getopt(argc, argv, "j:TWh")) != -1)
    {
        switch(option)
        {
            case 'j':
                json_file = optarg;
                break;
            case 'W':
                do_apply_weights = true;
                break;
            case 'T':
                do_apply_templates = true;
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

    TemplateWeights nd280_templates;
    TemplateWeights ingrid_templates;

    if(do_apply_templates)
    {
        std::cout << "Applying template weights to events." << std::endl;
        nd280_templates.ReadFile("nd280_templates_nuwro.root", 10);
        ingrid_templates.ReadFile("ingrid_templates_nuwro.root", 10);
    }

    TFile* out_file = TFile::Open(out_fname.c_str(), "RECREATE");
    TTree* out_seltree = new TTree(out_seltree_name.c_str(), out_seltree_name.c_str());
    TTree* out_trutree = new TTree(out_trutree_name.c_str(), out_trutree_name.c_str());

    const float mu_mass = 105.658374;
    int nutype, nutype_true;
    int reaction, reaction_true;
    int topology, topology_true;
    int target, target_true;
    int cut_branch;
    int beammode;
    float enu_true, enu_reco;
    float q2_true, q2_reco;
    float D1True, D1Reco;
    float D2True, D2Reco;
    float weight, weight_true;

    float selmu_mom_range;
    int pdg_true;

    out_seltree -> Branch("nutype", &nutype, "nutype/I");
    out_seltree -> Branch("reaction", &reaction, "reaction/I");
    out_seltree -> Branch("topology", &topology, "topology/I");
    out_seltree -> Branch("target", &target, "target/I");
    out_seltree -> Branch("cut_branch", &cut_branch, "cut_branch/I");
    out_seltree -> Branch("beammode", &beammode, "beammode/I");
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
    out_trutree -> Branch("beammode", &beammode, "beammode/I");
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
            f.beammode = file["beammode"];
            f.num_branches = file["num_branches"];
            f.cuts = file["cut_level"].get<std::vector<int>>();
            f.pot_norm = file["pot_norm"];

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
                  << TAG << "POT Norm: " << file.pot_norm << std::endl
                  << TAG << "Beam mode: " << file.beammode << std::endl
                  << TAG << "Num. Branches: " << file.num_branches << std::endl;

        std::cout << TAG << "Branch to Sample mapping:" << std::endl;
        for(const auto& kv : file.samples)
        {
            std::cout << TAG << "Sample " << kv.first << ": ";
            for(const auto& b : kv.second)
                std::cout << b << " ";
            std::cout << std::endl;
        }

        // Beam mode was previously read in from the .json config file:
        beammode = file.beammode;

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
        hl2_seltree -> SetBranchAddress("true_Q2", &q2_true);
        hl2_seltree -> SetBranchAddress("truelepton_pdg", &pdg_true);

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

            /*
            double emu_true = std::sqrt(selmu_mom_true * selmu_mom_true + mu_mass * mu_mass);
            q2_true = 2.0 * enu_true * (emu_true - selmu_mom_true * selmu_cos_true)
                - mu_mass * mu_mass;

            double emu_reco = std::sqrt(selmu_mom * selmu_mom + mu_mass * mu_mass);
            q2_reco = 2.0 * enu_reco * (emu_reco - selmu_mom * selmu_cos)
                - mu_mass * mu_mass;
            */

            q2_reco = q2_true;
            double q0_true = 0;
            double q3_true = 0;
            if(std::abs(pdg_true) == 14)
                q0_true = enu_true - selmu_mom_true;
            else
                q0_true = enu_true - std::sqrt(selmu_mom_true * selmu_mom_true + mu_mass * mu_mass);

            q3_true = std::sqrt(q2_true + q0_true*q0_true);

            if(do_apply_weights)
            {
                if(reaction == 0)
                    weight = weight * GetBeRPAWeight(q2_true / 1.0E6);
                //if(reaction == 1)
                //    weight = weight * GetCCResQ2Weight(q2_true / 1.0E6);
                //weight = weight * GetTestEnuWeight(enu_true);
                //if(topology == 0 || topology == 1 || topology == 2)
                //    weight *= 0.80;
                //if(topology == 0 || topology == 1 || topology == 2)
                //    weight = weight * GetCCZeroPiQ2Weight(q2_true / 1.0E6);
                //weight = weight * GetFluxWeightND5(enu_true);
            }

            if(do_apply_templates)
                weight *= nd280_templates(q0_true, q3_true, reaction);

            /*
            if(weight == 0)
            {
                std::cout << "w : " << weight << std::endl;
                std::cout << "r : " << reaction << std::endl;
                std::cout << "q0: " << q0_true << std::endl;
                std::cout << "q3: " << q3_true << std::endl;
            }
            */

            weight *= file.pot_norm;

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

        hl2_trutree -> SetBranchAddress("true_Q2", &q2_true);
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

            /*
            double emu_true = std::sqrt(selmu_mom_true * selmu_mom_true + mu_mass * mu_mass);
            q2_true = 2.0 * enu_true * (emu_true - selmu_mom_true * selmu_cos_true)
                - mu_mass * mu_mass;

            double q0_true = enu_true - emu_true;
            double q3_true = std::sqrt(q2_true + q0_true*q0_true);
            */

            double q0_true = 0;
            double q3_true = 0;
            if(std::abs(pdg_true) == 14)
                q0_true = enu_true - selmu_mom_true;
            else
                q0_true = enu_true - std::sqrt(selmu_mom_true * selmu_mom_true + mu_mass * mu_mass);

            q3_true = std::sqrt(q2_true + q0_true*q0_true);

            if(do_apply_weights)
            {
                if(reaction_true == 0)
                    weight_true = weight_true * GetBeRPAWeight(q2_true / 1.0E6);
                //if(reaction_true == 1)
                //    weight_true = weight_true * GetCCResQ2Weight(q2_true / 1.0E6);
                //weight_true = weight_true * GetTestEnuWeight(enu_true);
                //if(topology_true == 0 || topology_true == 1 || topology_true == 2)
                //    weight_true *= 0.80;
                //if(topology_true == 0 || topology_true == 1 || topology_true == 2)
                //    weight_true = weight_true * GetCCZeroPiQ2Weight(q2_true / 1.0E6);
                //weight_true = weight_true * GetFluxVarWeight(enu_true);
                //weight_true = weight_true * GetFluxWeightND5(enu_true);
            }

            if(do_apply_templates)
                weight_true *= nd280_templates(q0_true, q3_true, reaction_true);

            weight_true *= file.pot_norm;
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
            //D2True = angle_true;

            if(pmu_true == 0 && angle_true == 0)
            {
                D1True = 0;
                D2True = 0;
            }

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

            double q0_true = enu_true - emu_true;
            double q3_true = std::sqrt(q2_true + q0_true*q0_true);

            cut_branch = file.branch;
            weight = event_weight;

            if(do_apply_weights)
            {
                if(reaction == 0)
                    weight = weight * GetBeRPAWeight(q2_true / 1.0E6);
                //if(reaction == 1)
                //    weight = weight * GetCCResQ2Weight(q2_true / 1.0E6);
                //weight = weight * GetTestEnuWeight(enu_true);
                //if(topology == 0 || topology == 1 || topology == 2)
                //    weight *= 1.20;
                //if(topology == 0 || topology == 1 || topology == 2)
                //    weight = weight * GetCCZeroPiQ2Weight(q2_true / 1.0E6);
                //weight = weight * GetFluxWeightND2(enu_true);
            }

            if(do_apply_templates && reaction >= 0)
                weight *= ingrid_templates(q0_true, q3_true, reaction);

            if(weight == 0)
                weight = event_weight;

            /*
            if(weight == 0)
            {
                std::cout << "w : " << weight << std::endl;
                std::cout << "r : " << reaction << std::endl;
                std::cout << "q0: " << q0_true << std::endl;
                std::cout << "q3: " << q3_true << std::endl;
            }
            */

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
    else if(code > 30 && code < 53)
        reaction = 4;
    else if(code < 0)
        reaction = 5;
    else if(code == 17 || code == 22 || code == 23)
        reaction = 8;

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

double GetBeRPAWeight(double q2)
{
    const double A = 0.59;
    const double B = 1.05;
    const double D = 1.13;
    const double E = 0.88;
    const double U = 1.20;

    double weight = 1;
    if(q2 < U)
    {
        const double xp = q2 / U;
        const double mxp = 1.0 - xp;
        const double C = D + U * E * (D-1) / 3.0;

        weight = A * mxp * mxp * mxp
               + 3 * B * mxp * mxp * xp
               + 3 * C * mxp * xp * xp
               + D * xp * xp * xp;
    }
    else
        weight = 1 + (D-1) * std::exp(-E * (q2 - U));

    return weight;
}

double GetCCResQ2Weight(double q2)
{
    double weight = 1.0;

    if(q2 < 0.7)
        weight = 1.01 / (1 + std::exp(1 - (std::sqrt(q2) / 0.156)));

    return weight;
}

double GetCCZeroPiQ2Weight(double q2)
{
    std::vector<double> min_q2_bins = {0, 0.00625, 0.0125, 0.025, 0.0375, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 2, 4, 100};
    std::vector<double> min_q2_vals = {0.735455, 0.778532, 0.952489, 1.05089, 1.0566, 1.05761, 1.01586, 1.03872, 1.04451, 1.0605, 1.06971, 1.06523, 1.08958, 1.08486, 1.05823, 0.789291, 0.451755};

    int bin = -1;
    for(int i = 0; i < min_q2_bins.size(); ++i)
    {
        if(q2 < min_q2_bins.at(i))
        {
           bin = i-1;
           break;
        }
    }

    //std::cout << "Q2 : " << q2 << std::endl;
    //std::cout << "bin: " << bin << std::endl;
    //std::cout << "val: " << min_q2_vals[bin] << std::endl;

    return bin >= 0 ? min_q2_vals.at(bin) : 1.0;
}

double GetFluxWeightND5(double enu)
{
    std::vector<double> enu_bins = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0,
                                    1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 30.0};
    //std::vector<double> enu_vals = {1.01015, 1.006, 1.00655, 1.00926, 1.0102, 1.02011, 1.02012, 1.02127, 1.03202, 0.994345, 0.995292, 0.971514, 1.00326, 1.00371, 1.01274, 1.00728, 1.00433, 1.02093, 1.01379, 1.0};
    std::vector<double> enu_vals = {1.00829, 1.00595, 0.9982, 1.00807, 1.00048, 1.00893, 1.0062, 1.00645, 1.00593, 0.99711, 1.02508, 1.0165, 1.01261, 1.00304, 1.01192, 1.00151, 1.00258, 1.00053, 1.0116, 1.0};

    int bin = -1;
    double enu_gev = enu / 1000.0;
    for(int i = 0; i < enu_bins.size(); ++i)
    {
        if(enu_gev < enu_bins.at(i))
        {
           bin = i-1;
           break;
        }
    }

    //std::cout << "enu: " << enu << std::endl;
    //std::cout << "bin: " << bin << std::endl;
    //std::cout << "val: " << enu_vals[bin] << std::endl;

    return bin >= 0 ? 1.0 + 3 * (enu_vals.at(bin) - 1.0) : 1.0;
}

double GetFluxWeightND2(double enu)
{
    std::vector<double> enu_bins = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0,
                                    1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 30.0};
    //std::vector<double> enu_vals = {1.00741, 1.00906, 1.01141, 1.01169, 1.01883, 1.0151, 1.0113, 1.01153, 1.01486, 1.01775, 1.0179, 1.03152, 1.03926, 1.05449, 1.07131, 1.07828, 1.04941, 1.01617, 1.01962, 1.0};
    std::vector<double> enu_vals = {1.00962, 1.0062, 1.00624, 1.00513, 1.00824, 1.00179, 1.00215, 1.00755, 1.00688, 1.00393, 1.00644, 1.00959, 1.00436, 1.00538, 1.01224, 1.00602, 1.01365, 1.01655, 0.9999, 1.0};

    int bin = -1;
    double enu_gev = enu / 1000.0;
    for(int i = 0; i < enu_bins.size(); ++i)
    {
        if(enu_gev < enu_bins.at(i))
        {
           bin = i-1;
           break;
        }
    }

    return bin >= 0 ? 1.0 + 3 * (enu_vals.at(bin) - 1.0) : 1.0;
}
