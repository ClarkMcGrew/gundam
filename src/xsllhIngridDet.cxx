#include <algorithm>
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

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TStyle.h"
#include "TTree.h"

#include "json.hpp"
using json = nlohmann::json;

#include "BinManager.hh"
#include "ColorOutput.hh"
#include "ProgressBar.hh"
#include "ToyThrower.hh"

struct Event
{
    double pmu_true;
    double angle_true;
    double dist_reco;
    double angle_reco;
    double weight;
    int true_bin;
    int reco_bin;
};

int main(int argc, char** argv)
{
    const std::string TAG = color::GREEN_STR + "[xsIngridDet]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + color::BOLD_STR + "[ERROR]: " + color::RESET_STR;

    std::cout << "--------------------------------------------------------\n"
              << TAG << color::RainbowText("Welcome to the Super-xsLLh Detector Variation Interface.\n")
              << TAG << color::RainbowText("Initializing the variation machinery...") << std::endl;

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Reading Events "));

    ProgressBar pbar_toys(60, "#");
    pbar_toys.SetRainbow();
    pbar_toys.SetPrefix(std::string(TAG + "Throwing Toys "));

    std::string json_file;
    char option;
    while((option = getopt(argc, argv, "j:h")) != -1)
    {
        switch(option)
        {
            case 'j':
                json_file = optarg;
                break;
            case 'h':
                std::cout << "USAGE: " << argv[0] << "\nOPTIONS:\n"
                          << "-j : JSON input\n";
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

    json js;
    f >> js;

    std::string fname_input = js["fname_input"];
    std::string fname_output  = js["fname_output"];
    std::string cov_file_name = js["covariance_file"];
    std::string cov_mat_name  = js["covariance_name"];
    std::string true_binning = js["true_binning"];
    std::string reco_binning = js["reco_binning"];

    bool do_mc_stat = js["mc_stat_error"];
    unsigned int ntoys = js["num_toys"];
    double cov_offset = js["cov_offset"];

    std::cout << TAG << "Input event file: " << fname_input << std::endl
              << TAG << "Covariance file : " << cov_file_name << std::endl
              << TAG << "Covariance name : " << cov_mat_name << std::endl
              << TAG << "True binning file: " << true_binning << std::endl
              << TAG << "Reco binning file: " << reco_binning << std::endl
              << TAG << "Apply MC stat error: " << std::boolalpha << do_mc_stat << std::endl;
    std::cout << TAG << "Output ROOT file: " << fname_output << std::endl;

    TFile* ing_file = TFile::Open(fname_input.c_str(), "READ");
    TTree* ing_tree = (TTree*)ing_file -> Get("wtree");

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

    std::vector<Event> v_events;
    BinManager true_bm(true_binning);
    BinManager reco_bm(reco_binning);

    //true_bm.Print();
    //reco_bm.Print();

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

        if(selected_sample == 0)
        {
            Event e;
            e.pmu_true = pmu_true * 1000;
            e.angle_true = angle_true;
            //e.angle_true = TMath::Cos(angle_true * deg_to_rad);
            e.dist_reco = pmu_reco;
            e.angle_reco = angle_reco;
            e.weight = event_weight;
            e.true_bin = true_bm.GetBinIndex(std::vector<double>{e.angle_true, e.pmu_true});
            e.reco_bin = reco_bm.GetBinIndex(std::vector<double>{e.angle_reco, e.dist_reco});
            v_events.push_back(e);

            /*
            std::cout << "-----------" << std::endl;
            std::cout << "pmu_true : " << e.pmu_true << std::endl;
            std::cout << "pmu_reco : " << e.dist_reco << std::endl;
            std::cout << "angle_true : " << e.angle_true << std::endl;
            std::cout << "angle_reco : " << e.angle_reco << std::endl;
            std::cout << "true bin : " << e.true_bin << std::endl;
            std::cout << "reco bin : " << e.reco_bin << std::endl;
            */
        }

        if(i % 2000 == 0 || i == (nevents-1))
            pbar.Print(i, nevents-1);
    }

    ing_file->Close();

    std::cout << TAG << "Finished reading events." << std::endl;

    TFile* cov_file = TFile::Open(cov_file_name.c_str(), "READ");
    TH2D* cov_hist = nullptr;
    cov_file->GetObject(cov_mat_name.c_str(), cov_hist);

    const unsigned int nbins_pmu = 9;
    const unsigned int nbins_angle = 6;
    const unsigned int nbins_hist = nbins_pmu*nbins_angle;
    const unsigned int nbins_true = true_bm.GetNbins(); //(nbins_pmu-3) * (nbins_angle-2);
    const unsigned int nbins_reco = reco_bm.GetNbins();
    TMatrixTSym<double> cov_true(nbins_true);
    cov_true.Zero();
    TMatrixTSym<double> cov_reco(nbins_reco);
    cov_reco.Zero();
    TMatrixTSym<double> cor_reco(nbins_reco);
    cor_reco.Zero();

    std::vector<int> remove_bins = {5,6,11,12,17,18,23,24,29,30,
                                    31,32,33,34,35,36,37,38,39,40,
                                    41,42,43,44,45,46,47,48,49,50,
                                    51,52,53,54};

    std::vector<bool> row_mask_bool(nbins_hist, false);
    for(const auto& row : remove_bins)
        row_mask_bool.at(row-1) = true;

    std::cout << TAG << "Reading original covariance." << std::endl;

    unsigned int i_r{0}, j_r{0};
    for(int i = 0; i < nbins_hist; ++i)
    {
        if(!row_mask_bool[i])
        {
            for(int j = 0; j < nbins_hist; ++j)
            {
                if(!row_mask_bool[j])
                {
                    cov_true(i_r,j_r) = cov_hist->GetBinContent(i+1, j+1);
                    j_r++;
                }
            }
            i_r++;
            j_r = 0;
        }
    }

    for(int i = 0; i < nbins_true; ++i)
        cov_true(i,i) += cov_offset;

    std::cout << TAG << "Throwing toys..." << std::endl;

    std::vector<TH1D> h_toys;
    h_toys.reserve(ntoys);

    const unsigned int seed = 12257;
    ToyThrower toy_thrower(cov_true, seed, true, 1.0E-48);
    for(int i = 0; i < ntoys; ++i)
    {
        std::vector<double> toy(nbins_true, 0.0);
        toy_thrower.Throw(toy);

        //std::cout << "-------------" << std::endl;
        //for(const auto& val : toy)
        //    std::cout << val << std::endl;

        std::string temp = "toy" + std::to_string(i);
        TH1D temp_toy(temp.c_str(), temp.c_str(), nbins_reco, 0, nbins_reco);

        for(auto& event : v_events)
        {
            double reweight = event.weight * (1.0 + toy.at(event.true_bin));
            temp_toy.Fill(event.reco_bin + 0.5, reweight);
        }

        h_toys.emplace_back(temp_toy);

        if(i % 100 == 0 || i == (ntoys-1))
            pbar_toys.Print(i, ntoys-1);
    }

    TH1D mc_stat("mc_stat", "mc_stat", nbins_reco, 0, nbins_reco);
    for(auto& event : v_events)
        mc_stat.Fill(event.reco_bin + 0.5, event.weight);

    TH1D mc_stat_error("mc_stat_error", "mc_stat_error", nbins_reco, 0, nbins_reco);
    double* w = mc_stat.GetArray();
    double* w2 = mc_stat.GetSumw2()->GetArray();
    for(unsigned int i = 0; i < nbins_reco; ++i)
    {
        double rel_error = w2[i+1] / (w[i+1] * w[i+1]);
        if(std::isnan(rel_error))
            rel_error = 0;
        mc_stat_error.SetBinContent(i+1, rel_error);
    }

    std::cout << TAG << "Calculating mean and covariance." << std::endl;

    TH1D h_mean("", "", nbins_reco, 0, nbins_reco);
    for(const auto& hist : h_toys)
    {
        for(int i = 0; i < nbins_reco; ++i)
            h_mean.Fill(i + 0.5, hist.GetBinContent(i + 1));
    }
    h_mean.Scale(1.0 / (1.0 * ntoys));

    for(const auto& hist : h_toys)
    {
        for(int i = 0; i < nbins_reco; ++i)
        {
            for(int j = 0; j < nbins_reco; ++j)
            {
                const double x = 1.0 - hist.GetBinContent(i + 1) / h_mean.GetBinContent(i + 1);
                const double y = 1.0 - hist.GetBinContent(j + 1) / h_mean.GetBinContent(j + 1);
                cov_reco(i, j) += x * y / (1.0 * ntoys);
            }
        }
    }

    for(int i = 0; i < nbins_reco; ++i)
    {
        for(int j = 0; j < nbins_reco; ++j)
        {
            if(cov_reco(i,i) == 0 || std::isnan(cov_reco(i,i)))
                cov_reco(i,i) = 0.005;

            if(std::isnan(cov_reco(i,j)))
                cov_reco(i,j) = 0;

            cov_reco(i,i) += 1E-6;
        }
    }

    if(do_mc_stat)
    {
        for(int i = 0; i < nbins_reco; ++i)
            cov_reco(i,i) += mc_stat_error.GetBinContent(i+1);
    }

    for(int i = 0; i < nbins_reco; ++i)
    {
        for(int j = 0; j < nbins_reco; ++j)
        {
            const double x = cov_reco(i, i);
            const double y = cov_reco(j, j);
            const double z = cov_reco(i, j);
            cor_reco(i, j) = z / (std::sqrt(x * y));

            if(std::isnan(cov_reco(i, j)))
                cor_reco(i, j) = 0.0;
        }
    }

    std::cout << TAG << "Saving output." << std::endl;
    TFile* foutput = TFile::Open(fname_output.c_str(), "RECREATE");
    foutput->cd();

    h_mean.Write("h_mean");
    cov_true.Write("cov_true");
    cov_reco.Write("cov_reco");
    cor_reco.Write("cor_reco");

    if(do_mc_stat)
        mc_stat_error.Write("mc_stat_error");

    foutput->Close();

    std::cout << TAG << "Finished." << std::endl;
    std::cout << TAG << "\u3042\u308a\u304c\u3068\u3046\u3054\u3056\u3044\u307e\u3057\u305f\uff01"
              << std::endl;
    return 0;
}
