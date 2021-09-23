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
#include <TAxis.h>
#include <TH1D.h>
#include "TMatrixT.h"
#include "TMatrixDSym.h"

#include "json.hpp"
using json = nlohmann::json;

#include "ColorOutput.hh"
#include "ProgressBar.hh"


int main(int argc, char** argv)
{
    const std::string TAG = color::GREEN_STR + "[xsFluxCov]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + color::BOLD_STR
                            + "[ERROR]: " + color::RESET_STR;

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Reading Events "));

    std::cout << "-------------------------------------------------\n"
              << TAG << "Welcome to the Super-xsLLh Flux Covariance converter.\n"
              << TAG << "Initializing the tree machinery..." << std::endl;

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
                std::cout << "USAGE: "
                          << argv[0] << "\nOPTIONS:\n"
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

    json j;
    f >> j;

    std::cout << TAG << "Reading configuration options..." << std::endl;
    std::string in_fname  = j["input"]["fname"];
    std::string out_fname = j["output"]["fname"];
    // std::string out_seltree_name = j["output"]["sel_tree"];

    std::cout << TAG << "Input File: "  << in_fname << std::endl
              << TAG << "Output File: " << out_fname << std::endl;
    
    TFile* out_file = TFile::Open(out_fname.c_str(), "RECREATE");
    TFile* in_file  = TFile::Open(in_fname.c_str(),  "READ");
    
    if(!in_file)
    {
        std::cerr << ERR << "Failed to open " << in_fname << std::endl;
        return 1;
    }

    std::string input_mat_name  = j["input"]["matrix"];
    std::string output_mat_name = j["output"]["matrix"];
    std::string input_bin_name  = j["input"]["binning"];
    std::string output_bin_name = j["output"]["binning"];


    // Get binning from file
    //TAxis *axis_bins = (TAxis*)in_file -> Get(input_bin_name.c_str());
    //const unsigned int nbins = axis_bins -> GetNbins();
    const unsigned int nbins = 20;
    
    std::cout << TAG << "Number of bins : " << nbins << std::endl;

    // Get bin edges and store them into an array
    //TArrayD root_array = *(axis_bins -> GetXbins());
    //double* bin_array = root_array.GetArray();
    double bin_array[nbins];
    bin_array[0] = 0;
    bin_array[1] = 100;
    bin_array[2] = 200;
    bin_array[3] = 300;
    bin_array[4] = 400;
    bin_array[5] = 500;
    bin_array[6] = 600;
    bin_array[7] = 700;
    bin_array[8] = 800;
    bin_array[9] = 1000;
    bin_array[10] = 1200;
    bin_array[11] = 1500;
    bin_array[12] = 2000;
    bin_array[13] = 2500;
    bin_array[14] = 3000;
    bin_array[15] = 3500;
    bin_array[16] = 4000;
    bin_array[17] = 5000;
    bin_array[18] = 7000;
    bin_array[19] = 10000;
    bin_array[20] = 30000;
    
    for(int i = 0; i <= nbins; ++i)
        std::cout << "bin_array["<<i<<"] = " << bin_array[i] << std::endl;

    TH1D *h_binning = new TH1D(output_bin_name.c_str(), output_bin_name.c_str(), nbins, bin_array);

    // ND280 numu:
    h_binning->SetBinContent(1, 8.74355324e+10);
    h_binning->SetBinContent(2, 3.45058698e+11);
    h_binning->SetBinContent(3, 8.28720198e+11);
    h_binning->SetBinContent(4, 1.63632727e+12);
    h_binning->SetBinContent(5, 2.65118681e+12);
    h_binning->SetBinContent(6, 3.27433332e+12);
    h_binning->SetBinContent(7, 2.91284864e+12);
    h_binning->SetBinContent(8, 2.09249181e+12);
    h_binning->SetBinContent(9, 1.81540257e+12);
    h_binning->SetBinContent(10, 7.02215683e+11);
    h_binning->SetBinContent(11, 5.46322846e+11);
    h_binning->SetBinContent(12, 4.76502402e+11);
    h_binning->SetBinContent(13, 2.74465251e+11);
    h_binning->SetBinContent(14, 1.90793969e+11);
    h_binning->SetBinContent(15, 1.55897804e+11);
    h_binning->SetBinContent(16, 1.31488483e+11);
    h_binning->SetBinContent(17, 1.92737176e+11);
    h_binning->SetBinContent(18, 1.62935550e+11);
    h_binning->SetBinContent(19, 6.81140262e+10);
    h_binning->SetBinContent(20, 6.81140262e+10);
    
    // Get flux covariance matrix from file
    TMatrixDSym *flux_cov_tot = (TMatrixDSym*)in_file -> Get(input_mat_name.c_str());

    // ND280 Numu only:
    TMatrixDSym flux_cov(20);
    for(int i=160; i<180; i++)
        for(int j=160; j<180; j++)
            flux_cov(i-160,j-160) = (*flux_cov_tot)(i,j);

    std::cout << TAG << "Total covariance matrix dimension : " << flux_cov_tot -> GetNcols() <<"x" << flux_cov_tot -> GetNcols() << std::endl;
    std::cout << TAG << "Covariance submatrix dimension : " << flux_cov.GetNcols() <<"x" << flux_cov.GetNcols() << std::endl;



    // Write matrix and binning into file
    out_file -> cd();
    flux_cov.Write(output_mat_name.c_str());
    h_binning -> Write(output_bin_name.c_str());

    out_file -> Close();
    in_file  -> Close();

    std::cout << TAG << "Finished." << std::endl;
    return 0;
}
