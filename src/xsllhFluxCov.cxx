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
    TAxis *axis_bins = (TAxis*)in_file -> Get(input_bin_name.c_str());
    const unsigned int nbins = axis_bins -> GetNbins();
    
    std::cout << TAG << "Number of bins :" << nbins << std::endl;

    // Get bin edges and store them into an array
    TArrayD root_array = *(axis_bins -> GetXbins());
    double* bin_array = root_array.GetArray();
    
    // for(int i = 1; i <= nbins; i++)
    //     std::cout << "bin_array["<<i<<"] = " << bin_array[i] << std::endl;

    TH1D *h_binning = new TH1D(output_bin_name.c_str(), output_bin_name.c_str(), nbins, bin_array);


    // Get flux covariance matrix from file
    TMatrixDSym *flux_cov_tot = (TMatrixDSym*)in_file -> Get(input_mat_name.c_str());
    TMatrixDSym flux_cov(nbins);

    for(int i=0; i<nbins; i++)
        for(int j=0; j<nbins; j++)
            flux_cov(i,j) = (*flux_cov_tot)(i,j);
    
    std::cout << TAG << "Total covariance matrix dimension :" << flux_cov_tot -> GetNcols() <<"x" << flux_cov_tot -> GetNcols() << std::endl;
    std::cout << TAG << "Covariance submatrix dimension :" << flux_cov.GetNcols() <<"x" << flux_cov.GetNcols() << std::endl;



    // Write matrix and binning into file
    out_file -> cd();
    flux_cov.Write(output_mat_name.c_str());
    h_binning -> Write(output_bin_name.c_str());

    out_file -> Close();
    in_file  -> Close();

    std::cout << TAG << "Finished." << std::endl;
    return 0;
}