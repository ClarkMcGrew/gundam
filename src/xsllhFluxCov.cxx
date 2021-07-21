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

    /*
    bin_array.push_back(0);
    bin_array.push_back(100);
    bin_array.push_back(200);
    bin_array.push_back(300);
    bin_array.push_back(400);
    bin_array.push_back(500);
    bin_array.push_back(600);
    bin_array.push_back(700);
    bin_array.push_back(800);
    bin_array.push_back(1000);
    bin_array.push_back(1200);
    bin_array.push_back(1500);
    bin_array.push_back(2000);
    bin_array.push_back(2500);
    bin_array.push_back(3000);
    bin_array.push_back(3500);
    bin_array.push_back(4000);
    bin_array.push_back(5000);
    bin_array.push_back(7000);
    bin_array.push_back(10000);
    bin_array.push_back(30000);
    */
    
    for(int i = 0; i <= nbins; ++i)
        std::cout << "bin_array["<<i<<"] = " << bin_array[i] << std::endl;

    TH1D *h_binning = new TH1D(output_bin_name.c_str(), output_bin_name.c_str(), nbins, bin_array);


    // Get flux covariance matrix from file
    TMatrixDSym *flux_cov_tot = (TMatrixDSym*)in_file -> Get(input_mat_name.c_str());



    
    // INGRID only:
    TMatrixDSym flux_cov(20);
    for(int i=0; i<20; i++)
        for(int j=0; j<20; j++)
            flux_cov(i,j) = (*flux_cov_tot)(i,j);
    
    /*
    // ND280 only:
    TMatrixDSym flux_cov(80);
    for(int i=160; i<240; i++)
        for(int j=160; j<240; j++)
            flux_cov(i-160,j-160) = (*flux_cov_tot)(i,j);
    */
    /*
    // ND280 + INGRID:
    TMatrixDSym flux_cov(100);
    for(int i=0; i<20; i++)
        for(int j=0; j<20; j++)
            flux_cov(i,j) = (*flux_cov_tot)(i,j);
    for(int i=160; i<240; i++)
        for(int j=160; j<240; j++)
            flux_cov(i-140,j-140) = (*flux_cov_tot)(i,j);
    */





    /*
    for(int i=80; i<120; i++)
        for(int j=80; j<120; j++)
            flux_cov(i-40,j-40) = (*flux_cov_tot)(i,j);
    
    for(int i=0; i<40; i++)
        for(int j=80; j<120; j++)
            flux_cov(i,j-40) = (*flux_cov_tot)(i,j);
    */
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
