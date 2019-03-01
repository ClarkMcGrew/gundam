/*******************************************
* Author : Lucie Maret
* mail : lucie.maret@cern.ch
*******************************************/
#include <unistd.h>
#include <TFile.h>
#include "CommonStyle.h"


// Compute the pull for each parameter
void pull_study()
{
	
	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();


	std::cout << "------------------------" << std::endl;
	std::cout << "----- Begin script -----" << std::endl;

	plotStyle();

	int ntoys = 100;
	int nbins = 58;

	// string GenFilename = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toys/fit1_onlyFitPar_toy";
	string GenFilename = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toys/fit1_AllPar_toy";


	std::cout << "----- Plot the pull distribution of each toy -----" << std::endl;

	//for(int toyi = 1; toyi <= ntoys; toyi++)
	//{
	//	pull_of_one_toy(toyi, nbins, GenFilename);
	//}

	//std::cout << "----- Plot the parameter values of each toy -----" << std::endl;

	//parameters_of_each_toy(ntoys, nbins, GenFilename);


	std::vector< Double_t > pull_and_sigma_in_bin_j;
	TH1D *h_pull_mean   = new TH1D("h_pull_mean", "Template parameter pull", nbins, 0, nbins);
	TH1D *h_pull_sigma  = new TH1D("h_pull_sigma", "Pull width", nbins, 0, nbins);
	TH1D *h_pull_O        = new TH1D("h_pull_O", "Pull distribution for oxygen parameters", 20, -2.5, 2.5);
	TH1D *h_pull_C        = new TH1D("h_pull_C", "Pull distribution for carbon parameters", 20, -2.5, 2.5);


	std::cout << "----- Loop over the parameters to compute the mean among toys -----" << std::endl;

	for(int j = 1; j <= nbins; j++)
	{
		pull_and_sigma_in_bin_j = pull_of_one_bin(j, ntoys, GenFilename, false); // function defined below, returns the pull+error and its sigma+error

		h_pull_mean   -> SetBinContent(j, pull_and_sigma_in_bin_j[0]);
		h_pull_mean   -> SetBinError(j, pull_and_sigma_in_bin_j[1]);
		h_pull_sigma  -> SetBinContent(j, pull_and_sigma_in_bin_j[2]);
		h_pull_sigma  -> SetBinError(j, pull_and_sigma_in_bin_j[3]);

		if(j>nbins/2)
		{
			// h_pull_O -> Fill(pull_and_sigma_in_bin_j[0] / pull_and_sigma_in_bin_j[1]);
			h_pull_O -> Fill(pull_and_sigma_in_bin_j[0]);
		}
		else
		{
			// h_pull_C -> Fill(pull_and_sigma_in_bin_j[0] / pull_and_sigma_in_bin_j[1]);
			h_pull_C -> Fill(pull_and_sigma_in_bin_j[0]);
		}

	}


	std::cout << "----- Do a Gaussian fit -----" << std::endl;

	h_pull_O -> Fit("gaus");
	h_pull_C -> Fit("gaus");


	std::cout << "----- Plot results -----" << std::endl;

	TCanvas *c_pull_mean = new TCanvas("c_pull_mean","Parameter pull",900,500);
	h_pull_mean -> SetTitle("Pull; Parameter number; ");
	c_pull_mean -> SetGrid();
	h_pull_mean -> SetStats(kFALSE);
	h_pull_mean -> GetYaxis() -> SetRangeUser(-4, 5);

	h_pull_mean  -> SetMarkerColor(kBlue);
	h_pull_mean  -> SetMarkerStyle(kFullCircle);
	h_pull_sigma -> SetMarkerColor(kRed);
	h_pull_sigma -> SetMarkerStyle(kFullCircle);

	h_pull_mean  -> SetFillColor(kBlue-9);
	h_pull_mean  -> SetFillStyle(1001);
	h_pull_sigma -> SetFillColor(kRed-9);
	h_pull_sigma -> SetFillStyle(3144);

	TLine *line0 = new TLine(0,0,58,0);
	TLine *line1 = new TLine(0,1,58,1);
	TLine *line2 = new TLine(29,-4,29,5);

	h_pull_mean   -> Draw("P E2");
	h_pull_sigma  -> Draw("same P E2");
	line0 -> Draw("same");
	line1 -> Draw("same");
	line2 -> Draw("same");

	legend = new TLegend(0.6,0.7,0.85,0.9);
	legend->SetFillColor(0);
	legend->AddEntry(h_pull_mean,"Pull mean","p");
	legend->AddEntry(h_pull_sigma,"Pull sigma","p");
	legend->Draw(); 

	c_pull_mean -> Print("histos/pullstudy/pull_mean_and_sigma.pdf");


	TCanvas *c_pull_O = new TCanvas("c_pull_O","Pull distribution for oxygen parameters",800,500);
	h_pull_O -> SetTitle("; Pull distribution for oxygen parameters; ");
	c_pull_O -> SetGrid();

	h_pull_O -> Draw();


	c_pull_O -> Print("histos/pullstudy/pull_distribution_O.pdf");

	TCanvas *c_pull_C = new TCanvas("c_pull_C","Pull distribution for carbon parameters",800,500);
	h_pull_C -> SetTitle(" ; Pull distribution for carbon parameters; ");
	c_pull_C -> SetGrid();

	h_pull_C -> Draw();

	c_pull_C -> Print("histos/pullstudy/pull_distribution_C.pdf");


	std::cout << "----- End script -----" << std::endl;
	std::cout << "----------------------" << std::endl;

}








// Compute the pull mean for each bin
std::vector< Double_t > pull_of_one_bin(int bin, int Ntoys, string genFilename, bool draw_distr = false)
{
	string filename;
	TFile *fin;
	TH1D *h_paramhist_toyi;
	TH1D *h_paramerrhist_toyi;
	TH1D *h_pull_values_in_bin = new TH1D("h_pull_values_in_bin", "Distribution of the pull in bin i", 30, -10, 10);
	TH1D *h_param_values_in_bin = new TH1D("h_param_values_in_bin", "Distribution of the c_i", 20, 0, 2);
	TH1D *h_paramerr_values_in_bin = new TH1D("h_paramerr_values_in_bin", "Distribution of the error on c_i", 60, 0, 0.3);

	// Loop over the toys to fill the histogram with the distribution of the pull in bin "bin"
	for(int i = 1; i <= Ntoys; i++)
	{
		// Define the input file
		filename = Form("%s%d.root", genFilename.c_str(), i);
		TFile *fin = new TFile(filename.c_str());
		// std::cout << "File " << filename.c_str() << " is open." << std::endl;

		// Get the c_i and delta_c_i, then fill the pull distribution
		h_paramhist_toyi    = (TH1D*)fin->Get("hist_par_fit_result");
		h_paramerrhist_toyi = (TH1D*)fin->Get("hist_par_fit_error_final");
		
		h_pull_values_in_bin     -> Fill( (1 - h_paramhist_toyi -> GetBinContent(bin))/ (h_paramerrhist_toyi -> GetBinContent(bin)) );
		h_param_values_in_bin    -> Fill( h_paramhist_toyi    -> GetBinContent(bin) );
		h_paramerr_values_in_bin -> Fill( h_paramerrhist_toyi -> GetBinContent(bin) );

		fin->Close();
	}

	// Get the values of mean and RMS
	Double_t mean      = h_pull_values_in_bin -> GetMean();
	Double_t mean_err  = h_pull_values_in_bin -> GetMeanError();
	Double_t sigma     = h_pull_values_in_bin -> GetRMS();
	Double_t sigma_err = h_pull_values_in_bin -> GetRMSError();

	if(draw_distr)
	{ // Draw pull (optional)
		TCanvas *c_pull_distribution = new TCanvas("c_pull_distribution", Form("Pull distribution of parameter %d", bin-1),900,600);
		h_pull_values_in_bin -> SetTitle( Form("Pull distribution of parameter %d; Pull; ", bin-1) );
		h_pull_values_in_bin -> Draw("hist");
		// fitfct -> Draw("same");
		c_pull_distribution -> Print( Form("histos/pullstudy/pull_distribution_of_param_%d.pdf", bin-1) );

		TCanvas *c_param_distribution = new TCanvas("c_param_distribution", Form("Parameter c_%d distribution", bin-1),900,600);
		h_param_values_in_bin -> SetTitle( Form("Parameter c_%d distribution; Parameter value; ", bin-1) );
		h_param_values_in_bin -> Draw("hist");
		// fitfct -> Draw("same");
		c_param_distribution -> Print( Form("histos/pullstudy/param_distribution_%d.pdf", bin-1) );

		TCanvas *c_paramerr_distribution = new TCanvas("c_paramerr_distribution", Form("Error on parameter c_%d distribution", bin-1),900,600);
		h_paramerr_values_in_bin -> SetTitle( Form("Error on parameter c_%d distribution; Parameter error; ", bin-1) );
		h_paramerr_values_in_bin -> Draw("hist");
		// fitfct -> Draw("same");
		c_paramerr_distribution -> Print( Form("histos/pullstudy/paramerr_distribution_%d.pdf", bin-1) );

		// sleep(1);
	}

	std::vector< Double_t > pull_and_sigma;
	pull_and_sigma.push_back(mean);
	pull_and_sigma.push_back(mean_err);
	pull_and_sigma.push_back(sigma);
	pull_and_sigma.push_back(sigma_err);

	std::cout 	<< "Pull has been computed for bin " << bin
				<< ", " << std::endl
				<< " *** mean = " << mean << " +/- " << mean_err << std::endl
				<< " *** RMS = " << sigma << " +/- " << sigma_err << std::endl;

	delete h_pull_values_in_bin;
	delete h_param_values_in_bin;
	delete h_paramerr_values_in_bin;
	return pull_and_sigma;
}







// Plot the pull distribution for one toy
void pull_of_one_toy(int toy, int Nbins, string genFilename)
{
	TH1D *h_paramhist;
	TH1D *h_paramerrhist;
	TH1D *h_pull_values = new TH1D("h_pull_values", "Distribution of the pull in bin i", 30, -10, 10);

	// Define the input file
	string filename = Form("%s%d.root", genFilename.c_str(), toy);
	TFile *fin = new TFile(filename.c_str());
	
	// Get the c_i and delta_c_i, then fill the pull distribution
	h_paramhist    = (TH1D*)fin->Get("hist_par_fit_result");
	h_paramerrhist = (TH1D*)fin->Get("hist_par_fit_error_final");
	for(int i = 1; i<=Nbins; i++)
	{ 
		h_pull_values -> Fill( (1-h_paramhist -> GetBinContent(i)) / (h_paramerrhist -> GetBinContent(i)) );
	}

	fin->Close();
	
	// Get the values of mean and RMS
	Double_t mean      = h_pull_values -> GetMean();
	Double_t mean_err  = h_pull_values -> GetMeanError();
	Double_t sigma     = h_pull_values -> GetRMS();
	Double_t sigma_err = h_pull_values -> GetRMSError();

	// Draw pull
	TCanvas *c_pull_distribution = new TCanvas("c_pull_distribution", Form("Pull distribution for toy %d", toy),900,600);
	h_pull_values -> SetTitle( Form("Pull distribution for toy %d; Pull; ", toy) );
	h_pull_values -> Draw("hist");
	// fitfct -> Draw("same");
	c_pull_distribution -> Print( Form("histos/pullstudy/pull_distribution_for_toy_%d.pdf", toy) );

	sleep(1);

	std::cout 	<< "Pull has been computed for toy " << toy
				<< ", " << std::endl
				<< " *** mean = " << mean << " +/- " << mean_err << std::endl
				<< " *** RMS = " << sigma << " +/- " << sigma_err << std::endl;
}




void parameters_of_each_toy(int const Ntoys, int const Nbins, string genFilename)
{
	string filename;
	TFile *fin;
	TH1D *h_param_toyi;


	// Write the parameter histogram of each toy into one root file
	for(int toyi = 1; toyi <= Ntoys; toyi++)
	{
		// Define the input file
		filename = Form("%s%d.root", genFilename.c_str(), toyi);
		TFile *fin = new TFile(filename.c_str());
		if (fin->IsOpen() ) printf("Input file opened successfully\n");

		h_param_toyi = (TH1D*)fin->Get("hist_par_fit_result");
		// h_param_toyi -> SetTitle(Form("param_for_toy_%d", toyi));
		h_param_toyi -> SetLineColor(kRed+toyi/10);
		h_param_toyi -> SetDirectory(0);
		h_param_toyi -> GetYaxis() -> SetRangeUser(0.2, 1.8);

		fin->Close();

		TFile *outfile = new TFile("outputs/fitparam_for_each_toy.root", "UPDATE");
		if (outfile->IsOpen() ) printf("Output file opened successfully\n");
		
		h_param_toyi -> SetName(Form("param_for_toy_%d", toyi));
		h_param_toyi -> Write(Form("param_for_toy_%d", toyi), TObject::kWriteDelete);

		outfile->Print();
		outfile->Close();

		h_param_toyi -> Reset();

		std::cout << "End of toy i = " << toyi << std::endl;
	}
	
	TFile *param_file = new TFile("outputs/fitparam_for_each_toy.root", "READ");
	
	// Draw parameters
	TCanvas *c_parameters = new TCanvas("c_parameters", "Parameters for each toy",2000,1200);
	for(int toyi = 1; toyi <= Ntoys; toyi++)
	{
		(TH1D*)param_file->Get(Form("param_for_toy_%d", toyi)) -> Draw("hist same");	
	}
	c_parameters -> Print("histos/pullstudy/parameters_for_each_toy.pdf");

}


