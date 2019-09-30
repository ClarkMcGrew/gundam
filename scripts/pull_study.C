/*******************************************
* Author : Lucie Maret
* mail : lucie.maret@cern.ch
*******************************************/

/*

root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit1_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit1_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  26, "fit1_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1122, "fit1_statFluc")'

root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit3_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit3_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  26, "fit3_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1122, "fit3_statFluc")'

root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit3_asimov_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit3_asimov_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  26, "fit3_asimov_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1122, "fit3_asimov_statFluc")'

*/

#include "CommonStyle.h"
#include "CommonHeader.h"


// Compute the pull for each parameter
void pull_study(string par_name = "fit", string par_title = "Template parameters", int Nparam = 116, string fit_type = "fit3_statFluc")
{
	
	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gStyle->SetOptStat(1);
	gROOT->ForceStyle();


	std::cout << "------------------------" << std::endl;
	std::cout << "----- Begin script -----" << std::endl;

	// fit1 
	int toys_fit1[]   = {1, 10, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 11, 110, 111, 112, 113, 115, 116, 117, 118, 119, 12, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 13, 130, 131, 132, 133, 134, 135, 136, 137, 138, 14, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 15, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 16, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 17, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 18, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 19, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 2, 20, 200, 21, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 4, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 5, 50, 51, 52, 53, 54, 55, 56, 57, 58, 6, 60, 61, 62, 63, 64, 65, 66, 68, 69, 7, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 8, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99 };
	
	// fit3 
	int toys_fit3[] = {1, 10, 11, 12, 15, 16, 18, 19, 2, 20, 21, 22, 23, 25, 26, 27, 28, 29, 3, 30, 31, 32, 34, 35, 36, 37, 38, 39, 4, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 5, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 6, 60, 61, 62, 63, 64, 67, 69, 7, 70, 71, 72, 73, 74, 76, 77, 8, 84, 86, 9 };
	
	int ntoys;
	if(fit_type == "fit1_statFluc")             ntoys = (sizeof(toys_fit1)/sizeof(*toys_fit1));
	else if(fit_type == "fit3_statFluc")        ntoys = (sizeof(toys_fit3)/sizeof(*toys_fit3));
	
	std::cout << "----- Number of toys used = "<<ntoys<<" -----" << std::endl;

	int nbins = Nparam;

	string GenFilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toysStatFluc/%s_toy", fit_type.c_str());


	// std::cout << "----- Plot the pull distribution of each toy -----" << std::endl;

	//for(int toyi = 1; toyi <= ntoys; toyi++)
	//{
	//	pull_of_one_toy(toyi, nbins, GenFilename);
	//}

	//std::cout << "----- Plot the parameter values of each toy -----" << std::endl;

	//parameters_of_each_toy(ntoys, nbins, GenFilename);

	int binMin, binMax;
	if(fit_type == "fit1_statFluc")      { binMin = -5.0; binMax = 5.0; }
	else if(fit_type == "fit3_statFluc") { binMin = -1.0; binMax = 1.0; }

	std::vector< Double_t > pull_and_sigma_in_bin_j;
	TH1D *h_pull_mean   = new TH1D("h_pull_mean", "Template parameter pull", nbins, 0, nbins);
	TH1D *h_pull_sigma  = new TH1D("h_pull_sigma", "Pull width", nbins, 0, nbins);
	TH1D *h_pull_O      = new TH1D("h_pull_O", "Pull distribution for oxygen parameters", 40, binMin, binMax);
	TH1D *h_pull_C      = new TH1D("h_pull_C", "Pull distribution for carbon parameters", 40, binMin, binMax);
	TH1D *h_pull        = new TH1D("h_pull",   "Pull distribution for parameters",        40, binMin, binMax);


	std::cout << "----- Loop over the parameters to compute the mean among toys -----" << std::endl;

	for(int j = 1; j <= nbins; j++)
	{
		if(fit_type == "fit1_statFluc")
			pull_and_sigma_in_bin_j = pull_of_one_bin(par_name, j, toys_fit1, ntoys, GenFilename, false); // function defined below, returns the pull+error and its sigma+error
		else if(fit_type == "fit3_statFluc")
			pull_and_sigma_in_bin_j = pull_of_one_bin(par_name, j, toys_fit3, ntoys, GenFilename, false); // function defined below, returns the pull+error and its sigma+error

		h_pull_mean   -> SetBinContent(j, pull_and_sigma_in_bin_j[0]);
		h_pull_mean   -> SetBinError(j, pull_and_sigma_in_bin_j[1]);
		h_pull_sigma  -> SetBinContent(j, pull_and_sigma_in_bin_j[2]);
		h_pull_sigma  -> SetBinError(j, pull_and_sigma_in_bin_j[3]);

		if(par_name == "fit")
		{
			if(j>nbins/2)
				h_pull_O -> Fill(pull_and_sigma_in_bin_j[0]);
			else
				h_pull_C -> Fill(pull_and_sigma_in_bin_j[0]);
		}
		else
			h_pull -> Fill(pull_and_sigma_in_bin_j[0]);
	}


	std::cout << "----- Do a Gaussian fit -----" << std::endl;

	if(par_name == "fit")
	{
		h_pull_O -> Fit("gaus");
		h_pull_C -> Fit("gaus");
	}
	else
		h_pull -> Fit("gaus");


	std::cout << "----- Plot results -----" << std::endl;

	TCanvas *c_pull_mean = new TCanvas("c_pull_mean","Parameter pull",1000,600);
	h_pull_mean -> SetTitle(Form(";%s; ", par_title.c_str()));
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

	TLine *line0 = new TLine(0,0,nbins,0);
	TLine *line1 = new TLine(0,1,nbins,1);
	TLine *line2 = new TLine(nbins/2,-4,nbins/2,5);

	h_pull_mean   -> Draw("P E2");
	h_pull_sigma  -> Draw("same P E2");
	line0 -> Draw("same");
	line1 -> Draw("same");
	if(par_name=="fit") line2 -> Draw("same");

	legend = new TLegend(0.6,0.7,0.85,0.9);
	legend->SetFillColor(0);
	legend->AddEntry(h_pull_mean,"Pull mean","p");
	legend->AddEntry(h_pull_sigma,"Pull sigma","p");
	legend->Draw(); 

	c_pull_mean -> Print(Form("plots/pullstudy/%s_pull_mean_and_sigma_%s.pdf", par_name.c_str(), fit_type.c_str()));


	if(par_name == "fit")
	{

		TCanvas *c_pull_O = new TCanvas("c_pull_O","Pull distribution for oxygen parameters",900,700);
		h_pull_O -> SetTitle("; Pull for oxygen parameters; ");
		c_pull_O -> SetGrid();
		h_pull_O -> GetYaxis() -> SetRangeUser(0.0, 1.3 * h_pull_O -> GetMaximum());

		h_pull_O -> Draw();

		c_pull_O -> Print(Form("plots/pullstudy/%s_pull_distribution_O_%s.pdf", par_name.c_str(), fit_type.c_str()));


		TCanvas *c_pull_C = new TCanvas("c_pull_C","Pull distribution for carbon parameters",900,700);
		h_pull_C -> SetTitle(" ; Pull for carbon parameters; ");
		c_pull_C -> SetGrid();
		h_pull_C -> GetYaxis() -> SetRangeUser(0.0, 1.3 * h_pull_C -> GetMaximum());

		h_pull_C -> Draw();

		c_pull_C -> Print(Form("plots/pullstudy/%s_pull_distribution_C_%s.pdf", par_name.c_str(), fit_type.c_str()));
	}
	else
	{

		TCanvas *c_pull = new TCanvas("c_pull","Pull distribution for parameters",900,700);
		h_pull -> SetTitle(Form("; Pull for %s parameters; ", par_name.c_str()));
		c_pull -> SetGrid();
		h_pull -> GetYaxis() -> SetRangeUser(0.0, 1.3 * h_pull -> GetMaximum());

		h_pull -> Draw();

		c_pull -> Print(Form("plots/pullstudy/%s_pull_distribution_%s.pdf", par_name.c_str(), fit_type.c_str()));

	}


	std::cout << "----- End script -----" << std::endl;
	std::cout << "----------------------" << std::endl;

}








// Compute the pull mean for each bin
std::vector< Double_t > pull_of_one_bin(string par_name, int bin, int Toys[], int Ntoys, string genFilename, bool draw_distr = false)
{
	string filename;
	TFile *fin;
	TH1D *h_parampriorhist_toyi;
	TH1D *h_paramhist_toyi;
	TH1D *h_paramerrhist_toyi;
	TH1D *h_pull_values_in_bin = new TH1D("h_pull_values_in_bin", "Distribution of the pull in bin i", 40, -10, 10);
	TH1D *h_param_values_in_bin = new TH1D("h_param_values_in_bin", "Distribution of the c_i", 20, 0, 2);
	TH1D *h_paramerr_values_in_bin = new TH1D("h_paramerr_values_in_bin", "Distribution of the error on c_i", 60, 0, 0.3);

	// Loop over the toys to fill the histogram with the distribution of the pull in bin "bin"
	for(int i = 0; i < Ntoys; i++)
	{
		// Define the input file
		filename = Form("%s%d.root", genFilename.c_str(), Toys[i]);
		TFile *fin = new TFile(filename.c_str());
		// std::cout << "File " << filename.c_str() << " is open." << std::endl;

		// Get the c_i and delta_c_i, then fill the pull distribution
		h_parampriorhist_toyi = (TH1D*)fin->Get(Form("hist_par_%s_prior", par_name.c_str()));
		h_paramhist_toyi      = (TH1D*)fin->Get(Form("hist_par_%s_result", par_name.c_str()));
		h_paramerrhist_toyi   = (TH1D*)fin->Get(Form("hist_par_%s_error_final", par_name.c_str()));
		
		h_pull_values_in_bin     -> Fill( (h_parampriorhist_toyi -> GetBinContent(bin) - h_paramhist_toyi -> GetBinContent(bin))/ (h_paramerrhist_toyi -> GetBinContent(bin)) );
		h_param_values_in_bin    -> Fill( h_paramhist_toyi    -> GetBinContent(bin) );
		h_paramerr_values_in_bin -> Fill( h_paramerrhist_toyi -> GetBinContent(bin) );

		fin->Close();
		delete fin;
	}

	// Get the values of mean and RMS
	Double_t mean      = h_pull_values_in_bin -> GetMean();
	Double_t mean_err  = h_pull_values_in_bin -> GetMeanError();
	Double_t sigma     = h_pull_values_in_bin -> GetRMS();
	Double_t sigma_err = h_pull_values_in_bin -> GetRMSError();

	std::vector< Double_t > pull_and_sigma;
	pull_and_sigma.push_back(mean);
	pull_and_sigma.push_back(mean_err);
	pull_and_sigma.push_back(sigma);
	pull_and_sigma.push_back(sigma_err);

	if(draw_distr)
	{ // Draw pull (optional)
		TCanvas *c_pull_distribution = new TCanvas("c_pull_distribution", Form("Pull distribution of parameter %d", bin-1),900,600);
		h_pull_values_in_bin -> SetTitle( Form("Pull distribution of parameter %d; Pull; ", bin-1) );
		h_pull_values_in_bin -> Draw("hist");
		// fitfct -> Draw("same");
		c_pull_distribution -> Print( Form("plots/pullstudy/pull_distribution_of_param_%d.pdf", bin-1) );

		TCanvas *c_param_distribution = new TCanvas("c_param_distribution", Form("Parameter c_%d distribution", bin-1),900,600);
		h_param_values_in_bin -> SetTitle( Form("Parameter c_%d distribution; Parameter value; ", bin-1) );
		h_param_values_in_bin -> Draw("hist");
		// fitfct -> Draw("same");
		c_param_distribution -> Print( Form("plots/pullstudy/param_distribution_%d.pdf", bin-1) );

		TCanvas *c_paramerr_distribution = new TCanvas("c_paramerr_distribution", Form("Error on parameter c_%d distribution", bin-1),900,600);
		h_paramerr_values_in_bin -> SetTitle( Form("Error on parameter c_%d distribution; Parameter error; ", bin-1) );
		h_paramerr_values_in_bin -> Draw("hist");
		// fitfct -> Draw("same");
		c_paramerr_distribution -> Print( Form("plots/pullstudy/paramerr_distribution_%d.pdf", bin-1) );

		// sleep(1);
	}

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
	h_paramhist    = (TH1D*)fin->Get(Form("hist_par_%s_result", par_name.c_str()));
	h_paramerrhist = (TH1D*)fin->Get(Form("hist_par_%s_error_final", par_name.c_str()));
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
	c_pull_distribution -> Print( Form("plots/pullstudy/pull_distribution_for_toy_%d.pdf", toy) );

	sleep(1);

	std::cout 	<< "Pull has been computed for toy " << toy
				<< ", " << std::endl
				<< " *** mean = " << mean << " +/- " << mean_err << std::endl
				<< " *** RMS = " << sigma << " +/- " << sigma_err << std::endl;
}




void parameters_of_each_toy(int Toys[], int Ntoys, int const Nbins, string genFilename)
{
	string filename;
	TFile *fin;
	TH1D *h_param_toyi;


	// Write the parameter histogram of each toy into one root file
	for(int toyi = 1; toyi <= Ntoys; toyi++)
	{
		// Define the input file
		filename = Form("%s%d.root", genFilename.c_str(), Toys[toyi]);
		TFile *fin = new TFile(filename.c_str());
		if (fin->IsOpen() ) printf("Input file opened successfully\n");

		h_param_toyi = (TH1D*)fin->Get(Form("hist_par_%s_result", par_name.c_str()));
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
	for(int i = 1; i <= Ntoys; i++)
	{
		(TH1D*)param_file->Get(Form("param_for_toy_%d", Toys[toyi])) -> Draw("hist same");	
	}
	c_parameters -> Print("plots/pullstudy/parameters_for_each_toy.pdf");

}


