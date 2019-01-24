/*******************************************
* Author : Lucie Maret
* mail : lucie.maret@cern.ch
*******************************************/

#include <TFile.h>
// #include "/atlas/users/lmaret/h2_v2r21/highland2/highlandTools/v2r17p1/src/SetT2KStyle.H"


// Compute the statistical fluctuations in every bin
void stat_fluct()
{
	plotStyle();

	int ntoys = 100;
	int nbins = 58;


	// string GenFilename = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toys/fit1_onlyFitPar_toy";
	string GenFilename = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toys/fit1_AllPar_toy";


	Double_t fluct = 0;
	Double_t Mean, Mean_err, Sigma, Sigma_err;
	TH1D *h_fluct = new TH1D("h_fluct",  "Statistical fluctuations of parameters", nbins, 0, nbins);
	TH1D *h_par   = new TH1D("h_par",    "Parameter values", nbins, 0, nbins);

	// Loop over the bins
	for(int j = 1; j <= nbins; j++)
	{
		stat_fluct_of_one_bin(j, ntoys, GenFilename, Mean, Mean_err, Sigma, Sigma_err);
		h_fluct -> SetBinContent(j, Sigma);
		h_fluct -> SetBinError(j, Sigma_err);
		h_par   -> SetBinContent(j, Mean);
		h_par   -> SetBinError(j, Mean_err);
	}


	TCanvas *c_fluct = new TCanvas("c_fluct","Fluctuation of template parameters",1000,500);
	c_fluct -> SetGrid();
	h_fluct -> SetStats(kFALSE);
	h_fluct -> SetTitle("; Template parameters; Statistical fluctuation");

	h_fluct -> SetMarkerColor(kRed);
	h_fluct -> SetMarkerStyle(kFullCircle);
	h_fluct -> SetFillColor(kRed-9);
	h_fluct -> SetFillStyle(3144);

	h_fluct -> Draw("P E2");

	// legend = new TLegend(0.63, 0.76, 0.87, 0.9);
	// legend->SetFillColor(0);
	// legend->AddEntry(h_fluct_O,"FGD1 and FGD2","l");
	// legend->Draw();

	c_fluct -> Print("histos/statistical_fluctuations.eps");


	// ofstream outfile;
	TFile *outfile = new TFile("output_stat_fluct.root", "RECREATE");
	if (outfile->IsOpen() ) printf("Output file opened successfully\n");

	h_fluct -> Write();

	outfile->Print();

}



// Compute the statistical fluctuation in one bin
void stat_fluct_of_one_bin(int bin, int Ntoys, string genFilename, Double_t &mean, Double_t &mean_err, Double_t &sigma, Double_t &sigma_err, bool draw_distr = 0)
{
	string filename;
	TFile *fin;
	TH1D *h_paramhist_toyi;
	TH1D *h_param_values_in_bin = new TH1D("h_param_values_in_bin", "Distribution of the parameters in bin i", 30, -0.5, 2.5);

	// Loop over the toys to fill the histogram with the distribution of the parameter in bin "bin"
	for(int i = 1; i <= Ntoys; i++)
	{

		// Define the input file
		filename = Form("%s%d.root", genFilename.c_str(), i);
		TFile *fin = new TFile(filename.c_str());
		// std::cout << "File " << filename.c_str() << " is open." << std::endl;

		// if(h_paramhist_toyi != NULL)
		// {
		//     h_paramhist_toyi->Reset
		//     std::cout << "Histogram h_paramhist_toyi for i = " << i << " has been reset." << std::endl;
		// }

		h_paramhist_toyi = (TH1D*)fin->Get("hist_par_fit_result");
		h_param_values_in_bin -> Fill( h_paramhist_toyi -> GetBinContent(bin) );
		// std::cout << "Content of bin " << bin << " for toy i = " << i << " has been added." << std::endl;

		fin->Close();
	}

	// Do a Gaussian fit for the parameter distribution !!!!!!!!!! needs more than 100 toys for Gaussian fit!!!
	/////////////////////////////// 100 toys -> Use RMS of distribution.
	
	// Get mean and sigma values with Gaussian fit
	// h_param_values_in_bin->Fit("gaus");
	// TF1 *fitfct = h_param_values_in_bin->GetFunction("gaus");  
	// Double_t mean      = fitfct->GetParameter(1);
	// Double_t mean_err  = fitfct->GetParError(1);
	// Double_t sigma     = fitfct->GetParameter(2);
	// Double_t sigma_err = fitfct->GetParError(2);

	// Get mean and RMS values from histogram
	         mean      = h_param_values_in_bin -> GetMean();
	Double_t mean_err  = h_param_values_in_bin -> GetMeanError();
	         sigma     = h_param_values_in_bin -> GetRMS();
	Double_t sigma_err = h_param_values_in_bin -> GetRMSError();


	if(draw_distr)
	{// Draw parameter distribution
		TCanvas *c_param_distribution = new TCanvas("c_param_distribution", Form("Parameter distribution in bin %d", bin),900,600);
		h_param_values_in_bin -> SetTitle( Form("Parameter distribution in bin %d; Parameter; ", bin) );
		h_param_values_in_bin -> Draw("hist");
		fitfct -> Draw("same");
		c_param_distribution -> Print( Form("histos/param_distribution_in_bin_%d.eps", bin) );
	}

	std::cout 	<< "Statistical flucutation has been computed for bin " << bin
				<< ", " << std::endl
				<< " *** mean = " << mean << " +/- " << mean_err << std::endl
				<< " *** sigma = " << sigma << " +/- " << sigma_err << std::endl;

	delete h_param_values_in_bin;
}



// Set plot style
void plotStyle()
{
	// T2K style
	// -- WhichStyle -- // 1 = presentation large fonts // 2 = presentation small fonts // 3 = publication/paper 
	// TStyle* t2kstyle = SetT2KStyle(3, "T2K");
	// gROOT->SetStyle(t2kstyle->GetName());

	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);

	gStyle->SetHistLineWidth( Width_t(2.5) );

	std::cout << "Set plot style." << std::endl;
}
