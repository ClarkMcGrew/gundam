/*******************************************
* Author : Lucie Maret
* mail : lucie.maret@cern.ch
*******************************************/

// root -b -q 'regularisation_Lcurve.C+()'

#include "CommonHeader.h"
#include "CommonStyle.h"


// Compute the statistical fluctuations in every bin
void regularisation_Lcurve(string inputname = "fit3_statFluc", const std::string& dir_name = "fakedata/statFlucReg")
{
    //======================================================================================================  
    //=== Set common style
    CommonStyle();
    gROOT->ForceStyle();


	// if(inputname == "fit3_statFluc")
	// {

	// const int NTOYS = 8;
	// double      pReg[NTOYS]     = { 0.001,   0.01,   0.1,   0.3,   0.5,   1.0,   2.5,   5.0};
	// std::string pReg_str[NTOYS] = {"0.001", "0.01", "0.1", "0.3", "0.5", "1.0", "2.5", "5.0"};

	// }
	// else if(inputname == "fit1_asimov_statFluc")
	// {

	const int NTOYS = 9;
	double      pReg[NTOYS]     = { 0.01,   0.1,   0.5,   1.0,   1.5,   2.0,   2.5,   5.0,   10.0};
	std::string pReg_str[NTOYS] = {"0.01", "0.1", "0.5", "1.0", "1.5", "2.0", "2.5", "5.0", "10.0"};

	// }

	double chi2fit[NTOYS];
	double chi2reg[NTOYS];

	double chi2fit_best[1];
	double chi2reg_best[1];

	TFile *fin;
	TH1D *h_chi2_stat_periter;
	TH1D *h_chi2_sys_periter;
	TH1D *h_chi2_reg_periter;

	for(int i = 0; i < NTOYS; i++)
	{
		// Define the input file
		std::string filename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toysReg/%s_reg_toy%d.root", inputname.c_str(), i);
		cout << "Toy i = "<<i<<", open file " << filename << endl;
		fin = new TFile(filename.c_str());

		// Load histograms of the chi2
		h_chi2_stat_periter = (TH1D*)fin->Get("chi2_stat_periter");
		h_chi2_sys_periter  = (TH1D*)fin->Get("chi2_sys_periter");
		h_chi2_reg_periter  = (TH1D*)fin->Get("chi2_reg_periter");
		

		// Define the x and y of the plot
		chi2fit[i] =  h_chi2_stat_periter -> GetBinContent(h_chi2_stat_periter->GetNbinsX()-1);
		chi2fit[i] += h_chi2_sys_periter  -> GetBinContent(h_chi2_sys_periter->GetNbinsX()-1);
		chi2reg[i] = ( h_chi2_reg_periter -> GetBinContent(h_chi2_reg_periter->GetNbinsX()-1) ) / pReg[i];

		cout << "    chi2fit = " << chi2fit[i] << ", chi2reg = " << chi2reg[i] << endl;

		// if(i == i_best)
		// {
		// 	chi2fit_best[0] = chi2fit[i];
		// 	chi2reg_best[0] = chi2reg[i];
		// }

		fin->Close();
	}


	// Draw L-curve

	TCanvas *c_Lcurve = new TCanvas("c_fluct","L-curve", 1400, 900);

	TGraph* gr = new TGraph(NTOYS, chi2fit, chi2reg);
	TGraph* gr_best = new TGraph(2, chi2fit_best, chi2reg_best);
	gr->SetTitle("L-curve for p bins; chi2 stat + syst; chi2 reg / reg strength");
	gr_best->SetMarkerColor(2);

	int i=3;

	std::vector< TLatex* > label(NTOYS);

	for(int i=0; i<NTOYS; i++)
	{
		label[i] = new TLatex(gr->GetX()[i], gr->GetY()[i], pReg_str[i].c_str());
		gr->GetListOfFunctions()->Add(label[i]);
	}

	// gr->Draw("AC*");
	gr->Draw("alp*");

	c_Lcurve -> Print( Form("plots/fitteroutput/%s/Lcurve_%s.pdf", dir_name.c_str(), inputname.c_str()) );

}

