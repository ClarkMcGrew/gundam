//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  !!!!!  When running error propagation code : option -t to save toys !!!!!
//  Usage: root -b -q 'DrawXsecToys.C+("fit3_statFluc", 10000)'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


void DrawXsecToys(string inputname = "fit3_statFluc", const int Ntoys = 10000)
{

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());

	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store final output in a TFile =====" << std::endl;

	TFile* fin = new TFile(infilename.c_str());
	//======================================================================================================


	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with cross-section result =====" << std::endl;

	vector< TH1D* > h_sel_signal;
	vector< TH1D* > h_tru_signal;
	vector< TH1D* > h_eff_signal;
	vector< TH1D* > h_param;
	
	for(int itoy=0; itoy<Ntoys; itoy++)
	{
		h_sel_signal.push_back( (TH1D*)(fin -> Get( Form("sel_signal_toy%d", itoy) )) );
		h_tru_signal.push_back( (TH1D*)(fin -> Get( Form("tru_signal_toy%d", itoy) )) );
		h_eff_signal.push_back( (TH1D*)(fin -> Get( Form("eff_combined_toy%d", itoy) )) );
		h_param.push_back(      (TH1D*)(fin -> Get( Form("param_toy%d", itoy) )) );
	}

	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;

	// //======================================================================================================
	// std::cout << "===== Draw distribution for each toy separately =====" << std::endl;

	// TCanvas* c_sel_sign[Ntoys];
	// for(int itoy=0; itoy<Ntoys; itoy++)
	// {
	// 	c_sel_sign[itoy] = new TCanvas(Form("sel_signal_toy_%d", itoy), Form("sel_signal_toy_%d", itoy), 800, 600);
	
	// 	h_sel_signal[itoy] -> SetLineColor(kBlue+itoy*4.0/Ntoys);
	// 	// h_tru_signal[itoy] -> SetLineColor(kRed+itoy*4.0/Ntoys);

	// 	// h_tru_signal[itoy]->Draw("hist");
	// 	h_sel_signal[itoy]->Draw("hist");
	// 	c_sel_sign[itoy] -> Print( Form("plots/xsecResults/errpropToys/singleToys/sel_signal_toy%d.pdf", itoy) );
	// 	delete c_sel_sign[itoy];
	// }

	// TCanvas* c_tru_sign[Ntoys];
	// for(int itoy=0; itoy<Ntoys; itoy++)
	// {
	// 	c_tru_sign[itoy] = new TCanvas(Form("tru_signal_toy_%d", itoy), Form("tru_signal_toy_%d", itoy), 800, 600);
	
	// 	h_tru_signal[itoy] -> SetLineColor(kRed+itoy*4.0/Ntoys);

	// 	h_tru_signal[itoy]->Draw("hist");
	// 	c_tru_sign[itoy] -> Print( Form("plots/xsecResults/errpropToys/singleToys/tru_signal_toy%d.pdf", itoy) );
	// 	delete c_tru_sign[itoy];
	// }
	// //======================================================================================================



	//======================================================================================================
	std::cout << "===== Draw distribution for each toy on same plot =====" << std::endl;

	TCanvas* c_sel_sign_all = new TCanvas("sel_signal_all", "sel_signal_all",800,600);
	double max = 0;
	double min = 999;
	for(int itoy=0; itoy<Ntoys; itoy++)
	{
		max = MAX(max, h_sel_signal[itoy] -> GetMaximum());
		min = MIN(min, h_sel_signal[itoy] -> GetMinimum());
		h_sel_signal[itoy] -> Draw("same hist");

		if(h_sel_signal[itoy] -> GetMaximum() > 1E-10)
			std::cout << "WARNING : toys " << itoy << " has maximum value = " << h_sel_signal[itoy] -> GetMaximum() << std::endl;
	}

	h_sel_signal[0] -> GetYaxis() -> SetRangeUser(0.9*min, 1.1*max);
	c_sel_sign_all -> Print( Form("plots/xsecResults/errpropToys/sel_signal_toys_all.pdf") );
	gPad->SetLogy();
	c_sel_sign_all -> Print( Form("plots/xsecResults/errpropToys/sel_signal_toys_all_log.pdf") );
	
	//======================================================================================================

	TCanvas* c_tru_sign_all = new TCanvas("tru_signal_all", "tru_signal_all",800,600);
	max = 0;
	min = 999;
	for(int itoy=0; itoy<Ntoys; itoy++)
	{
		max = MAX(max, h_tru_signal[itoy] -> GetMaximum());
		min = MIN(min, h_tru_signal[itoy] -> GetMinimum());
		h_tru_signal[itoy] -> Draw("same hist");
	}

	h_tru_signal[0] -> GetYaxis() -> SetRangeUser(0.9*min, 1.1*max);
	c_tru_sign_all  -> Print( Form("plots/xsecResults/errpropToys/tru_signal_toys_all.pdf") );
	gPad->SetLogy();
	c_tru_sign_all -> Print( Form("plots/xsecResults/errpropToys/tru_signal_toys_all_log.pdf") );
	
	//======================================================================================================

	TCanvas* c_eff_sign_all = new TCanvas("eff_signal_all", "eff_signal_all",800,600);
	max = 0;
	min = 999;
	for(int itoy=0; itoy<Ntoys; itoy++)
	{
		max = MAX(max, h_eff_signal[itoy] -> GetMaximum());
		min = MIN(min, h_eff_signal[itoy] -> GetMinimum());
		h_eff_signal[itoy] -> Draw("same hist");
	}

	h_eff_signal[0] -> GetYaxis() -> SetRangeUser(0.9*min, 1.1*max);
	c_eff_sign_all -> Print( Form("plots/xsecResults/errpropToys/eff_signal_toys_all.pdf") );
	gPad->SetLogy();
	c_eff_sign_all -> Print( Form("plots/xsecResults/errpropToys/eff_signal_toys_all_log.pdf") );
	
	//======================================================================================================






	// //======================================================================================================
	// std::cout << "===== Draw parameters for each toy separately =====" << std::endl;

	// TCanvas* c_param[Ntoys];
	// for(int itoy=0; itoy<Ntoys; itoy++)
	// {
	// 	c_param[itoy] = new TCanvas(Form("param_toy_%d", itoy), Form("param_toy_%d", itoy), 800, 600);
		
	// 	h_param[itoy] -> SetLineColor(kBlue+itoy*4.0/Ntoys);

	// 	h_param[itoy]->Draw("hist");
	// 	c_param[itoy] -> Print( Form("plots/xsecResults/errpropToys/singleToys/param_toy%d.pdf", itoy) );
	// 	delete c_param[itoy];
	// }
	// // ======================================================================================================


	//======================================================================================================
	std::cout << "===== Draw parameters for each toy on same plot =====" << std::endl;

	TCanvas* c_param_all = new TCanvas("param_all", "param_all",800,600);
	max = 0;
	for(int itoy=0; itoy<Ntoys; itoy++)
	{
		max = MAX(max, h_param[itoy] -> GetMaximum());
		h_param[itoy] -> Draw("same hist");
	}

	h_param[0] -> GetYaxis() -> SetRangeUser(1E-6, 1.1*max);
	c_param_all -> Print( Form("plots/xsecResults/errpropToys/param_toys_all.pdf") );
	gPad->SetLogy();
	c_param_all -> Print( Form("plots/xsecResults/errpropToys/param_toys_all_log.pdf") );
	
	//======================================================================================================




	//======================================================================================================
	// delete fin;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}



