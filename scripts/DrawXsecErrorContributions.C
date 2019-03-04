//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root 'DrawXsecErrorContributions.C+()'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

void DrawXsecErrorContributions(string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename_fit  = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyFitPar.root";
	string infilename_flux = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyFluxPar.root";
	string infilename_det  = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyDetPar.root";
	string infilename_xsec = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyXsecPar.root";

	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set binning using BinningTools class =====" << std::endl;
	
	BinningTools bin; 
	bin.SetBinning(fbinning.c_str());

	//  int Nbins = bin.GetNbins();//Total number of bins
	int Nbins_costh = bin.GetNbins_costh();//Number of costheta bins

	int* Nmombins = new int[Nbins_costh];//Number of mom bin in each slice of costheta
	Nmombins = bin.GetMomNBinsInEachCosthSlice();

	double** mombins = new double*[Nbins_costh];//2-d array with momentum binning for each slice of cosine theta
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
			mombins[nbcth] = new double[nbm];

	mombins = bin.GetMomBins();

	TString h_title[] = { "-1  < cos#theta < 0.2",
	                     "0.2  < cos#theta < 0.6",
	                     "0.6  < cos#theta < 0.7",
	                     "0.7  < cos#theta < 0.8",
	                     "0.8  < cos#theta < 0.85",
	                     "0.85 < cos#theta < 0.9",
	                     "0.9  < cos#theta < 0.94",
	                     "0.94 < cos#theta < 0.98",
	                     "0.98 < cos#theta < 1.0"
	                 };
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store final output in a TFile =====" << std::endl;

	TFile* fin_fit  = new TFile(infilename_fit.c_str());
	TFile* fin_flux = new TFile(infilename_flux.c_str());
	TFile* fin_det  = new TFile(infilename_det.c_str());
	TFile* fin_xsec = new TFile(infilename_xsec.c_str());
	//======================================================================================================










///////////////////////////////////////////////////////////

// FROM HERE IT NEES TO BE MODIFIED 

///////////////////////////////////////////////////////////














	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with cross-section result =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec(Ntarget);
	vector< vector<TH1D*> > h_xsec_err(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec[itar].push_back(     (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_err[itar].push_back( (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
		}

	vector< TH1D* > h_ratio;
	vector< TH1D* > h_ratio_err;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio.push_back(     (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_err.push_back( (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_errs", nbcth)) );
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set histos with errors =====" << std::endl;
	
	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_err[itar][nbcth] -> Reset();
			for(int i=1; i<=h_xsec_err[itar][nbcth]->GetNbinsX(); i++)
			{
				double err  = h_xsec[itar][nbcth] -> GetBinError(i);
				double mean = h_xsec[itar][nbcth] -> GetBinContent(i);
				if(err/mean<100)
					h_xsec_err[itar][nbcth] -> SetBinContent(i, err/mean);
				else
					h_xsec_err[itar][nbcth] -> SetBinContent(i, 0.0);
			}
		}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_err[nbcth] -> Reset();
		for(int i=1; i<=h_ratio_err[nbcth]->GetNbinsX(); i++)
		{
			double err  = h_ratio[nbcth] -> GetBinError(i);
			double mean = h_ratio[nbcth] -> GetBinContent(i);
			if(err/mean<100)
				h_ratio_err[nbcth] -> SetBinContent(i, err/mean);
			else
				h_ratio_err[nbcth] -> SetBinContent(i, 0.0);
		}
	}
	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	TCanvas* c_xsec[Ntarget];
	TCanvas* c_xsec_err[Ntarget];
	// vector<TLegend*> leg(Ntarget);
	// vector<TLegend*> leg_err(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
	{
		c_xsec[itar] = new TCanvas(Form("Xsec_on_%s", targetlist[itar].c_str()),Form("Xsec_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec[itar][nbcth] -> SetLineColor(kBlue);

			c_xsec[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			h_xsec[itar][nbcth]->Draw("E"); // with error bars

		}
		c_xsec[itar]->Print(Form("plots/xsecResults/xsec_%s_%s.pdf", inputname.c_str(), targetlist[itar].c_str()));


		c_xsec_err[itar] = new TCanvas(Form("Xsec_error_on_%s", targetlist[itar].c_str()),Form("Xsec_error_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec_err[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_err[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_err[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec_err[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec_err[itar][nbcth] -> SetLineColor(kBlue);

			c_xsec_err[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			h_xsec_err[itar][nbcth]->Draw("hist"); // with error bars

		}
		c_xsec_err[itar]->Print(Form("plots/xsecResults/xsec_error_%s_%s.pdf", inputname.c_str(), targetlist[itar].c_str()));
	}


	TCanvas* c_ratio = new TCanvas("Xsec_ratio", "Xsec_ratios",1700,1000);
	c_ratio -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio[nbcth] -> SetLineColor(kBlue);

		c_ratio -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio[nbcth]->Draw("E"); // with error bars
	}
	c_ratio->Print( Form("plots/xsecResults/xsec_OCratio_%s.pdf", inputname.c_str()) );



	TCanvas* c_ratio_err = new TCanvas("Xsec_ratio", "Xsec_ratio",1700,1000);
	c_ratio_err -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_err[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_err[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_err[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio_err[nbcth] -> SetLineColor(kBlue);

		c_ratio_err -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio_err[nbcth]->Draw("hist"); // with error bars
	}
	c_ratio_err->Print( Form("plots/xsecResults/xsec_error_OCratio_%s.pdf", inputname.c_str()) );






	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with efficiency and plot it =====" << std::endl;
	
	TH1D* h_eff = (TH1D*)( fin->Get("eff_best_fit") );

	TCanvas* c_eff = new TCanvas("Efficiency", "Efficiency",800,600);
	
	h_eff -> SetTitle("Selection efficiency");
	h_eff -> GetXaxis() -> SetTitle("true analysis bins");
	h_eff -> GetYaxis() -> SetTitle("true events / selected events");
	h_eff -> SetLineColor(kBlue);

	// if(nbcth!=0) gPad->SetLogx();

	h_eff->Draw("hist"); // with error bars
	
	c_eff->Print( Form("plots/xsecResults/efficiency_%s.pdf", inputname.c_str()) );

	
	//======================================================================================================





	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

	delete fin;
}


