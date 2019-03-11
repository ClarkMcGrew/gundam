//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root -b -q 'DrawXsec.C+()'
//         root -b -q 'DrawXsec.C+("fit3_statFluc")'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

void DrawXsec(string inputname = "fit3_statFluc", const std::string& dir_name = "fakedata/statFluc", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());

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

	TFile* fin = new TFile(infilename.c_str());
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get pre-fit and post-fit chi2 =====" << std::endl;
	
	TH1D* chi2_tot = (TH1D*)(fin->Get("chi2_tot_periter"));
	Double_t chi2_tot_prefit  = chi2_tot -> GetBinContent(1);
	Double_t chi2_tot_postfit = chi2_tot -> GetBinContent(chi2_tot->GetNbinsX()-1);
	//======================================================================================================


	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with cross-section result and truth =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec_truth(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit_err(Ntarget);


	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_truth[itar].push_back(       (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d_truth", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_postfit[itar].push_back(     (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d",       targetlist[itar].c_str(), nbcth) )) );
			h_xsec_postfit_err[itar].push_back( (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d",       targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
		}

	vector< TH1D* > h_ratio_truth;
	vector< TH1D* > h_ratio_postfit;
	vector< TH1D* > h_ratio_postfit_err;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_truth.push_back(       (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d_truth", nbcth) )) );
		h_ratio_postfit.push_back(     (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_postfit_err.push_back( (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_errs", nbcth)) );
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set histos with errors =====" << std::endl;
	
	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_postfit_err[itar][nbcth] -> Reset();
			for(int i=1; i<=h_xsec_postfit_err[itar][nbcth]->GetNbinsX(); i++)
			{
				double err  = h_xsec_postfit[itar][nbcth] -> GetBinError(i);
				double mean = h_xsec_postfit[itar][nbcth] -> GetBinContent(i);
				if(err/mean<100)
					h_xsec_postfit_err[itar][nbcth] -> SetBinContent(i, err/mean);
				else
					h_xsec_postfit_err[itar][nbcth] -> SetBinContent(i, 0.0);
			}
		}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit_err[nbcth] -> Reset();
		for(int i=1; i<=h_ratio_postfit_err[nbcth]->GetNbinsX(); i++)
		{
			double err  = h_ratio_postfit[nbcth] -> GetBinError(i);
			double mean = h_ratio_postfit[nbcth] -> GetBinContent(i);
			if(err/mean<100)
				h_ratio_postfit_err[nbcth] -> SetBinContent(i, err/mean);
			else
				h_ratio_postfit_err[nbcth] -> SetBinContent(i, 0.0);
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
			h_xsec_truth[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_truth[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec_truth[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec_truth[itar][nbcth]   -> SetLineColor(kBlue+1);
			h_xsec_postfit[itar][nbcth] -> SetLineColor(kRed+1);

			c_xsec[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			h_xsec_truth[itar][nbcth]   -> Draw("hist"); // with error bars
			h_xsec_postfit[itar][nbcth] -> Draw("same E"); // with error bars

		}
		c_xsec[itar]->Print(Form("plots/xsecResults/%s/xsec_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));


		c_xsec_err[itar] = new TCanvas(Form("Xsec_error_on_%s", targetlist[itar].c_str()),Form("Xsec_error_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec_err[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_postfit_err[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_postfit_err[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec_postfit_err[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec_postfit_err[itar][nbcth] -> SetLineColor(kBlue);

			c_xsec_err[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			h_xsec_postfit_err[itar][nbcth] -> Draw("hist"); // with error bars

		}
		c_xsec_err[itar]->Print(Form("plots/xsecResults/%s/xsec_error_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));
	}


	TCanvas* c_ratio = new TCanvas("Xsec_ratio", "Xsec_ratios",1700,1000);
	c_ratio -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_truth[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_truth[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_truth[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");
		h_ratio_truth[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 2.0);
		
		h_ratio_truth[nbcth] -> SetLineColor(kBlue+1);
		h_ratio_postfit[nbcth] -> SetLineColor(kRed+1);

		c_ratio -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio_truth[nbcth]   -> Draw("hist"); // with error bars
		h_ratio_postfit[nbcth] -> Draw("same E"); // with error bars
	}
	c_ratio->Print( Form("plots/xsecResults/%s/xsec_OCratio_%s.pdf", dir_name.c_str(), inputname.c_str()) );



	TCanvas* c_ratio_err = new TCanvas("Xsec_ratio", "Xsec_ratio",1700,1000);
	c_ratio_err -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit_err[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_postfit_err[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_postfit_err[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio_postfit_err[nbcth] -> SetLineColor(kRed+1);

		c_ratio_err -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio_postfit_err[nbcth]->Draw("hist"); // with error bars
	}
	c_ratio_err->Print( Form("plots/xsecResults/%s/xsec_error_OCratio_%s.pdf", dir_name.c_str(), inputname.c_str()) );

	//======================================================================================================






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
	
	c_eff->Print( Form("plots/xsecResults/%s/efficiency_%s.pdf", dir_name.c_str(), inputname.c_str()) );
	//======================================================================================================




	//======================================================================================================
	delete fin;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}
