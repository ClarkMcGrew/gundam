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

void DrawXsec(string inputname = "fit3_statFluc", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
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
	std::cout << "===== Get histos with cross-section result =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec_postfit(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit_err(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_postfit[itar].push_back(     (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_postfit_err[itar].push_back( (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
		}

	vector< TH1D* > h_ratio_postfit;
	vector< TH1D* > h_ratio_postfit_err;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
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











	// //======================================================================================================
	// std::cout << "================================================" << std::endl;
	// std::cout << "===== Store pre-fit MC and data in a TFile =====" << std::endl;

	// TFile* fin_mc = new TFile("xsllh_nominal_noECal.root");
	// TFile* fin_dt = new TFile("xsllh_nominal_fakedata_noECal.root");

	// TTree* seltree_mc = (TTree*)fin_mc -> Get("selectedEvents");
	// TTree* seltree_dt = (TTree*)fin_dt -> Get("selectedEvents");

	// Int_t cut_branch_mc, cut_branch_dt;
	// seltree_mc -> Branch("cut_branch", &cut_branch_mc, "cut_branch/I");
	// seltree_dt -> Branch("cut_branch", &cut_branch_dt, "cut_branch/I");

	// Int_t fgdtarget_mc, fgdtarget_dt;
	// seltree_mc -> Branch("fgdtarget", &fgdtarget_mc, "fgdtarget/I");
	// seltree_dt -> Branch("fgdtarget", &fgdtarget_dt, "fgdtarget/I");

	// Float_t weight_mc, weight_dt;
	// seltree_mc -> Branch("weight", &weight_mc, "weight/F");
	// seltree_dt -> Branch("weight", &weight_dt, "weight/F");

	// Float_t D1True_mc, D1Reco_mc, D2True_mc, D2Reco_mc;
 //    seltree_mc -> Branch("D1True",    &D1True_mc,    "D1True/F");
 //    seltree_mc -> Branch("D1Reco",    &D1Reco_mc,    "D1Reco/F");
 //    seltree_mc -> Branch("D2True",    &D2True_mc,    "D2True/F");
 //    seltree_mc -> Branch("D2Reco",    &D2Reco_mc,    "D2Reco/F");
 //    Float_t D1Reco_dt, D2Reco_dt;
 //    seltree_dt -> Branch("D1Reco",    &D1Reco_dt,    "D1Reco/F");
 //    seltree_dt -> Branch("D2Reco",    &D2Reco_dt,    "D2Reco/F");
	// //======================================================================================================


	// //======================================================================================================
	// std::cout << "================================================" << std::endl;
	// std::cout << "===== Define empty histograms with p, theta binning =====" << std::endl;

	// for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	// {
	// 	h_tmp.push_back(new TH1D("","", Nmombins[nbcth], mombins[nbcth]));
	// 	// h_tmp[nbcth] -> SetTitle( Form("%s", h_title[nbcth].Data()) );
	// 	// h_tmp[nbcth] -> GetYaxis()->SetTitle("CC-0#pi events/100 MeV");
	// 	// h_tmp[nbcth] -> GetXaxis()->SetTitle("p^{#mu}_{true} [GeV/c]");
	// 	// h_tmp[nbcth] -> GetXaxis()->SetTitle("p^{#mu}_{reco} [GeV/c]");
	// }

	// vector< vector<TH1D*> > h_xsec_prefit(Ntarget);
	// vector< vector<TH1D*> > h_xsec_prefit_err(Ntarget);
	// vector< vector<TH1D*> > h_xsec_data(Ntarget);
	// vector< vector<TH1D*> > h_xsec_data_err(Ntarget);

	// for(int itar = 0; itar < Ntarget; itar++)
	// {
	// 	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	// 	{
	// 		h_xsec_prefit[itar].push_back(    (TH1D*)h_tmp[nbcth] -> Clone(Form("h_xsec_prefit_cth_%d",nbcth)));
	// 		h_xsec_prefit_err[itar].push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_xsec_prefit_err_cth_%d",nbcth)));
	// 		h_xsec_data[itar].push_back(      (TH1D*)h_tmp[nbcth] -> Clone(Form("h_xsec_data_cth_%d",nbcth)));
	// 		h_xsec_data_err[itar].push_back(  (TH1D*)h_tmp[nbcth] -> Clone(Form("h_xsec_data_err_cth_%d",nbcth)));
	// 	}
	// }
 //    //======================================================================================================


	// //======================================================================================================
	// std::cout << "================================================" << std::endl;
	// std::cout << "===== Fill pre-fit MC histograms =====" << std::endl;

	// for(int i = 0; i < seltree_mc -> GetEntries(); ++i)
 //    {
 //        seltree_mc -> GetEntry(i);

 //    }
 //    //======================================================================================================


	// //======================================================================================================
	// std::cout << "================================================" << std::endl;
	// std::cout << "===== Fill data histograms =====" << std::endl;

	// for(int i = 0; i < seltree_dt -> GetEntries(); ++i)
 //    {
 //        seltree_dt -> GetEntry(i);

 //    }
 //    //======================================================================================================







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
			h_xsec_postfit[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_postfit[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec_postfit[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec_postfit[itar][nbcth] -> SetLineColor(kBlue);

			c_xsec[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			h_xsec_postfit[itar][nbcth]->Draw("E"); // with error bars

		}
		c_xsec[itar]->Print(Form("plots/xsecResults/xsec_%s_%s.pdf", inputname.c_str(), targetlist[itar].c_str()));


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

			h_xsec_postfit_err[itar][nbcth]->Draw("hist"); // with error bars

		}
		c_xsec_err[itar]->Print(Form("plots/xsecResults/xsec_error_%s_%s.pdf", inputname.c_str(), targetlist[itar].c_str()));
	}


	TCanvas* c_ratio = new TCanvas("Xsec_ratio", "Xsec_ratios",1700,1000);
	c_ratio -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_postfit[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_postfit[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio_postfit[nbcth] -> SetLineColor(kBlue);

		c_ratio -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio_postfit[nbcth]->Draw("E"); // with error bars
	}
	c_ratio->Print( Form("plots/xsecResults/xsec_OCratio_%s.pdf", inputname.c_str()) );



	TCanvas* c_ratio_err = new TCanvas("Xsec_ratio", "Xsec_ratio",1700,1000);
	c_ratio_err -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit_err[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_postfit_err[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_postfit_err[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio_postfit_err[nbcth] -> SetLineColor(kBlue);

		c_ratio_err -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio_postfit_err[nbcth]->Draw("hist"); // with error bars
	}
	c_ratio_err->Print( Form("plots/xsecResults/xsec_error_OCratio_%s.pdf", inputname.c_str()) );

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
	
	c_eff->Print( Form("plots/xsecResults/efficiency_%s.pdf", inputname.c_str()) );
	//======================================================================================================




	//======================================================================================================
	delete fin;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}







// // To compute the cross section :
// // Normalise by the flux * nb of targets * efficiency * bin width
// void NormaliseXsec(TH1D* &hist, TH1D* efficiency, int target, BinningTools bintool)
// {
// 	double flux = 1.12E13;
// 	double Ntarg = 0;

// 	if(target==0)      Ntarg = 7.439E29;
// 	else if(target==0) Ntarg = 2.581E29;

// 	hist -> Scale(1.0 / (flux * Ntarg));
// 	hist -> Divide(efficiency);

// 	unit_scale = 0.1;

//     for(int i = 0; i < hist -> GetNbinsX(); ++i)
//     {
//         const double bin_width = bintool.GetMomBinWidth(i) * bintool.GetCosBinWidth(i) / unit_scale;
//         const double bin_value = hist.GetBinContent(i + 1);
//         hist.SetBinContent(i + 1, bin_value / bin_width);
//     }

// }