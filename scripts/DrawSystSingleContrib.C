//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root -b -q 'DrawSystSingleContrib.C+()'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


void DrawSystSingleContrib(const std::string& dir_name = "systSingleContrib", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename_TempPar = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyTemPar.root";
	string infilename_FluxPar = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyFluxPar.root";
	string infilename_XsecPar = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyXsecPar.root";
	string infilename_DetePar = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_fit1_asimov_onlyDetPar.root";

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
	std::cout << "===== Store all outputs in TFile =====" << std::endl;

	TFile* fin_TempPar = new TFile(infilename_TempPar.c_str());
	TFile* fin_FluxPar = new TFile(infilename_FluxPar.c_str());
	TFile* fin_XsecPar = new TFile(infilename_XsecPar.c_str());
	TFile* fin_DetePar = new TFile(infilename_DetePar.c_str());
	//======================================================================================================



	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with cross-section result =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec_TempPar(Ntarget);
	vector< vector<TH1D*> > h_xsec_FluxPar(Ntarget);
	vector< vector<TH1D*> > h_xsec_XsecPar(Ntarget);
	vector< vector<TH1D*> > h_xsec_DetePar(Ntarget);
	vector< vector<TH1D*> > h_xsec_TempPar_err(Ntarget);
	vector< vector<TH1D*> > h_xsec_FluxPar_err(Ntarget);
	vector< vector<TH1D*> > h_xsec_XsecPar_err(Ntarget);
	vector< vector<TH1D*> > h_xsec_DetePar_err(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_TempPar[itar].push_back(     (TH1D*)(fin_TempPar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_FluxPar[itar].push_back(     (TH1D*)(fin_FluxPar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_XsecPar[itar].push_back(     (TH1D*)(fin_XsecPar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_DetePar[itar].push_back(     (TH1D*)(fin_DetePar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_TempPar_err[itar].push_back( (TH1D*)(fin_TempPar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
			h_xsec_FluxPar_err[itar].push_back( (TH1D*)(fin_FluxPar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
			h_xsec_XsecPar_err[itar].push_back( (TH1D*)(fin_XsecPar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
			h_xsec_DetePar_err[itar].push_back( (TH1D*)(fin_DetePar->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
		}


	vector<TH1D*> h_ratio_TempPar;
	vector<TH1D*> h_ratio_FluxPar;
	vector<TH1D*> h_ratio_XsecPar;
	vector<TH1D*> h_ratio_DetePar;
	vector<TH1D*> h_ratio_TempPar_err;
	vector<TH1D*> h_ratio_FluxPar_err;
	vector<TH1D*> h_ratio_XsecPar_err;
	vector<TH1D*> h_ratio_DetePar_err;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_TempPar.push_back(     (TH1D*)(fin_TempPar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_FluxPar.push_back(     (TH1D*)(fin_FluxPar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_XsecPar.push_back(     (TH1D*)(fin_XsecPar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_DetePar.push_back(     (TH1D*)(fin_DetePar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_TempPar_err.push_back( (TH1D*)(fin_TempPar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_err", nbcth)) );
		h_ratio_FluxPar_err.push_back( (TH1D*)(fin_FluxPar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_err", nbcth)) );
		h_ratio_XsecPar_err.push_back( (TH1D*)(fin_XsecPar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_err", nbcth)) );
		h_ratio_DetePar_err.push_back( (TH1D*)(fin_DetePar->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_err", nbcth)) );
	}

	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set histos with errors =====" << std::endl;
	
	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_TempPar_err[itar][nbcth] -> Reset();
			h_xsec_FluxPar_err[itar][nbcth] -> Reset();
			h_xsec_XsecPar_err[itar][nbcth] -> Reset();
			h_xsec_DetePar_err[itar][nbcth] -> Reset();

			for(int i=1; i<=h_xsec_TempPar_err[itar][nbcth]->GetNbinsX(); i++)
			{
				double err_TempPar  = h_xsec_TempPar[itar][nbcth] -> GetBinError(i);
				double err_FluxPar  = h_xsec_FluxPar[itar][nbcth] -> GetBinError(i);
				double err_XsecPar  = h_xsec_XsecPar[itar][nbcth] -> GetBinError(i);
				double err_DetePar  = h_xsec_DetePar[itar][nbcth] -> GetBinError(i);

				double mean_TempPar = h_xsec_TempPar[itar][nbcth] -> GetBinContent(i);
				double mean_FluxPar = h_xsec_FluxPar[itar][nbcth] -> GetBinContent(i);
				double mean_XsecPar = h_xsec_XsecPar[itar][nbcth] -> GetBinContent(i);
				double mean_DetePar = h_xsec_DetePar[itar][nbcth] -> GetBinContent(i);
				
				h_xsec_TempPar_err[itar][nbcth] -> SetBinContent(i, err_TempPar/mean_TempPar);
				h_xsec_FluxPar_err[itar][nbcth] -> SetBinContent(i, err_FluxPar/mean_FluxPar);
				h_xsec_XsecPar_err[itar][nbcth] -> SetBinContent(i, err_XsecPar/mean_XsecPar);
				h_xsec_DetePar_err[itar][nbcth] -> SetBinContent(i, err_DetePar/mean_DetePar);
			}
		}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_TempPar_err[nbcth] -> Reset();
		h_ratio_FluxPar_err[nbcth] -> Reset();
		h_ratio_XsecPar_err[nbcth] -> Reset();
		h_ratio_DetePar_err[nbcth] -> Reset();

		for(int i=1; i<=h_ratio_TempPar_err[nbcth]->GetNbinsX(); i++)
		{
			double err_TempPar  = h_ratio_TempPar[nbcth] -> GetBinError(i);
			double err_FluxPar  = h_ratio_FluxPar[nbcth] -> GetBinError(i);
			double err_XsecPar  = h_ratio_XsecPar[nbcth] -> GetBinError(i);
			double err_DetePar  = h_ratio_DetePar[nbcth] -> GetBinError(i);

			double mean_TempPar = h_ratio_TempPar[nbcth] -> GetBinContent(i);
			double mean_FluxPar = h_ratio_FluxPar[nbcth] -> GetBinContent(i);
			double mean_XsecPar = h_ratio_XsecPar[nbcth] -> GetBinContent(i);
			double mean_DetePar = h_ratio_DetePar[nbcth] -> GetBinContent(i);
			
			h_ratio_TempPar_err[nbcth] -> SetBinContent(i, err_TempPar/mean_TempPar);
			h_ratio_FluxPar_err[nbcth] -> SetBinContent(i, err_FluxPar/mean_FluxPar);
			h_ratio_XsecPar_err[nbcth] -> SetBinContent(i, err_XsecPar/mean_XsecPar);
			h_ratio_DetePar_err[nbcth] -> SetBinContent(i, err_DetePar/mean_DetePar);
		}
	}
	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	TCanvas* c_xsec[Ntarget];
	TCanvas* c_xsec_err[Ntarget];
	vector<TLegend*> leg(Ntarget);
	double max;

	for(int itar = 0; itar < Ntarget; itar++)
	{
		c_xsec_err[itar] = new TCanvas(Form("Xsec_error_on_%s", targetlist[itar].c_str()),Form("Xsec_error_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec_err[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_DetePar_err[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_DetePar_err[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec_DetePar_err[itar][nbcth] -> GetYaxis() -> SetTitle("Relative fit error");

			// h_xsec_TempPar_err[itar][nbcth] -> SetLineWidth(2);

			h_xsec_TempPar_err[itar][nbcth] -> SetLineColor(kBlack);
			h_xsec_FluxPar_err[itar][nbcth] -> SetLineColor(kBlue);
			h_xsec_XsecPar_err[itar][nbcth] -> SetLineColor(kRed);
			h_xsec_DetePar_err[itar][nbcth] -> SetLineColor(kGreen+1);

			c_xsec_err[itar] -> cd(nbcth+1);

			if(nbcth!=0) gPad->SetLogx();

			h_xsec_DetePar_err[itar][nbcth] -> Draw("hist");
			h_xsec_FluxPar_err[itar][nbcth] -> Draw("same hist");
			h_xsec_XsecPar_err[itar][nbcth] -> Draw("same hist");
			h_xsec_TempPar_err[itar][nbcth] -> Draw("same hist");

			h_xsec_DetePar_err[itar][nbcth] -> GetYaxis() -> SetRangeUser(0.0, 0.2);

			if(nbcth==0)
			{
				leg[itar] = new TLegend(0.2,0.40,0.8,0.85);
				leg[itar] -> SetFillColor(0);
				leg[itar] -> SetBorderSize(1);
				leg[itar] -> SetFillStyle(0);
				//leg[itar]->SetTextSize(0.075);
				leg[itar] -> SetHeader("Systematic contributions");
				leg[itar] -> AddEntry(h_xsec_TempPar_err[itar][0], "Statistics", "l");
				leg[itar] -> AddEntry(h_xsec_FluxPar_err[itar][0], "Flux", "l");
				leg[itar] -> AddEntry(h_xsec_XsecPar_err[itar][0], "Cross-section model", "l");
				leg[itar] -> AddEntry(h_xsec_DetePar_err[itar][0], "Detector", "l");
				c_xsec_err[itar]->cd(1);
				leg[itar]->Draw();
			}

		}
		c_xsec_err[itar]->Print(Form("plots/%s/xsec_error_singleContrib_%s.pdf", dir_name.c_str(), targetlist[itar].c_str()));
	}





	TCanvas* c_ratio_err = new TCanvas("Xsec_ratio", "Xsec_ratio",1700,1000);
	c_ratio_err -> Divide(3,3);

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_DetePar_err[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_DetePar_err[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_DetePar_err[nbcth] -> GetYaxis() -> SetTitle("Relative fit error");

		// h_ratio_TempPar_err[nbcth] -> SetLineWidth(2);

		h_ratio_TempPar_err[nbcth] -> SetLineColor(kBlack);
		h_ratio_FluxPar_err[nbcth] -> SetLineColor(kBlue);
		h_ratio_XsecPar_err[nbcth] -> SetLineColor(kRed);
		h_ratio_DetePar_err[nbcth] -> SetLineColor(kGreen+1);

		c_ratio_err -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_ratio_DetePar_err[nbcth] -> Draw("hist");
		h_ratio_FluxPar_err[nbcth] -> Draw("same hist");
		h_ratio_XsecPar_err[nbcth] -> Draw("same hist");
		h_ratio_TempPar_err[nbcth] -> Draw("same hist");

		h_ratio_DetePar_err[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 0.2);

		if(nbcth==0)
		{
			TLegend* leg2 = new TLegend(0.2,0.40,0.8,0.85);
			leg2 -> SetFillColor(0);
			leg2 -> SetBorderSize(1);
			leg2 -> SetFillStyle(0);
			//leg2->SetTextSize(0.075);
			leg2 -> SetHeader("Systematic contributions");
			leg2 -> AddEntry(h_ratio_TempPar_err[0], "Statistics", "l");
			leg2 -> AddEntry(h_ratio_FluxPar_err[0], "Flux", "l");
			leg2 -> AddEntry(h_ratio_XsecPar_err[0], "Cross-section model", "l");
			leg2 -> AddEntry(h_ratio_DetePar_err[0], "Detector", "l");
			c_ratio_err->cd(1);
			leg2->Draw();
		}

	}
	c_ratio_err->Print(Form("plots/%s/xsec_error_singleContrib_OCratio.pdf", dir_name.c_str()));
	


	//======================================================================================================


	//======================================================================================================
	delete fin_TempPar;
	delete fin_FluxPar;
	delete fin_XsecPar;
	delete fin_DetePar;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================
}
