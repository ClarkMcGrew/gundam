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

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double calcChi2_M(TH1D*, TH1D*, TMatrixD);


void DrawXsec(string inputname = "fit3_statFluc", const std::string& dir_name = "fakedata/statFluc", bool modelComparison=false, string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt")
{

	bool drawNEUT = true;
	// if(inputname == "fit1_asimov" || inputname == "fit3_statFluc") drawNEUT = false;
	// if(inputname == "fit1_asimov") drawNEUT = false;

	bool drawFakeData = true;
	if(inputname == "fit2_dataCS_fakedataSignal" || inputname == "fit2_dataCS_fakedataSignal_noMichel" || inputname == "fit2_data" || inputname == "fit2_data_noMichel") drawFakeData = false;

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());
	string infilename_comp = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsecUpdatedNEUT.root";

	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set binning using BinningTools class =====" << std::endl;
	
	BinningTools bin; 
	bin.SetBinning(fbinning.c_str());

	int Nbins = bin.GetNbins();//Total number of bins
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

	TFile* fin            = new TFile(infilename.c_str());
	TFile* fin_otherModel = new TFile(infilename_comp.c_str());
	//======================================================================================================



	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with cross-section result and truth =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec_truth_data(Ntarget);
	vector< vector<TH1D*> > h_xsec_truth_nominal(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit(Ntarget);
	vector< vector<TH1D*> > h_xsec_otherModel(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit_err(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			if(inputname != "fit3_statFluc_reg" && inputname != "fit3_statFluc" && inputname != "fit3_asimov_statFluc_reg")
			{
				h_xsec_truth_data[itar].push_back(    (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_data",    targetlist[itar].c_str(), nbcth) )) );
				h_xsec_truth_nominal[itar].push_back( (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_nominal", targetlist[itar].c_str(), nbcth) )) );
				if (modelComparison)
					h_xsec_otherModel[itar].push_back((TH1D*)(fin_otherModel -> Get( Form("Hist_%d", nbcth) )) );
				// std::cout << "DEBUG : No stat+syst fluctuation fit." << std::endl;
			}
			else
			{
				h_xsec_truth_data[itar].push_back(    (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_nominal", targetlist[itar].c_str(), nbcth) )) );
				h_xsec_truth_nominal[itar].push_back( (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_data",    targetlist[itar].c_str(), nbcth) )) );
				if (modelComparison)
					h_xsec_otherModel[itar].push_back((TH1D*)(fin_otherModel -> Get( Form("Hist_%d", nbcth) )) );
				// std::cout << "DEBUG : Stat+syst fluctuation fit." << std::endl;
			}
			h_xsec_postfit[itar].push_back(       (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_postfit", targetlist[itar].c_str(), nbcth) )) );
			
			h_xsec_postfit_err[itar].push_back(   (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_postfit", targetlist[itar].c_str(), nbcth) )) -> Clone(Form("CC0pi%s_cos_bin%d_postfit_err", targetlist[itar].c_str(), nbcth)) );
			
		}

	vector< TH1D* > h_ratio_truth_data;
	vector< TH1D* > h_ratio_truth_nominal;
	vector< TH1D* > h_ratio_postfit;
	vector< TH1D* > h_ratio_postfit_err;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		if(inputname != "fit3_statFluc_reg" && inputname != "fit3_statFluc" && inputname != "fit3_asimov_statFluc_reg")
		{
			h_ratio_truth_data.push_back(    (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_data",    nbcth) )) );
			h_ratio_truth_nominal.push_back( (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_nominal", nbcth) )) );
			// std::cout << "DEBUG : No stat+syst fluctuation fit." << std::endl;
		}
		else
		{
			h_ratio_truth_data.push_back(    (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_nominal", nbcth) )) );
			h_ratio_truth_nominal.push_back( (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_data",    nbcth) )) );
			// std::cout << "DEBUG : Stat+syst fluctuation fit." << std::endl;
		}
		h_ratio_postfit.push_back(       (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_postfit", nbcth) )) );

		h_ratio_postfit_err.push_back(   (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_postfit", nbcth) )) -> Clone(Form("CC0piOCRatio_cos_bin%d_errs", nbcth)) );
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
	std::cout << "===== Get covariance matrix and compute chi2 =====" << std::endl;
	
	TMatrixDSym *cov_mat       = (TMatrixDSym*)(fin -> Get("xsec_cov"));
	TMatrixDSym *cov_mat_ratio = (TMatrixDSym*)(fin -> Get("ratio_cov"));

	TH1D* h_sel_postfit = (TH1D*)(fin -> Get("sel_best_fit"));
	TH1D* h_tru_nominal;
	TH1D* h_tru_data;
	TH1D* h_sel_postfit_ratio = (TH1D*)(fin -> Get("sel_best_fit_ratio"));
	TH1D* h_tru_nominal_ratio;
	TH1D* h_tru_data_ratio;
	if(inputname != "fit3_statFluc_reg" && inputname != "fit3_statFluc" && inputname != "fit3_asimov_statFluc_reg")
	{
		h_tru_nominal       = (TH1D*)(fin -> Get("tru_best_fit"));
		h_tru_data          = (TH1D*)(fin -> Get("fake_data_concat"));
		h_tru_nominal_ratio = (TH1D*)(fin -> Get("tru_best_fit_ratio"));
		h_tru_data_ratio    = (TH1D*)(fin -> Get("fake_data_ratio"));
		// std::cout << "DEBUG : No stat+syst fluctuation fit." << std::endl;
	}
	else
	{
		h_tru_nominal       = (TH1D*)(fin -> Get("fake_data_concat"));
		h_tru_data          = (TH1D*)(fin -> Get("tru_best_fit"));
		h_tru_nominal_ratio = (TH1D*)(fin -> Get("fake_data_ratio"));
		h_tru_data_ratio    = (TH1D*)(fin -> Get("tru_best_fit_ratio"));
		// std::cout << "DEBUG : Stat+syst fluctuation fit." << std::endl;
	}

	double chi2      = -999;
	double chi2_NEUT = -999;

	chi2      = calcChi2_M(h_sel_postfit, h_tru_data,    *cov_mat);
	chi2_NEUT = calcChi2_M(h_sel_postfit, h_tru_nominal, *cov_mat);

	std::cout << "===== Chi2 data = " << chi2      << " with probability (for nbins="<<2*Nbins<<") = " << TMath::Prob(chi2,      2*Nbins) << std::endl;
	std::cout << "===== Chi2 NEUT = " << chi2_NEUT << " with probability (for nbins="<<2*Nbins<<") = " << TMath::Prob(chi2_NEUT, 2*Nbins) << std::endl;

	double chi2_ratio      = -999;
	double chi2_NEUT_ratio = -999;

	chi2_ratio      = calcChi2_M(h_sel_postfit_ratio, h_tru_data_ratio,    *cov_mat_ratio);
	chi2_NEUT_ratio = calcChi2_M(h_sel_postfit_ratio, h_tru_nominal_ratio, *cov_mat_ratio);

	std::cout << "===== Chi2 data Ratio = " << chi2_ratio      << " with probability (for nbins="<<Nbins<<") = " << TMath::Prob(chi2_ratio,      Nbins) << std::endl;
	std::cout << "===== Chi2 NEUT Ratio = " << chi2_NEUT_ratio << " with probability (for nbins="<<Nbins<<") = " << TMath::Prob(chi2_NEUT_ratio, Nbins) << std::endl;

	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Draw cross-section in analysis binning =====" << std::endl;
	
	TCanvas* c_xsec_allBins = new TCanvas("c_xsec_allBins", "c_xsec_allBins",1800,800);
	h_tru_data    -> SetTitle("Cross-section in analysis bins");
	h_tru_data    -> GetXaxis() -> SetTitle("Analysis bins");
	h_tru_data    -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
	h_sel_postfit -> SetTitle("Cross-section in analysis bins");
	h_sel_postfit -> GetXaxis() -> SetTitle("Analysis bins");
	h_sel_postfit -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
	h_tru_nominal -> SetTitle("Cross-section in analysis bins");
	h_tru_nominal -> GetXaxis() -> SetTitle("Analysis bins");
	h_tru_nominal -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
	
	double max, min;
	max = MAX(h_tru_data->GetMaximum(), h_sel_postfit->GetMaximum());
	min = MIN(h_tru_data->GetMinimum(), h_sel_postfit->GetMinimum());
	h_tru_data    -> GetYaxis() -> SetRangeUser(0.1*min, 1.3*max);
	h_sel_postfit -> GetYaxis() -> SetRangeUser(0.1*min, 1.3*max);
	h_tru_nominal -> GetYaxis() -> SetRangeUser(0.1*min, 1.3*max);

	h_sel_postfit -> Draw(); // Need to do this in order for the lines to be drawn
	
	TLine *verline[2*Nbins_costh];
	int binSide = 0;
	for(int il = 0; il < Nbins_costh; il++)
	{
		binSide += Nmombins[il];
		verline[il]               = new TLine(binSide,         0, binSide,         1.3*max );
		verline[il + Nbins_costh] = new TLine(binSide + Nbins, 0, binSide + Nbins, 1.3*max );
	}

	for(int il = 0; il < 2*Nbins_costh-1; il++)
	{
		verline[il] -> SetLineWidth(1);
		verline[il] -> SetLineStyle(1);
		verline[il] -> SetLineColor(kGray);
		verline[il] -> Draw();
	}
	TLine *verlineTargets = new TLine(Nbins, 0.0, Nbins, 1.3*max );
	verlineTargets -> Draw();

	TLine *horlineTargets = new TLine(0.0, 0.0, 2*Nbins, 0.0 );
	horlineTargets -> Draw();

	if(drawFakeData)
		h_tru_data    -> SetLineColor(kGreen+2);
	if(drawNEUT)
		h_tru_nominal -> SetLineColor(kRed+2);
	h_sel_postfit     -> SetLineColor(kBlack);
	h_sel_postfit     -> SetMarkerColor(kBlue+2);
    h_sel_postfit     -> SetMarkerStyle(kFullCircle);
	
	if(drawFakeData)
		h_tru_data    -> Draw("same");
	if(drawNEUT)
		h_tru_nominal -> Draw("same");
	h_sel_postfit     -> Draw("same P E1");


	TLegend* legendAllBins = new TLegend(0.65,0.75,0.92,0.9);
	legendAllBins -> SetFillColor(0);
	legendAllBins -> SetBorderSize(1);
	legendAllBins -> AddEntry(h_sel_postfit,    "post-fit", "lep");
	if(drawNEUT)
		legendAllBins -> AddEntry(h_tru_nominal, Form("NEUT truth, #chi^{2} = %d",      (int)chi2_NEUT), "l");
	if(drawFakeData)
		legendAllBins ->     AddEntry(h_tru_data,    Form("Fake data truth, #chi^{2} = %d", (int)chi2),      "l");
	legendAllBins->Draw();

	c_xsec_allBins -> Print(Form("plots/xsecResults/%s/xsec_%s_allBins.pdf", dir_name.c_str(), inputname.c_str()));







	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots of cross-section =====" << std::endl;                                                                                                                      

	TCanvas* c_xsec[Ntarget];
	TCanvas* c_xsec_noLastBin[Ntarget];
	TCanvas* c_xsec_err[Ntarget];
	vector<TLegend*> leg(Ntarget+1);
	vector<TLegend*> leg_noLastBin(Ntarget+1);
	
	for(int itar = 0; itar < Ntarget; itar++)
	{
		c_xsec[itar] = new TCanvas(Form("Xsec_on_%s", targetlist[itar].c_str()),Form("Xsec_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_truth_data[itar][nbcth]    -> SetTitle(h_title[nbcth]);
			h_xsec_truth_data[itar][nbcth]    -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_truth_data[itar][nbcth]    -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
			h_xsec_postfit[itar][nbcth]       -> SetTitle(h_title[nbcth]);
			h_xsec_postfit[itar][nbcth]       -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_postfit[itar][nbcth]       -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
			h_xsec_truth_nominal[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_truth_nominal[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_truth_nominal[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");

			max = MAX(h_xsec_truth_data[itar][nbcth]->GetMaximum(), h_xsec_postfit[itar][nbcth]->GetMaximum());
			min = MIN(h_xsec_truth_data[itar][nbcth]->GetMinimum(), h_xsec_postfit[itar][nbcth]->GetMinimum());
			max = MAX(h_xsec_truth_nominal[itar][nbcth]->GetMaximum(), max);
			min = MIN(h_xsec_truth_nominal[itar][nbcth]->GetMinimum(), min);
			h_xsec_truth_data[itar][nbcth]    -> GetYaxis() -> SetRangeUser(1E-42, 1.3*max);
			h_xsec_truth_nominal[itar][nbcth] -> GetYaxis() -> SetRangeUser(1E-42, 1.3*max);
			h_xsec_postfit[itar][nbcth]       -> GetYaxis() -> SetRangeUser(1E-42, 1.3*max);

			if(drawFakeData)
				h_xsec_truth_data[itar][nbcth]    -> SetLineColor(kGreen+2);
			if(modelComparison)
				h_xsec_otherModel[itar][nbcth]    -> SetLineColor(kGreen+2);
			if(drawNEUT)
				h_xsec_truth_nominal[itar][nbcth] -> SetLineColor(kRed+2);
			h_xsec_postfit[itar][nbcth]           -> SetLineColor(kBlack);
			h_xsec_postfit[itar][nbcth]           -> SetMarkerColor(kBlue+2);
		    h_xsec_postfit[itar][nbcth]           -> SetMarkerStyle(kFullCircle);
	

			c_xsec[itar] -> cd(nbcth+1);

			if(nbcth!=0) gPad->SetLogx();
			gPad->SetLogy();

			if(min<0) std::cout << "ATTENTION : there is a negative cross-section value in costheta bin " << nbcth << std::endl;

			if(drawFakeData)
				h_xsec_truth_data[itar][nbcth]    -> Draw("hist");
			if(drawNEUT) 
				h_xsec_truth_nominal[itar][nbcth] -> Draw("same hist");
			if(modelComparison)
				h_xsec_otherModel[itar][nbcth]    -> Draw("same hist");
			h_xsec_postfit[itar][nbcth]           -> Draw("same P E1");

			if(nbcth==0)
			{
				leg[itar] = new TLegend(0.2,0.2,0.8,0.6);
				leg[itar] -> SetFillColor(0);
				leg[itar] -> SetBorderSize(0);
				leg[itar] -> SetFillStyle(0);
				leg[itar] ->     AddEntry(h_xsec_postfit[itar][0],      "post-fit", "lep");
				if(drawNEUT)
					leg[itar] -> AddEntry(h_xsec_truth_nominal[itar][0], Form("NEUT truth, #chi^{2} = %d",      (int)chi2_NEUT), "l");
				if(drawFakeData)
					leg[itar] -> AddEntry(h_xsec_truth_data[itar][0],    Form("Fake data truth, #chi^{2} = %d", (int)chi2),      "l");
				if(modelComparison)
					leg[itar] -> AddEntry(h_xsec_otherModel[itar][0],    Form("Updated NEUT"),                                   "l");
				c_xsec[itar]->cd(1);
				leg[itar]->Draw();
			}

		}
		c_xsec[itar]->Print(Form("plots/xsecResults/%s/xsec_%s_%s_withHighMomBin.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));







		c_xsec_noLastBin[itar] = new TCanvas(Form("Xsec_on_%s", targetlist[itar].c_str()),Form("Xsec_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec_noLastBin[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_truth_data[itar][nbcth]    -> SetTitle(h_title[nbcth]);
			h_xsec_truth_data[itar][nbcth]    -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_truth_data[itar][nbcth]    -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
			h_xsec_postfit[itar][nbcth]       -> SetTitle(h_title[nbcth]);
			h_xsec_postfit[itar][nbcth]       -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_postfit[itar][nbcth]       -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");
			h_xsec_truth_nominal[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_truth_nominal[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_truth_nominal[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");

			max = MAX(h_xsec_truth_data[itar][nbcth]->GetMaximum(), h_xsec_postfit[itar][nbcth]->GetMaximum());
			min = MIN(h_xsec_truth_data[itar][nbcth]->GetMinimum(), h_xsec_postfit[itar][nbcth]->GetMinimum());
			max = MAX(h_xsec_truth_nominal[itar][nbcth]->GetMaximum(), max);
			min = MIN(h_xsec_truth_nominal[itar][nbcth]->GetMinimum(), min);
			h_xsec_truth_data[itar][nbcth]    -> GetYaxis() -> SetRangeUser(0.5*min, 1.3*max);
			h_xsec_truth_nominal[itar][nbcth] -> GetYaxis() -> SetRangeUser(0.5*min, 1.3*max);
			h_xsec_postfit[itar][nbcth]       -> GetYaxis() -> SetRangeUser(0.5*min, 1.3*max);

			// Set range without high momentum bin
			double lastBinEdge = mombins[nbcth][(Nmombins[nbcth]-1)]; // Nmombins = number of mom bins in this costh slice
			// std::cout << "costh bin nb = " << nbcth << ", Nmombins = " << Nmombins[nbcth] << " and last bin edge = " << lastBinEdge << std::endl;
			h_xsec_truth_data[itar][nbcth]    -> GetXaxis() -> SetRangeUser(0.0, lastBinEdge);
			h_xsec_truth_nominal[itar][nbcth] -> GetXaxis() -> SetRangeUser(0.0, lastBinEdge);
			h_xsec_postfit[itar][nbcth]       -> GetXaxis() -> SetRangeUser(0.0, lastBinEdge);

			// if(drawFakeData)
			// 	h_xsec_truth_data[itar][nbcth]    -> SetLineColor(kGreen+2);
			// if(modelComparison)
			// 	h_xsec_otherModel[itar][nbcth]    -> SetLineColor(kGreen+2);
			// if(drawNEUT)
			// 	h_xsec_truth_nominal[itar][nbcth] -> SetLineColor(kRed+2);
			// h_xsec_postfit[itar][nbcth]           -> SetLineColor(kBlack);
			// h_xsec_postfit[itar][nbcth]           -> SetMarkerColor(kBlue+2);
		 //    h_xsec_postfit[itar][nbcth]           -> SetMarkerStyle(kFullCircle);
	

			c_xsec_noLastBin[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			if(drawFakeData)
				h_xsec_truth_data[itar][nbcth]    -> Draw("hist");
			if(drawNEUT) 
				h_xsec_truth_nominal[itar][nbcth] -> Draw("same hist");
			if(modelComparison)
				h_xsec_otherModel[itar][nbcth]    -> Draw("same hist");
			h_xsec_postfit[itar][nbcth]           -> Draw("same P E1");

			if(nbcth==0)
			{
				leg_noLastBin[itar] = new TLegend(0.2,0.2,0.8,0.6);
				leg_noLastBin[itar] -> SetFillColor(0);
				leg_noLastBin[itar] -> SetBorderSize(0);
				leg_noLastBin[itar] -> SetFillStyle(0);
				leg_noLastBin[itar] ->     AddEntry(h_xsec_postfit[itar][0],      "post-fit", "lep");
				if(drawNEUT)
					leg_noLastBin[itar] -> AddEntry(h_xsec_truth_nominal[itar][0], Form("NEUT truth, #chi^{2} = %d",      (int)chi2_NEUT), "l");
				if(drawFakeData)
					leg_noLastBin[itar] -> AddEntry(h_xsec_truth_data[itar][0],    Form("Fake data truth, #chi^{2} = %d", (int)chi2),      "l");
				if(modelComparison)
					leg_noLastBin[itar] -> AddEntry(h_xsec_otherModel[itar][0],    Form("Updated NEUT"),                                   "l");
				c_xsec_noLastBin[itar]->cd(1);
				leg_noLastBin[itar]->Draw();
			}

		}
		c_xsec_noLastBin[itar]->Print(Form("plots/xsecResults/%s/xsec_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));









		c_xsec_err[itar] = new TCanvas(Form("Xsec_error_on_%s", targetlist[itar].c_str()),Form("Xsec_error_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec_err[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_postfit_err[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_postfit_err[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_postfit_err[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^{2}}{nucleon GeV/c}]");

			h_xsec_postfit_err[itar][nbcth] -> SetLineColor(kBlue+2);

			c_xsec_err[itar] -> cd(nbcth+1);

			if(nbcth!=0) gPad->SetLogx();

			h_xsec_postfit_err[itar][nbcth] -> GetYaxis() -> SetRangeUser(0.0, 1.0);
			h_xsec_postfit_err[itar][nbcth] -> Draw("hist");

		}
		c_xsec_err[itar]->Print(Form("plots/xsecResults/%s/xsec_error_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));
	}



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots of O/C ratio =====" << std::endl;                                                                                                                      


	TCanvas* c_ratio = new TCanvas("Xsec_ratio", "Xsec_ratios",1700,1000);
	c_ratio -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_truth_data[nbcth]    -> SetTitle(h_title[nbcth]);
		h_ratio_truth_data[nbcth]    -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_ratio_truth_data[nbcth]    -> GetYaxis() -> SetTitle("#sigma^{O} / #sigma^{C}");
		h_ratio_truth_data[nbcth]    -> GetYaxis() -> SetRangeUser(0.0, 2.0);
		h_ratio_truth_nominal[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_truth_nominal[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_ratio_truth_nominal[nbcth] -> GetYaxis() -> SetTitle("#sigma^{O} / #sigma^{C}");
		h_ratio_truth_nominal[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 2.0);
		h_ratio_postfit[nbcth]       -> SetTitle(h_title[nbcth]);
		h_ratio_postfit[nbcth]       -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_ratio_postfit[nbcth]       -> GetYaxis() -> SetTitle("#sigma^{O} / #sigma^{C}");
		h_ratio_truth_data[nbcth]    -> GetYaxis() -> SetRangeUser(0.0, 3.0);
		h_ratio_truth_nominal[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 3.0);
		h_ratio_postfit[nbcth]       -> GetYaxis() -> SetRangeUser(0.0, 3.0);
		
		if(drawFakeData)
			h_ratio_truth_data[nbcth]    -> SetLineColor(kGreen+2);
		if(drawNEUT)
			h_ratio_truth_nominal[nbcth] -> SetLineColor(kRed+2);
		h_ratio_postfit[nbcth]           -> SetLineColor(kBlack);
		h_ratio_postfit[nbcth]           -> SetMarkerColor(kBlue+2);
		h_ratio_postfit[nbcth]           -> SetMarkerStyle(kFullCircle);

		c_ratio -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		if(drawFakeData)
			h_ratio_truth_data[nbcth]    -> Draw("hist");
		if(drawNEUT) 
			h_ratio_truth_nominal[nbcth] -> Draw("same hist");
		h_ratio_postfit[nbcth]           -> Draw("same P E1");

		if(nbcth==0)
		{
			TLegend* legratio = new TLegend(0.19,0.48,0.78,0.88);
			legratio -> SetFillColor(0);
			legratio -> SetBorderSize(0);
			legratio -> SetFillStyle(0);
			legratio ->     AddEntry(h_ratio_postfit[0],       "post-fit", "lep");
			if(drawNEUT)
				legratio -> AddEntry(h_ratio_truth_nominal[0], Form("NEUT truth, #chi^{2} = %d",      (int)chi2_NEUT_ratio), "l");
			if(drawFakeData)
				legratio -> AddEntry(h_ratio_truth_data[0],    Form("Fake data truth, #chi^{2} = %d", (int)chi2_ratio),      "l");
			c_ratio->cd(1);
			legratio->Draw();
		}
	}
	c_ratio->Print( Form("plots/xsecResults/%s/xsec_%s_OCratio.pdf", dir_name.c_str(), inputname.c_str()) );



	TCanvas* c_ratio_err = new TCanvas("Xsec_ratio", "Xsec_ratio",1700,1000);
	c_ratio_err -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit_err[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_postfit_err[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_ratio_postfit_err[nbcth] -> GetYaxis() -> SetTitle("#sigma^{O} / #sigma^{C}");

		h_ratio_postfit_err[nbcth] -> SetLineColor(kBlue+2);

		c_ratio_err -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_ratio_postfit_err[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 1.0);
		h_ratio_postfit_err[nbcth] -> Draw("hist");
	}
	c_ratio_err->Print( Form("plots/xsecResults/%s/xsec_error_%s_OCratio.pdf", dir_name.c_str(), inputname.c_str()) );

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

	h_eff->Draw("hist");
	
	c_eff->Print( Form("plots/xsecResults/%s/efficiency_%s.pdf", dir_name.c_str(), inputname.c_str()) );
	//======================================================================================================




	//======================================================================================================
	// delete fin;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}





double calcChi2_M(TH1D* h1, TH1D* h2, TMatrixD covar)
{
	double chi2=0;

	cout << "calcChi2_M::Starting chi2 calc" << endl;
	cout << "calcChi2_M::h1 (data) nbins: " << h1->GetXaxis()->GetNbins() << endl;
	cout << "calcChi2_M::h2 (MC) nbins: " << h2->GetXaxis()->GetNbins() << endl;
	cout << "calcChi2_M::matrix rows: " << covar.GetNrows() << endl;

	cout << "calcChi2_M::Printing data histo to be sure of content:" << endl;
	h1->Print("all");

	cout << "calcChi2_M::Determinant of covar is: " << covar.Determinant() << endl;

	//Matrix inversion section:
	bool inversionError = false;
	TMatrixDSym* covar_inv = new TMatrixDSym();
	TMatrixD covar_origin = covar;
	int count = 0;

	while(true)
	{
		//covar.Print();
		//covar.SetTol(1e-100);
		//covar.Invert();
		TDecompSVD LU = TDecompSVD(covar);
		covar_inv = new TMatrixDSym(covar.GetNrows(), LU.Invert().GetMatrixArray(), "");
		bool inversionError = false;

		//Check that the matrix inversion worked 
		TMatrixD test = covar*(*covar_inv);
		for(int i=0; i<covar.GetNrows(); i++)
		{
			for(int j=0; j<covar.GetNrows(); j++)
			{
				if(i==j && (test[i][j]>1.00001 || test[i][j]<0.99999))
				// if(i==j && (test[i][j]>1.001 || test[i][j]<0.999)) //LM
				{
					if(!inversionError) cout 	<< "****** WARNING: Issue with matrix inversion in chi2 calculation."
												<< " Element ("<<i<<","<<j<<") = " << test[i][j] << endl; 
					inversionError=true;
				}
				else if(i!=j && (test[i][j]>0.0000001))
				// else if(i!=j && (test[i][j]>0.00001)) //LM
				{
				if(!inversionError) cout 	<< "****** WARNING: Issue with matrix inversion in chi2 calculation."
											<< " Element ("<<i<<","<<j<<") = " << test[i][j] << endl; 
				inversionError=true;
				}
			}
		}
		if(inversionError)
		{ 
			cout << "DEBUG info after iteration " << count << endl;
			cout << "  Input matrix: determinant is: " << covar.Determinant() <<  endl;
			cout << "Will try adding 0.00000001pc to the diagonal ..." << endl;
			// cout << "Will try adding 0.0000001pc to the diagonal ..." << endl; //LM
			for(int i=0; i<covar.GetNrows(); i++)
			{
				covar[i][i]=covar[i][i]*1.00000001;
				// covar[i][i]=covar[i][i]*1.00000001; //LM
			}
			count++;
			//if(count>10)getchar();
		}
		else
		{
			cout << "chi2 calculated successfully after interation " << count << endl;
			//if(count>0) test.Print();
			break;
		}
	}

	//covar.Print();
	//cout << "Inverted covariance matrix" << endl;
	for(int i=0; i<covar.GetNrows(); i++)
	{
		for(int j=0; j<covar.GetNrows(); j++)
		{
			chi2+= ((h1->GetBinContent(i+1)) - (h2->GetBinContent(i+1)))*(*covar_inv)[i][j]*((h1->GetBinContent(j+1)) - (h2->GetBinContent(j+1)));
			// if(i==j)
			// {
			// 	cout << "calcChi2_M::bin " << i << ": " << h1->GetBinError(i+1) << " should be the same as " << sqrt(covar_origin[i][j]) << endl;
			// }
		}
	}
	cout << "calcChi2_M::chi2 is: " << chi2 << endl;
	return chi2;
}
