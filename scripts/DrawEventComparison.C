//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root -b -q 'DrawEventComparison.C+()'
//         root -b -q 'DrawEventComparison.C+("fit3_statFluc")'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

double calculateChi2(TH1D* hist1, TH1D* hist2);

void DrawEventComparison(string inputname = "fit3_statFluc", const std::string& dir_name = "fakedata/statFluc", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/%s.root", inputname.c_str());

	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store fitter output in a TFile =====" << std::endl;
	std::cout << "===== "<< infilename <<" =====" << std::endl;

	TFile* fin = new TFile(infilename.c_str());
	//======================================================================================================


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
	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store sample name in a vector of string =====" << std::endl;
	
	const int Nfgd = 3;
	string fgdNames[Nfgd] = {"FGD1", "FGD2x", "FGD2y"};

	const int Nsample = 8;
	vector< vector< string > > nameSamples(Nfgd);

	for(int is = 0; is < Nsample; is++)// add FGD1 distributions
		nameSamples[0].push_back(Form("evhist_sam%d", is));

	for(int is = 0; is < Nsample; is++)// add FGD2x distributions
		nameSamples[1].push_back(Form("evhist_sam%d", is+Nsample));
	
	for(int is = 0; is < Nsample; is++)// add FGD2y distributions
		nameSamples[2].push_back(Form("evhist_sam%d", is+2*Nsample));

	//======================================================================================================


	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TH1D* htemp_prefit;
	TH1D* htemp_postfit;
	TH1D* htemp_data;

	vector< vector< TH1D* > > hSample_prefit(Nfgd);
	vector< vector< TH1D* > > hSample_postfit(Nfgd);
	vector< vector< TH1D* > > hSample_data(Nfgd);
	
	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
		for(int ns=0; ns<Nsample; ns++)
		{
			htemp_prefit = (TH1D*)(fin->Get(Form("%s_iter0_pred", nameSamples[ifgd][ns].c_str())));
			hSample_prefit[ifgd].push_back(htemp_prefit);

			htemp_postfit = (TH1D*)(fin->Get(Form("%s_finaliter_pred", nameSamples[ifgd][ns].c_str())));
			hSample_postfit[ifgd].push_back(htemp_postfit);

			htemp_data = (TH1D*)(fin->Get(Form("%s_iter0_data", nameSamples[ifgd][ns].c_str())));
			hSample_data[ifgd].push_back(htemp_data);
		}

	// const int Nbins = hSample_prefit[0][0] -> GetNbinsX();
	string h_title[Nsample+1] = {"muTPC",
	                           "muTPCpTPC",
	                           "muTPCpFGD",
	                           "muFGDpTPC",
	                           "muFGD",
	                           "CC1pi",
	                           "CCOther",
	                           "CCMichel",
	                           "all muTPC"
	                           };


	vector< TH1D* > hSample_prefit_allSubDet;
	vector< TH1D* > hSample_postfit_allSubDet;
	vector< TH1D* > hSample_data_allSubDet;
	
	for(int ns=0; ns<Nsample; ns++)
	{
		htemp_prefit  = (TH1D*)(fin->Get(Form("%s_iter0_pred", nameSamples[0][ns].c_str())));
		htemp_postfit = (TH1D*)(fin->Get(Form("%s_finaliter_pred", nameSamples[0][ns].c_str())));
		htemp_data    = (TH1D*)(fin->Get(Form("%s_iter0_data", nameSamples[0][ns].c_str())));

		htemp_prefit  -> Reset();
		htemp_postfit -> Reset();
		htemp_data    -> Reset();

		for(int ifgd=0; ifgd<Nfgd; ifgd++)
		{
			htemp_prefit  -> Add(hSample_prefit[ifgd][ns]);
			htemp_postfit -> Add(hSample_postfit[ifgd][ns]);
			htemp_data    -> Add(hSample_data[ifgd][ns]);
		}

		hSample_prefit_allSubDet.push_back(htemp_prefit);
		hSample_postfit_allSubDet.push_back(htemp_postfit);
		hSample_data_allSubDet.push_back(htemp_data);
	}

	//======================================================================================================


	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Put all muTPC+... together =====" << std::endl;
	
	hSample_prefit_allSubDet.push_back(  (TH1D*)(hSample_prefit_allSubDet[0] -> Clone("hSample_prefit_allSubDet_allMuTPC")));
	hSample_postfit_allSubDet.push_back( (TH1D*)(hSample_postfit_allSubDet[0]-> Clone("hSample_postfit_allSubDet_allMuTPC")));
	hSample_data_allSubDet.push_back(    (TH1D*)(hSample_data_allSubDet[0]   -> Clone("hSample_data_allSubDet_allMuTPC")));

	for(int ns=1; ns<3; ns++)
	{
		hSample_prefit_allSubDet[Nsample] -> Add(hSample_prefit_allSubDet[ns]);
		hSample_postfit_allSubDet[Nsample]-> Add(hSample_postfit_allSubDet[ns]);
		hSample_data_allSubDet[Nsample]   -> Add(hSample_data_allSubDet[ns]);
	}
	//======================================================================================================




	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get pre-fit and post-fit chi2 =====" << std::endl;
	
	TH1D* chi2_tot  = (TH1D*)(fin->Get("chi2_tot_periter"));
	TH1D* chi2_stat = (TH1D*)(fin->Get("chi2_stat_periter"));
	TH1D* chi2_reg  = (TH1D*)(fin->Get("chi2_reg_periter"));
	Double_t chi2_tot_prefit   = chi2_tot  -> GetBinContent(1);
	Double_t chi2_tot_postfit  = chi2_tot  -> GetBinContent(chi2_tot  -> GetNbinsX()-1);
	Double_t chi2_stat_prefit  = chi2_stat -> GetBinContent(1);
	Double_t chi2_stat_postfit = chi2_stat -> GetBinContent(chi2_stat -> GetNbinsX()-1);
	Double_t chi2_reg_prefit   = chi2_reg  -> GetBinContent(1);
	Double_t chi2_reg_postfit  = chi2_reg  -> GetBinContent(chi2_reg  -> GetNbinsX()-1);
	//======================================================================================================




	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Draw comparison between data, nominal and post-fit distributions =====" << std::endl;
	std::cout << "===== For each subdetector (FGD1, FGD2x, FGD2y) =====" << std::endl;

	TCanvas* c_events[Nfgd];
	vector<TLegend*> leg(Nfgd);
	
	int ipad = 0;
	double ymax = 0;

	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		c_events[ifgd] = new TCanvas(Form("Nevents_in_%s",fgdNames[ifgd].c_str()),Form("Nevents_in_%s",fgdNames[ifgd].c_str()),1700,1000);
		c_events[ifgd] -> Divide(3,3);

		ipad = 0;
		for(int is=0; is<Nsample; is++)
		{
			ymax = TMath::Max(hSample_prefit[ifgd][is]->GetMaximum(), hSample_postfit[ifgd][is]->GetMaximum());
			ymax = TMath::Max(hSample_data[ifgd][is]->GetMaximum(), ymax);

			hSample_prefit[ifgd][is] -> GetYaxis()->SetRangeUser(0, ymax*1.3);
			hSample_prefit[ifgd][is] -> SetStats(false);

			hSample_prefit[ifgd][is] -> SetTitle(h_title[is].c_str());
			hSample_prefit[ifgd][is] -> GetXaxis() -> SetTitle("Analysis bins");
			hSample_prefit[ifgd][is] -> GetYaxis() -> SetTitle("# events");

			hSample_data[ifgd][is]    -> SetLineColor(kBlack);
			hSample_prefit[ifgd][is]  -> SetLineColor(kRed+2);
			hSample_postfit[ifgd][is] -> SetLineColor(kBlue+2);

			hSample_data[ifgd][is]    -> SetMarkerColor(kBlack);
			hSample_prefit[ifgd][is]  -> SetMarkerColor(kRed+2);
			hSample_postfit[ifgd][is] -> SetMarkerColor(kBlue+2);

			c_events[ifgd] -> cd(ipad+1);

			// TLine *verline[Nbins_costh];
			// int binSide = 0;
			// for(int il = 0; il < Nbins_costh; il++)
			// {
			// 	binSide += Nmombins[il];
			// 	verline[il]               = new TLine(binSide,         0, binSide,         1.3*ymax );
			// }

			hSample_prefit[ifgd][is] ->Draw();

			// for(int il = 0; il < Nbins_costh-1; il++)
			// {
			// 	verline[il] -> SetLineWidth(1);
			// 	verline[il] -> SetLineStyle(1);
			// 	verline[il] -> SetLineColor(kGray);
			// 	verline[il] -> Draw();
			// }

			hSample_postfit[ifgd][is]->Draw("SAME");
			hSample_data[ifgd][is]   ->Draw("E SAME"); // with error bars

			TPaveText *paveFGD = new TPaveText(0.6, 0.75, 0.9, 0.87, "NDC");
			paveFGD->SetFillColor(0);
			TText *t1=paveFGD->AddText(Form("#chi^{2}_{simpl.} = %d", (int)calculateChi2(hSample_data[ifgd][is], hSample_postfit[ifgd][is]) ));
			paveFGD -> Draw();
			c_events[ifgd] -> Update();

			if(ipad==4)
			{
				leg[ifgd] = new TLegend(0.15,0.24,0.95,0.64);
				leg[ifgd] -> SetFillColor(0);
				leg[ifgd] -> SetBorderSize(1);
				leg[ifgd] -> SetFillStyle(0);
				//leg[ifgd]->SetTextSize(0.075);
				if(inputname == "fit2_data" || inputname == "fit2_dataCS_fakedataSignal")
					leg[ifgd] -> AddEntry(hSample_data[ifgd][0],        "Data",        "lep");
				else
					leg[ifgd] -> AddEntry(hSample_data[ifgd][0],        "(Fake) data", "lep");
				leg[ifgd] -> AddEntry(hSample_prefit[ifgd][0], Form("Pre-fit  #chi^{2}_{tot} = %d",                                              (int)chi2_tot_prefit), "l");
				// leg[ifgd] -> AddEntry(hSample_postfit[ifgd][0],Form("Post-fit #chi^{2}_{tot} = %d,  #chi^{2}_{stat} = %d,  #chi^{2}_{reg} = %d", (int)chi2_tot_postfit, (int)chi2_stat_postfit, (int)chi2_reg_postfit),"l");
				leg[ifgd] -> AddEntry(hSample_postfit[ifgd][0],Form("Post-fit #chi^{2}_{tot} = %d,  #chi^{2}_{stat} = %d",                       (int)chi2_tot_postfit, (int)chi2_stat_postfit),"l");
				c_events[ifgd]->cd(6);
				leg[ifgd]->Draw();
			}

			if(ipad==4) ipad++;
			ipad++;
		}

		c_events[ifgd]->Print(Form("plots/fitteroutput/%s/eventDistr_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), fgdNames[ifgd].c_str()));
	}



	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Draw comparison between data, nominal and post-fit distributions =====" << std::endl;
	std::cout << "===== For all subdetector together (FGD1 + FGD2x + FGD2y) =====" << std::endl;

	TCanvas* c_events_allSubDet = new TCanvas("Nevents","Nevents",1700,1000);
	c_events_allSubDet -> Divide(3,3);
	TLegend* leg_allSubDet;
	
	ipad = 0;
	for(int is=0; is<Nsample; is++)
	{
		ymax = TMath::Max(hSample_prefit_allSubDet[is]->GetMaximum(), hSample_postfit_allSubDet[is]->GetMaximum());
		ymax = TMath::Max(hSample_data_allSubDet[is]->GetMaximum(), ymax);

		hSample_prefit_allSubDet[is] -> GetYaxis()->SetRangeUser(0, ymax*1.3);
		hSample_prefit_allSubDet[is] -> SetStats(false);

		hSample_prefit_allSubDet[is] -> SetTitle(h_title[is].c_str());
		hSample_prefit_allSubDet[is] -> GetXaxis() -> SetTitle("Analysis bins");
		hSample_prefit_allSubDet[is] -> GetYaxis() -> SetTitle("# events");

		hSample_data_allSubDet[is]    -> SetLineColor(kBlack);
		hSample_prefit_allSubDet[is]  -> SetLineColor(kRed+2);
		hSample_postfit_allSubDet[is] -> SetLineColor(kBlue+2);

		hSample_data_allSubDet[is]    -> SetMarkerColor(kBlack);
		hSample_prefit_allSubDet[is]  -> SetMarkerColor(kRed+2);
		hSample_postfit_allSubDet[is] -> SetMarkerColor(kBlue+2);

		c_events_allSubDet -> cd(ipad+1);

		// TLine *verline[Nbins_costh];
		// int binSide = 0;
		// for(int il = 0; il < Nbins_costh; il++)
		// {
		// 	binSide += Nmombins[il];
		// 	verline[il]               = new TLine(binSide,         0, binSide,         1.3*ymax );
		// }

		hSample_prefit_allSubDet[is] ->Draw();

		// for(int il = 0; il < Nbins_costh-1; il++)
		// {
		// 	verline[il] -> SetLineWidth(1);
		// 	verline[il] -> SetLineStyle(1);
		// 	verline[il] -> SetLineColor(kGray);
		// 	verline[il] -> Draw();
		// }
		
		hSample_postfit_allSubDet[is]->Draw("SAME");
		hSample_data_allSubDet[is]   ->Draw("E SAME"); // with error bars

		TPaveText *pave = new TPaveText(0.6, 0.75, 0.9, 0.87, "NDC");
		pave->SetFillColor(0);
		TText *t1=pave->AddText(Form("#chi^{2}_{simpl.} = %d", (int)calculateChi2(hSample_data_allSubDet[is], hSample_postfit_allSubDet[is]) ));
		pave->Draw();
		c_events_allSubDet->Update();

		if(ipad==4)
		{
			leg_allSubDet = new TLegend(0.15,0.24,0.95,0.64);
			leg_allSubDet -> SetFillColor(0);
			leg_allSubDet -> SetBorderSize(1);
			leg_allSubDet -> SetFillStyle(0);
			//leg_allSubDet->SetTextSize(0.075);
			if(inputname == "fit2_data" || inputname == "fit2_dataCS_fakedataSignal")
				leg_allSubDet -> AddEntry(hSample_data_allSubDet[0],   "Data",        "lep");
			else
				leg_allSubDet -> AddEntry(hSample_data_allSubDet[0],   "(Fake) data", "lep");
			leg_allSubDet -> AddEntry(hSample_prefit_allSubDet[0], Form("Pre-fit  #chi^{2}_{tot} = %d",                                              (int)chi2_tot_prefit), "l");
			// leg_allSubDet -> AddEntry(hSample_postfit_allSubDet[0],Form("Post-fit #chi^{2}_{tot} = %d,  #chi^{2}_{stat} = %d,  #chi^{2}_{reg} = %d", (int)chi2_tot_postfit, (int)chi2_stat_postfit, (int)chi2_reg_postfit),"l");
			leg_allSubDet -> AddEntry(hSample_postfit_allSubDet[0],Form("Post-fit #chi^{2}_{tot} = %d,  #chi^{2}_{stat} = %d",                       (int)chi2_tot_postfit, (int)chi2_stat_postfit),"l");
			c_events_allSubDet->cd(6);
			leg_allSubDet->Draw();
		}

		if(ipad==4) ipad++;
		ipad++;
	}

	c_events_allSubDet->Print(Form("plots/fitteroutput/%s/eventDistr_%s_allSubDet.pdf", dir_name.c_str(), inputname.c_str()));
	



	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Draw comparison between data, nominal and post-fit distributions =====" << std::endl;
	std::cout << "===== For all subdetector together (FGD1 + FGD2x + FGD2y) =====" << std::endl;
	std::cout << "===== With muTPC grouped regions =====" << std::endl;

	TCanvas* c_events_allSubDet_grouped = new TCanvas("NeventsGrouped","NeventsGrouped",1700,1000);
	c_events_allSubDet_grouped -> Divide(3,3);
	TLegend* leg_allSubDet_grouped;
	
	ipad = 0;

	int is=0;
	for(int ISam=0; ISam<Nsample; ISam++)
	{
		if(ISam==0) is = Nsample;
		else        is = ISam;

		ymax = TMath::Max(hSample_prefit_allSubDet[is]->GetMaximum(), hSample_postfit_allSubDet[is]->GetMaximum());
		ymax = TMath::Max(hSample_data_allSubDet[is]->GetMaximum(), ymax);

		hSample_prefit_allSubDet[is] -> GetYaxis()->SetRangeUser(0, ymax*1.3);
		hSample_prefit_allSubDet[is] -> SetStats(false);

		hSample_prefit_allSubDet[is] -> SetTitle(h_title[is].c_str());
		hSample_prefit_allSubDet[is] -> GetXaxis() -> SetTitle("Analysis bins");
		hSample_prefit_allSubDet[is] -> GetYaxis() -> SetTitle("# events");

		hSample_data_allSubDet[is]    -> SetLineColor(kBlack);
		hSample_prefit_allSubDet[is]  -> SetLineColor(kRed+2);
		hSample_postfit_allSubDet[is] -> SetLineColor(kBlue+2);

		hSample_data_allSubDet[is]    -> SetMarkerColor(kBlack);
		hSample_prefit_allSubDet[is]  -> SetMarkerColor(kRed+2);
		hSample_postfit_allSubDet[is] -> SetMarkerColor(kBlue+2);

		c_events_allSubDet_grouped -> cd(ipad+1);

		// TLine *verline[Nbins_costh];
		// int binSide = 0;
		// for(int il = 0; il < Nbins_costh; il++)
		// {
		// 	binSide += Nmombins[il];
		// 	verline[il]               = new TLine(binSide,         0, binSide,         1.3*ymax );
		// }

		hSample_prefit_allSubDet[is] ->Draw();

		// for(int il = 0; il < Nbins_costh-1; il++)
		// {
		// 	verline[il] -> SetLineWidth(1);
		// 	verline[il] -> SetLineStyle(1);
		// 	verline[il] -> SetLineColor(kGray);
		// 	verline[il] -> Draw();
		// }

		hSample_postfit_allSubDet[is]->Draw("SAME");
		hSample_data_allSubDet[is]   ->Draw("E SAME"); // with error bars

		TPaveText *pave = new TPaveText(0.6, 0.75, 0.9, 0.87, "NDC");
		pave->SetFillColor(0);
		TText *t1=pave->AddText(Form("#chi^{2}_{simpl.} = %d", (int)calculateChi2(hSample_data_allSubDet[is], hSample_postfit_allSubDet[is]) ));
		pave->Draw();
		c_events_allSubDet_grouped->Update();

		if(ISam==0) ISam = 2;
		ipad++;
	}

	leg_allSubDet_grouped = new TLegend(0.15,0.24,0.95,0.64);
	leg_allSubDet_grouped -> SetFillColor(0);
	leg_allSubDet_grouped -> SetBorderSize(1);
	leg_allSubDet_grouped -> SetFillStyle(0);
	//leg_allSubDet_grouped->SetTextSize(0.075);
	if(inputname == "fit2_data" || inputname == "fit2_dataCS_fakedataSignal")
		leg_allSubDet_grouped -> AddEntry(hSample_data_allSubDet[0],   "Data",    "lep");
	else
		leg_allSubDet_grouped -> AddEntry(hSample_data_allSubDet[0],   "(Fake) data",    "lep");
	leg_allSubDet_grouped -> AddEntry(hSample_prefit_allSubDet[0], Form("Pre-fit  #chi^{2}_{tot} = %d",                                              (int)chi2_tot_prefit), "l");
	// leg_allSubDet_grouped -> AddEntry(hSample_postfit_allSubDet[0],Form("Post-fit #chi^{2}_{tot} = %d,  #chi^{2}_{stat} = %d,  #chi^{2}_{reg} = %d", (int)chi2_tot_postfit, (int)chi2_stat_postfit, (int)chi2_reg_postfit),"l");
	leg_allSubDet_grouped -> AddEntry(hSample_postfit_allSubDet[0],Form("Post-fit #chi^{2}_{tot} = %d,  #chi^{2}_{stat} = %d",                       (int)chi2_tot_postfit, (int)chi2_stat_postfit),"l");
	c_events_allSubDet_grouped->cd(7);
	leg_allSubDet_grouped->Draw();

	c_events_allSubDet_grouped->Print(Form("plots/fitteroutput/%s/eventDistr_%s_allSubDet_muTPCgrouped.pdf", dir_name.c_str(), inputname.c_str()));
	

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}



double calculateChi2(TH1D* hist1, TH1D* hist2)
{
	double chi2 = 0.0;

	for(int i=1; i<=hist1->GetNbinsX(); i++)
		if(hist1->GetBinContent(i) > 0)
			chi2 = chi2 + (hist1->GetBinContent(i) - hist2->GetBinContent(i))*(hist1->GetBinContent(i) - hist2->GetBinContent(i))/ (hist1->GetBinContent(i));
		else
			std::cout << "WARNING : bin " << i << " has 0 content !!!" << std::endl;
	return chi2;
}