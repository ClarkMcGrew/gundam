//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root 'DrawEventComparison.C+()'
//         root 'DrawEventComparison.C+("fit1_statFluc")'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

void DrawEventComparison(string inputname = "fit1_asimov")
{

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/%s.root", inputname.c_str());

	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store fitter output in a TFile =====" << std::endl;

	TFile* fin = new TFile(infilename.c_str());
	
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store sample name in a vector of string =====" << std::endl;
	
	const int Nfgd = 3;
	string fgdNames[Nfgd] = {"FGD1", "FGD2x", "FGD2y"};

	const int Nsample = 8;
	vector< vector<string> > nameSamples(Nfgd);

	for(int is = 0; is < Nsample; is++)// add FGD1 distributions
		nameSamples[0].push_back(Form("evhist_sam%d", is));

	for(int is = 0; is < Nsample; is++)// add FGD2x distributions
		nameSamples[1].push_back(Form("evhist_sam%d", is+8));
	
	for(int is = 0; is < Nsample; is++)// add FGD2y distributions
		nameSamples[2].push_back(Form("evhist_sam%d", is+16));

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
	{
		for(int ns=0; ns<Nsample; ns++)
		{
			htemp_prefit = (TH1D*)(fin->Get(Form("%s_iter0_pred", nameSamples[ifgd][ns].c_str())));
			hSample_prefit[ifgd].push_back(htemp_prefit);

			htemp_postfit = (TH1D*)(fin->Get(Form("%s_finaliter_pred", nameSamples[ifgd][ns].c_str())));
			hSample_postfit[ifgd].push_back(htemp_postfit);

			htemp_data = (TH1D*)(fin->Get(Form("%s_iter0_data", nameSamples[ifgd][ns].c_str())));
			hSample_data[ifgd].push_back(htemp_data);
		}
	}

	//======================================================================================================





	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Draw comparison between data, nominal and post-fit distributions =====" << std::endl;

	const int Nbins = hSample_prefit[0][0] -> GetNbinsX();
	string h_title[Nsample] = {"muTPC",
	                           "muTPCpTPC",
	                           "muTPCpFGD",
	                           "muFGDpTPC",
	                           "muFGD",
	                           "CC1pi",
	                           "CCOther",
	                           "CCMichel"
	                           };

	TCanvas* c_events[Nfgd];
	vector<TLegend*> leg(Nfgd);

	double ymax = 0;

	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		c_events[ifgd] = new TCanvas(Form("Nevents_in_%s",fgdNames[ifgd].c_str()),Form("Nevents_in_%s",fgdNames[ifgd].c_str()),1700,1000);
		c_events[ifgd] -> Divide(3,3);

		int ipad = 0;
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
			hSample_prefit[ifgd][is]  -> SetLineColor(kRed);
			hSample_postfit[ifgd][is] -> SetLineColor(kBlue);

			hSample_data[ifgd][is]    -> SetMarkerColor(kBlack);
			hSample_prefit[ifgd][is]  -> SetMarkerColor(kRed);
			hSample_postfit[ifgd][is] -> SetMarkerColor(kBlue);

			c_events[ifgd] -> cd(ipad+1);

			hSample_prefit[ifgd][is] ->Draw();
			hSample_postfit[ifgd][is]->Draw("SAME");
			hSample_data[ifgd][is]   ->Draw("E SAME"); // with error bars

			if(ipad==4)
			{
				leg[ifgd] = new TLegend(0.2,0.24,0.8,0.64);
				leg[ifgd] -> SetFillColor(0);
				leg[ifgd] -> SetBorderSize(1);
				leg[ifgd] -> SetFillStyle(0);
				//leg[ifgd]->SetTextSize(0.075);
				leg[ifgd] -> AddEntry(hSample_data[ifgd][0],   "Data distribution","p");
				leg[ifgd] -> AddEntry(hSample_prefit[ifgd][0], "Pre-fit distribution","l");
				leg[ifgd] -> AddEntry(hSample_postfit[ifgd][0],"Post-fit distribution","l");
				c_events[ifgd]->cd(6);
				leg[ifgd]->Draw();
			}

			if(ipad==4) ipad++;
			ipad++;
		}

		c_events[ifgd]->Print(Form("plots/fitteroutput/asimov/eventDistr_%s_%s.pdf", inputname.c_str(), fgdNames[ifgd].c_str()));
	}

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}


