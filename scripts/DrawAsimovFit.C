//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root 'DrawAsimovFit.C+()'
//         root 'DrawAsimovFit.C+("fit1_statFluc")'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

void DrawAsimovFit(string inputname = "fit1_asimov", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/%s.root", inputname.c_str());

	//======================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
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
	{
		for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
		{
			mombins[nbcth] = new double[nbm];
		}
	}
	mombins = bin.GetMomBins();

	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store fitter output in a TFile =====" << std::endl;

	TFile* fin = new TFile(infilename.c_str());
	
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store sample name in a vector of string =====" << std::endl;
	
	const int Nfgd = 2;
	int Nsample[Nfgd] = {5, 10};
	vector< vector<string> > nameSamples(Nfgd);

	for(int is = 0; is < Nsample[0]; is++)
	{	// add FGD1 distributions
		nameSamples[0].push_back(Form("evhist_sam%d", is));
		std::cout << "add sample " << is << " to FGD1" << std::endl;
	}
	for(int is = 0; is < Nsample[0]; is++)
	{	// add FGD2x distributions
		nameSamples[1].push_back(Form("evhist_sam%d", is+8));
		std::cout << "add sample " << is+8 << " to FGD2" << std::endl;
	}
	for(int is = 0; is < Nsample[0]; is++)
	{	// add FGD2y distributions
		nameSamples[1].push_back(Form("evhist_sam%d", is+16));
		std::cout << "add sample " << is+16 << " to FGD2" << std::endl;
	}

	//======================================================================================================


	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TH1D* htemp_C;
	TH1D* htemp_O;

	//======================================================================================================  
	std::cout << "===== Get histos nominal and post-fit =====" << std::endl;
	
	vector< vector< TH1D* > > hSample_nom_C(Nfgd); 
	vector< vector< TH1D* > > hSample_nom_O(Nfgd); 
	vector< vector< TH1D* > > hSample_fit_C(Nfgd);
	vector< vector< TH1D* > > hSample_fit_O(Nfgd);

	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		for(int ns=0; ns<Nsample[ifgd]; ns++)
		{
			htemp_C = (TH1D*)(fin->Get(Form("%s_iter0_pred_C", nameSamples[ifgd][ns].c_str())));
			htemp_O = (TH1D*)(fin->Get(Form("%s_iter0_pred_O", nameSamples[ifgd][ns].c_str())));
			hSample_nom_C[ifgd].push_back(htemp_C);
			hSample_nom_O[ifgd].push_back(htemp_O);

			htemp_C = (TH1D*)(fin->Get(Form("%s_finaliter_pred_C", nameSamples[ifgd][ns].c_str())));
			htemp_O = (TH1D*)(fin->Get(Form("%s_finaliter_pred_O", nameSamples[ifgd][ns].c_str())));
			hSample_fit_C[ifgd].push_back(htemp_C);
			hSample_fit_O[ifgd].push_back(htemp_O);
		}
	}

	//=== Add all the histos post-fit
	TH1D* hSignal_nom_C[Nfgd];
	TH1D* hSignal_nom_O[Nfgd];
	TH1D* hSignal_fit_C[Nfgd];
	TH1D* hSignal_fit_O[Nfgd];
	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		hSignal_nom_C[ifgd] = (TH1D*)(hSample_nom_C[ifgd][0]->Clone( Form("hSignal_nom_C[%d]", ifgd) ));
		hSignal_nom_O[ifgd] = (TH1D*)(hSample_nom_O[ifgd][0]->Clone( Form("hSignal_nom_O[%d]", ifgd) ));
		hSignal_fit_C[ifgd] = (TH1D*)(hSample_fit_C[ifgd][0]->Clone( Form("hSignal_fit_C[%d]", ifgd) ));
		hSignal_fit_O[ifgd] = (TH1D*)(hSample_fit_O[ifgd][0]->Clone( Form("hSignal_fit_O[%d]", ifgd) ));
		
		for(int ns=1; ns<Nsample[ifgd]; ns++)
		{
			hSignal_nom_C[ifgd] -> Add(hSample_nom_C[ifgd][ns]);
			hSignal_nom_O[ifgd] -> Add(hSample_nom_O[ifgd][ns]);
			hSignal_fit_C[ifgd] -> Add(hSample_fit_C[ifgd][ns]);
			hSignal_fit_O[ifgd] -> Add(hSample_fit_O[ifgd][ns]);
		}
	}
	//======================================================================================================






	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Draw comparison between nominal and post-fit distribution =====" << std::endl;
	
	//======================================================================================================
	std::cout << "=== Set binning with a temp histogram" << std::endl;

	vector<TH1D*> h_tmp;
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

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_tmp.push_back(new TH1D("","", Nmombins[nbcth], mombins[nbcth]));
		h_tmp[nbcth] -> SetTitle( Form("%s", h_title[nbcth].Data()) );
		h_tmp[nbcth] -> GetYaxis()->SetTitle("CC-0#pi events/100 MeV");
		// h_tmp[nbcth] -> GetXaxis()->SetTitle("p^{#mu}_{true} [GeV/c]");
		h_tmp[nbcth] -> GetXaxis()->SetTitle("p^{#mu}_{reco} [GeV/c]");
	}

	vector< vector<TH1D*> > h_nom_C(Nfgd);
	vector< vector<TH1D*> > h_nom_O(Nfgd);
	vector< vector<TH1D*> > h_fit_C(Nfgd);
	vector< vector<TH1D*> > h_fit_O(Nfgd);
	
	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_nom_C[ifgd].push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_nom_C_cth_%d",nbcth)));
			h_nom_O[ifgd].push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_nom_O_cth_%d",nbcth)));
			h_fit_C[ifgd].push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_fit_C_cth_%d",nbcth)));
			h_fit_O[ifgd].push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_fit_O_cth_%d",nbcth)));
		}
	}

	//======================================================================================================
	std::cout << "=== Fill nominal and post-fit histos" << std::endl;

	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		int bincounter = 0;
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
			{
				// Set bin content
				double binwidth         = ((mombins[nbcth][nbm+1] - mombins[nbcth][nbm])/0.1);
				double binContent_nom_C = (hSignal_nom_C[ifgd] -> GetBinContent(bincounter+1)/binwidth);
				double binContent_nom_O = (hSignal_nom_O[ifgd] -> GetBinContent(bincounter+1)/binwidth);
				double binContent_fit_C = (hSignal_fit_C[ifgd] -> GetBinContent(bincounter+1)/binwidth);
				double binContent_fit_O = (hSignal_fit_O[ifgd] -> GetBinContent(bincounter+1)/binwidth);
				h_nom_C[ifgd][nbcth] -> SetBinContent(nbm+1, binContent_nom_C);
				h_nom_O[ifgd][nbcth] -> SetBinContent(nbm+1, binContent_nom_O);
				h_fit_C[ifgd][nbcth] -> SetBinContent(nbm+1, binContent_fit_C);
				h_fit_O[ifgd][nbcth] -> SetBinContent(nbm+1, binContent_fit_O);
				
				// Set bin error
				double binError_fit_C   = (hSignal_fit_C[ifgd] -> GetBinError(bincounter+1)/binwidth);
				double binError_fit_O   = (hSignal_fit_O[ifgd] -> GetBinError(bincounter+1)/binwidth);
				h_fit_C[ifgd][nbcth] -> SetBinError(nbm+1, binError_fit_C);
				h_fit_O[ifgd][nbcth] -> SetBinError(nbm+1, binError_fit_O);
				bincounter++;
			}
		}
	}
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	TCanvas* c_asimov[Nfgd];
	vector<TLegend*> leg(Nfgd);

	for(int ifgd = 0; ifgd < Nfgd; ifgd++)
	{
		c_asimov[ifgd] = new TCanvas(Form("AsimovFit_in_FGD%d",ifgd+1),Form("AsimovFit_in_FGD%d", ifgd+1),1700,1000);
		c_asimov[ifgd] -> Divide(3,3);

		double* ymax = new double[Nbins_costh];

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			ymax[nbcth] = TMath::Max(h_nom_C[ifgd][nbcth]->GetMaximum(), h_nom_O[ifgd][nbcth]->GetMaximum());

			h_nom_C[ifgd][nbcth] -> GetYaxis()->SetRangeUser(0, ymax[nbcth]*1.3);
			h_nom_C[ifgd][nbcth] -> SetStats(false);

			h_nom_C[ifgd][nbcth] -> SetLineColor(kBlue);
			h_fit_C[ifgd][nbcth] -> SetLineColor(kBlue+3);
			h_fit_C[ifgd][nbcth] -> SetMarkerColor(kBlue+3);
			h_nom_O[ifgd][nbcth] -> SetLineColor(kRed);
			h_fit_O[ifgd][nbcth] -> SetLineColor(kRed+3);
			h_fit_O[ifgd][nbcth] -> SetMarkerColor(kRed+3);

			c_asimov[ifgd] -> cd(nbcth+1);

			if(nbcth!=0) gPad->SetLogx();

			h_nom_C[ifgd][nbcth]->Draw();
			h_nom_O[ifgd][nbcth]->Draw("SAME");
			h_fit_C[ifgd][nbcth]->Draw("E SAME"); // with error bars
			h_fit_O[ifgd][nbcth]->Draw("E SAME"); // with error bars

			if(nbcth==0)
			{
				leg[ifgd] = new TLegend(0.2,0.24,0.8,0.64);
				leg[ifgd] -> SetFillColor(0);
				leg[ifgd] -> SetBorderSize(1);
				leg[ifgd] -> SetFillStyle(0);
				//leg[ifgd]->SetTextSize(0.075);
				leg[ifgd] -> AddEntry(h_nom_C[ifgd][0],"Nominal distribution, carbon","l");
				leg[ifgd] -> AddEntry(h_nom_O[ifgd][0],"Nominal distribution, oxygen","l");
				leg[ifgd] -> AddEntry(h_fit_C[ifgd][0],"Post fit distribution, carbon","p");
				leg[ifgd] -> AddEntry(h_fit_O[ifgd][0],"Post fit distribution, oxygen","p");
				c_asimov[ifgd]->cd(1);
				leg[ifgd]->Draw();
			}
		}

		c_asimov[ifgd]->Print(Form("plots/fitteroutput/asimov/eventDistr_%s_FGD%d.pdf", inputname.c_str(), ifgd+1));
	}

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}


