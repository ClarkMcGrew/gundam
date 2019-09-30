//  Usage:
// root -b -q 'signalEfficiencyRatio.C+()'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int GetTrueEvents(string, TTree*, std::vector< TH1D* >, double*, int);
int GetSelectedEvents(string, TTree*, std::vector< TH1D* >, double*, int);
int FindCosthBin(double, double*, int);


void signalEfficiencyRatio(string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt")
{

	string infilename_neut  = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_neut_prod6B/output_run2348_mc_500toys.root";
	// string infilename_genie = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_genie_prod6B/output_run2a_mc_500toys.root";
	// string infilename_nuwro = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_nuwro/output_nuwro_mc_all_500toys.root";


	//======================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
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
	{
		for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
		{
			mombins[nbcth] = new double[nbm];
		}
	}
	mombins = bin.GetMomBins();


	double* costhBins = new double[Nbins_costh];
	costhBins = bin.GetCosthBins();

	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store Highland2 outputs in a TFile =====" << std::endl;

	TFile* fin_neut  = new TFile(infilename_neut.c_str());
	// TFile* fin_genie = new TFile(infilename_genie.c_str());
	// TFile* fin_nuwro = new TFile(infilename_nuwro.c_str());
	//======================================================================================================


	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Initialise histograms with binning" << std::endl;

	std::vector< TH1D* > h_tmp;

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
		// h_tmp[nbcth] -> GetYaxis()->SetTitle("CC-0#pi events/100 MeV");
		// h_tmp[nbcth] -> GetXaxis()->SetTitle("p^{#mu}_{true} [GeV/c]");
		// h_tmp[nbcth] -> GetXaxis()->SetTitle("p^{#mu}_{reco} [GeV/c]");
	}

	std::vector< TH1D* > h_sel_neut_C;
	std::vector< TH1D* > h_sel_neut_O;
	std::vector< TH1D* > h_tru_neut_C;
	std::vector< TH1D* > h_tru_neut_O;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_sel_neut_C.push_back((TH1D*)h_tmp[nbcth]  -> Clone(Form("h_sel_neut_C_%d",nbcth)));
		h_sel_neut_O.push_back((TH1D*)h_tmp[nbcth]  -> Clone(Form("h_sel_neut_O_%d",nbcth)));
		h_tru_neut_C.push_back((TH1D*)h_tmp[nbcth]  -> Clone(Form("h_tru_neut_C_%d",nbcth)));
		h_tru_neut_O.push_back((TH1D*)h_tmp[nbcth]  -> Clone(Form("h_tru_neut_O_%d",nbcth)));
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TTree * truth_tree_neut  = (TTree*) fin_neut  -> Get("truth");
	TTree * selec_tree_neut  = (TTree*) fin_neut  -> Get("default");

	int NselEvents_neut_C;
	int NselEvents_neut_O;
	int NtruEvents_neut_C;
	int NtruEvents_neut_O;

	std::cout << "Get NEUT events..." << std::endl;
	NtruEvents_neut_C = GetTrueEvents("C", truth_tree_neut,     h_tru_neut_C, costhBins, Nbins_costh);
	NtruEvents_neut_O = GetTrueEvents("O", truth_tree_neut,     h_tru_neut_O, costhBins, Nbins_costh);
	NselEvents_neut_C = GetSelectedEvents("C", selec_tree_neut, h_sel_neut_C, costhBins, Nbins_costh);
	NselEvents_neut_O = GetSelectedEvents("O", selec_tree_neut, h_sel_neut_O, costhBins, Nbins_costh);
	std::cout << "--> NEUT signal efficiency in C = " << (1.0*NselEvents_neut_C)/(1.0*NtruEvents_neut_C) << std::endl;
	std::cout << "--> NEUT signal efficiency in O = " << (1.0*NselEvents_neut_O)/(1.0*NtruEvents_neut_O) << std::endl;

	//======================================================================================================



	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute the efficiency" << std::endl;

	std::vector< TH1D* > h_eff_neut_C;
	std::vector< TH1D* > h_eff_neut_O;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut_C.push_back((TH1D*) h_sel_neut_C[nbcth]  -> Clone(Form("h_eff_neut_C_%d",nbcth)));
		h_eff_neut_O.push_back((TH1D*) h_sel_neut_O[nbcth]  -> Clone(Form("h_eff_neut_O_%d",nbcth)));
	}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut_C[nbcth]  -> Divide(h_tru_neut_C[nbcth]);
		h_eff_neut_O[nbcth]  -> Divide(h_tru_neut_O[nbcth]);
	}
	//======================================================================================================

	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute the efficiency RATIO" << std::endl;

	std::vector< TH1D* > h_eff_neut_ratio;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut_ratio.push_back((TH1D*) h_eff_neut_O[nbcth]  -> Clone(Form("h_eff_neut_ratio_%d",nbcth)));
	}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut_ratio[nbcth]  -> Divide(h_eff_neut_C[nbcth]);
	}
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      


	TCanvas* c_eff = new TCanvas("efficiency", "efficiency",1700,1000);
	c_eff -> Divide(3,3);

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut_ratio[nbcth] -> SetTitle(h_title[nbcth]);
		h_eff_neut_ratio[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_eff_neut_ratio[nbcth] -> GetYaxis() -> SetTitle("Signal efficiency ratio O/C");
		h_eff_neut_ratio[nbcth] -> GetYaxis() -> SetRangeUser(0, 2.0);

		h_eff_neut_ratio[nbcth]  -> SetLineColor(kBlue+1);
		h_eff_neut_ratio[nbcth]  -> SetLineWidth(2);

		c_eff -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_eff_neut_ratio[nbcth]  -> Draw("hist");

		// if(nbcth==0)
		// {
		// 	TLegend* leg = new TLegend(0.2,0.40,0.8,0.85);
		// 	leg -> SetFillColor(0);
		// 	leg -> SetBorderSize(1);
		// 	leg -> SetFillStyle(0);
		// 	leg -> AddEntry(h_eff_neut_ratio[0],  "NEUT",  "l");
		// 	c_eff->cd(1);
		// 	leg->Draw();
		// }
	}
	c_eff -> Print("plots/efficiencies/efficiency_ratio_NEUT.pdf");



	//======================================================================================================

}



// Get true signal events from truth tree
int GetTrueEvents(string target, TTree* tree, std::vector< TH1D* > h_events, double* CosthBins, int Nbins_costh)
{
	int fgdtarget, fgd2fgdtarget;
	float mom, costh;
	int costhbin;
	int signal=999;

	tree -> SetBranchAddress("truelepton_mom",       &mom);
	tree -> SetBranchAddress("truelepton_costheta",  &costh);
	tree -> SetBranchAddress("fgdtargetCCZeroPi",    &fgdtarget);
	tree -> SetBranchAddress("fgd2fgdtargetCCZeroPi",&fgd2fgdtarget);

	if(target=="C")      signal = 3;
	else if(target=="O") signal = 0;

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);
		int costhbin = FindCosthBin(costh, CosthBins, Nbins_costh);

		if(fgdtarget == signal || fgd2fgdtarget == signal ) // if numu CC0pi event on carbon or oxygen
		{
			h_events[costhbin] -> Fill(mom);
			Nevents++;
		}

	}
	std::cout << "Number of true signal events : " << Nevents << std::endl;
	return Nevents;
}

// Get selected events from default tree
int GetSelectedEvents(string target, TTree* tree, std::vector< TH1D* > h_events, double* CosthBins, int Nbins_costh)
{
	int sample, fgdtarget, fgd2fgdtarget;
	float mom, costh, weight;
	int costhbin;
	int signal=999;

	tree -> SetBranchAddress("sample_clst_fgd2layer_xsec", &sample);
	tree -> SetBranchAddress("truelepton_mom",             &mom);
	tree -> SetBranchAddress("truelepton_costheta",        &costh);
	tree -> SetBranchAddress("weight_corr_total",          &weight);
	tree -> SetBranchAddress("fgdtargetCCZeroPi",          &fgdtarget);
	tree -> SetBranchAddress("fgd2fgdtargetCCZeroPi",      &fgd2fgdtarget);

	if(target=="C")      signal = 3;
	else if(target=="O") signal = 0;

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		int costhbin = FindCosthBin(costh, CosthBins, Nbins_costh);

		if( (sample==1  || sample==2  || sample==3  || sample==4  || sample==5  || sample==9  || sample==10 || sample==11 ||
		     sample==12 || sample==13 || sample==14 || sample==15 || sample==16 || sample==20 || sample==21 || sample==22 ||
		     sample==23 || sample==24 || sample== 25|| sample==26 || sample==27 || sample==31 || sample==32 || sample==33 ) // if selected event
		    && (fgdtarget == signal || fgd2fgdtarget == signal) // if signal on carbon or oxygen
		    && weight>0.0 && weight<10.0 )
		{
			h_events[costhbin] -> Fill(mom, weight);
			Nevents++;
		}
	}
	std::cout << "Number of selected events : " << Nevents << std::endl;
	return Nevents;
}


// Find costheta bin index of a costh value
int FindCosthBin(double Costh, double* CosthBins, int Nbins_costh)
{
	int finalBin=999;
	for(int i=0; i<Nbins_costh; i++)
		if(Costh > CosthBins[i] && Costh <= CosthBins[i+1])
			finalBin = i;

	return finalBin;
}
