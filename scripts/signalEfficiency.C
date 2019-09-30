//  Usage:
// root -b -q 'signalEfficiency.C+("C")'
// root -b -q 'signalEfficiency.C+("O")'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int GetTrueEvents(string, TTree*, std::vector< TH1D* >, double*, int);
int GetSelectedEvents(string, TTree*, std::vector< TH1D* >, double*, int);
int FindCosthBin(double, double*, int);


void signalEfficiency(string target = "C", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt")
{

	if(target != "C" && target != "O")
		std::cout << "ERROR : First argument must be C (for carbon) or O (for oxygen) !!!" << std::endl;

	string infilename_neut  = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_neut_prod6B/output_run2348_mc_500toys.root";
	string infilename_genie = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_genie_prod6B/output_run2a_mc_500toys.root";
	// string infilename_nuwro = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_nuwro/output_nuwro_mc_all_500toys.root";
	string infilename_nuwro = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_nuwro/output_nuwro_mc_new.root";


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
	TFile* fin_genie = new TFile(infilename_genie.c_str());
	TFile* fin_nuwro = new TFile(infilename_nuwro.c_str());
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

	std::vector< TH1D* > h_sel_neut;
	std::vector< TH1D* > h_sel_genie;
	std::vector< TH1D* > h_sel_nuwro;
	std::vector< TH1D* > h_tru_neut;
	std::vector< TH1D* > h_tru_genie;
	std::vector< TH1D* > h_tru_nuwro;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_sel_neut.push_back((TH1D*)h_tmp[nbcth]  -> Clone(Form("h_sel_neut_%d",nbcth)));
		h_sel_genie.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_sel_genie_%d",nbcth)));
		h_sel_nuwro.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_sel_nuwro_%d",nbcth)));
		h_tru_neut.push_back((TH1D*)h_tmp[nbcth]  -> Clone(Form("h_tru_neut_%d",nbcth)));
		h_tru_genie.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_tru_genie_%d",nbcth)));
		h_tru_nuwro.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_tru_nuwro_%d",nbcth)));
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TTree * truth_tree_neut  = (TTree*) fin_neut  -> Get("truth");
	TTree * truth_tree_genie = (TTree*) fin_genie -> Get("truth");
	TTree * truth_tree_nuwro = (TTree*) fin_nuwro -> Get("truth");
	TTree * selec_tree_neut  = (TTree*) fin_neut  -> Get("default");
	TTree * selec_tree_genie = (TTree*) fin_genie -> Get("default");
	TTree * selec_tree_nuwro = (TTree*) fin_nuwro -> Get("default");

	// neut
	// genie
	// nuwro

	int NselEvents_neut, NtruEvents_neut;
	int NselEvents_genie, NtruEvents_genie;
	int NselEvents_nuwro, NtruEvents_nuwro;

	std::cout << "Get NEUT events..." << std::endl;
	NtruEvents_neut = GetTrueEvents(target, truth_tree_neut,     h_tru_neut, costhBins, Nbins_costh);
	NselEvents_neut = GetSelectedEvents(target, selec_tree_neut, h_sel_neut, costhBins, Nbins_costh);
	std::cout << "--> NEUT signal efficiency in "<<target<<" = " << (1.0*NselEvents_neut)/(1.0*NtruEvents_neut) << std::endl;

	std::cout << "Get GENIE events..." << std::endl;
	NtruEvents_genie = GetTrueEvents(target, truth_tree_genie,     h_tru_genie, costhBins, Nbins_costh);
	NselEvents_genie = GetSelectedEvents(target, selec_tree_genie, h_sel_genie, costhBins, Nbins_costh);
	std::cout << "--> GENIE signal efficiency in "<<target<<" = " << (1.0*NselEvents_genie)/(1.0*NtruEvents_genie) << std::endl;

	std::cout << "Get NuWro events..." << std::endl;
	NtruEvents_nuwro = GetTrueEvents(target, truth_tree_nuwro,     h_tru_nuwro, costhBins, Nbins_costh);
	NselEvents_nuwro = GetSelectedEvents(target, selec_tree_nuwro, h_sel_nuwro, costhBins, Nbins_costh);
	std::cout << "--> NuWro signal efficiency in "<<target<<" = " << (1.0*NselEvents_nuwro)/(1.0*NtruEvents_nuwro) << std::endl;

	//======================================================================================================



	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute the efficiency" << std::endl;

	std::vector< TH1D* > h_eff_neut;
	std::vector< TH1D* > h_eff_genie;
	std::vector< TH1D* > h_eff_nuwro;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut.push_back((TH1D*) h_sel_neut[nbcth]  -> Clone(Form("h_eff_neut_%d",nbcth)));
		h_eff_genie.push_back((TH1D*)h_sel_genie[nbcth] -> Clone(Form("h_eff_genie_%d",nbcth)));
		h_eff_nuwro.push_back((TH1D*)h_sel_nuwro[nbcth] -> Clone(Form("h_eff_nuwro_%d",nbcth)));
	}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut[nbcth]  -> Divide(h_tru_neut[nbcth]);
		h_eff_genie[nbcth] -> Divide(h_tru_genie[nbcth]);
		h_eff_nuwro[nbcth] -> Divide(h_tru_nuwro[nbcth]);
	}
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      


	TCanvas* c_eff = new TCanvas("efficiency", "efficiency",1700,1000);
	c_eff -> Divide(3,3);

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_neut[nbcth] -> SetTitle(h_title[nbcth]);
		h_eff_neut[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_eff_neut[nbcth] -> GetYaxis() -> SetTitle("Signal efficiency");
		h_eff_neut[nbcth] -> GetYaxis() -> SetRangeUser(0, 1.0);

		h_eff_neut[nbcth]  -> SetLineColor(kBlue+1);
		h_eff_genie[nbcth] -> SetLineColor(kRed+1);
		h_eff_nuwro[nbcth] -> SetLineColor(kGreen+2);

		h_eff_neut[nbcth]  -> SetLineWidth(2);
		h_eff_genie[nbcth] -> SetLineWidth(2);
		h_eff_nuwro[nbcth] -> SetLineWidth(2);

		c_eff -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_eff_neut[nbcth]  -> Draw("hist");
		h_eff_genie[nbcth] -> Draw("same hist");
		h_eff_nuwro[nbcth] -> Draw("same hist");

		if(nbcth==0)
		{
			TLegend* leg = new TLegend(0.2,0.40,0.8,0.85);
			leg -> SetFillColor(0);
			leg -> SetBorderSize(1);
			leg -> SetFillStyle(0);
			leg -> AddEntry(h_eff_neut[0],  "NEUT",  "l");
			leg -> AddEntry(h_eff_genie[0], "GENIE", "l");
			leg -> AddEntry(h_eff_nuwro[0], "NuWro", "l");
			c_eff->cd(1);
			leg->Draw();
		}
	}
	c_eff -> Print(Form("plots/efficiencies/efficiency_generators_%s.pdf", target.c_str()));



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
