//  Usage:
// source /usr/local/root/pro/bin/thisroot.sh
// root -b -q 'signalEfficiency2D_OandC.C+("C")'
// root -b -q 'signalEfficiency2D_OandC.C+("O")'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int GetTrueEvents(string, TTree*, TH2D*);
int GetSelectedEvents(string, TTree*, TH2D*);


void signalEfficiency2D_OandC(string target)
{

	if(target != "C" && target != "O")
		std::cout << "ERROR : First argument must be C (for carbon) or O (for oxygen) !!!" << std::endl;

	string infilename_neut  = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_neut_prod6B/output_run2348_mc_500toys.root";
	

	//======================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();
	gROOT->ForceStyle();

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set binning  =====" << std::endl;

	int Nbins = 50;
	double pmu_min =  0.0;
	double pmu_max =  5.0;
	double cth_min = -1.0;
	double cth_max =  1.0;
	
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store Highland2 outputs in a TFile =====" << std::endl;

	TFile* fin_neut  = new TFile(infilename_neut.c_str());
	//======================================================================================================


	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Initialise histograms with binning" << std::endl;

	TH2D* h_sel_neut = new TH2D("h_sel_neut", "h_sel_neut", Nbins, cth_min, cth_max, Nbins, pmu_min, pmu_max);
	TH2D* h_tru_neut = new TH2D("h_tru_neut", "h_tru_neut", Nbins, cth_min, cth_max, Nbins, pmu_min, pmu_max);

	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TTree * truth_tree_neut  = (TTree*) fin_neut  -> Get("truth");
	TTree * selec_tree_neut  = (TTree*) fin_neut  -> Get("default");

	int NselEvents_neut, NtruEvents_neut;

	std::cout << "Get NEUT events..." << std::endl;
	NtruEvents_neut = GetTrueEvents(target, truth_tree_neut,     h_tru_neut);
	NselEvents_neut = GetSelectedEvents(target, selec_tree_neut, h_sel_neut);
	std::cout << "--> NEUT signal efficiency = " << (1.0*NselEvents_neut)/(1.0*NtruEvents_neut) << std::endl;

	//======================================================================================================



	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute the efficiency" << std::endl;

	h_sel_neut -> Divide(h_tru_neut);

	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      


	TCanvas* c_eff = new TCanvas("efficiency", "efficiency",1400,1200);
	
	gStyle->SetPalette(kCherry);
	TColor::InvertPalette();

	gROOT->ForceStyle();
	
	h_sel_neut -> GetXaxis() -> SetNdivisions(5);
	h_sel_neut -> GetYaxis() -> SetNdivisions(5);
	h_sel_neut -> SetTitle(Form("Selection efficiency in %s; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]", target.c_str()));
	
	h_sel_neut -> Draw("colz");
	
	c_eff -> Print(Form("plots/efficiencies/selection_eff_2D_%s.pdf", target.c_str()));

	//======================================================================================================

}



// Get true signal events from truth tree
int GetTrueEvents(string target, TTree* tree, TH2D* h_events)
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
		if(fgdtarget == signal || fgd2fgdtarget == signal) // if numu CC0pi event on carbon or oxygen
		{
			h_events -> Fill(costh, mom/1000.0);
			Nevents++;
		}

	}
	std::cout << "Number of true signal events : " << Nevents << std::endl;
	return Nevents;
}

// Get selected events from default tree
int GetSelectedEvents(string target, TTree* tree, TH2D* h_events)
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

		if( (sample==1  || sample==2  || sample==3  || sample==4  || sample==5  || sample==9  || sample==10 || sample==11 ||
		     sample==12 || sample==13 || sample==14 || sample==15 || sample==16 || sample==20 || sample==21 || sample==22 ||
		     sample==23 || sample==24 || sample== 25|| sample==26 || sample==27 || sample==31 || sample==32 || sample==33 ) // if selected event
		    && (fgdtarget == signal || fgd2fgdtarget == signal) // if signal on carbon or oxygen
		    && weight>0.0 && weight<10.0 )
		{
			h_events -> Fill(costh, mom/1000.0, weight);
			Nevents++;
		}
	}
	std::cout << "Number of selected events : " << Nevents << std::endl;
	return Nevents;
}
