//  Usage:
// source /usr/local/root/pro/bin/thisroot.sh
// root -b -q 'GetQ2.C+("C")'
// root -b -q 'GetQ2.C+("O")'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


void GetQ2(string, TTree*, std::vector< TH2D* >, double*, int);
int FindCosthBin(double, double*, int);


void trueQ2(string target = "C", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt")
{

	if(target != "C" && target != "O")
		std::cout << "ERROR : First argument must be C (for carbon) or O (for oxygen) !!!" << std::endl;

	string infilename  = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_neut_prod6B/output_run2348_mc_500toys.root";


	//======================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.25);
	// gROOT->ForceStyle();

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

	const int Nq2bins = 50;
	double q2bin_high = 1500E3;
	// double q2bins[Nq2bins+1] = {0.0, 100E3, 200E3, 300E3, 400E3, 500E3, 600E3, 700E3, 800E3, 900E3, 1000E3, 2000E3, 3000E3, 6000E3};

	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store Highland2 outputs in a TFile =====" << std::endl;

	TFile* fin  = new TFile(infilename.c_str());
	//======================================================================================================


	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Initialise histograms with binning" << std::endl;

	std::vector< TH2D* > h_tmp;

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
		h_tmp.push_back(new TH2D("","", Nmombins[nbcth], mombins[nbcth], Nq2bins, 0.0, q2bin_high));
		h_tmp[nbcth] -> SetTitle( Form("%s", h_title[nbcth].Data()) );
	}

	std::vector< TH2D* > h_mom_q2;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_mom_q2.push_back((TH2D*)h_tmp[nbcth]  -> Clone(Form("h_mom_q2_%d",nbcth)));
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TTree * truth_tree_neut  = (TTree*) fin -> Get("truth");

	std::cout << "Get momentum VS Q^2 distributions..." << std::endl;
	GetQ2(target, truth_tree_neut, h_mom_q2, costhBins, Nbins_costh);

	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      


	TCanvas* c_q2 = new TCanvas("efficiency", "efficiency",1700,1000);
	c_q2 -> Divide(3,3);

	gStyle->SetPalette(kCherry);
	TColor::InvertPalette();

	gROOT->ForceStyle();

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_mom_q2[nbcth] -> SetTitle(h_title[nbcth]);
		h_mom_q2[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_mom_q2[nbcth] -> GetYaxis() -> SetTitle("q^{2}_{true} [(MeV/c)^{2}]");

		// h_mom_q2[nbcth] -> SetMaximum(45000);

		c_q2 -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		gPad -> SetLogz();

		// Set range without high momentum bin
		double lastBinEdge = mombins[nbcth][(Nmombins[nbcth]-1)]; // Nmombins = number of mom bins in this costh slice
		// std::cout << "costh bin nb = " << nbcth << ", Nmombins = " << Nmombins[nbcth] << " and last bin edge = " << lastBinEdge << std::endl;
		h_mom_q2[nbcth] -> GetXaxis() -> SetRangeUser(0.0, lastBinEdge);
		

		h_mom_q2[nbcth] -> Draw("colz");

	}
	c_q2 -> Print(Form("plots/momentumTransfer/momentum_vs_q2_distribution_%s.pdf", target.c_str()));


	delete fin;
	delete c_q2;

	//======================================================================================================

}



// Get true signal events from truth tree
void GetQ2(string target, TTree* tree, std::vector< TH2D* > h_events, double* CosthBins, int Nbins_costh)
{
	float true_Q2;
	int fgdtarget, fgd2fgdtarget;
	float mom, costh;
	int costhbin;
	int signal=999;

	
	tree -> SetBranchAddress("true_Q2",       &true_Q2);
	tree -> SetBranchAddress("truelepton_mom",       &mom);
	tree -> SetBranchAddress("truelepton_costheta",  &costh);
	tree -> SetBranchAddress("fgdtargetCCZeroPi",    &fgdtarget);
	tree -> SetBranchAddress("fgd2fgdtargetCCZeroPi",&fgd2fgdtarget);

	if(target=="C")      signal = 3;
	else if(target=="O") signal = 0;

	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree -> GetEntry(i);
		int costhbin = FindCosthBin(costh, CosthBins, Nbins_costh);

		if(fgdtarget == signal || fgd2fgdtarget == signal ) // if numu CC0pi event on carbon or oxygen
		{
			h_events[costhbin] -> Fill(mom, true_Q2);
		}

	}
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
