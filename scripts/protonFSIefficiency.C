//  Usage: root -b -q 'protonFSIefficiency.C+()'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void GetTrueEvents(TTree*, std::vector< TH1D* >, double*, int);
void GetSelectedEvents(TTree*, std::vector< TH1D* >, double*, int);
int FindCosthBin(double, double*, int);


void protonFSIefficiency(string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt")
{

	string infilename_noFSI = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_nuwro/output_nuwro_mc_fgd1_noFSI_500toys.root";
	string infilename_wiFSI = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_nuwro/output_nuwro_mc_fgd1_500toys.root";

	string fitpath = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/";

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

	TFile* fin_noFSI = new TFile(infilename_noFSI.c_str());
	TFile* fin_wiFSI = new TFile(infilename_wiFSI.c_str());
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

	std::vector< TH1D* > h_sel_noFSI;
	std::vector< TH1D* > h_sel_wiFSI;
	std::vector< TH1D* > h_tru_noFSI;
	std::vector< TH1D* > h_tru_wiFSI;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_sel_noFSI.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_sel_noFSI_%d",nbcth)));
		h_sel_wiFSI.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_sel_wiFSI_%d",nbcth)));
		h_tru_noFSI.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_tru_noFSI_%d",nbcth)));
		h_tru_wiFSI.push_back((TH1D*)h_tmp[nbcth] -> Clone(Form("h_tru_wiFSI_%d",nbcth)));
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TTree * truth_tree_noFSI = (TTree*) fin_noFSI -> Get("truth");
	TTree * truth_tree_wiFSI = (TTree*) fin_wiFSI -> Get("truth");
	TTree * selec_tree_noFSI = (TTree*) fin_noFSI -> Get("default");
	TTree * selec_tree_wiFSI = (TTree*) fin_wiFSI -> Get("default");

	std::cout << "Get events without FSI." << std::endl;
	GetTrueEvents(truth_tree_noFSI, h_tru_noFSI, costhBins, Nbins_costh);
	GetSelectedEvents(selec_tree_noFSI, h_sel_noFSI, costhBins, Nbins_costh);

	std::cout << "Get events with FSI." << std::endl;
	GetTrueEvents(truth_tree_wiFSI, h_tru_wiFSI, costhBins, Nbins_costh);
	GetSelectedEvents(selec_tree_wiFSI, h_sel_wiFSI, costhBins, Nbins_costh);


	//======================================================================================================



	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute the efficiency" << std::endl;

	std::vector< TH1D* > h_eff_noFSI;
	std::vector< TH1D* > h_eff_wiFSI;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_noFSI.push_back((TH1D*)h_sel_noFSI[nbcth] -> Clone(Form("h_eff_noFSI_%d",nbcth)));
		h_eff_wiFSI.push_back((TH1D*)h_sel_wiFSI[nbcth] -> Clone(Form("h_eff_wiFSI_%d",nbcth)));
	}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_noFSI[nbcth] -> Divide(h_tru_noFSI[nbcth]);
		h_eff_wiFSI[nbcth] -> Divide(h_tru_wiFSI[nbcth]);
	}
	//======================================================================================================


	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute relative error on the efficiency" << std::endl;
	
	std::vector< TH1D* > h_eff_error;
	std::vector< TH1D* > h_eff_errorTEST;

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_error.push_back((TH1D*)h_eff_wiFSI[nbcth] -> Clone(Form("h_eff_error_%d",nbcth)));
	}

	double tmp_effwiFSI, tmp_effnoFSI;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
		{
			tmp_effnoFSI = h_eff_noFSI[nbcth] -> GetBinContent(nbm+1);
			tmp_effwiFSI = h_eff_wiFSI[nbcth] -> GetBinContent(nbm+1);

			h_eff_error[nbcth] -> SetBinContent(nbm+1, fabs(tmp_effwiFSI - tmp_effnoFSI)/tmp_effwiFSI );
			// h_eff_error[nbcth] -> SetBinContent(nbm+1, (tmp_effwiFSI - tmp_effnoFSI)/tmp_effwiFSI );
		}
	//======================================================================================================


	//===================================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "=== Compute covariance matrix" << std::endl;
	
	TH1D* effWiFSI = new TH1D("effWiFSI", "effWiFSI", 2*Nbins, 0, 2*Nbins);
	TH1D* effNoFSI = new TH1D("effNoFSI", "effNoFSI", 2*Nbins, 0, 2*Nbins);

	int nbinindex = 1;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
		{
			effWiFSI -> SetBinContent(nbinindex, h_eff_wiFSI[nbcth] -> GetBinContent(nbm+1) );
			effNoFSI -> SetBinContent(nbinindex, h_eff_noFSI[nbcth] -> GetBinContent(nbm+1) );
			
			// Since we don't have FGD2, we set the oxygen bins to the values computed for carbon.
			effWiFSI -> SetBinContent(nbinindex + Nbins, h_eff_wiFSI[nbcth] -> GetBinContent(nbm+1) );
			effNoFSI -> SetBinContent(nbinindex + Nbins, h_eff_noFSI[nbcth] -> GetBinContent(nbm+1) );
			
			nbinindex++;
		}

	TMatrixDSym *cov_mat = new TMatrixDSym(2*Nbins);
	TMatrixDSym *cor_mat = new TMatrixDSym(2*Nbins);

	// // Compute matrix with co-variances

	for(int i=0; i<2*Nbins; i++)
		for(int j=0; j<2*Nbins; j++)
			(*cov_mat)(i,j) = (effWiFSI->GetBinContent(i+1) - effNoFSI->GetBinContent(i+1)) * (effWiFSI->GetBinContent(j+1) - effNoFSI->GetBinContent(j+1)) / (effWiFSI->GetBinContent(i+1)*effWiFSI->GetBinContent(j+1));

	for(int i=0; i<2*Nbins; i++)
		for(int j=0; j<2*Nbins; j++)
			(*cor_mat)(i,j) = (*cov_mat)(i,j)/(TMath::Sqrt((*cov_mat)(i,i))*TMath::Sqrt((*cov_mat)(j,j)));


	// Compute matrix without co-variances

	// for(int i=0; i<2*Nbins; i++)
	// 	for(int j=0; j<2*Nbins; j++)
	// 	{
	// 		if(i==j) (*cov_mat)(i,j) = (effWiFSI->GetBinContent(i+1) - effNoFSI->GetBinContent(i+1)) * (effWiFSI->GetBinContent(j+1) - effNoFSI->GetBinContent(j+1)) / (effWiFSI->GetBinContent(i+1)*effWiFSI->GetBinContent(j+1));
	// 		else     (*cov_mat)(i,j) = 0;
	// 	}

	TFile* file_output = TFile::Open(Form("%s/inputs/fgd1fgd2Fit/xsllh_covarProtonFSI.root", fitpath.c_str()), "RECREATE");
    file_output->cd();

    cov_mat -> Write("cov_mat");
    cor_mat -> Write("cor_mat");

	//======================================================================================================




	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	TCanvas *c_matrix = new TCanvas("c_matrix","c_matrix",1000,900);
	cor_mat -> Draw("colz");
	c_matrix -> Print("plots/protonFSI/protonFSI_correlation.pdf");


	TCanvas* c_eff = new TCanvas("nuwro_efficiency", "nuwro_efficiency",1700,1000);
	c_eff -> Divide(3,3);

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_noFSI[nbcth] -> SetTitle(h_title[nbcth]);
		h_eff_noFSI[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_eff_noFSI[nbcth] -> GetYaxis() -> SetTitle("Signal efficiency");
		h_eff_noFSI[nbcth] -> GetYaxis() -> SetRangeUser(0, 1.0);

		h_eff_noFSI[nbcth] -> SetLineColor(kRed+1);
		h_eff_wiFSI[nbcth] -> SetLineColor(kBlue+1);

		c_eff -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_eff_noFSI[nbcth] -> Draw("hist");
		h_eff_wiFSI[nbcth] -> Draw("same hist");

		if(nbcth==0)
		{
			TLegend* leg = new TLegend(0.2,0.40,0.8,0.85);
			leg -> SetFillColor(0);
			leg -> SetBorderSize(1);
			leg -> SetFillStyle(0);
			leg -> SetHeader("NuWro");
			leg -> AddEntry(h_eff_noFSI[0], "no FSI",   "l");
			leg -> AddEntry(h_eff_wiFSI[0], "with FSI", "l");
			c_eff->cd(1);
			leg->Draw();
		}
	}
	c_eff -> Print("plots/protonFSI/protonFSI_efficiency.pdf");




	TCanvas* c_error = new TCanvas("protonFSIsyst", "protonFSIsyst",1700,1000);
	c_error -> Divide(3,3);

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_eff_error[nbcth] -> SetTitle(h_title[nbcth]);
		h_eff_error[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
		h_eff_error[nbcth] -> GetYaxis() -> SetTitle("Proton FSI systematics");
		// h_eff_error[nbcth] -> GetYaxis() -> SetRangeUser(0, 1.0);
		h_eff_error[nbcth] -> SetLineColor(kBlue+1);

		c_error -> cd(nbcth+1);
		if(nbcth!=0) gPad->SetLogx();

		h_eff_error[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 0.4);
		// h_eff_error[nbcth] -> GetYaxis() -> SetRangeUser(-0.15, 0.15);
		h_eff_error[nbcth] -> Draw("hist"); 
	}
	c_error -> Print("plots/protonFSI/protonFSI_syst.pdf");



	//======================================================================================================





}



// Get true signal events from truth tree
void GetTrueEvents(TTree* tree, std::vector< TH1D* > h_events, double* CosthBins, int Nbins_costh)
{
	int fgdtarget;
	float mom, costh;
	int costhbin;

	tree -> SetBranchAddress("fgdtargetCCZeroPi",   &fgdtarget);
	tree -> SetBranchAddress("truelepton_mom",      &mom);
	tree -> SetBranchAddress("truelepton_costheta", &costh);

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);
		int costhbin = FindCosthBin(costh, CosthBins, Nbins_costh);

		if(fgdtarget==3) // if numu CC0pi event on carbon
		{
			h_events[costhbin] -> Fill(mom);
			Nevents++;
		}

	}
	std::cout << "Number of true signal events : " << Nevents << std::endl;
}

// Get selected events from default tree
void GetSelectedEvents(TTree* tree, std::vector< TH1D* > h_events, double* CosthBins, int Nbins_costh)
{
	int sample, fgdtarget;
	float mom, costh, weight;
	int costhbin;

	tree -> SetBranchAddress("sample_clst_fgd2layer_xsec", &sample);
	tree -> SetBranchAddress("fgdtargetCCZeroPi",          &fgdtarget);
	tree -> SetBranchAddress("truelepton_mom",             &mom);
	tree -> SetBranchAddress("truelepton_costheta",        &costh);
	tree -> SetBranchAddress("weight_corr_total",          &weight);

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		int costhbin = FindCosthBin(costh, CosthBins, Nbins_costh);

		if( (sample==1 || sample==2 || sample==3 || sample==4 || sample==5 || sample==9 || sample==10 || sample==11 ) // if selected event
		    && fgdtarget==3 && weight>0.0 && weight<10.0 ) // if signal
		{
			h_events[costhbin] -> Fill(mom, weight);
			Nevents++;
		}
	}
	std::cout << "Number of selected events : " << Nevents << std::endl;
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
