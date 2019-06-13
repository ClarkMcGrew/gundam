//  Usage:
// source /usr/local/root/pro/bin/thisroot.sh
// root -b -q 'signalEfficiency2D_regions.C+()'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int GetTrueEvents(string, TTree*, TH2D*);
int GetSelectedEvents(string, TTree*, TH2D*);
int GetSelectedEventsRegion(string, int, TTree*, TH2D*);
TH2D* effByRegions(string, int);

void signalEfficiency2D_regions()
{

	TH2D* h_eff_fgd1_allregions = (TH2D*)effByRegions("FGD1", 0);
	TH2D* h_eff_fgd2_allregions = (TH2D*)effByRegions("FGD2", 0);

	std::vector< TH2D* > h_eff_fgd1_region;
	std::vector< TH2D* > h_eff_fgd2_region;

	for(int i=0; i<5; i++)
	{
		h_eff_fgd1_region.push_back( (TH2D*)effByRegions("FGD1", i+1) );
		h_eff_fgd2_region.push_back( (TH2D*)effByRegions("FGD2", i+1) );
	}


	// FGD1-FGD2 efficiency difference

	TH2D* h_eff_diff_allregions = (TH2D*)h_eff_fgd2_allregions -> Clone("h_eff_diff_allregions");
	h_eff_diff_allregions -> Scale(-1.0);
	h_eff_diff_allregions -> Add(h_eff_fgd1_allregions);

	std::vector< TH2D* > h_eff_diff_region;
	for(int i=0; i<5; i++)
	{
		h_eff_diff_region.push_back( (TH2D*)h_eff_fgd2_region[i] -> Clone(Form("h_eff_diff_region%d",i)) );

		h_eff_diff_region[i] -> Scale(-1.0);
		h_eff_diff_region[i] -> Add(h_eff_fgd1_region[i]);
	}
	

	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	// DRAW OPTIONS
	gStyle->SetPalette(kCherry);
	TColor::InvertPalette();
	gROOT->ForceStyle();

	// Draw efficiency for all signal regions

	TCanvas* c_eff_allregions = new TCanvas("c_eff_allregions", "c_eff_allregions",1900,800);
	c_eff_allregions -> Divide(2,1);

	c_eff_allregions -> cd(1);
	
	h_eff_fgd1_allregions -> SetMinimum(0.0);
	h_eff_fgd1_allregions -> SetMaximum(1.0);

	h_eff_fgd1_allregions -> GetXaxis() -> SetNdivisions(5);
	h_eff_fgd1_allregions -> GetYaxis() -> SetNdivisions(5);
	h_eff_fgd1_allregions -> SetTitle("Selection efficiency in FGD1; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]");
	
	h_eff_fgd1_allregions -> Draw("colz");
	gPad->SetRightMargin(0.2);
	gPad->Update();

	c_eff_allregions -> cd(2);
	
	h_eff_fgd2_allregions -> SetMinimum(0.0);
	h_eff_fgd2_allregions -> SetMaximum(1.0);

	h_eff_fgd2_allregions -> GetXaxis() -> SetNdivisions(5);
	h_eff_fgd2_allregions -> GetYaxis() -> SetNdivisions(5);
	h_eff_fgd2_allregions -> SetTitle("Selection efficiency in FGD2; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]");
	
	h_eff_fgd2_allregions -> Draw("colz");
	gPad->SetRightMargin(0.2);
	gPad->Update();

	c_eff_allregions -> Print("plots/efficiencies/selection_eff_2D_signal_allregions.pdf");



	// Draw efficiency for each signal regions

	//FGD1
	TCanvas* c_eff_regions_FGD1 = new TCanvas("c_eff_regions_FGD1", "c_eff_regions_FGD1",1800,1000);
	c_eff_regions_FGD1 -> Divide(3,2);

	for(int i=0; i<5; i++)
	{
		c_eff_regions_FGD1 -> cd(i+1);
		
		h_eff_fgd1_region[i] -> SetMinimum(0.0);
		h_eff_fgd1_region[i] -> SetMaximum(1.0);

		h_eff_fgd1_region[i] -> GetXaxis() -> SetNdivisions(5);
		h_eff_fgd1_region[i] -> GetYaxis() -> SetNdivisions(5);
		h_eff_fgd1_region[i] -> SetTitle(Form("FGD1, signal region %d; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]", i+1));
		
		h_eff_fgd1_region[i] -> Draw("colz");
		gPad->SetRightMargin(0.2);
		gPad->Update();
	}
	
	c_eff_regions_FGD1 -> Print("plots/efficiencies/selection_eff_2D_signal_each_region_FGD1.pdf");

	//FGD2
	TCanvas* c_eff_regions_FGD2 = new TCanvas("c_eff_regions_FGD2", "c_eff_regions_FGD2",1800,1000);
	c_eff_regions_FGD2 -> Divide(3,2);

	for(int i=0; i<5; i++)
	{
		c_eff_regions_FGD2 -> cd(i+1);
		
		h_eff_fgd2_region[i] -> SetMinimum(0.0);
		h_eff_fgd2_region[i] -> SetMaximum(1.0);

		h_eff_fgd2_region[i] -> GetXaxis() -> SetNdivisions(5);
		h_eff_fgd2_region[i] -> GetYaxis() -> SetNdivisions(5);
		h_eff_fgd2_region[i] -> SetTitle(Form("FGD2, signal region %d; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]", i+1));
		
		h_eff_fgd2_region[i] -> Draw("colz");
		gPad->SetRightMargin(0.2);
		gPad->Update();
	}
	
	c_eff_regions_FGD2 -> Print("plots/efficiencies/selection_eff_2D_signal_each_region_FGD2.pdf");




	//======================================================================================================


	const Int_t NRGBs = 3;
	const Int_t NCont = 20;
	Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 1.00, 1.00 };
	Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
	Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

	// Draw FGD1-FGD2 efficiency DIFFERENCE for all signal regions

	TCanvas* c_diff_allregions = new TCanvas("c_diff_allregions", "c_diff_allregions",900,800);
	
	h_eff_diff_allregions -> SetMinimum(-1.0);
	h_eff_diff_allregions -> SetMaximum(1.0);

	h_eff_diff_allregions -> GetXaxis() -> SetNdivisions(5);
	h_eff_diff_allregions -> GetYaxis() -> SetNdivisions(5);
	h_eff_diff_allregions -> SetTitle("FGD1-FGD2 difference; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]");
	
	h_eff_diff_allregions -> Draw("colz");
	gPad->SetRightMargin(0.2);
	gPad->Update();

	c_diff_allregions -> Print("plots/efficiencies/selection_eff_2D_signal_allregions_FGD1FGD2diff.pdf");



	// Draw efficiency for each signal regions

	TCanvas* c_diff_regions = new TCanvas("c_diff_regions", "c_diff_regions",1800,1000);
	c_diff_regions -> Divide(3,2);

	for(int i=0; i<5; i++)
	{
		c_diff_regions -> cd(i+1);
		
		h_eff_diff_region[i] -> SetMinimum(-1.0);
		h_eff_diff_region[i] -> SetMaximum(1.0);

		h_eff_diff_region[i] -> GetXaxis() -> SetNdivisions(5);
		h_eff_diff_region[i] -> GetYaxis() -> SetNdivisions(5);
		h_eff_diff_region[i] -> SetTitle(Form("FGD1-FGD2 diff., signal region %d; True cos(#theta_{#mu}); True p_{#mu} [GeV/c]", i+1));
		
		h_eff_diff_region[i] -> Draw("colz");
		gPad->SetRightMargin(0.2);
		gPad->Update();
	}
	
	c_diff_regions -> Print("plots/efficiencies/selection_eff_2D_signal_each_region_FGD1FGD2diff.pdf");



	std::cout << "================================================" << std::endl;
	//======================================================================================================
}




//======================================================================================================//======================================================================================================
// Compute the efficiency for one signal region (1,2,3,4,5) or all signal regions together (0)
//======================================================================================================
TH2D* effByRegions(string fgd = "FGD1", int region = 0)
{

	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Start to compute efficiency for " << fgd << " in signal region " << region << " =====" << std::endl;

	if(fgd != "FGD1" && fgd != "FGD2")
		std::cout << "ERROR : First argument must be FGD1 or FGD2 !!!" << std::endl;

	string infilename_neut  = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_neut_prod6B/output_run2348_mc_500toys.root";
	// string infilename_neut  = "/sps/t2k/lmaret/outputs_highland2/numuCCZeroPi/outputs_neut_prod6B/output_run2w_mc_500toys.root";
	

	//======================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();
	gROOT->ForceStyle();

	int Nbins = 50;
	double pmu_min =  0.0;
	double pmu_max =  5.0;
	double cth_min = -1.0;
	double cth_max =  1.0;
	
	//======================================================================================================


	//======================================================================================================
	std::cout << "===== Store Highland2 file and initialise histograms =====" << std::endl;

	TFile* fin_neut  = new TFile(infilename_neut.c_str());

	TH2D* h_sel_neut = new TH2D("h_sel_neut", "h_sel_neut", Nbins, cth_min, cth_max, Nbins, pmu_min, pmu_max);
	TH2D* h_tru_neut = new TH2D("h_tru_neut", "h_tru_neut", Nbins, cth_min, cth_max, Nbins, pmu_min, pmu_max);

	//======================================================================================================

	//======================================================================================================  
	std::cout << "===== Get histograms containing number of events =====" << std::endl;
	
	TTree * truth_tree_neut  = (TTree*) fin_neut  -> Get("truth");
	TTree * selec_tree_neut  = (TTree*) fin_neut  -> Get("default");

	int NselEvents_neut, NtruEvents_neut;

	std::cout << "Get NEUT events..." << std::endl;
	NtruEvents_neut = GetTrueEvents(fgd, truth_tree_neut,     h_tru_neut);
	if(region==0) NselEvents_neut = GetSelectedEvents(fgd, selec_tree_neut, h_sel_neut);
	else          NselEvents_neut = GetSelectedEventsRegion(fgd, region, selec_tree_neut, h_sel_neut);
	std::cout << "--> NEUT signal efficiency = " << (1.0*NselEvents_neut)/(1.0*NtruEvents_neut) << std::endl;

	//======================================================================================================


	//===================================================================================================================
	std::cout << "=== Compute the efficiency" << std::endl;

	h_sel_neut -> Divide(h_tru_neut);

	TH2D *h_eff_2D = new TH2D(*h_sel_neut);
	//======================================================================================================

	// delete fin_neut;
	// delete h_sel_neut;
	// delete h_tru_neut;

	//======================================================================================================  
	std::cout << "=== Done." << std::endl;

	return h_eff_2D;
}




//======================================================================================================
// Get true signal events from truth tree
//======================================================================================================
int GetTrueEvents(string fgd, TTree* tree, TH2D* h_events)
{
	int fgdtarget, fgd2fgdtarget, detector;
	int det = -999;
	float mom, costh;
	int costhbin;
	int signal=999;

	tree -> SetBranchAddress("truelepton_mom",       &mom);
	tree -> SetBranchAddress("truelepton_costheta",  &costh);
	tree -> SetBranchAddress("fgdtargetCCZeroPi",    &fgdtarget);
	tree -> SetBranchAddress("fgd2fgdtargetCCZeroPi",&fgd2fgdtarget);
	tree -> SetBranchAddress("detector",             &detector);

	if(fgd=="FGD1")      det = 0;
	else if(fgd=="FGD2") det = 1;

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);
		if( (fgdtarget == 0 || fgd2fgdtarget == 0 ||fgdtarget == 3 || fgd2fgdtarget == 3) && detector == det ) // if numu CC0pi event on carbon or oxygen
		{
			h_events -> Fill(costh, mom/1000.0);
			Nevents++;
		}

	}
	std::cout << "Number of true signal events : " << Nevents << std::endl;
	return Nevents;
}



//======================================================================================================
// Get selected events of all signal regions from default tree
//======================================================================================================
int GetSelectedEvents(string fgd, TTree* tree, TH2D* h_events)
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

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		if(fgd=="FGD1") // FGD1
		{
			if( (sample==1  || sample==2  || sample==3  || sample==4  || sample==5  || sample==9  || sample==10 || sample==11 ) // if selected event
			    && (fgdtarget == 0 || fgdtarget == 3) // if signal on carbon or oxygen
			    && weight>0.0 && weight<10.0 )
			{
				h_events -> Fill(costh, mom/1000.0, weight);
				Nevents++;
			}
		}
		else if(fgd=="FGD2") // FGD2
		{
			if( (sample==12 || sample==13 || sample==14 || sample==15 || sample==16 || sample==20 || sample==21 || sample==22 ||
			     sample==23 || sample==24 || sample== 25|| sample==26 || sample==27 || sample==31 || sample==32 || sample==33 ) // if selected event
			    && ( fgd2fgdtarget == 0 || fgd2fgdtarget == 3) // if signal on carbon or oxygen
			    && weight>0.0 && weight<10.0 )
			{
				h_events -> Fill(costh, mom/1000.0, weight);
				Nevents++;
			}
		}
	}
	std::cout << "Number of selected events : " << Nevents << std::endl;
	return Nevents;
}



//======================================================================================================
// Get selected events of one signal region from default tree
//======================================================================================================
int GetSelectedEventsRegion(string fgd, int region, TTree* tree, TH2D* h_events)
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

	int Nevents = 0;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		if(fgd == "FGD1") // FGD1
		{
			if( ( sample == region ) // if selected event
			    && (fgdtarget == 0 || fgdtarget == 3) // if signal on carbon or oxygen
			    && weight>0.0 && weight<10.0 )
			{
				h_events -> Fill(costh, mom/1000.0, weight);
				Nevents++;
			}
		}
		else if(fgd == "FGD2") // FGD2
		{
			if( ( sample == region+11 || sample == region+22 ) // if selected event
			    && ( fgd2fgdtarget == 0 || fgd2fgdtarget == 3) // if signal on carbon or oxygen
			    && weight>0.0 && weight<10.0 )
			{
				h_events -> Fill(costh, mom/1000.0, weight);
				Nevents++;
			}
		}
	}
	std::cout << "Number of selected events : " << Nevents << std::endl;
	return Nevents;
}
