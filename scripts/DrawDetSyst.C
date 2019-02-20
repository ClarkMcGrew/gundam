// Usage:
// root 'DrawDetSyst.C+()'

#include "CommonHeader.h"
#include "CommonStyle.h"
#include "BinningTools.cc"

using namespace std;

void DrawDetSyst()
{
	//====================================================================================================                      
	
	string infile = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/xsllh_detcovmat_4samples.root";
	// string infile = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/xsllh_detcovmat.root";
	string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt";

	//=== Set binning using BinningTools class                                                                                
	BinningTools bin;
	bin.SetBinning(fbinning.c_str());
	int Nbins = bin.GetNbins();//Total number of bins

	int NsamplesTot = 12; // To draw whole signal region + control regions
	// int NsamplesTot = 24; // To draw each signal region + control regions

	//====================================================================================================                      


	//====================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();

	//======================================================================================================

	// Open file and get covariance matrix
	TFile *findetcov = TFile::Open(infile.c_str());
	TMatrixDSym *cov_mat = (TMatrixDSym*)findetcov -> Get("cov_mat");
	// TMatrixDSym *cor_mat = (TMatrixDSym*)findetcov -> Get("cor_mat");

	cout << "Number of samples = "      << NsamplesTot << endl;
	cout << "Number of bins = "         << Nbins << endl;
	cout << "Number matrix elements = " << cov_mat->GetNcols() << endl;
	if(cov_mat->GetNcols() != NsamplesTot*Nbins) cout << "ERROR : matrix dimension is not equal to nb of bins * nb of samples !!!" << endl;

	//== Get matrix diagonal (= error^2)

	TH1F* h_err[NsamplesTot];
	int index=0;
	float cont=0;

	for(int is = 0; is<NsamplesTot; is++)
	{
		h_err[is] = new TH1F(Form("h_err_%d",is), Form("h_err_%d",is), Nbins, 0, Nbins);

		// cout << "--------- SAMPLE " << is << " ---------" << endl;
		for(int ib = 0; ib<Nbins; ib++)
		{
			index = is*Nbins + ib;
			cont  = TMath::Sqrt( (*cov_mat)(index, index) );
			h_err[is] -> SetBinContent(ib+1, cont);

			// cout << "bin = " << ib+1 << " and index = " << index << " --> error = " << cont << endl;
		}
	}


	//=== Draw errors
	
	if(NsamplesTot%3 != 0) std::cout << "WARNING : number of samples " << NsamplesTot << " not a multiple of 3 !!!" << std::endl;
	int NsamplesDraw = NsamplesTot/3;

	string sampleNames[NsamplesDraw];

	sampleNames[0] = "CC-0Pion";
	sampleNames[1] = "CC-1PiPlus";
	sampleNames[2] = "CC-Other";
	sampleNames[3] = "CC-1PiMichel";

	// sampleNames[0] = "muTPC";
	// sampleNames[1] = "muTPCpTPC";
	// sampleNames[2] = "muTPCpFGD";
	// sampleNames[3] = "muFGDpTPC";
	// sampleNames[4] = "muFGD";
	// sampleNames[5] = "CC-1PiPlus";
	// sampleNames[6] = "CC-Other";
	// sampleNames[7] = "CC-1PiMichel";

	TCanvas *c4[NsamplesDraw];

	for(int is=0; is < NsamplesDraw; is++)
	{
		c4[is] = new TCanvas(Form("c4_%d", is),Form("c4_%d", is),800,600);
		c4[is] -> Draw();
		
		h_err[is] -> SetTitle(Form("%s;(p,cos#theta) bin number;Relative detector uncertainty", sampleNames[is].c_str()));

		double Max = std::max( h_err[is]->GetMaximum(), h_err[is+NsamplesDraw]->GetMaximum() );
		Max        = std::max( h_err[is+2*NsamplesDraw]->GetMaximum(), Max );
			
		h_err[is] -> GetYaxis() -> SetRangeUser(0.0, Max * 1.3);

		h_err[is]                -> SetLineColor(kRed);
		h_err[is+NsamplesDraw]   -> SetLineColor(kBlue);
		h_err[is+2*NsamplesDraw] -> SetLineColor(kGreen+2);

		h_err[is]                -> Draw("same hist");
		h_err[is+NsamplesDraw]   -> Draw("same hist");
		h_err[is+2*NsamplesDraw] -> Draw("same hist");

		TLegend* leg = new TLegend(0.68,0.7,0.9,0.87);
		leg -> SetFillColor(0);
		leg -> AddEntry(h_err[is]               , "FGD1" ,  "L");
		leg -> AddEntry(h_err[is+NsamplesDraw]  , "FGD2x",  "L");
		leg -> AddEntry(h_err[is+2*NsamplesDraw], "FGD2y",  "L");
		leg -> Draw();


		c4[is] -> Print(Form("plots/detSyst/detSyst_cov_%s.pdf", sampleNames[is].c_str()));
	}
}
