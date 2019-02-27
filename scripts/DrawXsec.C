//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root 'DrawXsec.C+()'
//         root 'DrawXsec.C+("fit1_statFluc")'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

void DrawXsec(string inputname = "fit1_statFluc", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());

	//======================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();

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
		for(int nbm=0; nbm<Nmombins[nbcth]; nbm++)
			mombins[nbcth] = new double[nbm];

	mombins = bin.GetMomBins();

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
	//======================================================================================================


	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Store final output in a TFile =====" << std::endl;

	TFile* fin = new TFile(infilename.c_str());
	//======================================================================================================


	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with cross-section result =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
			h_xsec[itar].push_back( (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d", targetlist[itar].c_str(), nbcth) )) );

	vector< TH1D* > h_ratio;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		h_ratio.push_back( (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );

	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	TCanvas* c_xsec[Ntarget];
	vector<TLegend*> leg(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
	{
		c_xsec[itar] = new TCanvas(Form("Xsec_on_%s", targetlist[itar].c_str()),Form("Xsec_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec[itar][nbcth] -> SetStats(false);

			h_xsec[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec[itar][nbcth] -> SetLineColor(kBlue);

			c_xsec[itar] -> cd(nbcth+1);

			// if(nbcth!=0) gPad->SetLogx();

			h_xsec[itar][nbcth]->Draw("E"); // with error bars

			// if(nbcth==0)
			// {
			// 	leg[ifgd] = new TLegend(0.2,0.24,0.8,0.64);
			// 	leg[ifgd] -> SetFillColor(0);
			// 	leg[ifgd] -> SetBorderSize(1);
			// 	leg[ifgd] -> SetFillStyle(0);
			// 	//leg[ifgd]->SetTextSize(0.075);
			// 	leg[ifgd] -> AddEntry(h_nom_C[ifgd][0],"Nominal distribution, carbon","l");
			// 	leg[ifgd] -> AddEntry(h_nom_O[ifgd][0],"Nominal distribution, oxygen","l");
			// 	c_xsec[ifgd]->cd(1);
			// 	leg[ifgd]->Draw();
			// }
		}

		c_xsec[itar]->Print(Form("plots/xsecResults/xsec_%s_%s.pdf", inputname.c_str(), targetlist[itar].c_str()));
	}


	TCanvas* c_ratio = new TCanvas("Xsec_ratio", "Xsec_ratios",1700,1000);
	c_ratio -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio[nbcth] -> SetStats(false);

		h_ratio[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio[nbcth] -> SetLineColor(kBlue);

		c_ratio -> cd(nbcth+1);

		// if(nbcth!=0) gPad->SetLogx();

		h_ratio[nbcth]->Draw("E"); // with error bars

		// if(nbcth==0)
		// {
		// 	leg[ifgd] = new TLegend(0.2,0.24,0.8,0.64);
		// 	leg[ifgd] -> SetFillColor(0);
		// 	leg[ifgd] -> SetBorderSize(1);
		// 	leg[ifgd] -> SetFillStyle(0);
		// 	//leg[ifgd]->SetTextSize(0.075);
		// 	leg[ifgd] -> AddEntry(h_nom_C[ifgd][0],"Nominal distribution, carbon","l");
		// 	leg[ifgd] -> AddEntry(h_nom_O[ifgd][0],"Nominal distribution, oxygen","l");
		// 	c_xsec[ifgd]->cd(1);
		// 	leg[ifgd]->Draw();
		// }
	}

	c_ratio->Print( Form("plots/xsecResults/xsec_OCratio_%s.pdf", inputname.c_str()) );

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

	delete fin;
}


