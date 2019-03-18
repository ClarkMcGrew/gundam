//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root -b -q 'DrawXsec.C+()'
//         root -b -q 'DrawXsec.C+("fit3_statFluc")'
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double calcChi2_M(TH1D*, TH1D*, TMatrixD);


void DrawXsec(string inputname = "fit3_statFluc", const std::string& dir_name = "fakedata/statFluc", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());

	//======================================================================================================  
	//=== Set common style
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
	std::cout << "===== Get histos with cross-section result and truth =====" << std::endl;
	
	vector< vector<TH1D*> > h_xsec_truth(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit(Ntarget);
	vector< vector<TH1D*> > h_xsec_postfit_err(Ntarget);


	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_truth[itar].push_back(       (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d_truth", targetlist[itar].c_str(), nbcth) )) );
			h_xsec_postfit[itar].push_back(     (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d",       targetlist[itar].c_str(), nbcth) )) );
			h_xsec_postfit_err[itar].push_back( (TH1D*)(fin->Get( Form("CC0pi%s_cos_bin%d",       targetlist[itar].c_str(), nbcth) ))->Clone(Form("CC0pi%s_cos_bin%d_err", targetlist[itar].c_str(), nbcth)) );
		}

	vector< TH1D* > h_ratio_truth;
	vector< TH1D* > h_ratio_postfit;
	vector< TH1D* > h_ratio_postfit_err;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_truth.push_back(       (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d_truth", nbcth) )) );
		h_ratio_postfit.push_back(     (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) )) );
		h_ratio_postfit_err.push_back( (TH1D*)(fin->Get( Form("CC0piOCRatio_cos_bin%d", nbcth) ))->Clone(Form("CC0piOCRatio_cos_bin%d_errs", nbcth)) );
	}
	//======================================================================================================

	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Set histos with errors =====" << std::endl;
	
	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_postfit_err[itar][nbcth] -> Reset();
			for(int i=1; i<=h_xsec_postfit_err[itar][nbcth]->GetNbinsX(); i++)
			{
				double err  = h_xsec_postfit[itar][nbcth] -> GetBinError(i);
				double mean = h_xsec_postfit[itar][nbcth] -> GetBinContent(i);
				if(err/mean<100)
					h_xsec_postfit_err[itar][nbcth] -> SetBinContent(i, err/mean);
				else
					h_xsec_postfit_err[itar][nbcth] -> SetBinContent(i, 0.0);
			}
		}

	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit_err[nbcth] -> Reset();
		for(int i=1; i<=h_ratio_postfit_err[nbcth]->GetNbinsX(); i++)
		{
			double err  = h_ratio_postfit[nbcth] -> GetBinError(i);
			double mean = h_ratio_postfit[nbcth] -> GetBinContent(i);
			if(err/mean<100)
				h_ratio_postfit_err[nbcth] -> SetBinContent(i, err/mean);
			else
				h_ratio_postfit_err[nbcth] -> SetBinContent(i, 0.0);
		}
	}
	//======================================================================================================






	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get covariance matrix and compute chi2 =====" << std::endl;
	
	TMatrixDSym *cov_mat = (TMatrixDSym*)(fin -> Get("xsec_cov"));
	TH1D* h_sel_best_fit = (TH1D*)(fin->Get("sel_best_fit"));
	TH1D* h_tru_best_fit = (TH1D*)(fin->Get("tru_best_fit"));

	double chi2 = calcChi2_M(h_sel_best_fit, h_tru_best_fit, *cov_mat);

	// h_sel_best_fit -> Draw();
	// h_tru_best_fit -> Draw("same");

	//======================================================================================================








	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Make final plots =====" << std::endl;                                                                                                                      

	TCanvas* c_xsec[Ntarget];
	TCanvas* c_xsec_err[Ntarget];
	vector<TLegend*> leg(Ntarget+1);
	double max;

	for(int itar = 0; itar < Ntarget; itar++)
	{
		c_xsec[itar] = new TCanvas(Form("Xsec_on_%s", targetlist[itar].c_str()),Form("Xsec_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_truth[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_truth[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [MeV/c]");
			h_xsec_truth[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon MeV/c}]");

			max = MAX(h_xsec_truth[itar][nbcth]->GetMaximum(), h_xsec_postfit[itar][nbcth]->GetMaximum());
			h_xsec_truth[itar][nbcth] -> GetYaxis() -> SetRangeUser(0, 1.3*max);

			h_xsec_truth[itar][nbcth]   -> SetLineColor(kGreen+2);
			h_xsec_postfit[itar][nbcth] -> SetLineColor(kBlue+1);

			c_xsec[itar] -> cd(nbcth+1);

			if(nbcth!=0) gPad->SetLogx();

			h_xsec_truth[itar][nbcth]   -> Draw("hist");
			h_xsec_postfit[itar][nbcth] -> Draw("same E"); // with error bars

			if(nbcth==0)
			{
				leg[itar] = new TLegend(0.2,0.24,0.8,0.64);
				leg[itar] -> SetFillColor(0);
				leg[itar] -> SetBorderSize(1);
				leg[itar] -> SetFillStyle(0);
				//leg[itar]->SetTextSize(0.075);
				leg[itar] -> SetHeader(Form("#chi^{2}  = %d", (int)chi2));
				leg[itar] -> AddEntry(h_xsec_postfit[itar][0], "post-fit", "lep");
				leg[itar] -> AddEntry(h_xsec_truth[itar][0],   "truth",       "l");
				c_xsec[itar]->cd(1);
				leg[itar]->Draw();
			}

		}
		c_xsec[itar]->Print(Form("plots/xsecResults/%s/xsec_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));


		c_xsec_err[itar] = new TCanvas(Form("Xsec_error_on_%s", targetlist[itar].c_str()),Form("Xsec_error_on_%s", targetlist[itar].c_str()),1700,1000);
		c_xsec_err[itar] -> Divide(3,3);

		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			h_xsec_postfit_err[itar][nbcth] -> SetTitle(h_title[nbcth]);
			h_xsec_postfit_err[itar][nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
			h_xsec_postfit_err[itar][nbcth] -> GetYaxis() -> SetTitle("#frac{d#sigma}{dp_{#mu} dcos#theta_{#mu}} [#frac{cm^2}{nucleon GeV/c}]");

			h_xsec_postfit_err[itar][nbcth] -> SetLineColor(kBlue+1);

			c_xsec_err[itar] -> cd(nbcth+1);

			if(nbcth!=0) gPad->SetLogx();

			h_xsec_postfit_err[itar][nbcth] -> GetYaxis() -> SetRangeUser(0.0, 0.25);
			h_xsec_postfit_err[itar][nbcth] -> Draw("hist");

		}
		c_xsec_err[itar]->Print(Form("plots/xsecResults/%s/xsec_error_%s_%s.pdf", dir_name.c_str(), inputname.c_str(), targetlist[itar].c_str()));
	}


	TCanvas* c_ratio = new TCanvas("Xsec_ratio", "Xsec_ratios",1700,1000);
	c_ratio -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_truth[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_truth[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_truth[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");
		h_ratio_truth[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 2.0);
		
		h_ratio_truth[nbcth] -> SetLineColor(kGreen+2);
		h_ratio_postfit[nbcth] -> SetLineColor(kBlue+1);

		c_ratio -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_ratio_truth[nbcth]   -> Draw("hist");
		h_ratio_postfit[nbcth] -> Draw("same E"); // with error bars

		// if(nbcth==0)
		// {
		// 	leg[Ntarget] = new TLegend(0.2,0.24,0.8,0.64);
		// 	leg[Ntarget] -> SetFillColor(0);
		// 	leg[Ntarget] -> SetBorderSize(1);
		// 	leg[Ntarget] -> SetFillStyle(0);
		// 	//leg[Ntarget]->SetTextSize(0.075);
		// 	leg[Ntarget] -> SetHeader(Form("#chi^{2}  = %d", (int)chi2));
		// 	leg[Ntarget] -> AddEntry(h_ratio_postfit[0], "(fake) data", "lep");
		// 	leg[Ntarget] -> AddEntry(h_ratio_truth[0],   "truth",       "l");
		// 	c_ratio->cd(1);
		// 	leg[Ntarget]->Draw();
		// }
	}
	c_ratio->Print( Form("plots/xsecResults/%s/xsec_%s_OCratio.pdf", dir_name.c_str(), inputname.c_str()) );



	TCanvas* c_ratio_err = new TCanvas("Xsec_ratio", "Xsec_ratio",1700,1000);
	c_ratio_err -> Divide(3,3);
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		h_ratio_postfit_err[nbcth] -> SetTitle(h_title[nbcth]);
		h_ratio_postfit_err[nbcth] -> GetXaxis() -> SetTitle("p^{#mu}_{true} [GeV/c]");
		h_ratio_postfit_err[nbcth] -> GetYaxis() -> SetTitle("#frac{#sigma^{O}}{#sigma^{C}}");

		h_ratio_postfit_err[nbcth] -> SetLineColor(kBlue+1);

		c_ratio_err -> cd(nbcth+1);

		if(nbcth!=0) gPad->SetLogx();

		h_ratio_postfit_err[nbcth] -> GetYaxis() -> SetRangeUser(0.0, 0.5);
		h_ratio_postfit_err[nbcth] -> Draw("hist");
	}
	c_ratio_err->Print( Form("plots/xsecResults/%s/xsec_error_%s_OCratio.pdf", dir_name.c_str(), inputname.c_str()) );

	//======================================================================================================






	//======================================================================================================  
	std::cout << "================================================" << std::endl;
	std::cout << "===== Get histos with efficiency and plot it =====" << std::endl;
	
	TH1D* h_eff = (TH1D*)( fin->Get("eff_best_fit") );

	TCanvas* c_eff = new TCanvas("Efficiency", "Efficiency",800,600);
	
	h_eff -> SetTitle("Selection efficiency");
	h_eff -> GetXaxis() -> SetTitle("true analysis bins");
	h_eff -> GetYaxis() -> SetTitle("true events / selected events");
	h_eff -> SetLineColor(kBlue);

	// if(nbcth!=0) gPad->SetLogx();

	h_eff->Draw("hist");
	
	c_eff->Print( Form("plots/xsecResults/%s/efficiency_%s.pdf", dir_name.c_str(), inputname.c_str()) );
	//======================================================================================================




	//======================================================================================================
	delete fin;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}





double calcChi2_M(TH1D* h1, TH1D* h2, TMatrixD covar)
{
	double chi2=0;

	cout << "calcChi2_M::Starting chi2 calc" << endl;
	cout << "calcChi2_M::h1 (data) nbins: " << h1->GetXaxis()->GetNbins() << endl;
	cout << "calcChi2_M::h2 (MC) nbins: " << h2->GetXaxis()->GetNbins() << endl;
	cout << "calcChi2_M::matrix rows: " << covar.GetNrows() << endl;

	cout << "calcChi2_M::Printing data histo to be sure of content:" << endl;
	h1->Print("all");

	cout << "calcChi2_M::Determinant of covar is: " << covar.Determinant() << endl;

	//Matrix inversion section:
	bool inversionError = false;
	TMatrixDSym* covar_inv = new TMatrixDSym();
	TMatrixD covar_origin = covar;
	int count = 0;

	while(true)
	{
		//covar.Print();
		//covar.SetTol(1e-100);
		//covar.Invert();
		TDecompSVD LU = TDecompSVD(covar);
		covar_inv = new TMatrixDSym(covar.GetNrows(), LU.Invert().GetMatrixArray(), "");
		bool inversionError = false;

		//Check that the matrix inversion worked 
		TMatrixD test = covar*(*covar_inv);
		for(int i=0; i<covar.GetNrows(); i++)
		{
			for(int j=0; j<covar.GetNrows(); j++)
			{
				if(i==j && (test[i][j]>1.00001 || test[i][j]<0.99999))
				{
					if(!inversionError) cout 	<< "****** WARNING: Issue with matrix inversion in chi2 calculation."
												<< " Element ("<<i<<","<<j<<") = " << test[i][j] << endl; 
					inversionError=true;
				}
				else if(i!=j && (test[i][j]>0.0000001))
				{
				if(!inversionError) cout 	<< "****** WARNING: Issue with matrix inversion in chi2 calculation."
											<< " Element ("<<i<<","<<j<<") = " << test[i][j] << endl; 
				inversionError=true;
				}
			}
		}
		if(inversionError)
		{ 
			cout << "DEBUG info after iteration " << count << endl;
			cout << "  Input matrix: determinent is: " << covar.Determinant() <<  endl;
			// cout << "Will try adding 0.00000001pc to the diagonal ..." << endl;
			cout << "Will try adding 0.001pc to the diagonal ..." << endl;
			for(int i=0; i<covar.GetNrows(); i++)
			{
				// covar[i][i]=covar[i][i]*1.00000001;
				covar[i][i]=covar[i][i]*1.001;
			}
			count++;
			//if(count>10)getchar();
		}
		else
		{
			cout << "chi2 calculated successfully after interaction " << count << endl;
			//if(count>0) test.Print();
			break;
		}
	}

	//covar.Print();
	//cout << "Inverted covariance matrix" << endl;
	for(int i=0; i<covar.GetNrows(); i++)
	{
		for(int j=0; j<covar.GetNrows(); j++)
		{
			chi2+= ((h1->GetBinContent(i+1)) - (h2->GetBinContent(i+1)))*(*covar_inv)[i][j]*((h1->GetBinContent(j+1)) - (h2->GetBinContent(j+1)));
			if(i==j)
			{
				cout << "calcChi2_M::bin " << i << ": " << h1->GetBinError(i+1) << " should be the same as " << sqrt(covar_origin[i][j]) << endl;
			}
		}
	}
	cout << "calcChi2_M::chi2 is: " << chi2 << endl;
	return chi2;
}