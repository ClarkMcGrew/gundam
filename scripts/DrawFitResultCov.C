// Usage:
// root 'DrawFitResultCov.C+()'

#include "CommonHeader.h"
#include "CommonStyle.h"
#include "BinningTools.cc"

void DrawFitResultCov(string inputname = "fit1")
{
	string infile = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/%s.root", inputname.c_str());

	int NbinsAna = 58;
	int NbinsFlux = 20;  
	int NbinsXsec = 22;
	int NbinsDet = 3*8*NbinsAna;                                                   

	int Nbins = 2*NbinsAna + NbinsXsec + NbinsDet + NbinsFlux;

	

	//====================================================================================================                      
	//=== Set boundary lines between samples   

	const int Nlines = 4;
	double boundary[Nlines+1];
	boundary[0] = 0;
	boundary[1] = 2*NbinsAna;
	boundary[2] = 2*NbinsAna + NbinsFlux;
	boundary[3] = 2*NbinsAna + NbinsFlux + NbinsXsec;
	boundary[4] = 2*NbinsAna + NbinsFlux + NbinsXsec + NbinsDet;
	

	//=== Draw an horizontal line for parameter type
	TLine *verline[2*Nlines];
	for(int il = 0; il < Nlines; il++)
	{
		verline[il] = new TLine(boundary[il], boundary[il], boundary[il], boundary[il+1] );
		verline[il] -> SetLineWidth(2);
		verline[il] -> SetLineStyle(2);
		verline[il+Nlines] = new TLine(boundary[il+1], boundary[il], boundary[il+1], boundary[il+1] );
		verline[il+Nlines] -> SetLineWidth(2);
		verline[il+Nlines] -> SetLineStyle(2);
	}

	//=== Draw a vertical line for parameter type
	TLine *orline[2*Nlines];
	for(int il = 0; il < Nlines; il++)
	{
		orline[il] = new TLine(boundary[il], boundary[il], boundary[il+1], boundary[il] );
		orline[il] -> SetLineWidth(2);
		orline[il] -> SetLineStyle(2);
		orline[il+Nlines] = new TLine(boundary[il], boundary[il+1], boundary[il+1], boundary[il+1] );
		orline[il+Nlines] -> SetLineWidth(2);
		orline[il+Nlines] -> SetLineStyle(2);
	}


	const int Nsamples = 23;
	double boundarySample[Nsamples+1];
	for(int il = 0; il < Nsamples; il++)
	{
		boundarySample[il] = 2*NbinsAna+NbinsFlux+NbinsXsec + NbinsAna*(il+1);
	}

	//=== Draw an horizontal line for det samples
	TLine *verlineSample[Nsamples];
	for(int il = 0; il < Nsamples; il++)
	{
		verlineSample[il] = new TLine(boundarySample[il], 2*NbinsAna+NbinsFlux+NbinsXsec, boundarySample[il], Nbins );
		verlineSample[il] -> SetLineWidth(1);
		verlineSample[il] -> SetLineStyle(2);
	}

	//=== Draw a vertical line for det samples
	TLine *orlineSample[Nsamples];
	for(int il = 0; il < Nsamples; il++)
	{
		orlineSample[il] = new TLine(2*NbinsAna+NbinsFlux+NbinsXsec, boundarySample[il], Nbins, boundarySample[il] );
		orlineSample[il] -> SetLineWidth(1);
		orlineSample[il] -> SetLineStyle(2);
	}
	//====================================================================================================                      


	//====================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();

	//======================================================================================================

	TFile *findetcov = TFile::Open(infile.c_str());

	//==== setup enu bins and covm for flux 
	TMatrixDSym *cov_mat = (TMatrixDSym*)findetcov -> Get("res_cov_matrix");
	TMatrixDSym *cor_mat = (TMatrixDSym*)findetcov -> Get("res_cor_matrix");

	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetPadTopMargin(0.02);
	gStyle->SetPadBottomMargin(0.08);
	// gStyle->SetHistMinimumZero(kFALSE);
	// gStyle->SetPaintTextFormat("5.3f");

	// Set a custom color palette
	const Int_t NRGBs = 3;
	const Int_t NCont = 20;
	Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 1.00, 1.00 };
	Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
	Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);


	//=== Draw covariance matrix
	TCanvas *c2 = new TCanvas("c2","c2",1000,900);
	c2 -> Draw();
	cov_mat -> Draw("colz");

	for(int il=0; il<2*Nlines; il++) orline[il]  -> Draw();
	for(int il=0; il<2*Nlines; il++) verline[il] -> Draw();

	for(int il=0; il<Nsamples; il++) orlineSample[il]  -> Draw();
	for(int il=0; il<Nsamples; il++) verlineSample[il] -> Draw();


	//=== Draw correlation matrix
	TCanvas *c3 = new TCanvas("c3","c3",1000,900);
	c3 -> Draw();
	cor_mat -> Draw("colz");

	for(int il=0; il<2*Nlines; il++) orline[il]  -> Draw();
	for(int il=0; il<2*Nlines; il++) verline[il] -> Draw();

	for(int il=0; il<Nsamples; il++) orlineSample[il]  -> Draw();
	for(int il=0; il<Nsamples; il++) verlineSample[il] -> Draw();


	// c2->Print("plots/covariancematrices/FinalCovMatrix.pdf");
	c3->Print(Form("plots/covariancematrices/FinalCorrMatrix_%s.pdf", inputname.c_str()));
	c3->Print(Form("plots/covariancematrices/FinalCorrMatrix_%s.png", inputname.c_str()));

}
