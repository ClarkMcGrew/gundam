// Usage:
// root -b -q 'DrawCovariancePostfit_noMichel.C+("fit2_dataCS_fakedataSignal_noMichel", "data")'

#include "CommonHeader.h"
#include "CommonStyle.h"
#include "BinningTools.cc"

void DrawCovariancePostfit_noMichel(string inputname = "fit2_dataCS_fakedataSignal_noMichel", const std::string& dir_name = "data")
{
	
	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();
	//======================================================================================================  


	string infile = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/%s.root", inputname.c_str());

	TFile *findetcov = TFile::Open(infile.c_str());

	//==== setup enu bins and covm for flux 
	TMatrixDSym *cov_mat = (TMatrixDSym*)findetcov -> Get("res_cov_matrix");
	TMatrixDSym *cor_mat = (TMatrixDSym*)findetcov -> Get("res_cor_matrix");

	
	int NbinsAna  = 58;
	int NbinsFlux = 20;  
	int NbinsXsec = 26;
	int NbinsDet[] = {58, 58, 58, 10, 16, 58, 58,
		              58, 58, 58, 10, 16, 58, 58,
		              58, 58, 58, 10, 16, 58, 58 };

	int Nbins = cov_mat -> GetNrows();

	

	//====================================================================================================                      
	//=== Set boundary lines between different types of parameters   

	const int Nlines = 4;
	double boundary[Nlines+1];
	boundary[0] = 0;
	boundary[1] = 2*NbinsAna;
	boundary[2] = 2*NbinsAna + NbinsFlux;
	boundary[3] = 2*NbinsAna + NbinsFlux + NbinsXsec;
	boundary[4] = Nbins;
	
	//=== Draw an horizontal line for parameter type
	TLine *verline[2*Nlines];
	for(int il = 0; il < Nlines; il++)
	{
		verline[il] = new TLine(boundary[il], boundary[il], boundary[il], boundary[il+1] );
		verline[il] -> SetLineWidth(2);
		verline[il] -> SetLineStyle(1);
		verline[il+Nlines] = new TLine(boundary[il+1], boundary[il], boundary[il+1], boundary[il+1] );
		verline[il+Nlines] -> SetLineWidth(2);
		verline[il+Nlines] -> SetLineStyle(1);
	}

	//=== Draw a vertical line for parameter type
	TLine *orline[2*Nlines];
	for(int il = 0; il < Nlines; il++)
	{
		orline[il] = new TLine(boundary[il], boundary[il], boundary[il+1], boundary[il] );
		orline[il] -> SetLineWidth(2);
		orline[il] -> SetLineStyle(1);
		orline[il+Nlines] = new TLine(boundary[il], boundary[il+1], boundary[il+1], boundary[il+1] );
		orline[il+Nlines] -> SetLineWidth(2);
		orline[il+Nlines] -> SetLineStyle(1);
	}


	//====================================================================================================                      
	//=== Set boundary lines between samples in detector covariance sub-matrix   

	const int Nsamples = 20;
	double boundarySample[Nsamples+1];
	boundarySample[0] = 2*NbinsAna + NbinsFlux + NbinsXsec + NbinsDet[0];
	for(int il = 1; il < Nsamples; il++)
	{
		boundarySample[il] = boundarySample[il-1] + NbinsDet[il];
	}

	//=== Draw an horizontal line for det samples
	TLine *verlineSample[Nsamples];
	for(int il = 0; il < Nsamples; il++)
	{
		verlineSample[il] = new TLine(boundarySample[il], 2*NbinsAna+NbinsFlux+NbinsXsec, boundarySample[il], Nbins );
		verlineSample[il] -> SetLineWidth(1);
		verlineSample[il] -> SetLineStyle(1);
		verlineSample[il] -> SetLineColor(kGray+2);
	}

	//=== Draw a vertical line for det samples
	TLine *orlineSample[Nsamples];
	for(int il = 0; il < Nsamples; il++)
	{
		orlineSample[il] = new TLine(2*NbinsAna+NbinsFlux+NbinsXsec, boundarySample[il], Nbins, boundarySample[il] );
		orlineSample[il] -> SetLineWidth(1);
		orlineSample[il] -> SetLineStyle(1);
		orlineSample[il] -> SetLineColor(kGray+2);
	}
	//====================================================================================================                      


	//====================================================================================================  
	//=== Draw matrices

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

	gStyle->SetLabelSize(0.04, "x"); 
	gStyle->SetLabelSize(0.04, "y"); 
	gROOT->ForceStyle();

	//=== Draw covariance matrix
	TCanvas *c2 = new TCanvas("c2","c2",1050,900);
	c2 -> Draw();
	cov_mat -> Draw("colz");

	for(int il=0; il<2*Nlines; il++) orline[il]  -> Draw();
	for(int il=0; il<2*Nlines; il++) verline[il] -> Draw();

	for(int il=0; il<Nsamples; il++) orlineSample[il]  -> Draw();
	for(int il=0; il<Nsamples; il++) verlineSample[il] -> Draw();


	//=== Draw correlation matrix
	TCanvas *c3 = new TCanvas("c3","c3",1050,900);
	c3 -> Draw();

	TH2D *cor_mat_th2 = new TH2D(*cor_mat);
	cor_mat_th2 -> SetMinimum(-1.0);
	cor_mat_th2 -> SetMaximum(1.0);
	
	cor_mat_th2 -> Draw("colz");

	if(cor_mat -> GetNrows() > 1000)
	{
		for(int il=0; il<2*Nlines; il++) orline[il]  -> Draw();
		for(int il=0; il<2*Nlines; il++) verline[il] -> Draw();

		for(int il=0; il<Nsamples; il++) orlineSample[il]  -> Draw();
		for(int il=0; il<Nsamples; il++) verlineSample[il] -> Draw();
	}

	c3->Print(Form("plots/fitteroutput/%s/PostFitCorrMatrix_%s.pdf", dir_name.c_str(), inputname.c_str()));
	c3->Print(Form("plots/fitteroutput/%s/PostFitCorrMatrix_%s.png", dir_name.c_str(), inputname.c_str()));


	// Draw matrix diagonal elements
	// for(int i=0; i<Nbins; i++)
	// 	std::cout << "Diag element "<<i<<" --> " << (*cov_mat)(i,i) << std::endl;

}
