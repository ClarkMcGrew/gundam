// Usage:
// root 'DrawCovarianceFinal.C+("fit3_statFluc")'

#include "CommonHeader.h"
#include "CommonStyle.h"
#include "BinningTools.cc"

void DrawCovarianceFinal(string inputname = "fit3_statFluc", const std::string& dir_name = "fakedata/statFluc")
{
	
	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gROOT->ForceStyle();
	//======================================================================================================  


	string infile = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());

	TFile *findetcov = TFile::Open(infile.c_str());

	//==== setup enu bins and covm for flux 
	TMatrixDSym *cov_mat = (TMatrixDSym*)findetcov -> Get("xsec_cov");
	TMatrixDSym *cor_mat = (TMatrixDSym*)findetcov -> Get("xsec_cor");

	int NbinsCos[] = {1, 5, 6, 6, 7, 8, 7, 10, 8,
					  1, 5, 6, 6, 7, 8, 7, 10, 8 };

	int Nbins = cov_mat -> GetNrows();

	

	//====================================================================================================                      
	//=== Set boundary lines between different types of parameters   

	const int Nlines = 18;
	double boundary[Nlines+1];
	boundary[0] = NbinsCos[0];
	for(int il = 1; il < Nlines; il++)
	{
		boundary[il] = boundary[il-1] + NbinsCos[il];
	}

	//=== Draw an horizontal line for det samples
	TLine *verline[Nlines];
	for(int il = 0; il < Nlines; il++)
	{
		verline[il] = new TLine(boundary[il], 0, boundary[il], Nbins );
		verline[il] -> SetLineWidth(1);
		verline[il] -> SetLineStyle(2);
	}

	//=== Draw a vertical line for det samples
	TLine *orline[Nlines];
	for(int il = 0; il < Nlines; il++)
	{
		orline[il] = new TLine(0, boundary[il], Nbins, boundary[il] );
		orline[il] -> SetLineWidth(1);
		orline[il] -> SetLineStyle(2);
	}

	TLine *verlineTar = new TLine(Nbins/2, 0, Nbins/2, Nbins);
	TLine *orlineTar  = new TLine(0, Nbins/2, Nbins, Nbins/2);
	verlineTar -> SetLineWidth(3);
	orlineTar  -> SetLineWidth(3);
	verlineTar -> SetLineStyle(2);
	orlineTar  -> SetLineStyle(2);
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

	for(int il=0; il<Nlines; il++) orline[il]  -> Draw();
	for(int il=0; il<Nlines; il++) verline[il] -> Draw();
	verlineTar ->Draw();
	orlineTar  ->Draw();

	c2->Print(Form("plots/xsecResults/%s/FinalXsecCovMatrix_%s.png", dir_name.c_str(), inputname.c_str()));


	// //=== Draw correlation matrix
	TCanvas *c3 = new TCanvas("c3","c3",1050,900);
	c3 -> Draw();

	TH2D *cor_mat_th2= new TH2D(*cor_mat);
	cor_mat_th2 -> SetMinimum(-1.0);
	cor_mat_th2 -> SetMaximum(1.0);

	cor_mat_th2 -> Draw("colz");

	for(int il=0; il<Nlines; il++) orline[il]  -> Draw();
	for(int il=0; il<Nlines; il++) verline[il] -> Draw();
	verlineTar ->Draw();
	orlineTar  ->Draw();

	c3->Print(Form("plots/xsecResults/%s/FinalXsecCorrMatrix_%s.pdf", dir_name.c_str(), inputname.c_str()));
	c3->Print(Form("plots/xsecResults/%s/FinalXsecCorrMatrix_%s.png", dir_name.c_str(), inputname.c_str()));
	

	// Draw matrix diagonal elements
	// for(int i=0; i<Nbins; i++)
	// 	std::cout << "Diag element "<<i<<" --> " << (*cov_mat)(i,i) << std::endl;


}
