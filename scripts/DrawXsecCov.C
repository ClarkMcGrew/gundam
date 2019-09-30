// Usage:
// root 'DrawXsecCov.C+()'

#include "CommonHeader.h"
#include "CommonStyle.h"
#include "BinningTools.cc"

void DrawXsecCov(	int Nbins = 26,
						string infile = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/xsllh_xseccovmat.root",
						string matrixname = "xsec_cov")
{

	//====================================================================================================                      
	//=== Set boundary lines between samples                                                                               
	int Nboundary = 2;

	double boundary[Nboundary];

	boundary[0] = 9;
	boundary[1] = 20;
	

	//=== Draw an horizontal line for each angular boundaries                                                                   
	TLine *orline[Nboundary];
	for(int il = 0; il < Nboundary; il++)
	{
		orline[il] = new TLine(boundary[il], 0, boundary[il], Nbins );
		orline[il] -> SetLineWidth(2);
		orline[il] -> SetLineStyle(2);
	}

	//=== Draw a vertical line for each angular boundaries                                                                     
	TLine *verline[Nboundary];
	for(int il = 0; il < Nboundary; il++)
	{
		verline[il] = new TLine(0, boundary[il], Nbins, boundary[il] );
		verline[il] -> SetLineWidth(2);
		verline[il] -> SetLineStyle(2);
	}



	//====================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	// CommonStyle();
	gStyle->SetOptStat(0);

	//======================================================================================================

	TFile *findetcov = TFile::Open(infile.c_str());

	//==== setup enu bins and covm for flux 
	TMatrixDSym *cov_mat = (TMatrixDSym*)findetcov -> Get(matrixname.c_str());

	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetPadTopMargin(0.02);
	gStyle->SetPadBottomMargin(0.1);
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
	TCanvas *c2 = new TCanvas("c2","c2",1800,800);
	c2 -> Draw();

	gROOT->ForceStyle();

	TH2D *cov_mat_th2= new TH2D(*cov_mat);
	cov_mat_th2 -> SetMinimum(-1.0);
	cov_mat_th2 -> SetMaximum(1.0);

	cov_mat_th2 -> Draw("text colz");

	for(int il=0; il<Nboundary; il++) orline[il]  -> Draw();
	for(int il=0; il<Nboundary; il++) verline[il] -> Draw();

	// //=== Compute correlation matrix
	// int dim_det = cov_mat->GetNcols();

	// TMatrixDSym cor_mat(dim_det);
	// for(int i=0; i<dim_det; i++)
	// {
	//   for(int j=0; j<dim_det; j++)
	//   {
	//     cor_mat(i,j) = (*cov_mat)(i,j)/(TMath::Sqrt((*cov_mat)(i,i))*TMath::Sqrt((*cov_mat)(j,j)));
	//   }
	// }

	//=== Draw correlation matrix
	// TCanvas *c3 = new TCanvas("c3","c3",1000,900);
	// c3 -> Draw();
	// cor_mat -> Draw("text colz");

	// for(int il=0; il<Nboundary; il++) orline[il]  -> Draw();
	// for(int il=0; il<Nboundary; il++) verline[il] -> Draw();

	c2 -> Print("plots/covariancematrices/XsecCovMatrix.pdf");
	// c3 -> Print("plots/covariancematrices/DetCorrMatrix.pdf");

}
