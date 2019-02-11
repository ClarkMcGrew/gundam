// Usage:
// root 'DrawDetCov.C+()'

#include "CommonHeader.h"
#include "CommonStyle.h"
#include "BinningTools.cc"

void DrawDetCov(	int NsamplesTot = 24,
					string infile = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/xsllh_detcovmat.root",
					string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/tn338_binning_format.txt")
{

	//====================================================================================================                      
	//=== Set binning using BinningTools class                                                                                
	BinningTools bin;
	bin.SetBinning(fbinning.c_str());
	int Nbins = bin.GetNbins();//Total number of bins
	int Nsamples = NsamplesTot - 1;
	int Nfgds = 3;                                                             

	//====================================================================================================                      
	//=== Set boundary lines between samples                                                                               
	double boundary[Nsamples];
	for(int il=0; il<Nsamples; il++)
	{
		boundary[il] = Nbins * (il+1);
	}

	//=== Draw an horizontal line for each angular boundaries                                                                   
	TLine *orline[Nsamples];
	for(int il = 0; il < Nsamples; il++)
	{
		orline[il] = new TLine(boundary[il], 0, boundary[il], Nbins*(Nsamples+1) );
		orline[il] -> SetLineWidth(1);
		orline[il] -> SetLineStyle(2);
	}

	//=== Draw a vertical line for each angular boundaries                                                                     
	TLine *verline[Nsamples];
	for(int il = 0; il < Nsamples; il++)
	{
		verline[il] = new TLine(0, boundary[il], Nbins*(Nsamples+1), boundary[il] );
		verline[il] -> SetLineWidth(1);
		verline[il] -> SetLineStyle(2);
	}


	//====================================================================================================                      
	//=== Set boundary lines between FGDs                                                                               
	double boundaryfgd[Nfgds-1];
	int Nsamplesperfgd = NsamplesTot/Nfgds;
	for(int il=0; il<Nfgds-1; il++)
	{
		boundaryfgd[il] = Nsamplesperfgd * Nbins * (il+1);
	}

	//=== Draw an horizontal line for each angular boundaries                                                                   
	TLine *orlinefgd[Nfgds-1];
	for(int il = 0; il < Nfgds-1; il++)
	{
		orlinefgd[il] = new TLine(boundaryfgd[il], 0, boundaryfgd[il], Nbins*(Nsamples+1) );
		orlinefgd[il] -> SetLineWidth(3);
		orlinefgd[il] -> SetLineStyle(2);
	}

	//=== Draw a vertical line for each angular boundaries                                                                     
	TLine *verlinefgd[Nfgds-1];
	for(int il = 0; il < Nfgds-1; il++)
	{
		verlinefgd[il] = new TLine(0, boundaryfgd[il], Nbins*(Nsamples+1), boundaryfgd[il] );
		verlinefgd[il] -> SetLineWidth(3);
		verlinefgd[il] -> SetLineStyle(2);
	}

	//====================================================================================================                      


	//====================================================================================================  
	//=== Set T2K style. In CommonStyle.h option 1 has been setted
	CommonStyle();

	//======================================================================================================

	TFile *findetcov = TFile::Open(infile.c_str());

	//==== setup enu bins and covm for flux 
	TMatrixDSym *cov_mat = (TMatrixDSym*)findetcov -> Get("cov_mat");
	TMatrixDSym *cor_mat = (TMatrixDSym*)findetcov -> Get("cor_mat");

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

	for(int il=0; il<Nsamples; il++) orline[il]  -> Draw();
	for(int il=0; il<Nsamples; il++) verline[il] -> Draw();

	for(int il=0; il<Nfgds-1; il++) orlinefgd[il]  -> Draw();
	for(int il=0; il<Nfgds-1; il++) verlinefgd[il] -> Draw();

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
	TCanvas *c3 = new TCanvas("c3","c3",700,600);
	c3 -> Draw();
	cor_mat -> Draw("colz");

	// cor_mat -> SetMinimum(-1);

	for(int il=0; il<Nsamples; il++) orline[il]  -> Draw();
	for(int il=0; il<Nsamples; il++) verline[il] -> Draw();

	for(int il=0; il<Nfgds-1; il++) orlinefgd[il]  -> Draw();
	for(int il=0; il<Nfgds-1; il++) verlinefgd[il] -> Draw();


	c2->Print("plots/covariancematrices/DetCovMatrix.pdf");
	c3->Print("plots/covariancematrices/DetCorrMatrix.pdf");
	c3->Print("plots/covariancematrices/DetCorrMatrix.png");

}
