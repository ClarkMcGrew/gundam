//////////////////////////////////////////////////////////////////////////////////////
// Macro to plot the fractional error for the xsec systematics
//
//
//  Created: 28 November 2016 
//
//  Modified: 
//
//
//  Usage: root -b -q 'EvaluateXsec.C+()'
//         root -b -q 'EvaluateXsec.C+("fit2_data", 1, 0)'
// 
// 
// 			firstBin = 0; // include the first costh bin
// 			firstBin = 1; // does NOT include the first costh bin for a restricted phase space
// 
// 			isData = 0; // post-fit values
// 			isData = 1; // data truth values
// 			isData = 2; // nominal truth values
// 
//
///////////////////////////////////////////////////////////////////////////////////////

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"


void EvaluateXsec(string inputname = "fit2_data", int firstBin = 0, int isData = 0, string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_format.txt")
{

	const int Ntarget = 2;
	string targetlist[Ntarget] = {"Carbon", "Oxygen"};

	string infilename = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/xsec_%s.root", inputname.c_str());

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
	
	vector< vector<TH1D*> > h_xsec(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
		{
			if(isData == 0)
				h_xsec[itar].push_back( (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_postfit", targetlist[itar].c_str(), nbcth) )) );
			else if(isData == 1)
				h_xsec[itar].push_back( (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_data", targetlist[itar].c_str(), nbcth) )) );
			else if(isData == 2)
				h_xsec[itar].push_back( (TH1D*)(fin -> Get( Form("CC0pi%s_cos_bin%d_nominal", targetlist[itar].c_str(), nbcth) )) );
			
		}

	vector< TH1D* > h_ratio_postfit;
	for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
	{
		if(isData == 0)
			h_ratio_postfit.push_back( (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_postfit", nbcth) )) );
		else if(isData == 1)
			h_ratio_postfit.push_back( (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_data", nbcth) )) );
		else if(isData == 2)
			h_ratio_postfit.push_back( (TH1D*)(fin -> Get( Form("CC0piCarbon_cos_bin%d_ratio_nominal", nbcth) )) );
	}
	//======================================================================================================



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Evaluate integrated cross-section for "<< inputname ;
		if     (isData == 0) std::cout << " post-fit value";
		else if(isData == 1) std::cout << " fakedata truth value";
		else if(isData == 2) std::cout << " nominal MC value";
		if     (firstBin==0) std::cout << " in full phase space";
		else if(firstBin==1) std::cout << " in restricted phase space";
	std::cout <<" =====" << std::endl;
	
	// Calculate the integrated cross-section per costheta slice
	vector< vector< double > > xsec_per_costh(Ntarget);

	for(int itar = 0; itar < Ntarget; itar++)
		for(int nbcth=0; nbcth<Nbins_costh; nbcth++)
			xsec_per_costh[itar].push_back(  h_xsec[itar][nbcth] -> Integral("WIDTH") );

	std::cout << "===== Get covariance matrices =====" << std::endl;
	TMatrixDSym *cov_mat       = (TMatrixDSym*)(fin -> Get("xsec_cov"));
	TMatrixDSym *cov_mat_ratio = (TMatrixDSym*)(fin -> Get("ratio_cov"));


	// Normalise each slice by costheta width
	// and calculate the total integrated cross section + error
	vector< double > xsec_int;
	vector< double > xsec_err;

	// TH1D* err_diag_C = new TH1D("err_diag_C", "err_diag_C", Nbins, 0, 58);
	// TH1D* err_diag_O = new TH1D("err_diag_O", "err_diag_O", Nbins, 0, 58);
	// TH1D* err_diag_ratio = new TH1D("err_diag_ratio", "err_diag_ratio", Nbins, 0, 58);

	for(int itar = 0; itar < Ntarget; itar++)
	{
		// Calculate cross section
		xsec_int.push_back(0.0);
		int nb = firstBin;
		for(int nbcth=firstBin; nbcth<Nbins_costh; nbcth++)
		{
			double binWidthCosth = bin.GetCosBinWidth_i(nb);
			xsec_int[itar] += xsec_per_costh[itar][nbcth] * binWidthCosth;
			nb += Nmombins[nbcth];
		}

		// Calculate error
		double err = 0.0;
		for(int ib = firstBin; ib < Nbins; ib++)
		{
			for(int jb = firstBin; jb < Nbins; jb++)
			{
				// if(ib==jb) // not include correlations
				// {
					{
						double binWidthMom   = bin.GetMomBinWidth_i(ib) * bin.GetMomBinWidth_i(jb);
						double binWidthCosth = bin.GetCosBinWidth_i(ib) * bin.GetCosBinWidth_i(jb);
						double err_ij = binWidthMom * binWidthCosth * (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins);
						err += err_ij;

						// if(ib==jb && itar==0)
						// {
						// 	err_diag_C     -> SetBinContent(ib+1, binWidthMom * binWidthCosth * (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins) );
						// 	err_diag_ratio -> SetBinContent(ib+1, binWidthMom * binWidthCosth * (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins) );
						// 	// err_diag_C     -> SetBinContent(ib+1, (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins) );
						// 	// err_diag_ratio -> SetBinContent(ib+1, (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins) );
						// }
						// else if(ib==jb && itar==1)
						// {
						// 	err_diag_O     -> SetBinContent(ib+1, binWidthMom * binWidthCosth * (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins) );
						// 	// err_diag_O     -> SetBinContent(ib+1, (*cov_mat)(ib + itar*Nbins, jb + itar*Nbins) );
						// err_diag_ratio -> SetBinContent(ib+1, err_diag_ratio->GetBinContent(ib+1) / err_diag_O->GetBinContent(ib+1) );
						// }
					}
				// }
			}
		}
		xsec_err.push_back( sqrt(err) );
	}


	
	// NOTE : divide by 1000 because the binning is in MeV but the bin content is in per GeV
	xsec_int[0] /= 1000;
	xsec_int[1] /= 1000;
	xsec_err[0] /= 1000;
	xsec_err[1] /= 1000;
	
	std::cout << "===== Integrated cross-section on carbon = (" << xsec_int[0]/1e-38 << " +/- " << xsec_err[0]/1e-38 << ") x 1e-38 cm^2/nucleon ("<< 100*(xsec_err[0])/(xsec_int[0]) <<" % error)" << std::endl;
	std::cout << "===== Integrated cross-section on oxygen = (" << xsec_int[1]/1e-38 << " +/- " << xsec_err[1]/1e-38 << ") x 1e-38 cm^2/nucleon ("<< 100*(xsec_err[1])/(xsec_int[1]) <<" % error)" << std::endl;
	



	//======================================================================================================
	std::cout << "================================================" << std::endl;
	std::cout << "===== Evaluate cross-section ratio =====" << std::endl;
	
	double ratio_int = xsec_int[1]/xsec_int[0];

	// without O / C correlation :
	double ratio_err_noCorr = sqrt( xsec_err[1]*xsec_err[1]/(xsec_int[0]*xsec_int[0]) + xsec_err[0]*xsec_err[0] * ( xsec_int[1]*xsec_int[1]/(xsec_int[0]*xsec_int[0]*xsec_int[0]*xsec_int[0]) ) );

	vector< double > derivative;

	for(int i = 0; i < Nbins; i++)
	{
		double binWidthMom   = bin.GetMomBinWidth_i(i) / 1000.0;
		double binWidthCosth = bin.GetCosBinWidth_i(i);
		derivative.push_back( 1.0 * binWidthMom*binWidthCosth / xsec_int[0] );
		// std::cout << "i="<<i<<" binWidthMom   : " << binWidthMom <<", binWidthCosth : " << binWidthCosth << std::endl;
	}

	for(int i = 0; i < Nbins; i++)
	{
		double binWidthMom   = bin.GetMomBinWidth_i(i) / 1000.0;
		double binWidthCosth = bin.GetCosBinWidth_i(i);
		derivative.push_back( -1.0 * binWidthMom*binWidthCosth * xsec_int[1] / (xsec_int[0]*xsec_int[0]) );
	}

	// Calculate error
	double ratio_err;

	double err = 0.0;
	for(int i = 0; i < 2*Nbins; i++)
		for(int j = 0; j < 2*Nbins; j++)
			if(firstBin==0 || (i!=0 && i!=Nbins && j!=0 && j!=Nbins) )
			{
				err += derivative[i] * (*cov_mat)(i, j) * derivative[j];
				// std::cout << "CONTRIBUTE i = " << i << ", j = " << j << std::endl;
			}

	ratio_err = sqrt(err);

	std::cout << "===== Integrated cross-section ratio = " << ratio_int << " +/- " << ratio_err << " ("<< 100*ratio_err/ratio_int <<" % error)" << std::endl;
	          std::cout << "===== (no correlations: +/- " << ratio_err_noCorr << ", "<< 100*ratio_err_noCorr/ratio_int <<" % error)" << std::endl;
	

	//======================================================================================================



	//======================================================================================================

	// std::cout << "================================================" << std::endl;
	// for(int ib = 0; ib < Nbins; ib++)
	// {
	// 	std::cout << "Errors in bin i="<<ib<<" are :  " << std::endl
	// 	                                        << "      C   -> " << sqrt((*cov_mat)(ib, ib)) << std::endl
	// 	                                        << "      O   -> " << sqrt((*cov_mat)(ib+Nbins, ib+Nbins)) << std::endl
	// 	                                        << "      O/C -> " << sqrt((*cov_mat_ratio)(ib, ib)) << std::endl;
	// }
	// std::cout << "================================================" << std::endl;
	
	// for(int ib = 0; ib < Nbins; ib++)
	// {
	// 	std::cout << "Relative errors in bin i="<<ib<<" are : " << std::endl
	// 	                                         << "      C   -> " << sqrt((*cov_mat)(ib, ib)) / xsec_int[0] << std::endl
	// 	                                         << "      O   -> " << sqrt((*cov_mat)(ib+Nbins, ib+Nbins)) / xsec_int[1] << std::endl
	// 	                                         << "      O/C -> " << sqrt((*cov_mat_ratio)(ib, ib)) / ratio_int << std::endl;
	// }
	// std::cout << "================================================" << std::endl;
	

	// //======================================================================================================
	// std::cout << "================================================" << std::endl;
	
	//======================================================================================================

	// TCanvas *c_err = new TCanvas("c_err", "c_err", 900, 600);

	// err_diag_C -> SetLineColor(kRed);

	// err_diag_C -> Draw("hist");
	// err_diag_O -> Draw("same hist");



	// TCanvas *c_ratio = new TCanvas("c_ratio", "c_ratio", 900, 600);

	// err_diag_ratio -> Draw("hist");


	//======================================================================================================
	delete fin;

	std::cout << "================================================" << std::endl;
	std::cout << "===== End of script =====" << std::endl;
	std::cout << "================================================" << std::endl;                                                                                               
	//======================================================================================================

}
