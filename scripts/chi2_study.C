/*******************************************
* Author : Lucie Maret
* mail : lucie.maret@cern.ch
*******************************************/
/*

root -b -q 'chi2_study.C("fit1_statFluc")'
root -b -q 'chi2_study.C("fit3_statFluc")'

*/

#include "CommonStyle.h"
#include "CommonHeader.h"


// Compute the pull for each parameter
void chi2_study(string fit_type = "fit3_statFluc")
{
	
	//======================================================================================================  
	//=== Set common style
	CommonStyle();
	gStyle->SetOptStat(1);
	gROOT->ForceStyle();


	std::cout << "------------------------" << std::endl;
	std::cout << "----- Begin script -----" << std::endl;

	//ALL int toys[] = {1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200};
	// int toys_fit1[] = {1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200};
	// int toys_fit3[] = {1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200};
	
	int toys_fit1[] = {10 , 100, 101, 102, 103, 104, 105, 106, 107, 11 , 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 12 , 120, 121, 122, 123, 124, 126, 127, 128, 129, 13 , 130, 131, 132, 134, 135, 136, 137, 138, 139, 14 , 140, 141, 142, 143, 145, 146, 147, 148, 149, 15 , 150, 151, 152, 153, 154, 156, 157, 158, 159, 16 , 161, 162, 164, 165, 166, 167, 168, 169, 17 , 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 18 , 180, 181, 183, 185, 187, 189, 19 , 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 2  , 20 , 200, 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 3  , 30 , 31 , 32 , 33 , 34 , 35 , 36 , 37 , 38 , 39 , 4  , 40 , 41 , 42 , 43 , 45 , 47 , 48 , 49 , 5  , 50 , 51 , 52 , 53 , 54 , 55 , 56 , 57 , 58 , 59 , 6  , 61 , 64 , 65 , 66 , 67 , 68 , 69 , 7  , 70 , 71 , 72 , 73 , 74 , 75 , 76 , 78 , 79 , 8  , 80 , 81 , 82 , 84 , 85 , 86 , 87 , 88 , 89 , 9  , 90 , 91 , 93 , 94 , 95 , 96 , 97 , 98 , 99 };
	int toys_fit3[] = {1, 10, 100, 101, 102, 103, 104, 105, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 4, 40, 43, 44, 45, 46, 47, 48, 49, 5, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 6, 60, 61, 62, 63, 64, 66, 67, 68, 69, 7, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 8, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 9, 90, 91, 92, 93, 94, 95, 97, 98, 99 };
	

	int ntoys;
	if(fit_type == "fit1_statFluc")      ntoys = (sizeof(toys_fit1)/sizeof(*toys_fit1));
	else if(fit_type == "fit3_statFluc") ntoys = (sizeof(toys_fit3)/sizeof(*toys_fit3));

	std::cout << "----- Number of toys used = "<<ntoys<<" -----" << std::endl;

	string GenFilename_fits = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toys/%s_toy", fit_type.c_str());
	string GenFilename_xsec = Form("/sps/t2k/lmaret/softwares/xsLLhFitterLM/outputs/toys/xsec_%s_toy", fit_type.c_str());

	double bin_min_minuit = 500;
	double bin_max_minuit = 1500;
	double bin_min_xsec   = 50;
	double bin_max_xsec   = 150;

	std::vector< double > chi2_xsec;
	std::vector< double > chi2_minuit;

	TH1D *h_chi2_minuit = new TH1D("h_chi2_minuit", "MINUIT chi2 distribution of toys", 50, bin_min_minuit, bin_max_minuit);
	TH1D *h_chi2_xsec   = new TH1D("h_chi2_xsec",   "Final chi2 distribution of toys",  50, bin_min_xsec,   bin_max_xsec);


	// =================================================================
	std::cout << "----- Compute the MINUIT fit chi2 -----" << std::endl;
	
	for(int itoy = 0; itoy<ntoys; itoy++)
	{
		string filename;
		if(fit_type == "fit1_statFluc")      filename = Form("%s%d.root", GenFilename_fits.c_str(), toys_fit1[itoy]);
		else if(fit_type == "fit3_statFluc") filename = Form("%s%d.root", GenFilename_fits.c_str(), toys_fit3[itoy]);
		TFile* fin = new TFile(filename.c_str());
		
		TH1D* chi2_per_iter    = (TH1D*)(fin->Get("chi2_tot_periter"));
		double chi2_minuit_toy = chi2_per_iter -> GetBinContent(chi2_per_iter->GetNbinsX()-1);

		chi2_minuit.push_back(chi2_minuit_toy);
		h_chi2_minuit -> Fill(chi2_minuit_toy);

		std::cout << "===== Chi2 for toy "<<itoy<<" = " << chi2_minuit[itoy] << std::endl;

		delete fin;
	}



	// ========================================================================================
	std::cout << "----- Compute the final xsec chi2 w.r.t. fake data truth -----" << std::endl;

	for(int itoy = 0; itoy<ntoys; itoy++)
	{
		string filename;
		if(fit_type == "fit1_statFluc")      filename = Form("%s%d.root", GenFilename_xsec.c_str(), toys_fit1[itoy]);
		else if(fit_type == "fit3_statFluc") filename = Form("%s%d.root", GenFilename_xsec.c_str(), toys_fit3[itoy]);
		// TFile* fin = new TFile(filename.c_str());

		// TMatrixDSym *cov_mat = (TMatrixDSym*)(fin -> Get("xsec_cov"));
		// TH1D* h_sel_best_fit = (TH1D*)(fin->Get("sel_best_fit"));
		// TH1D* h_tru_best_fit = (TH1D*)(fin->Get("fake_data_concat"));

		// double chi2_xsec_toy = calcChi2_M(h_sel_best_fit, h_tru_best_fit, *cov_mat);
		double chi2_xsec_toy = calc_chisq(filename);

		chi2_xsec.push_back(chi2_xsec_toy);
		h_chi2_xsec -> Fill(chi2_xsec_toy);

		std::cout << "===== Chi2 for toy "<<itoy<<" = " << chi2_xsec[itoy] << std::endl;

		// delete fin;
	}


	// ==================================================
	std::cout << "----- Plot results -----" << std::endl;

	int dof_fit = 1322;
	int dof_xsec = 116;


	TCanvas *c_chi2_minuit = new TCanvas("c_chi2_minuit","Pull distribution for parameters",900,700);
	h_chi2_minuit -> SetTitle(Form(";MINUIT fit #chi^{2}; "));
	c_chi2_minuit -> SetGrid();
	h_chi2_minuit -> GetYaxis() -> SetRangeUser(0.0, 1.3 * h_chi2_minuit -> GetMaximum());
	h_chi2_minuit -> GetXaxis() -> SetNdivisions(5);

	TF1* f_chisq = new TF1("f_chisq", "ROOT::Math::chisquared_pdf(x,[0],0)",0,1500);
	f_chisq -> SetParameter(0, dof_fit);

	// double norm_factor_minuit = f_chisq->GetMaximum() / h_chi2_minuit->GetMaximum();
	double norm_factor_minuit = 1.0 / (h_chi2_minuit->Integral());

	h_chi2_minuit -> Scale( norm_factor_minuit );

	h_chi2_minuit -> Draw();
	f_chisq -> Draw("same");

	c_chi2_minuit -> Print(Form("plots/pullstudy/chi2_postfit_distribution_%s.pdf", fit_type.c_str()));




	TCanvas *c_chi2_xsec = new TCanvas("c_chi2_xsec","Pull distribution for parameters",900,700);
	h_chi2_xsec -> SetTitle(Form(";Final xsec #chi^{2}; "));
	c_chi2_xsec -> SetGrid();
	h_chi2_xsec -> GetYaxis() -> SetRangeUser(0.0, 1.3 * h_chi2_xsec -> GetMaximum());
	h_chi2_xsec -> GetXaxis() -> SetNdivisions(5);

	h_chi2_xsec -> Draw();

	c_chi2_xsec -> Print(Form("plots/pullstudy/chi2_xsec_distribution_%s.pdf", fit_type.c_str()));


	std::cout << "----- End script -----" << std::endl;
	std::cout << "----------------------" << std::endl;

}




double calcChi2_M(TH1D* h1, TH1D* h2, TMatrixD covar)
{
	double chi2=0;

	cout << "calcChi2_M::Starting chi2 calc" << endl;
	cout << "calcChi2_M::h1 (data) nbins: " << h1->GetXaxis()->GetNbins() << endl;
	cout << "calcChi2_M::h2 (MC) nbins: "   << h2->GetXaxis()->GetNbins() << endl;
	cout << "calcChi2_M::matrix rows: "     << covar.GetNrows() << endl;

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
			cout << "  Input matrix: determinant is: " << covar.Determinant() <<  endl;
			cout << "Will try adding 0.00000001pc to the diagonal ..." << endl;
			for(int i=0; i<covar.GetNrows(); i++)
			{
				covar[i][i]=covar[i][i]*1.00000001;
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
			// if(i==j)
			// {
			// 	cout << "calcChi2_M::bin " << i << ": " << h1->GetBinError(i+1) << " should be the same as " << sqrt(covar_origin[i][j]) << endl;
			// }
		}
	}
	cout << "calcChi2_M::chi2 is: " << chi2 << endl;
	return chi2;
}





double calc_chisq(std::string filename)
{
    // if(ROOT_VERSION_CODE < ROOT_VERSION(6,0,0))
    // {
    //     std::cout << "[ERROR]: ROOT version too old! Need ROOT 6 or higher." << std::endl;
    //     return -1;
    // }

    std::string h1_name("sel_best_fit");
    std::string h2_name("fake_data_concat");
    std::string cov_name("xsec_cov");

    unsigned int dof = 0;

    // std::cout << "Reading " << filename << std::endl;
    TFile* fin = TFile::Open(filename.c_str(), "READ");

    TH1D* h1 = (TH1D*)(fin -> Get(h1_name.c_str()));
    TH1D* h2 = (TH1D*)(fin -> Get(h2_name.c_str()));
    TMatrixDSym* cov_sym = (TMatrixDSym*)(fin -> Get(cov_name.c_str()));

    
    if(h1->GetNbinsX() != h2->GetNbinsX())
    {
        std::cout << "[ERROR]: Histograms bin numbers do not match!" << std::endl;
        return -1;
    }

    TMatrixD cov(*cov_sym);
    TMatrixD inv(*cov_sym);

    double det = 0;
    if(!TDecompLU::InvertLU(inv, 1E-100, &det))
    {
        std::cout << "[ERROR]: Failed to invert matrix." << std::endl;
        return -1;
    }

    double chisq = 0;
    const unsigned int nbins = cov.GetNrows();
    for(unsigned int i = 0; i < nbins; ++i)
    {
        for(unsigned int j = 0; j < nbins; ++j)
        {
            double x = h1->GetBinContent(i+1) - h2->GetBinContent(i+1);
            double y = h1->GetBinContent(j+1) - h2->GetBinContent(j+1);
            chisq += x * y * inv[i][j];
        }
    }

    std::cout << "Chisq = " << chisq << std::endl;
    delete fin;
    return chisq;
}