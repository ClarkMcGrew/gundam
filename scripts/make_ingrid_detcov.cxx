void make_ingrid_detcov(const std::string& nd280_cov_name, const std::string& ingrid_cov_name = "")
{
    TFile* nd280_file = TFile::Open(nd280_cov_name.c_str(), "READ");
    TMatrixTSym<double>* nd280_cov = (TMatrixTSym<double>*)nd280_file -> Get("cov_mat");

    TFile* ingrid_file = TFile::Open(ingrid_cov_name.c_str(), "READ");
    TMatrixTSym<double>* ingrid_cov = (TMatrixTSym<double>*)ingrid_file -> Get("cov_reco");

    unsigned int ingrid_bins = ingrid_cov -> GetNrows();
    unsigned int nd280_bins = nd280_cov -> GetNrows();
    unsigned int total_bins = nd280_bins + ingrid_bins;

    TMatrixTSym<double> ingrid_cov_full(ingrid_bins * 2);
    ingrid_cov_full.Zero();

    total_bins += ingrid_bins;
    for(int i = 0; i < ingrid_bins; ++i)
    {
        for(int j = 0; j < ingrid_bins; ++j)
        {
            ingrid_cov_full(i,j) = (*ingrid_cov)(i,j);
            ingrid_cov_full(i+ingrid_bins,j+ingrid_bins) = (*ingrid_cov)(i,j);
        }
    }

    //TMatrixTSym<double> ingrid_cov(ingrid_bins);
    //ingrid_cov.Zero();
    TMatrixTSym<double> combined_cov(total_bins);
    combined_cov.Zero();
    TMatrixTSym<double> combined_cor(total_bins);
    combined_cor.Zero();

    //const double ingrid_err = 0.079;
    //for(int i = 0; i < ingrid_bins; ++i)
    //    ingrid_cov(i,i) = ingrid_err * ingrid_err;

    for(int i = 0; i < total_bins; ++i)
    {
        for(int j = 0; j < total_bins; ++j)
        {
            if(i < nd280_bins && j < nd280_bins)
                combined_cov(i,j) = (*nd280_cov)(i,j);

            if(i >= nd280_bins && j >= nd280_bins)
                combined_cov(i,j) = ingrid_cov_full(i-nd280_bins, j-nd280_bins);

            if(i == j)
                combined_cov(i,j) += 1E-6;
        }
    }

    for(int i = 0; i < total_bins; ++i)
    {
        for(int j = 0; j < total_bins; ++j)
        {
            const double x = combined_cov(i, i);
            const double y = combined_cov(j, j);
            const double z = combined_cov(i, j);
            combined_cor(i, j) = z / (sqrt(x * y));

            if(std::isnan(combined_cor(i, j)))
                combined_cor(i, j) = 0.0;
        }
    }

    TFile* output_file = TFile::Open("combined_detcov.root", "RECREATE");
    output_file -> cd();

    combined_cov.Write("det_cov");
    combined_cor.Write("det_cor");
    ingrid_cov->Write("ingrid_cov");
    nd280_cov->Write("nd280_cov");

    output_file -> Close();
    return;
}
