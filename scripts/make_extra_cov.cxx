void make_extra_cov(const std::string& fname_output, unsigned int num_samples)
{
    TMatrixDSym extra_cov(num_samples);
    extra_cov.Zero();

    for(unsigned int i = 0; i < num_samples; ++i)
    {
        extra_cov(i,i) = 0.25;
    }

    TFile* f = TFile::Open(fname_output.c_str(), "RECREATE");
    f->cd();

    extra_cov.Write("extra_cov");
    f->Close();
}
