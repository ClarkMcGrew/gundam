// Usage :
// root 'plot_fit_output.C("fit1_statFluc")'

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void plot_fit_output(const std::string& file_name = "fit1_statFluc")
{
    gStyle -> SetOptStat(0);

    const std::string output_dir = "histos/fitteroutput/asimov/";
    const std::string file_name_input = Form("../outputs/%s.root", file_name.c_str());

    std::string name;
    std::stringstream ss;


    TFile* file = TFile::Open(file_name_input.c_str(), "READ");

    const int Npar = 4;
    std::string par_name[Npar] = {"par_fit", "par_flux", "par_xsec", "par_det"};
    //std::string par_name[Npar] = {"par_fit", "par_xsec"};
    // const int Npar = 1;
    // std::string par_name[Npar] = {"par_fit"};



    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "-------- Drawing result and errors --------" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // for(const auto& name : par_name)
    for(int j = 0; j<Npar; j++)
    {
        name = par_name[j];
        std::cout << "----- for " << name << " -----" << std::endl;

        std::cout << "----- get histograms -----" << std::endl;

        ss.str("");
        ss << "hist_" << name << "_prior";
        TH1D* h_prior = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_result";
        TH1D* h_final = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_error_prior";
        TH1D* h_err_prior = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_error_final";
        TH1D* h_err_final = (TH1D*)file -> Get(ss.str().c_str());

        for(unsigned int i = 1; i <= h_prior->GetNbinsX(); ++i)
        {
            h_prior -> SetBinError(i, h_err_prior -> GetBinContent(i));
            h_final -> SetBinError(i, h_err_final -> GetBinContent(i));
        }


        std::cout << "----- Draw parameters and errors -----" << std::endl;

        TCanvas *c = new TCanvas("c", "c", 1400, 900);
        gStyle -> SetOptTitle(0);

        h_prior -> SetMarkerColor(kBlue);
        h_prior -> SetMarkerStyle(kFullCircle);
        h_final -> SetMarkerColor(kRed);
        h_final -> SetMarkerStyle(kFullCircle);

        h_prior -> SetFillColor(kBlue-9);
        h_prior -> SetFillStyle(1001);
        h_final -> SetFillColor(kRed-9);
        h_final -> SetFillStyle(3144);

        h_prior -> Draw("P E2");
        h_final -> Draw("P E2 same");

        if(name == "par_xsec") // for xsec parameters
            gPad->SetLogy();

        TLegend* legend = new TLegend(0.7,0.75,0.9,0.9);
        legend -> SetFillColor(0);
        legend -> AddEntry(h_prior, "Prior","p");
        legend -> AddEntry(h_final, "Final","p");
        legend -> Draw();

        ss.str("");
        ss << output_dir << name << "_overlay_" << file_name << ".pdf";
        c -> Print(ss.str().c_str());
        delete c;




        std::cout << "----- Draw errors only -----" << std::endl;

        // Divide by parameter value to have the relative error
        if(name == "par_xsec") // for xsec parameters
        {
	        h_err_prior -> Divide(h_prior);
	        h_err_final -> Divide(h_final);
	    }

        TCanvas *c_err = new TCanvas("c_err", "c_err", 1400, 900);
        gStyle -> SetOptTitle(0);

        double Max = MAX(h_err_prior->GetMaximum(), h_err_final->GetMaximum());
        h_err_final -> GetYaxis() -> SetRangeUser(0.0, 1.3*Max );

        h_err_prior -> SetLineColor(kBlue);
        h_err_prior -> SetLineWidth(2);
        h_err_final -> SetLineColor(kRed);
        h_err_final -> SetLineWidth(2);

        h_err_final -> Draw("hist");
        h_err_prior -> Draw("hist same");

        TLegend* legend = new TLegend(0.55,0.75,0.9,0.9);
        legend -> SetFillColor(0);
        legend -> AddEntry(h_err_prior, "Prior error","l");
        legend -> AddEntry(h_err_final, "Final error","l");
        legend -> Draw();

        ss.str("");
        ss << output_dir << name << "_errors_" << file_name << ".pdf";
        c_err -> Print(ss.str().c_str());
        delete c_err;
    }








    std::cout << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "------------- Drawing chi2 -------------" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::endl;
    
    
    TH1D* h_chi2_stat = (TH1D*)file -> Get("chi2_stat_periter");
    TH1D* h_chi2_sys  = (TH1D*)file -> Get("chi2_sys_periter");
    TH1D* h_chi2_reg  = (TH1D*)file -> Get("chi2_reg_periter");
    TH1D* h_chi2_tot  = (TH1D*)file -> Get("chi2_tot_periter");

    TCanvas *c = new TCanvas("c", "c", 1400, 900);
    gStyle -> SetOptTitle(0);

    // gPad->SetLogy();

    h_chi2_stat -> SetLineColor(kBlue+1);
    h_chi2_sys  -> SetLineColor(kRed+1);
    h_chi2_reg  -> SetLineColor(kGreen+2);
    h_chi2_tot  -> SetLineColor(kBlack);
    h_chi2_tot  -> SetLineWidth(2);

    double max = h_chi2_stat->GetMaximum();
    double min = h_chi2_sys->GetMinimum();
    std::cout << "min = " << min << ", max = " << max << std::endl;

    h_chi2_stat-> GetYaxis() -> SetRangeUser(1.1*min, 1.1*max);

    h_chi2_stat -> Draw("hist");
    h_chi2_sys  -> Draw("hist same");
    // h_chi2_reg  -> Draw("hist same");
    h_chi2_tot  -> Draw("hist same");

    TLegend* legend = new TLegend(0.55,0.75,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> AddEntry(h_chi2_stat, "chi2 stat","l");
    legend -> AddEntry(h_chi2_sys , "chi2 sys","l");
    // legend -> AddEntry(h_chi2_reg , "chi2 reg","l");
    legend -> AddEntry(h_chi2_tot , "chi2 total","l");
    legend -> Draw();

    ss.str("");
    ss << output_dir << "chi2_" << file_name << ".pdf";
    c -> Print(ss.str().c_str());
    // delete c;
    

    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::endl;

}
