void plot_fit_output(const std::string& file_name_input)
{
    gStyle -> SetOptStat(0);

    TFile* file = TFile::Open(file_name_input.c_str(), "READ");

    // std::vector<std::string> par_name = {"par_fit", "par_flux", "par_xsec", "par_det"};
    // const int Npar = 4;
    // std::string par_name[Npar] = {"par_fit", "par_flux", "par_xsec", "par_det"};
    // const int Npar = 3;
    // std::string par_name[Npar] = {"par_fit", "par_flux", "par_det"};
    const int Npar = 1;
    std::string par_name[Npar] = {"par_fit"};

    std::string name;
    std::stringstream ss;

    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "------------- Drawing results -------------" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;

    // for(const auto& name : par_name)
    for(int j = 0; j<Npar; j++)
    {
        name = par_name[j];
        std::cout << "----- for " << name << " -----" << std::endl;
        
        ss.str("");
        ss << "hist_" << name << "_prior";
        TH1D* h_prior = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_result";
        TH1D* h_final = (TH1D*)file -> Get(ss.str().c_str());

        TCanvas *c = new TCanvas("c", "c", 1400, 900);
        h_prior -> SetLineColor(kBlue);
        h_prior -> SetLineWidth(2);
        h_final -> SetLineColor(kRed);
        h_final -> SetLineWidth(2);

        h_final -> Draw("hist");
        h_prior -> Draw("hist same");

        TLegend* legend = new TLegend(0.7,0.75,0.9,0.9);
        legend -> SetFillColor(0);
        legend -> AddEntry(h_prior, "Prior","l");
        legend -> AddEntry(h_final, "Final","l");
        legend -> Draw();

        ss.str("");
        ss << "../outputs/" << name << "_results.pdf";
        c -> Print(ss.str().c_str());
        delete c;
    }

    
    std::cout << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "------------- Drawing errors -------------" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // for(const auto& name : par_name)
    for(int j = 0; j<Npar; j++)
    {
        name = par_name[j];
        std::cout << "----- for " << name << " -----" << std::endl;

        ss.str("");
        ss << "hist_" << name << "_error_prior";
        TH1D* h_err_prior = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_error_final";
        TH1D* h_err_final = (TH1D*)file -> Get(ss.str().c_str());

        TCanvas *c = new TCanvas("c", "c", 1400, 900);
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
        ss << "../outputs/" << name << "_errors.pdf";
        c -> Print(ss.str().c_str());
        delete c;
    }


    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "------------- Drawing overlay -------------" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // for(const auto& name : par_name)
    for(int j = 0; j<Npar; j++)
    {
        name = par_name[j];
        std::cout << "----- for " << name << " -----" << std::endl;

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

        TCanvas *c = new TCanvas("c", "c", 1400, 900);
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

        TLegend* legend = new TLegend(0.7,0.75,0.9,0.9);
        legend -> SetFillColor(0);
        legend -> AddEntry(h_prior, "Prior","p");
        legend -> AddEntry(h_final, "Final","p");
        legend -> Draw();

        ss.str("");
        ss << "../outputs/" << name << "_overlay.pdf";
        c -> Print(ss.str().c_str());
        delete c;
    }

}
