// Usage :
// root 'DrawParameters.C("fit3_statFluc")'

#include "BinningTools.cc"
#include "CommonHeader.h"
#include "CommonStyle.h"


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void DrawParameters(const std::string& file_name = "fit3_statFluc", const std::string& dir_name = "fakedata/statFluc", string fbinning = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/binning/tn337_binning_GeV_format.txt")
{
    //======================================================================================================  
    //=== Set common style
    CommonStyle();
    gROOT->ForceStyle();

    bool drawParThrow = false;
    if(file_name == "fit3_statFluc") drawParThrow = true;

    const std::string output_dir = Form("plots/fitteroutput/%s/", dir_name.c_str());
    const std::string file_name_input = Form("../outputs/%s.root", file_name.c_str());

    std::string name;
    std::string title;
    std::stringstream ss;

    TFile* file = TFile::Open(file_name_input.c_str(), "READ");

    const int Npar = 4;
    std::string par_name[Npar]  = {"par_fit", "par_flux", "par_xsec", "par_det"};
    std::string par_title[Npar] = {"Template parameters", "Flux parameters", "Xsec model parameters", "Detector parameters"};


    // Get binning info
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



    std::cout << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "-------- Drawing parameter results --------" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    // for(const auto& name : par_name)
    for(int j = 0; j<Npar; j++)
    {
        name = par_name[j];
        title = par_title[j];

        std::cout << "----- for " << title << " (" << name << ") -----" << std::endl;

        // std::cout << "----- get histograms -----" << std::endl;

        ss.str("");
        ss << "hist_" << name << "_prior";
        TH1D* h_prior = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_iter0";
        TH1D* h_iter0 = new TH1D();
        if(drawParThrow) h_iter0 = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_result";
        TH1D* h_final = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_error_prior";
        TH1D* h_err_prior = (TH1D*)file -> Get(ss.str().c_str());

        ss.str("");
        ss << "hist_" << name << "_error_final";
        TH1D* h_err_final = (TH1D*)file -> Get(ss.str().c_str());

        for(int i = 1; i <= h_prior->GetNbinsX(); ++i)
        {
            h_prior -> SetBinError(i, h_err_prior -> GetBinContent(i));
            h_final -> SetBinError(i, h_err_final -> GetBinContent(i));
        }


        int division = 10;
        for(int i = 0; i < h_prior->GetNbinsX(); ++i)
        {
            if(par_name[j] == "par_det") division = 100;

            if(i%division == 0)
            {
                h_err_final -> GetXaxis() -> SetBinLabel(i+1, Form("%d", i));
                h_prior     -> GetXaxis() -> SetBinLabel(i+1, Form("%d", i));
            }
            else if(i%division != 0)
            {
                h_err_final -> GetXaxis() -> SetBinLabel(i+1, "");
                h_prior     -> GetXaxis() -> SetBinLabel(i+1, "");
            }
        }


        // std::cout << "----- Draw parameters and errors -----" << std::endl;

        TCanvas *c = new TCanvas("c", "c", 1400, 900);
        // gStyle -> SetOptTitle(0);

        h_prior -> SetTitle(Form(";%s;Parameter value", par_title[j].c_str() ));

        h_prior -> SetMarkerColor(kBlue+1);
        h_prior -> SetMarkerStyle(kFullSquare);
        if(drawParThrow) h_iter0 -> SetMarkerColor(kGreen+3);
        if(drawParThrow) h_iter0 -> SetMarkerStyle(kFullTriangleUp);
        h_final -> SetMarkerColor(kRed+1);
        h_final -> SetMarkerStyle(kFullCircle);

        h_prior -> SetMarkerSize(2);
        if(drawParThrow) h_iter0 -> SetMarkerSize(2);
        h_final -> SetMarkerSize(2);

        h_prior -> SetFillColor(kBlue-9);
        h_prior -> SetFillStyle(1001);
        h_final -> SetFillColor(kRed-9);
        h_final -> SetFillStyle(3144);

        if(name == "par_xsec") // for xsec parameters
            gPad->SetLogy();
        if(name == "par_flux") // for flux parameters
            h_prior -> GetYaxis() -> SetRangeUser(0.7,1.3);

        TLine *verline[2*Nbins_costh];
        int binSide = 0;
        for(int il = 0; il < Nbins_costh; il++)
        {
            binSide += Nmombins[il];
            verline[il]               = new TLine(binSide,         0, binSide,         2.1 );
            verline[il + Nbins_costh] = new TLine(binSide + Nbins, 0, binSide + Nbins, 2.1 );
        }

        h_prior -> Draw("P E2");

        if(name == "par_fit")
        {
            for(int il = 0; il < 2*Nbins_costh-1; il++)
            {
                verline[il] -> SetLineWidth(1);
                verline[il] -> SetLineStyle(1);
                verline[il] -> SetLineColor(kGray);
                verline[il] -> Draw();
            }
        }
        TLine *verlineTargets = new TLine(Nbins, 0, Nbins, 2.1 );
        if(name == "par_fit")
            verlineTargets -> Draw();


        if(drawParThrow) h_iter0 -> Draw("P same");
        h_final -> Draw("P E2 same");

        // TLegend* legend = new TLegend(0.55,0.67,0.9,0.9);
        TLegend* legend = new TLegend(0.75,0.75,0.9,0.9);
        legend -> SetFillColor(0);
        // legend -> SetHeader(Form("%s", title.c_str()));
        legend -> AddEntry(h_prior, "Prior","p");
        if(drawParThrow) legend -> AddEntry(h_iter0, "Initial","p");
        legend -> AddEntry(h_final, "Final","p");
        legend -> Draw();

        ss.str("");
        ss << output_dir << name << "_overlay_" << file_name << ".pdf";
        c -> Print(ss.str().c_str());
        delete c;




        // std::cout << "----- Draw errors only -----" << std::endl;

        // Divide by parameter value to have the relative error
        if(name == "par_xsec") // for xsec parameters
        {
	        h_err_prior -> Divide(h_prior);
	        h_err_final -> Divide(h_final);
	    }

        TCanvas *c_err = new TCanvas("c_err", "c_err", 1400, 900);
        // gStyle -> SetOptTitle(0);
         h_err_final -> SetTitle(Form(";%s;Parameter error", par_title[j].c_str() ));

        double Max = MAX(h_err_prior->GetMaximum(), h_err_final->GetMaximum());
        h_err_final -> GetYaxis() -> SetRangeUser(0.0, 1.3*Max );

        h_err_prior -> SetLineColor(kBlue+1);
        h_err_prior -> SetLineWidth(2);
        h_err_final -> SetLineColor(kRed+1);
        h_err_final -> SetLineWidth(2);

        h_err_final -> Draw("hist");
        h_err_prior -> Draw("hist same");

        // TLegend* legend = new TLegend(0.55,0.65,0.9,0.9);
        TLegend* legend2 = new TLegend(0.65,0.75,0.9,0.9);
        legend2 -> SetFillColor(0);
        // legend2 -> SetHeader(Form("%s", title.c_str()));
        legend2 -> AddEntry(h_err_prior, "Prior error","l");
        legend2 -> AddEntry(h_err_final, "Final error","l");
        legend2 -> Draw();

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
    ss << output_dir << "chi2_" << file_name << ".png";
    c -> Print(ss.str().c_str());
    // delete c;
    

    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::endl;

}
