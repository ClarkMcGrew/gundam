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
    if(file_name == "fit3_statFluc" || file_name == "fit3_statFluc_reg" ) drawParThrow = true;

    const std::string output_dir = Form("plots/fitteroutput/%s/", dir_name.c_str());
    const std::string file_name_input = Form("../outputs/%s.root", file_name.c_str());

    std::string name;
    std::string title;
    std::stringstream ss;

    TFile* file = TFile::Open(file_name_input.c_str(), "READ");

    const int Npar = 4;
    std::string par_name[Npar]  = {"par_fit", "par_flux", "par_xsec", "par_det"};
    std::string par_title[Npar] = {"Template parameters", "Flux parameters", "Xsec model parameters", "Detector parameters"};

    std::string xsec_par_name[26] = {"M_{A}^{QE}", "p_{F}^{C}", "MEC C", "E_{B}^{C}", "p_{F}^{O}", "MEC O", "E_{B}^{O}", "MEC O shape", "MEC C shape", 
                                        "C_{A}^{5}", "M_{A}^{RES}", "Bkg Res", "CC-#nu_{e}", "CC-DIS shape", "CC coh", "NC-Other", "CC-1#pi norm. low E", "CC-1#pi norm. high E", "CC-multi-#pi norm.", "CC-DIS norm.",
                                        "FSI inel. low E", "FSI #pi abs.", "FSI charge exch. low E", "FSI #pi prod.", "FSI inel. high E", "FSI charge exch. high E" };

    const int lowlim = -2.0;
    const int highlim = 4.0;

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
    // for(int j = 2; j<3; j++)
    for(int j = 0; j<Npar; j++)
    {
        name = par_name[j];
        title = par_title[j];

        std::cout << "----- for " << title << " (" << name << ") -----" << std::endl;

        if( file_name == "fit3_statFluc_reg_onlyXsec" && (j==1 || j==3) )
            continue;

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


        if(name == "par_xsec")
        {
            for(int i = 1; i <= h_prior->GetNbinsX(); ++i)
            {
                h_final -> SetBinError(i, h_err_final -> GetBinContent(i) / h_prior->GetBinContent(i));
                h_prior -> SetBinError(i, h_err_prior -> GetBinContent(i) / h_prior->GetBinContent(i));
            }

            h_err_prior -> Divide(h_prior);
            h_err_final -> Divide(h_final);

            for(int i = 1; i <= h_prior->GetNbinsX(); ++i)
            {
                h_final                  -> SetBinContent(i, h_final -> GetBinContent(i) / h_prior->GetBinContent(i));
                if(drawParThrow) h_iter0 -> SetBinContent(i, h_iter0 -> GetBinContent(i) / h_prior->GetBinContent(i));
                h_prior                  -> SetBinContent(i, h_prior -> GetBinContent(i) / h_prior->GetBinContent(i));
            }
        }




        // Define bin labels on the parameter axis

        int division = 10;
        for(int i = 0; i < h_prior->GetNbinsX(); ++i)
        {
            if(name == "par_det")  division = 100;
            if(name == "par_xsec") division = 1;

            if(i%division == 0)
            {
                if(name == "par_xsec")
                {
                    h_err_final -> GetXaxis() -> SetBinLabel(i+1, xsec_par_name[i].c_str());
                    h_prior     -> GetXaxis() -> SetBinLabel(i+1, xsec_par_name[i].c_str());
                }
                else
                {
                    h_err_final -> GetXaxis() -> SetBinLabel(i+1, Form("%d", i));
                    h_prior     -> GetXaxis() -> SetBinLabel(i+1, Form("%d", i));
                }
            }
            else if(i%division != 0)
            {
                h_err_final -> GetXaxis() -> SetBinLabel(i+1, "");
                h_prior     -> GetXaxis() -> SetBinLabel(i+1, "");
            }
        }

        
        int fontsize = 65;

        // X axis title
        h_err_final -> GetXaxis() -> SetLabelFont(43);
        if(name == "par_xsec") h_err_final -> GetXaxis() -> SetLabelSize(25);
        else                   h_err_final -> GetXaxis() -> SetLabelSize(fontsize-10);
        h_err_final -> GetXaxis() -> SetTitleFont(43);
        h_err_final -> GetXaxis() -> SetTitleSize(fontsize);

        if(name == "par_xsec") h_err_final -> GetXaxis() -> SetLabelOffset(0.02);
        h_err_final -> GetXaxis() -> SetTitleOffset(4.3);

        // Y axis title and labels for the uncertainty
        h_err_final -> GetYaxis() -> SetLabelFont(43);
        h_err_final -> GetYaxis() -> SetLabelSize(fontsize-10);
        h_err_final -> GetYaxis() -> SetTitleFont(43);
        h_err_final -> GetYaxis() -> SetTitleSize(fontsize);

        h_err_final -> GetYaxis() -> SetNdivisions(7);

        h_err_final -> GetYaxis() -> SetTitleOffset(1.5);

        // Y axis title and labels for the parameters
        h_prior -> GetYaxis() -> SetLabelFont(43);
        h_prior -> GetYaxis() -> SetLabelSize(fontsize-10);
        h_prior -> GetYaxis() -> SetTitleFont(43);
        h_prior -> GetYaxis() -> SetTitleSize(fontsize);

        h_prior -> GetYaxis() -> SetTitleOffset(1.5);

        


        // std::cout << "----- Draw parameters and errors -----" << std::endl;

        TCanvas *c = new TCanvas("c", "c", 1800, 1400);
        // gStyle -> SetOptTitle(0);
   
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        // pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        

        h_prior -> SetTitle(Form(";%s;Parameter value", par_title[j].c_str() ));

        h_prior -> SetMarkerColor(kRed+2);
        h_prior -> SetMarkerStyle(kFullSquare);
        if(drawParThrow) h_iter0 -> SetMarkerColor(kGreen+3);
        if(drawParThrow) h_iter0 -> SetMarkerStyle(kFullTriangleUp);
        h_final -> SetMarkerColor(kBlue+2);
        h_final -> SetMarkerStyle(kFullCircle);
        
        // Define error bar style
        h_final -> SetLineColor(kBlack);
        h_final -> SetBarWidth(0.8);

        if(name == "par_det")
        {
            h_prior -> SetMarkerSize(1);
            if(drawParThrow) h_iter0 -> SetMarkerSize(1);
            h_final -> SetMarkerSize(1);
            h_final -> SetLineWidth(1);
        }
        else
        {
            h_prior -> SetMarkerSize(2);
            if(drawParThrow) h_iter0 -> SetMarkerSize(2);
            h_final -> SetMarkerSize(2);
            h_final -> SetLineWidth(2);
        }

        h_prior -> SetFillColor(kRed-9);
        // h_final -> SetFillColor(kBlue-9);
        // h_prior -> SetFillStyle(1001);
        h_prior -> SetFillStyle(3144);
        // h_final -> SetFillStyle(3144);

        if(name == "par_xsec") // for xsec parameters
        {
            // gPad->SetLogy();
            // double max_xsec = h_final->GetMaximum();
            // double min_xsec = h_final->GetMinimum();
            // h_prior -> GetYaxis() -> SetRangeUser(0.6*min_xsec, 1.4*max_xsec);
            h_prior -> GetYaxis() -> SetRangeUser(0.0, 3.0);
        }
        else if(name == "par_flux") // for flux parameters
            h_prior -> GetYaxis() -> SetRangeUser(0.7,1.3);
        else if(name == "par_fit")
            h_prior -> GetYaxis() -> SetRangeUser(lowlim, highlim);
        else if(name == "par_det") // for flux parameters
            h_prior -> GetYaxis() -> SetRangeUser(0.3,1.7);


        h_prior -> Draw("P E2");


        // Draw vertical lines for fit parameters

        TLine *verline[2*Nbins_costh];
        int binSide = 0;
        for(int il = 0; il < Nbins_costh; il++)
        {
            binSide += Nmombins[il];
            verline[il]               = new TLine(binSide,         lowlim, binSide,         highlim );
            verline[il + Nbins_costh] = new TLine(binSide + Nbins, lowlim, binSide + Nbins, highlim );
        }

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
        TLine *verlineTargets = new TLine(Nbins, lowlim, Nbins, highlim );
        TLine *horlineTargets = new TLine(0.0, 0.0, 2*Nbins, 0.0 );
        if(name == "par_fit")
        {
            verlineTargets -> SetLineColor(kGray+1);
            horlineTargets -> SetLineColor(kGray+1);
            verlineTargets -> Draw();
            horlineTargets -> Draw();
        }



        // Draw vertical lines for detector parameters

        int NbinsDetTot = h_prior -> GetNbinsX();
        int NbinsDet[] = {58, 58, 58, 10, 16, 58, 58, 58,
                          58, 58, 58, 10, 16, 58, 58, 58,
                          58, 58, 58, 10, 16, 58, 58, 58 };
        int Nsamples = 3*7 -1;
        int Nfgds = 3;
    
        double boundary[Nsamples];
        boundary[0] = NbinsDet[0];
        for(int il=1; il<Nsamples; il++)
        {
            boundary[il] = boundary[il-1] + NbinsDet[il];
        }                    

        TLine *verlineDet[Nsamples];
        if(name == "par_det")
        {
            for(int il = 0; il < Nsamples; il++)
            {
                verlineDet[il] = new TLine(boundary[il], 0.3, boundary[il], 1.7 );
                verlineDet[il] -> SetLineWidth(1);
                verlineDet[il] -> SetLineStyle(1);
                verlineDet[il] -> SetLineColor(kGray);
                verlineDet[il] -> Draw();
            }
        }
        TLine *verlineFGD1 = new TLine(NbinsDetTot/3,   0.3, NbinsDetTot/3,   1.7);
        TLine *verlineFGD2 = new TLine(NbinsDetTot*2/3, 0.3, NbinsDetTot*2/3, 1.7);
        if(name == "par_det")
        {
            verlineFGD1 -> SetLineColor(kGray+1);
            verlineFGD2 -> SetLineColor(kGray+1);
            verlineFGD1 -> Draw();
            verlineFGD2 -> Draw();
        }



        // Draw vertical lines for xsec parameters

        TLine *verlineXsec1 = new TLine(9.0,  0.0, 9.0,  3.0 );
        TLine *verlineXsec2 = new TLine(20.0, 0.0, 20.0, 3.0 );
        if(name == "par_xsec")
        {
            verlineXsec1 -> SetLineColor(kGray+1);
            verlineXsec2 -> SetLineColor(kGray+1);
            verlineXsec1 -> Draw();
            verlineXsec2 -> Draw();
        }



        // Draw parameters

        if(drawParThrow) h_iter0 -> Draw("P same");
        h_final -> Draw("P E1 same"); // option X0 to remove error bar on X axis

        // h_prior -> GetYaxis() -> SetLabelSize(0.);

        // TLegend* legend = new TLegend(0.55,0.67,0.9,0.9);
        TLegend* legend = new TLegend(0.75,0.75,0.9,0.9);
        legend -> SetFillColor(0);
        // legend -> SetHeader(Form("%s", title.c_str()));
        legend -> AddEntry(h_prior, "Prior","p");
        if(drawParThrow) legend -> AddEntry(h_iter0, "Initial","p");
        legend -> AddEntry(h_final, "Final","p");
        legend -> Draw();

        // ss.str("");
        // ss << output_dir << name << "_overlay_" << file_name << ".pdf";
        // c -> Print(ss.str().c_str());
        // delete c;










        // std::cout << "----- Draw errors only -----" << std::endl;

        // Divide by parameter value to have the relative error
     //    if(name == "par_xsec") // for xsec parameters
     //    {
	    //     h_err_prior -> Divide(h_prior);
	    //     h_err_final -> Divide(h_final);
	    // }

        // TCanvas *c_err = new TCanvas("c_err", "c_err", 1400, 900);
   
        // lower plot will be in pad
        c->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.03, 1, 0.4);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.4);

        pad2->SetBorderMode(0);
        pad2->SetBorderSize(2);
        pad2->SetFrameBorderMode(0);

        // pad2->SetGridx(); // vertical grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad




        gStyle->SetPadBottomMargin(0.1);
        gROOT->ForceStyle();

        // gStyle -> SetOptTitle(0);
        h_err_final -> SetTitle(Form(";%s;Relative error", par_title[j].c_str() ));

        double Max = MAX(h_err_prior->GetMaximum(), h_err_final->GetMaximum());
        h_err_final -> GetYaxis() -> SetRangeUser(0.0, 1.3*Max );

        h_err_prior -> SetLineColor(kRed+2);
        h_err_prior -> SetLineWidth(2);
        h_err_final -> SetLineColor(kBlue+2);
        h_err_final -> SetLineWidth(2);

        TLine *verlineErr[2*Nbins_costh];
        binSide = 0;
        for(int il = 0; il < Nbins_costh; il++)
        {
            binSide += Nmombins[il];
            verlineErr[il]               = new TLine(binSide,         0, binSide,         1.3*Max );
            verlineErr[il + Nbins_costh] = new TLine(binSide + Nbins, 0, binSide + Nbins, 1.3*Max );
        }

        h_err_final -> Draw("hist");

        if(name == "par_fit")
        {
            for(int il = 0; il < 2*Nbins_costh-1; il++)
            {
                verlineErr[il] -> SetLineWidth(1);
                verlineErr[il] -> SetLineStyle(1);
                verlineErr[il] -> SetLineColor(kGray);
                verlineErr[il] -> Draw();
            }
        }
        TLine *verlineErrTargets = new TLine(Nbins, 0, Nbins, 1.3*Max );
        TLine *horlineErrTargets = new TLine(0.0, 0.0, 2*Nbins, 0.0 );
        if(name == "par_fit")
        {
            verlineErrTargets -> SetLineColor(kGray+1);
            horlineErrTargets -> SetLineColor(kGray+1);
            verlineErrTargets -> Draw();
            horlineErrTargets -> Draw();
        }

        TLine *verlineDetErr[Nsamples];
        if(name == "par_det")
        {
            for(int il = 0; il < Nsamples; il++)
            {
                verlineDetErr[il] = new TLine(boundary[il], 0.0, boundary[il], 1.3*Max );
                verlineDetErr[il] -> SetLineWidth(1);
                verlineDetErr[il] -> SetLineStyle(1);
                verlineDetErr[il] -> SetLineColor(kGray);
                verlineDetErr[il] -> Draw();
            }
        }
        TLine *verlineFGD1Err = new TLine(NbinsDetTot/3,   0.0, NbinsDetTot/3,   1.3*Max);
        TLine *verlineFGD2Err = new TLine(NbinsDetTot*2/3, 0.0, NbinsDetTot*2/3, 1.3*Max);
        if(name == "par_det")
        {
            verlineFGD1Err -> SetLineColor(kGray+1);
            verlineFGD2Err -> SetLineColor(kGray+1);
            verlineFGD1Err -> Draw();
            verlineFGD2Err -> Draw();
        }


        // Draw vertical lines for xsec parameters

        TLine *verlineXsec1Err = new TLine(9.0,  0.0, 9.0,  1.3*Max );
        TLine *verlineXsec2Err = new TLine(20.0, 0.0, 20.0, 1.3*Max );
        if(name == "par_xsec")
        {
            verlineXsec1Err -> SetLineColor(kGray+1);
            verlineXsec2Err -> SetLineColor(kGray+1);
            verlineXsec1Err -> Draw();
            verlineXsec2Err -> Draw();
        }



        h_err_final -> Draw("hist same");
        h_err_prior -> Draw("hist same");

        // // TLegend* legend = new TLegend(0.55,0.65,0.9,0.9);
        // TLegend* legend2 = new TLegend(0.65,0.75,0.9,0.9);
        // legend2 -> SetFillColor(0);
        // // legend2 -> SetHeader(Form("%s", title.c_str()));
        // legend2 -> AddEntry(h_err_prior, "Prior error","l");
        // legend2 -> AddEntry(h_err_final, "Final error","l");
        // legend2 -> Draw();

        // ss.str("");
        // ss << output_dir << name << "_errors_" << file_name << ".pdf";
        // c_err -> Print(ss.str().c_str());
        // delete c_err;

        ss.str("");
        ss << output_dir << name << "_overlay_" << file_name << ".pdf";
        c -> Print(ss.str().c_str());
        delete c;

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

    h_chi2_stat -> SetLineColor(kBlue+2);
    h_chi2_sys  -> SetLineColor(kRed+2);
    h_chi2_reg  -> SetLineColor(kGreen+2);
    h_chi2_tot  -> SetLineColor(kBlack);
    
    h_chi2_tot  -> SetLineWidth(2);
    h_chi2_stat -> SetLineWidth(2);
    h_chi2_sys  -> SetLineWidth(2);
    h_chi2_reg  -> SetLineWidth(2);

    double max = h_chi2_stat->GetMaximum();
    double min = h_chi2_sys->GetMinimum();
    std::cout << "min = " << min << ", max = " << max << std::endl;

    h_chi2_stat -> GetYaxis() -> SetRangeUser(1.1*min, 1.1*max);

    h_chi2_stat -> GetXaxis() -> SetNdivisions(5);

    h_chi2_stat -> Draw("hist");
    h_chi2_sys  -> Draw("hist same");
    h_chi2_reg  -> Draw("hist same");
    h_chi2_tot  -> Draw("hist same");

    TLegend* legend = new TLegend(0.55,0.75,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> AddEntry(h_chi2_stat, "chi2 stat","l");
    legend -> AddEntry(h_chi2_sys , "chi2 sys","l");
    legend -> AddEntry(h_chi2_reg , "chi2 reg","l");
    legend -> AddEntry(h_chi2_tot , "chi2 total","l");
    legend -> Draw();

    ss.str("");
    ss << output_dir << "chi2_" << file_name << ".png";
    c -> Print(ss.str().c_str());
    delete c;
    

    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::endl;

}
