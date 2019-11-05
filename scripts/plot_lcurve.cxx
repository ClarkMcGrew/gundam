void plot_lcurve(const std::string& filename, const std::string& saveas = "lcurve")
{
    std::cout << "Reading " << filename << " for data." << std::endl;

    std::vector<double> v_x;
    std::vector<double> v_y;
    std::vector<double> v_r;
    std::ifstream fin(filename, std::ios::in);
    if(!fin.is_open())
    {
        std::cout << "[ERROR]: Failed to open " << filename << std::endl;
        return -1;
    }
    else
    {
        std::string line;
        double r, x, y;
        while(std::getline(fin, line))
        {
            std::stringstream s(line);
            s >> r; s >> x; s >> y;

            std::cout << "Adding (" << x << "," << y << ") to graph." << std::endl;
            v_r.push_back(r);
            v_x.push_back(x);
            v_y.push_back(y);
        }
    }

    if(v_x.size() != v_y.size())
    {
        std::cout << "X and Y vector sizes do not match! Exiting!" << std::endl;
        return;
    }

    const unsigned int npoints = v_x.size();
    TCanvas* c = new TCanvas("c", "c", 1200, 900);
    TGraph* g = new TGraph(npoints);
    TText* t = nullptr;

    for(unsigned int i = 0; i < npoints; ++i)
    {
        g->SetPoint(i, v_x[i], v_y[i]);
    }

    g->SetTitle("");
    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerColor(kBlue);
    g->SetMarkerSize(1);
    g->GetXaxis()->SetDecimals(true);
    g->GetYaxis()->SetDecimals(true);
    g->GetXaxis()->SetTitle("#chi^{2}_{total}");
    g->GetYaxis()->SetTitle("#chi^{2}_{reg}/#lambda");
    g->Draw("AP");

    for(unsigned int i = 0; i < npoints; ++i)
    {
        // For some strange reason, the iostream format options are not
        // working. Instead fall back and use sprintf, which does work.
        // std::stringstream ss;
        // ss << std::setprecision(2) << std::fixed << std::to_string(v_r[i]);

        char buf[16];
        sprintf(buf, "%2.2f", v_r[i]);
        t = new TText(v_x[i]+0.1, v_y[i]+0.05, buf);
        t->SetTextSize(0.025);
        t->Draw("same");
    }

    std::string save_str;
    save_str = saveas + ".pdf";
    c->Print(save_str.c_str());

    return;
}
