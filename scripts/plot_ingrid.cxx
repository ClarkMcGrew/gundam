void plot_ingrid(const std::string fname)
{
    TFile* input_file = TFile::Open(fname.c_str(), "READ");
    TTree* ing_tree = (TTree*)input_file -> Get("wtree");

    int interaction_type, fsi_int;
    int muon_track, selected_sample;
    int track_sample;
    float enu, event_weight;
    float pmu_true, angle_true;
    float pmu_reco, angle_reco;
    bool is_anti, is_nue, is_fv;
    bool new_event;

    ing_tree -> SetBranchAddress("InteractionType", &interaction_type);
    ing_tree -> SetBranchAddress("FSIInt", &fsi_int);
    ing_tree -> SetBranchAddress("SelectedSample", &selected_sample);
    ing_tree -> SetBranchAddress("MuonCandidateTrack", &muon_track);
    ing_tree -> SetBranchAddress("Enu", &enu);
    ing_tree -> SetBranchAddress("IsAnti", &is_anti);
    ing_tree -> SetBranchAddress("IsNuE", &is_nue);
    ing_tree -> SetBranchAddress("IsFV", &is_fv);
    ing_tree -> SetBranchAddress("NewEvent", &new_event);
    ing_tree -> SetBranchAddress("weight", &event_weight);
    ing_tree -> SetBranchAddress("TrueMomentumMuon", &pmu_true);
    ing_tree -> SetBranchAddress("TrueAngleMuon", &angle_true);

    const float deg_to_rad = TMath::Pi() / 180;
    const float iron_carbon_ratio = 7.640777;
    int sample_track0, sample_track1, sample_track2;
    float iron_track0, plastic_track0, angle_track0;
    float iron_track1, plastic_track1, angle_track1;
    float iron_track2, plastic_track2, angle_track2;

    ing_tree -> SetBranchAddress("Sample_track0", &sample_track0);
    ing_tree -> SetBranchAddress("TrackAngle_track0", &angle_track0);
    ing_tree -> SetBranchAddress("IronDistance_track0", &iron_track0);
    ing_tree -> SetBranchAddress("PlasticDistance_track0", &plastic_track0);
    ing_tree -> SetBranchAddress("Sample_track1", &sample_track1);
    ing_tree -> SetBranchAddress("TrackAngle_track1", &angle_track1);
    ing_tree -> SetBranchAddress("IronDistance_track1", &iron_track1);
    ing_tree -> SetBranchAddress("PlasticDistance_track1", &plastic_track1);
    ing_tree -> SetBranchAddress("Sample_track2", &sample_track2);
    ing_tree -> SetBranchAddress("TrackAngle_track2", &angle_track2);
    ing_tree -> SetBranchAddress("IronDistance_track2", &iron_track2);
    ing_tree -> SetBranchAddress("PlasticDistance_track2", &plastic_track2);

    TH1D* h1d_angle = new TH1D("h1d_angle", "", 30, 0, 90);
    TH1D* h1d_dist = new TH1D("h1d_dist", "", 50, 0, 100);
    TH2D* h2d_angle = new TH2D("h2d_angle", "", 30, 0, 90, 30, 0, 90);
    TH2D* h2d_dist = new TH2D("h2d_dist", "", 50, 0, 100, 30, 0, 1.5);
    //TH2D* h2d_dist = new TH2D("h2d_dist", "", 15, 0, 1.5, 15, 0, 1.5);

    TH1D* h1d_angle_fsi0 = new TH1D("angle_fsi0", "", 30, 0, 90);
    h1d_angle_fsi0 -> SetFillColor(kOrange+10);
    TH1D* h1d_angle_fsi1 = new TH1D("angle_fsi1", "", 30, 0, 90);
    h1d_angle_fsi1 -> SetFillColor(kRed);
    TH1D* h1d_angle_fsi2 = new TH1D("angle_fsi2", "", 30, 0, 90);
    h1d_angle_fsi2 -> SetFillColor(kOrange+8);
    TH1D* h1d_angle_fsi3 = new TH1D("angle_fsi3", "", 30, 0, 90);
    h1d_angle_fsi3 -> SetFillColor(kAzure+10);
    TH1D* h1d_angle_fsi4 = new TH1D("angle_fsi4", "", 30, 0, 90);
    h1d_angle_fsi4 -> SetFillColor(kAzure+7);
    TH1D* h1d_angle_fsi5 = new TH1D("angle_fsi5", "", 30, 0, 90);
    h1d_angle_fsi5 -> SetFillColor(kBlue+2);
    TH1D* h1d_angle_fsi6 = new TH1D("angle_fsi6", "", 30, 0, 90);
    h1d_angle_fsi6 -> SetFillColor(kGray);
    TH1D* h1d_angle_fsi7 = new TH1D("angle_fsi7", "", 30, 0, 90);
    h1d_angle_fsi7 -> SetFillColor(kYellow);

    TH1D* h1d_dist_fsi0 = new TH1D("dist_fsi0", "", 20, 0, 100);
    h1d_dist_fsi0 -> SetFillColor(kOrange+10);
    TH1D* h1d_dist_fsi1 = new TH1D("dist_fsi1", "", 20, 0, 100);
    h1d_dist_fsi1 -> SetFillColor(kRed);
    TH1D* h1d_dist_fsi2 = new TH1D("dist_fsi2", "", 20, 0, 100);
    h1d_dist_fsi2 -> SetFillColor(kOrange+8);
    TH1D* h1d_dist_fsi3 = new TH1D("dist_fsi3", "", 20, 0, 100);
    h1d_dist_fsi3 -> SetFillColor(kAzure+10);
    TH1D* h1d_dist_fsi4 = new TH1D("dist_fsi4", "", 20, 0, 100);
    h1d_dist_fsi4 -> SetFillColor(kAzure+7);
    TH1D* h1d_dist_fsi5 = new TH1D("dist_fsi5", "", 20, 0, 100);
    h1d_dist_fsi5 -> SetFillColor(kBlue+2);
    TH1D* h1d_dist_fsi6 = new TH1D("dist_fsi6", "", 20, 0, 100);
    h1d_dist_fsi6 -> SetFillColor(kGray);
    TH1D* h1d_dist_fsi7 = new TH1D("dist_fsi7", "", 20, 0, 100);
    h1d_dist_fsi7 -> SetFillColor(kYellow);

    double sel = 0;
    double gen = 0;
    double sig = 0;

    const unsigned int nevents = ing_tree -> GetEntries();
    for(unsigned int i = 0; i < nevents; ++i)
    {
        ing_tree -> GetEntry(i);
        //if(pmu_true == 0 || angle_true == 0)
        //    continue;

        switch(muon_track)
        {
            case 0:
                pmu_reco = iron_track0 + (plastic_track0 / iron_carbon_ratio);
                angle_reco = angle_track0;
                track_sample = sample_track0;
                break;
            case 1:
                pmu_reco = iron_track1 + (plastic_track1 / iron_carbon_ratio);
                angle_reco = angle_track1;
                track_sample = sample_track1;
                break;
            case 2:
                pmu_reco = iron_track2 + (plastic_track2 / iron_carbon_ratio);
                angle_reco = angle_track2;
                track_sample = sample_track2;
                break;
            default:
                pmu_reco = -999;
                angle_reco = -999;
                track_sample = 999;
                break;
        }

        if(selected_sample == 0 && track_sample < 4)
            sel += event_weight;

        if(is_fv && !is_anti && !is_nue && fsi_int < 3 && new_event == 1)
            gen += event_weight;

        if(selected_sample == 0 && fsi_int < 3 && !is_anti && !is_nue && track_sample < 4)
            sig += event_weight;

        if(selected_sample == 0 && track_sample < 4)
        {
            h1d_angle->Fill(angle_reco, event_weight);
            h1d_dist->Fill(pmu_reco, event_weight);

            //pmu_reco = (pmu_reco * 0.0114127 + 0.230608);
            h2d_angle->Fill(angle_reco, angle_true, event_weight);
            h2d_dist->Fill(pmu_reco, pmu_true, event_weight);

            if(fsi_int == 0)
            {
                h1d_angle_fsi0->Fill(angle_reco, event_weight);
                h1d_dist_fsi0->Fill(pmu_reco, event_weight);
            }
            else if(fsi_int == 1)
            {
                h1d_angle_fsi1->Fill(angle_reco, event_weight);
                h1d_dist_fsi1->Fill(pmu_reco, event_weight);
            }
            else if(fsi_int == 2)
            {
                h1d_angle_fsi2->Fill(angle_reco, event_weight);
                h1d_dist_fsi2->Fill(pmu_reco, event_weight);
            }
            else if(fsi_int == 3)
            {
                h1d_angle_fsi3->Fill(angle_reco, event_weight);
                h1d_dist_fsi3->Fill(pmu_reco, event_weight);
            }
            else if(fsi_int == 4)
            {
                h1d_angle_fsi4->Fill(angle_reco, event_weight);
                h1d_dist_fsi4->Fill(pmu_reco, event_weight);
            }
            else if(fsi_int == 5)
            {
                h1d_angle_fsi5->Fill(angle_reco, event_weight);
                h1d_dist_fsi5->Fill(pmu_reco, event_weight);
            }
            else if(fsi_int == 6)
            {
                h1d_angle_fsi6->Fill(angle_reco, event_weight);
                h1d_dist_fsi6->Fill(pmu_reco, event_weight);
            }
            else
            {
                h1d_angle_fsi7->Fill(angle_reco, event_weight);
                h1d_dist_fsi7->Fill(pmu_reco, event_weight);
            }
        }
    }

    std::cout << "True Signal Events: " << sig << std::endl;
    std::cout << "Sel  Signal Events: " << sel << std::endl;
    std::cout << "Gen  Signal Events: " << gen << std::endl;
    std::cout << "Eff: " << sel / gen << std::endl;
    std::cout << "Pur: " << sig / sel << std::endl;

    TCanvas* c = new TCanvas("c","c",1200,900);
    c->Divide(2,2);
    c->cd(1);
    h1d_angle->Draw("hist");
    c->cd(2);
    h1d_dist->Draw("hist");
    c->cd(3);
    h2d_angle->Draw("colz");
    c->cd(4);
    h2d_dist->Draw("colz");

    TCanvas* d = new TCanvas("d","d",1200,900);
    THStack* hs_angle_fsi = new THStack("hs_angle_fsi", "");
    hs_angle_fsi->Add(h1d_angle_fsi0);
    hs_angle_fsi->Add(h1d_angle_fsi1);
    hs_angle_fsi->Add(h1d_angle_fsi2);
    hs_angle_fsi->Add(h1d_angle_fsi3);
    hs_angle_fsi->Add(h1d_angle_fsi4);
    hs_angle_fsi->Add(h1d_angle_fsi5);
    hs_angle_fsi->Add(h1d_angle_fsi6);
    hs_angle_fsi->Add(h1d_angle_fsi7);
    hs_angle_fsi->Draw("hist");
    d->Print("ingrid_angle.png");

    TCanvas* e = new TCanvas("e","e",1200,900);
    THStack* hs_dist_fsi = new THStack("hs_dist_fsi", "");
    hs_dist_fsi->Add(h1d_dist_fsi0);
    hs_dist_fsi->Add(h1d_dist_fsi1);
    hs_dist_fsi->Add(h1d_dist_fsi2);
    hs_dist_fsi->Add(h1d_dist_fsi3);
    hs_dist_fsi->Add(h1d_dist_fsi4);
    hs_dist_fsi->Add(h1d_dist_fsi5);
    hs_dist_fsi->Add(h1d_dist_fsi6);
    hs_dist_fsi->Add(h1d_dist_fsi7);
    hs_dist_fsi->Draw("hist");
    e->Print("ingrid_distance.png");
}

