///////////////////////////////////////////////////////////////////////////////////////
// 
//
//
//  Created: 16 November 2016 by Ciro Riccio
//
//  Modified: 19 February 2019 by Lucie Maret
//
//
//  Usage: root -b -q 'DrawSelectionEfficiency.C+()'
//         
//
//
///////////////////////////////////////////////////////////////////////////////////////

#include "CommonHeader.h"
#include "CommonStyle.h"

TString infilename = "/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/xsllh_nominal.root";
void DrawEffBySample(Int_t sampleID, TString sampleName);



void DrawSelectionEfficiency()
{

  const Int_t Nsample = 24;
  Int_t sampleID_tab[Nsample]   = { 0,  1,  2,  3,  4,  5,  6,  7,
                                    8,  9,  10, 11, 12, 13, 14, 15,
                                    16, 17, 18, 19, 20, 21, 22, 23 };

  TString samplename_tab[Nsample] = { "FGD1MuTPC",  "FGD1MuTPCpTPC",  "FGD1MuTPCpFGD",  "FGD1MuFGDpTPC",  "FGD1MuFGD",  "FGD1CC1pi",  "FGD1CCDIS",  "FGD1CC1piMichelEle",
                                      "FGD2xMuTPC", "FGD2xMuTPCpTPC", "FGD2xMuTPCpFGD", "FGD2xMuFGDpTPC", "FGD2xMuFGD", "FGD2xCC1pi", "FGD2xCCDIS", "FGD2xCC1piMichelEle",
                                      "FGD2yMuTPC", "FGD2yMuTPCpTPC", "FGD2yMuTPCpFGD", "FGD2yMuFGDpTPC", "FGD2yMuFGD", "FGD2yCC1pi", "FGD2yCCDIS", "FGD2yCC1piMichelEle" };

  for(int i=0; i<Nsample; i++)
  {
    DrawEffBySample(sampleID_tab[i], samplename_tab[i]);
  }


}


void DrawEffBySample(Int_t sampleID, TString sampleName)
{

  //======================================================================================================  
  //=== Set T2K style. In CommonStyle.h option 1 has been setted
  CommonStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.25);
  gStyle->SetOptTitle(0);
  //======================================================================================================

  //======================================================================================================
  //======= Read Tree
  
  //=== Define the variable of the tree
  Float_t Pmutrue, CTHmutrue, reco_weight;
  Int_t    reco_reaction, topology, cut_branch;

  TFile *infile = new TFile(infilename.Data(), "OPEN");

  //=== Open tree
  TTree *recotree    = (TTree*) infile->Get("selectedEvents");

  //=== SetBranchAddress
  recotree->SetBranchAddress( "D1True"    , &Pmutrue      );
  recotree->SetBranchAddress( "D2True"    , &CTHmutrue    );
  // recotree->SetBranchAddress( "reaction"  , &reco_reaction);
  recotree->SetBranchAddress( "topology"  , &topology     );  
  recotree->SetBranchAddress( "cut_branch", &cut_branch   );  
  recotree->SetBranchAddress( "weight"    , &reco_weight  );  


  //=== Define the variable of the tree
  Float_t mumom_true, mucostheta_true, true_weight;
  Int_t    true_topology;

  //=== Open tree
  TTree *truthtree  = (TTree*) infile->Get("trueEvents");

  //=== SetBranchAddress
  truthtree->SetBranchAddress( "D1True"   , &mumom_true      );
  truthtree->SetBranchAddress( "D2True"   , &mucostheta_true );
  truthtree->SetBranchAddress( "topology" , &true_topology     );  
  // truthtree->SetBranchAddress( "reaction" , &true_reaction   );
  truthtree->SetBranchAddress( "weight"   , &true_weight     );  

  //======================================================================================================

  //======================================================================================================
  //======= Set binning

  //=== Muon momentum finner binning
  Int_t    mom_nbin   = 50;
  Double_t mom_bins[mom_nbin+1]; 
  for(int i = 0; i<mom_nbin+1; i++){
    mom_bins[i] = 0.1*i;
  }
  
  //=== Muon costheta binning
  // Linear fine binning
  // Double_t costheta_bins[] = { 0 , 0.05 , 0.1 , 0.15 , 0.2 , 0.25 , 0.3 , 0.35 , 0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.95 , 1 };
  // Non-linear analysis binning for costh
  Double_t costheta_bins[] = { -1.0, 0.2, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 0.98, 1.0 };
  Int_t    costheta_nbin   = sizeof(costheta_bins)/sizeof(Double_t) - 1;

  //======================================================================================================

  //======================================================================================================

  TH2D *hGen = new TH2D("Generated","Generated", mom_nbin, mom_bins, costheta_nbin, costheta_bins);
  TH2D *hSel = new TH2D("Selected", "Selected",  mom_nbin, mom_bins, costheta_nbin, costheta_bins);
  TH2D *hAllSel = new TH2D("All Selected", "ALl Selected",  mom_nbin, mom_bins, costheta_nbin, costheta_bins);

  TH2D *hEff = new TH2D("Efficiency","Efficiency", mom_nbin, mom_bins, costheta_nbin, costheta_bins);

  Long64_t true_nentries = truthtree->GetEntries();//Number of entries in seltree

  for (Int_t i=0; i<true_nentries; i++)
  {
    truthtree->GetEntry(i);
    if(true_topology==0) 
        hGen -> Fill(mumom_true/1000, mucostheta_true, true_weight);
  }

  Long64_t reco_nentries = recotree->GetEntries();//Number of entries in seltree

  for (Int_t i=0; i<reco_nentries; i++){
    recotree->GetEntry(i);
      
      if(sampleID==cut_branch)
      {
        if(topology==0) hSel    -> Fill(Pmutrue/1000, CTHmutrue, reco_weight);
                        
                        hAllSel -> Fill(Pmutrue/1000, CTHmutrue, reco_weight);
      }
  }

  hEff -> Divide(hSel,hGen);
  
  TLine *orline[14];
  orline[0] = new TLine(0.3,0.2,0.3,0.9);
  orline[1] = new TLine(0.4,0.2,0.4,0.98);
  orline[2] = new TLine(0.5,0.2,0.5,1.);
  orline[3] = new TLine(0.6,0.2,0.6,0.98);
  orline[4] = new TLine(0.7,0.98,0.7,1.);
  orline[5] = new TLine(0.8,0.6,0.8,0.98);
  orline[6] = new TLine(0.9,0.98,0.9,1.);
  orline[7] = new TLine(1.,0.8,1.,0.9);
  orline[8] = new TLine(1.25,0.9,1.25,1.);
  orline[9] = new TLine(1.5,0.85,1.5,0.9);
  orline[10] = new TLine(1.5,0.94,1.5,0.98);
  orline[11] = new TLine(2.,0.9,2.,1.);
  orline[12] = new TLine(3.,0.94,3.,1.);
  orline[13] = new TLine(4.,0.98,4.,1.);
  
  for(int il=0; il<14; il++)
  {
    orline[il]->SetLineWidth(2);
    orline[il]->SetLineColor(kRed+2);
  }
  
  //=== Draw a vertical line for each angular boundaries                                                                     
  TLine *verline[8];
  verline[0] = new TLine(0,0.2,5,0.2);
  verline[1] = new TLine(0,0.6,5,0.6);
  verline[2] = new TLine(0,0.7,5,0.7);
  verline[3] = new TLine(0,0.8,5,0.8);
  verline[4] = new TLine(0,0.85,5,0.85);
  verline[5] = new TLine(0,0.9,5,0.9);
  verline[6] = new TLine(0,0.94,5,0.94);
  verline[7] = new TLine(0,0.98,5,0.98);
  
  for(int il=0; il<8; il++)
  {
    verline[il]->SetLineWidth(2);
    verline[il]->SetLineColor(kRed+2);
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Draw();
  hEff->SetMaximum(1.);
  hEff->GetYaxis()->SetTitle("True Muon cos#theta");
  hEff->GetXaxis()->SetTitle("True Muon Momentum [GeV/c]");
  hEff->Draw("colz");
  for(int il=0; il<8; il++) verline[il]->Draw();
  for(int il=0; il<14; il++) orline[il]->Draw();

  c1 -> Print( Form("plots/2Dplots/Eff2D_%s.pdf", sampleName.Data()) );


  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  c2->Draw();
  // hAllSel->SetMaximum(1.);
  hAllSel->GetYaxis()->SetTitle("True Muon cos#theta");
  hAllSel->GetXaxis()->SetTitle("True Muon Momentum [GeV/c]");
  hAllSel->Draw("colz");
  for(int il=0; il<8; il++) verline[il]->Draw();
  for(int il=0; il<14; il++) orline[il]->Draw();

  c2 -> Print( Form("plots/2Dplots/sampleDistr2D_%s.pdf", sampleName.Data()) );

} 
 
