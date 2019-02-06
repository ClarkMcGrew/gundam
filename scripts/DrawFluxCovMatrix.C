#include "CommonHeader.h"
#include "CommonStyle.h"

void DrawFluxCovMatrix(){

  //======================================================================================================  
  //== Set T2K style. In CommonStyle.h option 1 has been setted
  CommonStyle();

  //======================================================================================================
  
  TFile *finfluxcov = TFile::Open("/sps/t2k/lmaret/softwares/xsLLhFitterLM/inputs/fgd1fgd2Fit/flux_covariance_banff_13av2.root"); //contains flux and det. systematics info
   
  //==== setup enu bins and covm for flux 
  TMatrixDSym *cov_flux_in = (TMatrixDSym*)finfluxcov->Get("total_flux_cov");
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  c1->Draw();
  finfluxcov->Get("total_flux_cov")->Draw("colz");
 
  TAxis *nd_numu_bins       = (TAxis*)finfluxcov->Get("nd5_numode_numu_bins");//numu in nu beam
  TAxis *nd_antinumuws_bins = (TAxis*)finfluxcov->Get("nd5_numode_numub_bins");//anti-numu in nu beam
  TAxis *nd_nue_bins        = (TAxis*)finfluxcov->Get("nd5_numode_nue_bins");//nue in nu beam
  TAxis *nd_antinue_bins    = (TAxis*)finfluxcov->Get("nd5_numode_nueb_bins");//anti-nue in nu beam
  TAxis *nd_antinumu_bins   = (TAxis*)finfluxcov->Get("nd5_anumode_numub_bins");//anti-numu in anti-nu beam
  TAxis *nd_numuws_bins     = (TAxis*)finfluxcov->Get("nd5_anumode_numu_bins");//numu in anti-nu beam
 
  vector<double> enubins;
  enubins.push_back(nd_numu_bins->GetBinLowEdge(1));
  for(int i=0;i<nd_numu_bins->GetNbins();i++){
    enubins.push_back(nd_numu_bins->GetBinUpEdge(i+1));}

  vector<double> eantinuwsbins;
  eantinuwsbins.push_back(nd_antinumuws_bins->GetBinLowEdge(1));  
  for(int i=0;i<nd_antinumuws_bins->GetNbins();i++){
    eantinuwsbins.push_back(nd_antinumuws_bins->GetBinUpEdge(i+1));}
   
  vector<double> enuwsbins;
  enuwsbins.push_back(nd_numuws_bins->GetBinLowEdge(1));
  for(int i=0;i<nd_numuws_bins->GetNbins();i++){
    enuwsbins.push_back(nd_numuws_bins->GetBinUpEdge(i+1));}

  vector<double> eantinubins;
  eantinubins.push_back(nd_antinumu_bins->GetBinLowEdge(1));  
  for(int i=0;i<nd_antinumu_bins->GetNbins();i++){
    eantinubins.push_back(nd_antinumu_bins->GetBinUpEdge(i+1));}

  const int dimfluxcov=nd_numu_bins->GetNbins()+nd_antinumuws_bins->GetNbins()+nd_antinumu_bins->GetNbins()+nd_numuws_bins->GetNbins();
  const int dimfluxcov_nu=nd_numu_bins->GetNbins()+nd_antinumuws_bins->GetNbins();
  const int nuenueb = nd_nue_bins->GetNbins() + nd_antinue_bins->GetNbins();
  TMatrixDSym cov_flux(dimfluxcov);
  for(int i=0;i<dimfluxcov;i++){
    for(int j=0;j<dimfluxcov;j++){
      if (i<=dimfluxcov_nu && j<=dimfluxcov_nu){
	cov_flux(i, j) = (*cov_flux_in)(i,j);
      }
      else if (i>=dimfluxcov_nu && j<=dimfluxcov_nu) {
        cov_flux(i, j) = (*cov_flux_in)(i+9,j);
      }
      else if (i<=dimfluxcov_nu && j>=dimfluxcov_nu) {
	cov_flux(i, j) = (*cov_flux_in)(i,j+9);
      }
      else if (i>=dimfluxcov_nu && j>=dimfluxcov_nu) {
        cov_flux(i, j) = (*cov_flux_in)(i+9,j+9);
      }
    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.25);
  gStyle->SetPadLeftMargin(0.25);

  TCanvas *c2=new TCanvas("c2","c2",1200,900);
  c2->Draw();
  // TFile cov("final_fluxcov","recreate");
  // cov_flux.Write("finalfluxcov");
  //TH2D* hcov = new TH2D("","",dimfluxcov,0,dimfluxcov,dimfluxcov,0,dimfluxcov); 
  //hcov = (TH2D*)cov_flux;
  //hcov->Draw("colz");
  cov_flux.Draw("colz");
  TLine* olbin[3];
  olbin[0] = new TLine(11,0,11,32);
  olbin[1] = new TLine(16,0,16,32);
  olbin[2] = new TLine(21,0,21,32);
  for(int l=0;l<3;l++){
    olbin[l]->SetLineWidth(5);
    olbin[l]->Draw();
  }
  TLine* vlbin[3];
  vlbin[0] = new TLine(0,11,32,11);
  vlbin[1] = new TLine(0,16,32,16);
  vlbin[2] = new TLine(0,21,32,21);
  for(int l=0;l<3;l++){
    vlbin[l]->SetLineWidth(5);
    vlbin[l]->Draw();
  }
  c2->Print("plots/fluxCovMatrix.pdf");
  
  TMatrixDSym corr_flux(dimfluxcov);
  for(int i=0;i<dimfluxcov;i++){
    for(int j=0;j<dimfluxcov;j++){
      corr_flux(i,j) = cov_flux(i,j)/(TMath::Sqrt(cov_flux(i,i))*TMath::Sqrt(cov_flux(j,j)));
    }
  }

  TCanvas *c3=new TCanvas("c3","c3",1200,900);
  c3->Draw();
  corr_flux.Draw("colz");
  
  for(int l=0;l<3;l++){
    olbin[l]->SetLineWidth(5);
    olbin[l]->Draw();
  }
  for(int l=0;l<3;l++){
    vlbin[l]->SetLineWidth(5);
    vlbin[l]->Draw();
  }

  c3->Print("plots/covariancematrices/fluxCorrMatrix.pdf");

  //cov_flux.Close(); 
  finfluxcov->Close(); 
}
