#ifndef __BinningTools_hh__
#define __BinningTools_hh__

#include "CommonHeader.h"

using namespace std;

class BinningTools{
  
 public:
  BinningTools();
  ~BinningTools();
  
  void SetBinning(const string fname);

  vector<pair<double, double> > Getpair_costh();
  vector<pair<double, double> > Getpair_mom();

  int GetNbins() { return m_nbins; }
  int GetNbins_costh() { return m_nbins_costh; }
  int GetNbins_mom(int cthidx) { return m_nbins_mom[cthidx]; }
  int* GetMomNBinsInEachCosthSlice();

  double* GetCosthBins();		
  double** GetMomBins();		

  double* GetMomBinWidth();
  double* GetCosBinWidth();

 private:

  vector<pair<double, double> > m_momedges;
  vector<pair<double, double> > m_costhedges;
  vector< pair<double, double> > m_costhslices;
  
  int m_nbins;
  int m_nbins_costh;
  vector <int> m_nbins_mom;

  double *m_bins_costh;
  vector < vector<double> > m_bins_mom_vec;
  double **m_bins_mom_arr;

};

#endif
