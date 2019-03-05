#include "BinningTools.hh"

//====================================================================================================== 
BinningTools::BinningTools(){
//====================================================================================================== 

  m_nbins_costh = 0;
  m_nbins = 0;

}

//====================================================================================================== 
BinningTools::~BinningTools(){
//======================================================================================================

  m_costhedges.clear();
  m_momedges.clear();
  m_costhslices.clear();
  m_nbins_mom.clear();
  delete[] m_bins_costh;

}

//======================================================================================================
void BinningTools::SetBinning(const string fname){
//======================================================================================================
  
  ifstream fin(fname.c_str());
  assert(fin.is_open());
  double p1, p2, cth1, cth2;
  int nlines = 0;
  while (1) {
    fin >> cth1 >> cth2 >> p1 >> p2;
    if(!fin.good()) break;
    m_costhedges.push_back(make_pair(cth1,cth2));
    m_momedges.push_back(make_pair(p1,p2));
    nlines++;
  }
  fin.close();
  
  m_nbins = m_momedges.size();

  int counter=0;

  for(size_t i = 0; i < m_costhedges.size(); i++){
    counter++;
    if ((m_costhedges[i].first != m_costhedges[i+1].first) && (m_costhedges[i].second != m_costhedges[i+1].second)){
      m_costhslices.push_back(make_pair(m_costhedges[i].first,m_costhedges[i].second));
      m_nbins_mom.push_back(counter);
      counter=0;
    }
  }

  m_nbins_costh = m_costhslices.size();

  m_bins_costh = new double[m_nbins_costh+1];

  for(int i = 0; i < m_nbins_costh+1; i++){
    if (i!=m_nbins_costh) m_bins_costh[i] = m_costhslices[i].first;
    else  m_bins_costh[i] = m_costhslices[m_nbins_costh-1].second;    
  }

  //mom bins stored 2-d vector
  vector <double> tmp;//temporary vector
  counter=0;
  int c=0;

  for(size_t i = 0; i < m_momedges.size(); i++){
    tmp.push_back(m_momedges[i].first);
    if ((m_costhedges[i].first != m_costhedges[i+1].first)) tmp.push_back(m_momedges[i].second);
    counter++;
    if(counter==m_nbins_mom[c]){
      m_bins_mom_vec.push_back(tmp);
      tmp.clear();
      c++;
      counter=0;
    }
  }

} 

//========================================================================================================
vector<std::pair<double, double> > BinningTools::Getpair_costh(){
//======================================================================================================
  
  return m_costhedges;

}

//========================================================================================================
vector<std::pair<double, double> > BinningTools::Getpair_mom(){
//======================================================================================================
  
  return m_momedges;

}

//========================================================================================================
double* BinningTools::GetCosthBins() {
//======================================================================================================

  return m_bins_costh;

}

//========================================================================================================
double** BinningTools::GetMomBins() {
//======================================================================================================
  
  m_bins_mom_arr = new double*[m_nbins_costh];//mom bins stored as 2-d array
  
  for(size_t i = 0; i < m_bins_mom_vec.size(); i++){
    for(size_t j=0; j<m_bins_mom_vec[i].size(); j++){
      m_bins_mom_arr[i] = new double[j];
    }
  }
  
  for(size_t i = 0; i < m_bins_mom_vec.size(); i++){
    for(size_t j=0; j<m_bins_mom_vec[i].size(); j++){
      m_bins_mom_arr[i][j] = m_bins_mom_vec[i][j];
    }
  }
  
  return m_bins_mom_arr;

}

//========================================================================================================
int* BinningTools::GetMomNBinsInEachCosthSlice() {
//======================================================================================================

  int* a = new int[m_bins_mom_vec.size()];

  for(size_t i=0; i<m_bins_mom_vec.size();i++) a[i] = m_bins_mom_vec[i].size()-1;
  return a;

}

//========================================================================================================
double* BinningTools::GetMomBinWidth() {
//======================================================================================================

  double* binwidth = new double[m_nbins];

  for(int nb=0; nb<m_nbins; nb++)
    binwidth[nb] = m_momedges[nb].second - m_momedges[nb].first;  
 
  return binwidth;
 
}


//========================================================================================================
double* BinningTools::GetCosBinWidth() {
//======================================================================================================

  double* binwidth = new double[m_nbins];

  for(int nb=0; nb<m_nbins; nb++)
    binwidth[nb] = m_costhedges[nb].second - m_costhedges[nb].first;  
 
  return binwidth;
 
}

