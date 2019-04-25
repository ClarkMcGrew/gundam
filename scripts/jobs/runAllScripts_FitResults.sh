
cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/scripts

############################################################################################

# # Draw pull study
# root -b -q 'pull_study.C("fit",  "Template parameters",   116)'
# root -b -q 'pull_study.C("flux", "Flux parameters",        20)'
# root -b -q 'pull_study.C("xsec", "Xsec model parameters",  22)'
# root -b -q 'pull_study.C("det",  "Detector parameters",  1164)'

# root -b -q 'chi2_study.C()'

############################################################################################

# Draw parameter results
root -b -q 'DrawParameters.C("fit1_asimov",   "asimov")'
root -b -q 'DrawParameters.C("fit3_statFluc", "fakedata/statFluc")'

root -b -q 'DrawParameters.C("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawParameters.C("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawParameters.C("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawParameters.C("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawParameters.C("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawParameters.C("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawParameters.C("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawParameters.C("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawParameters.C("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawParameters.C("fit1_fakedata_nuwro", "fakedata/nuwro")'



# Draw event distributions
root -b -q 'DrawEventComparison.C("fit1_asimov",   "asimov")'
root -b -q 'DrawEventComparison.C("fit3_statFluc", "fakedata/statFluc")'

root -b -q 'DrawEventComparison.C("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawEventComparison.C("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawEventComparison.C("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawEventComparison.C("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawEventComparison.C("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawEventComparison.C("fit1_fakedata_genie", "fakedata/genie")' 
root -b -q 'DrawEventComparison.C("fit1_fakedata_nuwro", "fakedata/nuwro")' 



# Draw post-fit covariance matrix
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov",   "asimov")'
root -b -q 'DrawCovariancePostfit.C+("fit3_statFluc", "fakedata/statFluc")'

root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'



############################################################################################



# Draw final xsec covariance matrix
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov",   "asimov")'
root -b -q 'DrawCovarianceFinal.C+("fit3_statFluc", "fakedata/statFluc")'

root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyDetPar",  "singleContrib")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyXsecPar", "singleContrib")'



# Draw cross sections
root -b -q 'DrawXsec.C+("fit1_asimov",   "asimov")'
root -b -q 'DrawXsec.C+("fit3_statFluc", "fakedata/statFluc")'

root -b -q 'DrawXsec.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawXsec.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawXsec.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawXsec.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawXsec.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawXsec.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
root -b -q 'DrawXsec.C+("fit1_asimov_onlyDetPar",  "singleContrib")'
root -b -q 'DrawXsec.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
root -b -q 'DrawXsec.C+("fit1_asimov_onlyXsecPar", "singleContrib")'



# Draw statistical and systematic uncertainty from each type of parameters
root -b -q 'DrawSystSingleContrib.C+()'


cd -
