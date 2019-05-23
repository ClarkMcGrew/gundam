
cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/scripts

############################################################################################

# Draw pull study plots
root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit1_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit1_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  22, "fit1_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1164, "fit1_statFluc")'

root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit3_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit3_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  22, "fit3_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1164, "fit3_statFluc")'

# Draw chi2 study plots
root -b -q 'chi2_study.C("fit1_statFluc")'
root -b -q 'chi2_study.C("fit3_statFluc")'

###########################################################################################

# Draw parameter results
root -b -q 'DrawParameters.C+("fit1_asimov",          "asimov")'
root -b -q 'DrawParameters.C+("fit1_asimov_statFluc", "fakedata/statFluc")'
root -b -q 'DrawParameters.C+("fit1_statFluc",        "fakedata/statFluc")'
root -b -q 'DrawParameters.C+("fit3_statFluc",        "fakedata/statFluc")'

root -b -q 'DrawParameters.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawParameters.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawParameters.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawParameters.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawParameters.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawParameters.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawParameters.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawParameters.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawParameters.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawParameters.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'



# Draw event distributions
root -b -q 'DrawEventComparison.C+("fit1_asimov",          "asimov")'
root -b -q 'DrawEventComparison.C+("fit1_asimov_statFluc", "fakedata/statFluc")'
root -b -q 'DrawEventComparison.C+("fit1_statFluc",        "fakedata/statFluc")'
root -b -q 'DrawEventComparison.C+("fit3_statFluc",        "fakedata/statFluc")'

root -b -q 'DrawEventComparison.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'



# Draw post-fit covariance matrix
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov",          "asimov")'
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_statFluc", "fakedata/statFluc")'
root -b -q 'DrawCovariancePostfit.C+("fit1_statFluc",        "fakedata/statFluc")'
root -b -q 'DrawCovariancePostfit.C+("fit3_statFluc",        "fakedata/statFluc")'

root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'



############################################################################################



# Draw final xsec covariance matrix
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov",          "asimov")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_statFluc", "fakedata/statFluc")'
root -b -q 'DrawCovarianceFinal.C+("fit1_statFluc",        "fakedata/statFluc")'
root -b -q 'DrawCovarianceFinal.C+("fit3_statFluc",        "fakedata/statFluc")'

root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyDetPar",  "singleContrib")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_onlyXsecPar", "singleContrib")'



# Draw cross sections
root -b -q 'DrawXsec.C+("fit1_asimov",          "asimov")'
root -b -q 'DrawXsec.C+("fit1_asimov_statFluc", "fakedata/statFluc")'
root -b -q 'DrawXsec.C+("fit1_statFluc",        "fakedata/statFluc")'
root -b -q 'DrawXsec.C+("fit3_statFluc",        "fakedata/statFluc")'

root -b -q 'DrawXsec.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawXsec.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

root -b -q 'DrawXsec.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawXsec.C+("fit1_fakedata_m30Res", "fakedata/varyBKG")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Res", "fakedata/varyBKG")'
root -b -q 'DrawXsec.C+("fit1_fakedata_m30DIS", "fakedata/varyBKG")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30DIS", "fakedata/varyBKG")'

root -b -q 'DrawXsec.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
root -b -q 'DrawXsec.C+("fit1_asimov_onlyDetPar",  "singleContrib")'
root -b -q 'DrawXsec.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
root -b -q 'DrawXsec.C+("fit1_asimov_onlyXsecPar", "singleContrib")'






############################################################################################

# Draw statistical and systematic uncertainty from each type of parameters


root -b -q 'DrawParameters.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
root -b -q 'DrawParameters.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
root -b -q 'DrawParameters.C+("fit1_asimov_onlyXsecPar", "singleContrib")'
root -b -q 'DrawParameters.C+("fit1_asimov_onlyDetPar",  "singleContrib")'

root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyXsecPar", "singleContrib")'
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyDetPar",  "singleContrib")'


root -b -q 'DrawSystSingleContrib.C+()'


cd -
