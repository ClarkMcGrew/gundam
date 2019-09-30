
cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/scripts

############################################################################################

# Draw pull study plots
root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit1_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit1_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  26, "fit1_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1122, "fit1_statFluc")'

root -b -q 'pull_study.C("fit",  "Template parameters",   116, "fit3_statFluc")'
root -b -q 'pull_study.C("flux", "Flux parameters",        20, "fit3_statFluc")'
root -b -q 'pull_study.C("xsec", "Xsec model parameters",  26, "fit3_statFluc")'
root -b -q 'pull_study.C("det",  "Detector parameters",  1122, "fit3_statFluc")'

# Draw chi2 study plots
root -b -q 'chi2_study.C("fit1_statFluc")'
root -b -q 'chi2_study.C("fit3_statFluc")'

# ###########################################################################################

# # Draw L-curves
# root -b -q 'regularisation_Lcurve.C+("fit1_asimov_statFluc",     "fakedata/statFlucReg")'
# root -b -q 'regularisation_Lcurve.C+("fit3_statFluc",            "fakedata/statFlucReg)'

############################################################################################


# Draw parameter results
root -b -q 'DrawParameters.C+("fit1_asimov",               "asimov")'
root -b -q 'DrawParameters.C+("fit1_asimov_statFluc",      "fakedata/statFluc")'
root -b -q 'DrawParameters.C+("fit3_statFluc",             "fakedata/statFluc")'

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

root -b -q 'DrawParameters.C+("fit1_asimov_statFluc_onlyFGD2", "onlyFGD2")'

root -b -q 'DrawParameters.C+("fit2_dataCS_fakedataSignal", "data")'
root -b -q 'DrawParameters.C+("fit2_dataCS_fakedataSignal_noMichel", "data")'

root -b -q 'DrawParameters.C+("fit2_data",           "data")'
root -b -q 'DrawParameters.C+("fit2_data_noMichel",  "data")'


# Draw event distributions
root -b -q 'DrawEventComparison.C+("fit1_asimov",               "asimov")'
root -b -q 'DrawEventComparison.C+("fit1_asimov_statFluc",      "fakedata/statFluc")'
root -b -q 'DrawEventComparison.C+("fit3_statFluc",             "fakedata/statFluc")'

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

root -b -q 'DrawEventComparison.C+("fit2_dataCS_fakedataSignal",                   "data")'
root -b -q 'DrawEventComparison_noMichel.C+("fit2_dataCS_fakedataSignal_noMichel", "data")'

root -b -q 'DrawEventComparison.C+(          "fit2_data",           "data")'
root -b -q 'DrawEventComparison_noMichel.C+( "fit2_data_noMichel",  "data")'


# Draw post-fit covariance matrix
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov",               "asimov")'
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_statFluc",      "fakedata/statFluc")'
root -b -q 'DrawCovariancePostfit.C+("fit3_statFluc",             "fakedata/statFluc")'

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

root -b -q 'DrawCovariancePostfit.C+("fit2_dataCS_fakedataSignal",                   "data")'
root -b -q 'DrawCovariancePostfit_noMichel.C+("fit2_dataCS_fakedataSignal_noMichel", "data")'

root -b -q 'DrawCovariancePostfit.C+("fit2_data",                   "data")'
root -b -q 'DrawCovariancePostfit_noMichel.C+("fit2_data_noMichel", "data")'


############################################################################################




# Draw cross sections
root -b -q 'DrawXsec.C+("fit1_asimov",               "asimov")'
root -b -q 'DrawXsec.C+("fit1_asimov_statFluc",      "fakedata/statFluc")'
root -b -q 'DrawXsec.C+("fit3_statFluc",             "fakedata/statFluc")'

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

root -b -q 'DrawXsec.C+("fit1_asimov_statFluc_onlyFGD2", "onlyFGD2")'

root -b -q 'DrawXsec.C+("fit2_dataCS_fakedataSignal",          "data")'
root -b -q 'DrawXsec.C+("fit2_dataCS_fakedataSignal_noMichel", "data")'

root -b -q 'DrawXsec.C+("fit2_data",           "data")'
root -b -q 'DrawXsec.C+("fit2_data_noMichel",  "data")'




# Draw final xsec covariance matrix
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov",               "asimov")'
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_statFluc",      "fakedata/statFluc")'
root -b -q 'DrawCovarianceFinal.C+("fit3_statFluc",             "fakedata/statFluc")'

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

root -b -q 'DrawCovarianceFinal.C+("fit1_asimov_statFluc_onlyFGD2", "onlyFGD2")'

root -b -q 'DrawCovarianceFinal.C+("fit2_dataCS_fakedataSignal",          "data")'
root -b -q 'DrawCovarianceFinal.C+("fit2_dataCS_fakedataSignal_noMichel", "data")'

root -b -q 'DrawCovarianceFinal.C+("fit2_data",           "data")'
root -b -q 'DrawCovarianceFinal.C+("fit2_data_noMichel",  "data")'




############################################################################################

# Draw statistical and systematic uncertainty from each type of parameters


# root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyTemPar",  "singleContrib")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyFluxPar", "singleContrib")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyXsecPar", "singleContrib")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_asimov_onlyDetPar",  "singleContrib")'


# root -b -q 'DrawSystSingleContrib.C+()'



############################################################################################

# Draw selection efficiencies

# source /usr/local/root/pro/bin/thisroot.sh
# root -b -q 'signalEfficiency2D_OandC.C+("C")'
# root -b -q 'signalEfficiency2D_OandC.C+("O")'
# root -b -q 'signalEfficiency2D_regions.C+()'  #by region
# root -b -q 'signalEfficiency.C+("C")'         #by generator
# root -b -q 'signalEfficiency.C+("O")'         #by generator
# root -b -q 'protonFSIefficiency.C+()'         #for proton FSI (nuwro)

############################################################################################


xsecfitter
cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/scripts/jobs
