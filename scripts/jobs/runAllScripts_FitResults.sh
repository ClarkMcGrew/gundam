
cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/scripts

# # Draw parameter results
# root -b -q 'DrawParameters.C("fit1_asimov",   "asimov")'
# root -b -q 'DrawParameters.C("fit3_statFluc", "fakedata/statFluc")'

# root -b -q 'DrawParameters.C("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawParameters.C("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawParameters.C("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
# root -b -q 'DrawParameters.C("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

# root -b -q 'DrawParameters.C("fit1_fakedata_genie", "fakedata/genie")'
# root -b -q 'DrawParameters.C("fit1_fakedata_nuwro", "fakedata/nuwro")'


# # Draw event distributions
# root -b -q 'DrawEventComparison.C+("fit1_asimov",   "asimov")'
# root -b -q 'DrawEventComparison.C+("fit3_statFluc", "fakedata/statFluc")'

# root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
# root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

# root -b -q 'DrawEventComparison.C("fit1_fakedata_genie", "fakedata/genie")' 
# root -b -q 'DrawEventComparison.C("fit1_fakedata_nuwro", "fakedata/nuwro")' 



# # Draw post-fit covariance matrix
# root -b -q 'DrawCovariancePostfit.C+("fit1_asimov",   "asimov")'
# root -b -q 'DrawCovariancePostfit.C+("fit3_statFluc", "fakedata/statFluc")'

# root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

# root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_genie", "fakedata/genie")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'


# ############################################################################################


# # Draw statistical and systematic uncertainty from each type of parameters
# root -b -q 'DrawSystSingleContrib.C+()'

# # Draw final xsec covariance matrix
# root -b -q 'DrawCovarianceFinal.C+("fit1_asimov",   "asimov")'
# root -b -q 'DrawCovarianceFinal.C+("fit3_statFluc", "fakedata/statFluc")'

# root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
# root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
# root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

# root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_genie", "fakedata/genie")'
# root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'


# Draw cross sections
# root -b -q 'DrawXsec.C+("fit1_asimov",   "asimov")'
root -b -q 'DrawXsec.C+("fit3_statFluc", "fakedata/statFluc")'

root -b -q 'DrawXsec.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

root -b -q 'DrawXsec.C+("fit1_fakedata_genie", "fakedata/genie")'
root -b -q 'DrawXsec.C+("fit1_fakedata_nuwro", "fakedata/nuwro")'

cd -