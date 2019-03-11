
# Draw parameter results
root -b -q 'DrawParameters.C("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawParameters.C("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawParameters.C("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawParameters.C("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'


# Draw event distributions
root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawEventComparison.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'


# Draw post-fit covariance matrix
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawCovariancePostfit.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'



# Draw cross sections
root -b -q 'DrawXsec.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawXsec.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'

# Draw final xsec covariance matrix
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Carbon", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_m30Oxygen", "fakedata/varyOC")'
root -b -q 'DrawCovarianceFinal.C+("fit1_fakedata_p30Oxygen", "fakedata/varyOC")'
