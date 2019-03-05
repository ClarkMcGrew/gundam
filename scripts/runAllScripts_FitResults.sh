
# Draw parameter results
root -b -q 'DrawParameters.C("fit1_asimov")'
# root -b -q 'DrawParameters.C("fit1_statFluc")'
root -b -q 'DrawParameters.C("fit3_statFluc")'


# Draw event distributions
root -b -q 'DrawEventComparison.C+("fit1_asimov")'
# root -b -q 'DrawEventComparison.C+("fit1_statFluc")'
root -b -q 'DrawEventComparison.C+("fit3_statFluc")'


# Draw cross sections
root -b -q 'DrawXsec.C+("fit1_asimov")'
# root -b -q 'DrawXsec.C+("fit1_statFluc")'
root -b -q 'DrawXsec.C+("fit3_statFluc")'


# Draw pre-fit detector covariance matrix
root 'DrawDetCov.C+()'


# Draw post-fit covariance matrix
root -b -q 'DrawCovariancePostfit.C+("fit1_asimov")'
# root -b -q 'DrawCovariancePostfit.C+("fit1_statFluc")'
root -b -q 'DrawCovariancePostfit.C+("fit3_statFluc")'


# Draw final xsec covariance matrix
root -b -q 'DrawCovarianceFinal.C+("fit1_asimov")'
# root -b -q 'DrawCovarianceFinal.C+("fit1_statFluc")'
root -b -q 'DrawCovarianceFinal.C+("fit3_statFluc")'
