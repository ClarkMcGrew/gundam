

# Compute detector covariance matrix
xsllhDetVar -j ../inputs/fgd1fgd2Fit/detvar.json

# Prepare flux covariance matrix in correct format
# xsllhFluxCov -j ../inputs/fgd1fgd2Fit/fluxcov.json

# Prepare xsec covariance matrix in correct format
# xsllhXsecCov -i ../inputs/fgd1fgd2Fit/xsec_cov.txt -o ../inputs/fgd1fgd2Fit/xsllh_xseccovmat.root


# DON'T FORGET : generate splines with T2KReWeight

# Get event selection from Highland2 outputs
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree.json