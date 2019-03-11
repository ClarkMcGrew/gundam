
# Script containing all commands to get fitter inputs ready.
# 
# To submit job:
# qsub -l sps=1 -M maretlucie@gmail.com -m be runAllScripts_FitInputs.sh


cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/
source setup.sh
cd scripts/

# Compute detector covariance matrix
# xsllhDetVar -j ../inputs/fgd1fgd2Fit/detvar.json

# Draw pre-fit detector covariance matrix
# root -b -q 'DrawDetCov.C+()'

# Prepare flux covariance matrix in correct format
# xsllhFluxCov -j ../inputs/fgd1fgd2Fit/fluxcov.json

# Prepare xsec covariance matrix in correct format
# xsllhXsecCov -i ../inputs/fgd1fgd2Fit/xsec_cov.txt -o ../inputs/fgd1fgd2Fit/xsllh_xseccovmat.root


# DON'T FORGET : generate splines with T2KReWeight

# Get event selection from Highland2 outputs
#xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree.json
#xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_p30Carbon.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_m30Carbon.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_p30Oxygen.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_m30Oxygen.json
