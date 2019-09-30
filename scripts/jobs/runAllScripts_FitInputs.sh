# Script containing all commands to get fitter inputs ready.
# 
# To submit job:
# qsub -l sps=1 -M maretlucie@gmail.com -m be runAllScripts_FitInputs.sh


cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/
source setup.sh
cd /sps/t2k/lmaret/softwares/xsLLhFitterLM/scripts/



# Compute detector covariance matrix
xsllhDetVar -j ../inputs/fgd1fgd2Fit/detvar.json
xsllhDetVar -j ../inputs/fgd1fgd2Fit/detvar_onlyFGD2.json

# Draw pre-fit detector covariance matrix
root -b -q 'DrawDetCov.C+()'

# Prepare flux covariance matrix in correct format
# xsllhFluxCov -j ../inputs/fgd1fgd2Fit/fluxcov.json

# Prepare flux tuning used to compute the flux integrated
# root -b -q 'weight_flux.C()'

# Prepare xsec covariance matrix in correct format
# xsllhXsecCov -i ../inputs/fgd1fgd2Fit/xsec_cov.txt -o ../inputs/fgd1fgd2Fit/xsllh_xseccovmat.root


# DON'T FORGET : generate splines with T2KReWeight

# Get event selection from Highland2 outputs

xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_nominal_normalised.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_nominal_normalisedToNuWro_new.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_nominal_normalised_onlyFGD2.json

xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_genie.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_nuwro_new.json
# xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_nuwro.json

xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_p30Carbon.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_m30Carbon.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_p30Oxygen.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_m30Oxygen.json

xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_p30Res.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_m30Res.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_p30DIS.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_fakedata_m30DIS.json


xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_data.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_dataCS_fakedataSignal.json
xsllhTreeConvert -j ../inputs/fgd1fgd2Fit/tree_dataCS_fakedataSignal_noMichel.json