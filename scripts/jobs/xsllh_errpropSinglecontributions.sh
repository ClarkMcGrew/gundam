#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit
# rm job_xsllh_errprop_singleContr_*.sh

# Write scripts with jobs to submit

echo -e "cd $WORKDIR; source setup.sh; xsllhCalcXsec -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/errprop_onlyTemPar.json"  > job_xsllh_errprop_singleContr_onlyTemPar.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhCalcXsec -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/errprop_onlyFluxPar.json" > job_xsllh_errprop_singleContr_onlyFluxPar.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhCalcXsec -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/errprop_onlyXsecPar.json" > job_xsllh_errprop_singleContr_onlyXsecPar.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhCalcXsec -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/errprop_onlyDetPar.json"  > job_xsllh_errprop_singleContr_onlyDetPar.sh

#source /usr/local/shared/bin/openmpi_env.sh

# Submit jobs

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_errprop_singleContr_onlyTemPar.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_errprop_singleContr_onlyFluxPar.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_errprop_singleContr_onlyXsecPar.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_errprop_singleContr_onlyDetPar.sh

cd -

###########################################################
