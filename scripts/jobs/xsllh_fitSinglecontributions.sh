#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit
# rm job_xsllh_singleContr*.sh

# Write scripts with jobs to submit

echo -e "lscpu; cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/config_onlyTemPar.json"  > job_xsllh_singleContr_onlyTemPar.sh
echo -e "lscpu; cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/config_onlyFluxPar.json" > job_xsllh_singleContr_onlyFluxPar.sh
echo -e "lscpu; cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/config_onlyXsecPar.json" > job_xsllh_singleContr_onlyXsecPar.sh
echo -e "lscpu; cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/singleContrib/config_onlyDetPar.json"  > job_xsllh_singleContr_onlyDetPar.sh

source /usr/local/shared/bin/openmpi_env.sh

# Submit jobs

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_singleContr_onlyTemPar.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_singleContr_onlyFluxPar.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_singleContr_onlyXsecPar.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_singleContr_onlyDetPar.sh

cd - 

###########################################################
