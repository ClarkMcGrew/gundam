#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit
rm job_xsllh_fit*.sh

# Write scripts with jobs to submit

# Asimov fit
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_asimov.json" >> job_xsllh_fit1_asimov.sh

# Fake data with parameter thrown and statistical fluctuations
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit3_statFluc.json" >> job_xsllh_fit3_statFluc.sh

# Fake data with altered carbon or oxygen content
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30Carbon.json" >> job_xsllh_fit1_fakedata_m30Carbon.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30Carbon.json" >> job_xsllh_fit1_fakedata_p30Carbon.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30Oxygen.json" >> job_xsllh_fit1_fakedata_m30Oxygen.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30Oxygen.json" >> job_xsllh_fit1_fakedata_p30Oxygen.sh

# Fake data with altered background content
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30Res.json" >> job_xsllh_fit1_fakedata_m30Res.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30Res.json" >> job_xsllh_fit1_fakedata_p30Res.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30DIS.json" >> job_xsllh_fit1_fakedata_m30DIS.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30DIS.json" >> job_xsllh_fit1_fakedata_p30DIS.sh

# Fake data with GENIE and NuWro inputs
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_genie.json" >> job_xsllh_fit1_fakedata_genie.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_nuwro.json" >> job_xsllh_fit1_fakedata_nuwro.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_neutProd6D.json" >> job_xsllh_fit1_fakedata_neutProd6D.sh

source /usr/local/shared/bin/openmpi_env.sh

# Submit jobs

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_asimov.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit3_statFluc.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_genie.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_nuwro.sh
# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_neutProd6D.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30Carbon.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30Carbon.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30Oxygen.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30Oxygen.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30Res.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30Res.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30DIS.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30DIS.sh


cd - 

###########################################################
