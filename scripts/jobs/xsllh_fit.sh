#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit
# rm job_xsllh_fit*.sh

# Write scripts with jobs to submit

# Asimov fit
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_asimov.json" > job_xsllh_fit1_asimov.sh

# Fake data with parameter thrown and statistical fluctuations
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_asimov_statFluc.json"     > job_xsllh_fit1_asimov_statFluc.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit3_statFluc.json"            > job_xsllh_fit3_statFluc.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_asimov_statFluc_reg.json" > job_xsllh_fit1_asimov_statFluc_reg.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit3_statFluc_reg.json"        > job_xsllh_fit3_statFluc_reg.sh

echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit3_statFluc_reg_onlyXsec.json" > job_xsllh_fit3_statFluc_reg_onlyXsec.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit3_asimov_statFluc_reg.json"   > job_xsllh_fit3_asimov_statFluc_reg.sh

# Fake data with altered carbon or oxygen content
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30Carbon.json" > job_xsllh_fit1_fakedata_m30Carbon.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30Carbon.json" > job_xsllh_fit1_fakedata_p30Carbon.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30Oxygen.json" > job_xsllh_fit1_fakedata_m30Oxygen.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30Oxygen.json" > job_xsllh_fit1_fakedata_p30Oxygen.sh

# Fake data with altered background content
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30Res.json" > job_xsllh_fit1_fakedata_m30Res.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30Res.json" > job_xsllh_fit1_fakedata_p30Res.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_m30DIS.json" > job_xsllh_fit1_fakedata_m30DIS.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_p30DIS.json" > job_xsllh_fit1_fakedata_p30DIS.sh

# Fake data with GENIE and NuWro inputs
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_genie.json"     > job_xsllh_fit1_fakedata_genie.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_nuwro.json"     > job_xsllh_fit1_fakedata_nuwro.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_fakedata_nuwro_reg.json" > job_xsllh_fit1_fakedata_nuwro_reg.sh

# Fit only FGD2
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_asimov_onlyFGD2.json"          > job_xsllh_fit1_asimov_onlyFGD2.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit1_asimov_statFluc_onlyFGD2.json" > job_xsllh_fit1_asimov_statFluc_onlyFGD2.sh



# Data
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit2_dataCS_fakedataSignal.json"          > job_xsllh_fit2_dataCS_fakedataSignal.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit2_dataCS_fakedataSignal_noMichel.json" > job_xsllh_fit2_dataCS_fakedataSignal_noMichel.sh

echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit2_data.json"           > job_xsllh_fit2_data.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit2_data_noMichel.json"  > job_xsllh_fit2_data_noMichel.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit2_data_onlyCC1pi.json" > job_xsllh_fit2_data_onlyCC1pi.sh
echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/config_fit2_data_noFGD2ySB.json" > job_xsllh_fit2_data_noFGD2ySB.sh



source /pbs/software/centos-7-x86_64/openmpi/ccenv.sh

# Submit jobs

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_asimov.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_asimov_statFluc.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit3_statFluc.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_genie.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_nuwro.sh

# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_asimov_statFluc_reg.sh
# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit3_statFluc_reg.sh
# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit3_statFluc_reg_onlyXsec.sh
# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit3_asimov_statFluc_reg.sh
# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_nuwro_reg.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30Carbon.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30Carbon.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30Oxygen.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30Oxygen.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30Res.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30Res.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_m30DIS.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_fakedata_p30DIS.sh

# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_asimov_onlyFGD2.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit1_asimov_statFluc_onlyFGD2.sh



qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit2_dataCS_fakedataSignal.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit2_dataCS_fakedataSignal_noMichel.sh

qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit2_data.sh
qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit2_data_noMichel.sh


# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit2_data_onlyCC1pi.sh
# qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long -M maretlucie@gmail.com -m be job_xsllh_fit2_data_noFGD2ySB.sh


cd -

###########################################################
