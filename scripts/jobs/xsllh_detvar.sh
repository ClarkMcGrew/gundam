#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit

# Write scripts with jobs to submit

echo -e "cd $WORKDIR; source setup.sh; xsllhDetVar -j $WORKDIR/inputs/fgd1fgd2Fit/detvar.json" > job_xsllh_detvar.sh


# Submit jobs

qsub -l sps=1 -M maretlucie@gmail.com -m be job_xsllh_detvar.sh


cd -

###########################################################
