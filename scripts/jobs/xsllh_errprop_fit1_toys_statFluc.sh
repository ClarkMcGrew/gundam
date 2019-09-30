#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit

source /pbs/software/centos-7-x86_64/openmpi/ccenv.sh


for N in {1..100}; # loop over toys
do
	# Write a file containing the inputs for the fitter
	echo -e "{
    \"input_fit_file\" : \"/outputs/toysStatFluc/fit1_statFluc_toy${N}.root\",
    \"output_file\" : \"/outputs/toysStatFluc/xsec_fit1_statFluc_toy${N}.root\",
    \"extra_hists\" : \"\",
    \"proton_fsi_cov\" : \"/inputs/fgd1fgd2Fit/xsllh_covarProtonFSI.root\",
    \"num_toys\" : 10000,
    \"rng_seed\" : $((12258 + N)),
    \"sel_config\" : \"/inputs/fgd1fgd2Fit/toysStatFluc/config_fit1_toy${N}.json\",
    \"tru_config\" : \"/inputs/fgd1fgd2Fit/toysStatFluc/config_fit1_toy${N}.json\",
    \"read_data_events\" : true,
    \"decomposition\" : {
        \"do_force_posdef\" : true,
        \"force_posdef_val\" : 1E-6,
        \"incomplete_chol\" : false,
        \"drop_tolerance\" : 1E-9
    },
    \"sig_norm\" : {
        \"CC0piCarbon\" : {
            \"flux_file\" : \"/inputs/flux/weighted_flux13av4_run2-3-4-8.root\",
            \"flux_hist\" : \"flux_rebin\",
            \"flux_int\"  : 2.13159E+13,
            \"flux_err\"  : 0.085,
            \"use_flux_fit\" : true,
            \"num_targets_val\" : 7.439E29,
            \"num_targets_err\" : 0.0048,
            \"relative_err\" : true
        },
        \"CC0piOxygen\" : {
            \"flux_file\" : \"/inputs/flux/weighted_flux13av4_run2-3-4-8.root\",
            \"flux_hist\" : \"flux_rebin\",
            \"flux_int\"  : 2.13159E+13,
            \"flux_err\"  : 0.085,
            \"use_flux_fit\" : true,
            \"num_targets_val\" : 2.581E29,
            \"num_targets_err\" : 0.0073,
            \"relative_err\" : true
        }
    }
}
    
" > $WORKDIR/inputs/fgd1fgd2Fit/toysStatFluc/errprop_fit1_toy${N}.json

	# Write a file containing the job
	echo -e "cd $WORKDIR; source setup.sh; xsllhCalcXsec -j $WORKDIR/inputs/fgd1fgd2Fit/toysStatFluc/errprop_fit1_toy${N}.json" > job_xsllh_errprop_fit1_toy${N}.sh

	# Submit job
	qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long job_xsllh_errprop_fit1_toy${N}.sh

done

cd -

###########################################################
