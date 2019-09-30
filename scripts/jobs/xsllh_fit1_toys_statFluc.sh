#! /bin/bash  

export WORKDIR='/sps/t2k/lmaret/softwares/xsLLhFitterLM/'

cd $WORKDIR/scripts/jobs/toSubmit

source /pbs/software/centos-7-x86_64/openmpi/ccenv.sh

for N in {1..100}; # loop over toys
do
	# Write a file containing the inputs for the fitter
	echo -e "{
    \"data_file\" : \"xsllh_nominal.root\",
    \"mc_file\"   : \"xsllh_nominal.root\",
    \"output_file\" : \"/outputs/toysStatFluc/fit1_statFluc_toy${N}.root\",
    \"input_dir\" : \"/inputs/fgd1fgd2Fit/\",
    \"fit_type\"  : 1,
    \"stat_fluc\" : true,
    \"zero_syst\" : false,
    \"data_POT\"  : 1.0,
    \"mc_POT\"    : 1.0,
    \"rng_seed\"  : $((12258 + N)),
    \"num_threads\" : 16,
    \"sample_topology\" : [\"cc0pi\", \"cc1pi\", \"ccother\", \"bkg\", \"oofv\"],
    \"min_settings\" : {
        \"minimizer\" : \"Minuit2\",
        \"algorithm\" : \"Migrad\",
        \"print_level\" : 2,
        \"tolerance\" : 1E-1,
        \"strategy\" : 1,
        \"max_iter\" : 1E6,
        \"max_fcn\" : 1E9
    },
    \"flux_cov\" : {
        \"file\" : \"xsllh_fluxcovmat.root\",
        \"matrix\" : \"flux_numu_cov\",
        \"binning\" : \"flux_binning\",
        \"throw\"  : false,
        \"decomp\" : true,
        \"variance\" : 0.99,
        \"fit_par\" : true
    },
    \"det_cov\" : {
        \"file\" : \"xsllh_detcovmat.root\",
        \"matrix\" : \"cov_mat\",
        \"throw\"  : false,
        \"decomp\" : true,
        \"variance\" : 0.99,
        \"fit_par\" : true
    },
    \"xsec_cov\" : {
        \"file\" : \"xsllh_xseccovmat.root\",
        \"matrix\" : \"xsec_cov\",
        \"throw\"  : false,
        \"decomp\" : false,
        \"variance\" : 1.00,
        \"fit_par\" : true
    },
    \"regularisation\" : {
        \"enable\" : false,
        \"strength\" : 1.5,
        \"method\" : \"kL2Reg\"
    },
    \"detectors\" : [
        {
            \"name\" : \"ND280\",
            \"xsec_config\" : \"ptheta_dials.json\",
            \"xsec_truth_config\" : \"ptheta_dials_truth.json\",
            \"binning\" : \"binning/tn337_binning.txt\",
            \"template_par\" : [
                {
                    \"name\" : \"CC0piCarbon\",
                    \"binning\" : \"binning/tn337_binning.txt\",
                    \"signal\" : {\"fgdtarget\" : [3]},
                    \"use\" : true
                },
                {
                    \"name\" : \"CC0piOxygen\",
                    \"binning\" : \"binning/tn337_binning.txt\",
                    \"signal\" : {\"fgdtarget\" : [0]},
                    \"use\" : true
                }
            ],
            \"use_detector\" : true
        }
    ],
    \"samples\" : [
        {
            \"cut_branch\" : -1,
            \"name\" : \"truth\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 0,
            \"name\" : \"FGD1muTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1muTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 1,
            \"name\" : \"FGD1muTPCpTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1muTPCpTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 2,
            \"name\" : \"FGD1muTPCpFGD\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1muTPCpFGD.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 3,
            \"name\" : \"FGD1muFGDpTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1muFGDpTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 4,
            \"name\" : \"FGD1muFGD\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1muFGD.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 5,
            \"name\" : \"FGD1CC1pi\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1CC1pi.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 6,
            \"name\" : \"FGD1CCDIS\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1CCDIS.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 7,
            \"name\" : \"FGD1CC1piMichelEle\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD1CC1piMichelEle.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 8,
            \"name\" : \"FGD2xmuTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xmuTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 9,
            \"name\" : \"FGD2xmuTPCpTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xmuTPCpTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 10,
            \"name\" : \"FGD2xmuTPCpFGD\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xmuTPCpFGD.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 11,
            \"name\" : \"FGD2xmuFGDpTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xmuFGDpTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 12,
            \"name\" : \"FGD2xmuFGD\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xmuFGD.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 13,
            \"name\" : \"FGD2xCC1pi\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xCC1pi.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 14,
            \"name\" : \"FGD2xCCDIS\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xCCDIS.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 15,
            \"name\" : \"FGD2xCC1piMichelEle\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2xCC1piMichelEle.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 16,
            \"name\" : \"FGD2ymuTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2ymuTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 17,
            \"name\" : \"FGD2ymuTPCpTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2ymuTPCpTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 18,
            \"name\" : \"FGD2ymuTPCpFGD\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2ymuTPCpFGD.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 19,
            \"name\" : \"FGD2ymuFGDpTPC\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2ymuFGDpTPC.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 20,
            \"name\" : \"FGD2ymuFGD\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2ymuFGD.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 21,
            \"name\" : \"FGD2yCC1pi\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2yCC1pi.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 22,
            \"name\" : \"FGD2yCCDIS\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2yCCDIS.txt\",
            \"use_sample\" : true
        },
        {
            \"cut_branch\" : 23,
            \"name\" : \"FGD2yCC1piMichelEle\",
            \"detector\" : \"ND280\",
            \"binning\" : \"binning/tn337_binning_FGD2yCC1piMichelEle.txt\",
            \"use_sample\" : true
        }
    ]
}
" > $WORKDIR/inputs/fgd1fgd2Fit/toysStatFluc/config_fit1_toy${N}.json

	# Write a file containing the job
	echo -e "cd $WORKDIR; source setup.sh; xsllhFit -j $WORKDIR/inputs/fgd1fgd2Fit/toysStatFluc/config_fit1_toy${N}.json" > job_xsllh_fit1_toy${N}.sh

	# Submit job
	qsub -l os=cl7,sps=1 -pe openmpi 16 -q pa_long job_xsllh_fit1_toy${N}.sh

done

cd -

###########################################################