#!/bin/bash

echo "ls before download:"
ls -lh

echo "Now downloading json file:"
cp /afs/cern.ch/work/c/cschloes/T2K_final_for_real/Fitter/xsLLhFitter/inputs/Numu_AntiNumu_Fit/Fit/joint/template_flux_detector_xsec/HTCondor_syst_throws/config_Fit.json .

echo "ls after download:"
ls -lh

#seed=$1
#output_filename=$2

echo "Now running fit:"
xsllhFit -j config_Fit.json -s $1
echo "Done running fit"

echo "ls after fit:"
ls -lh

echo "copy results back:"
cp fit_result_joint.root /afs/cern.ch/work/c/cschloes/T2K_final_for_real/Fitter/xsLLhFitter/inputs/Numu_AntiNumu_Fit/Fit/joint/template_flux_detector_xsec/HTCondor_syst_throws/output/fit_result_joint_$1.root

#echo "remove json file and fit_result file:"
#rm config_Fit.json
#rm fit_result_Numu_AntiNumu.root

#echo "ls after removal:"
#ls -lh


