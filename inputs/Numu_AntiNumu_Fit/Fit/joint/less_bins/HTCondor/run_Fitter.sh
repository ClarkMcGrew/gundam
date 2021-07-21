#!/bin/bash

echo "ls before download:"
ls -lh

echo "Now downloading json file:"
cp /afs/cern.ch/work/c/cschloes/T2K_final_for_real/Fitter/xsLLhFitter/inputs/Numu_AntiNumu_Fit/Fit/joint/less_bins/$1.json .

echo "ls after download:"
ls -lh

#seed=$1
#output_filename=$2

echo "Now running fit:"
xsllhFit -j $1.json -o fitter_output_$1.root
echo "Done running fit"

echo "ls after fit:"
ls -lh

echo "copy results back:"
cp fitter_output_$1.root /afs/cern.ch/work/c/cschloes/T2K_final_for_real/Fitter/xsLLhFitter/inputs/Numu_AntiNumu_Fit/Fit/joint/less_bins/HTCondor/output/

#echo "remove json file and fit_result file:"
#rm config_Fit.json
#rm fit_result_Numu_AntiNumu.root

#echo "ls after removal:"
#ls -lh


