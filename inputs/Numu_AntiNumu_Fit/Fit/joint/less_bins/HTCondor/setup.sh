#!/bin/bash

# Creates a list with all the oaAnalysis files which the analysis then runs over:
#python make_filelist.py

# Delete all previous HTCondor log, error, out files. Creates the /htcondor_OutErrLog_files directory if it does not exist yet:
mkdir -p HTCondor_OutErrLog_files
rm HTCondor_OutErrLog_files/*

mkdir -p output

# Submits the jobs to HTCondor and queues the files from the generated file list:
condor_submit submit.sub
