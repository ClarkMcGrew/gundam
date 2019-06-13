#!/bin/sh

# Set cmake and gcc version required for xsec llh fitter
#ccenv cmake 3.11.0
#ccenv gcc 4.8.5
export CC=/pbs/software/scientific-6-x86_64/gcc/4.8.5/bin/gcc
export CXX=/pbs/software/scientific-6-x86_64/gcc/4.8.5/bin/g++
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/pbs/software/scientific-6-x86_64/gcc/4.8.5/lib64

#We will set up to build in a subdir of the source tree
#If it was sourced as . setup.sh then you can't scrub off the end... assume that
#we are in the correct directory.
if ! echo "${BASH_SOURCE}" | grep --silent "/"; then
  SETUPDIR=$(readlink -f $PWD)
else
  SETUPDIR=$(readlink -f ${BASH_SOURCE%/*})
fi
export XSLLHFITTER=${SETUPDIR}

BUILDSETUP="${XSLLHFITTER}/build/$(uname)/setup.sh"

echo "[INFO]: XSLLHFITTER root expected at: ${XSLLHFITTER}"

#source /usr/local/root/pro/bin/thisroot.sh

source ${XSLLHFITTER}/cmake/CMakeSetup.sh

if [ ! -e ${BUILDSETUP} ]; then
  echo "[INFO]: Cannot find build setup script where expected: ${BUILDSETUP}"
else
  echo "[INFO]: Sourcing build setup script."
  source ${BUILDSETUP}
  echo "[INFO]: \$ which xsllhFit: $(which xsllhFit)"
fi

echo "[INFO]: gcc Version `gcc --version` from: ${CXX}"
echo "[INFO]: ROOT Version `root-config --version` from: ${ROOTSYS}"



unset SETUPDIR
unset BUILDSETUP
