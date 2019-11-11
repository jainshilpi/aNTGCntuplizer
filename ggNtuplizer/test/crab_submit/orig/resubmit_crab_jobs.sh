#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh;
source /cvmfs/cms.cern.ch/crab3/crab.sh;
eval `scramv1 runtime -sh`;
voms-proxy-init --voms cms



commandval=$1
inputFile=$2

while IFS= read -r line
do
	if [[ $line == *"Please use 'crab status"* ]]; then
		crab_dir=${line//"Please use 'crab status -d "/}
		crab_dir=${crab_dir//"' to check how the submission process proceeds."/}

		if [ ! -d "$crab_dir" ]; then
			echo "Error! Directory "${crab_dir} "does not exist!"
			continue;
		fi
		echo "Resubmitting "${crab_dir}
		crab $commandval -d ${crab_dir}
	fi
done < ${inputFile}
