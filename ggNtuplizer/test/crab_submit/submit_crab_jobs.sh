#!/bin/bash

input_datasets="METbackgroundsV1.txt"
# input_datasets="2016_mc_samples.txt"
#writedir="/afs/cern.ch/work/m/mwadud/private/naTGC/CMSSW_9_4_13/src/ggAnalysis/ggNtuplizer/test/crab_submit/jobsMETv2/"
writedir="/afs/hep.wisc.edu/home/shilpi/private/forMETStudy/ntuplesfor2017_METv3/CMSSW_9_4_13/src/ggAnalysis/ggNtuplizer/test/crab_submit/MCv3"

psetname="/afs/hep.wisc.edu/home/shilpi/private/forMETStudy/ntuplesfor2017_METv3/CMSSW_9_4_13/src/ggAnalysis/ggNtuplizer/test/run_mc2017_94X.py"

writeSite="T2_US_Wisconsin"
mainOutputDir='/store/user/shilpi/aNTGC/ggNtuplizerSkim/METv3/mc'

lumiMaskFile="https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"

#maxFiles=5000
inputDBS=global
# inputDBS=phys03


################################################################################################################
source /cvmfs/cms.cern.ch/cmsset_default.sh;
source /cvmfs/cms.cern.ch/crab3/crab.sh;
eval `scramv1 runtime -sh`;
voms-proxy-init --voms cms


crab_cfg_template=crab_submit.py
units_perjob=2
splitting='FileBased'
#blacklist="'T2_US_Florida','T2_US_Vanderbilt','T3_US_PuertoRico'"
whitelist="'T3_US_UCR','T3_US_FNALLPC','T2_US_Purdue','T3_US_Rice','T3_US_Rutgers','T3_US_FIU','T3_US_FIT','T3_US_PSC','T3_US_OSU','T3_US_TAMU','T3_US_UMD','T3_US_VC3_NotreDame','T3_US_SDSC','T3_US_Colorado','T3_US_OSG','T3_US_Princeton_ICSE','T3_US_NERSC','T3_US_Baylor','T2_US_Nebraska','T2_US_UCSD','T2_US_Wisconsin','T2_US_MIT','T3_US_TACC','T3_US_TTU','T3_US_UMiss','T2_US_Caltech'"
blacklist=""
#whitelist="'T3_US_FNALLPC'"

mkdir -p ${writedir}
submit_log_file=${writedir}/crab_submitted_datasets_$(date '+%d_%m_%Y_%H_%M_%S').log

echo "*************************************************************************" | tee --append ${submit_log_file}
echo "Git Info" | tee --append ${submit_log_file}
git remote -v | tee --append ${submit_log_file}
git branch -a | tee --append ${submit_log_file}
git rev-parse HEAD | tee --append ${submit_log_file}
echo "*************************************************************************" | tee --append ${submit_log_file}


for dataset in `sed '/^$/d' ${input_datasets}`;
do
	 jobName=$(echo ${dataset} | cut -f1,2 -d'/')
	 jobName=${jobName#"/"}

#	jobName=$(echo ${dataset////_})
#	jobName=${jobName#"_"}
	jobName=$(echo ${jobName} | sed 's/[^a-zA-Z0-9]//g')

	jobDir=${writedir}/${jobName}/

	# rm -rf $jobDir

	# rm -rf crab_${jobName}/

	echo "Submitting " ${jobName}
	echo "Job directory: "${jobDir}
	echo "Write site: " ${writeSite}

	if [ -d "$jobDir" ]; then
		echo "Error! Directory "$jobDir "exists!"
		echo -e "\t\t\tSkipping!\n\n"
		# exit
		continue
	fi

	mkdir -p ${jobDir}

	crab_cfg_file=${jobDir}/crab_${jobName}.py
	cp ${crab_cfg_template} ${crab_cfg_file}

	sed -i 's|#psetname|'$psetname'|g' ${crab_cfg_file}
	sed -i 's|#workarea|'$jobDir'|g' ${crab_cfg_file}
	sed -i 's|#jobname|'$jobName'|g' ${crab_cfg_file}
	sed -i 's|#inputDataset|'$dataset'|g' ${crab_cfg_file}
	sed -i 's|#storageSite|'$writeSite'|g' ${crab_cfg_file}
	sed -i 's|#whitelist|'$whitelist'|g' ${crab_cfg_file}
	sed -i 's|#blacklist|'$blacklist'|g' ${crab_cfg_file}
	sed -i 's|#splitting|'$splitting'|g' ${crab_cfg_file}
	sed -i 's|#unitsperjob|'$units_perjob'|g' ${crab_cfg_file}
	sed -i 's|#mainOutputDir|'$mainOutputDir'|g' ${crab_cfg_file}

	sed -i 's|#inputDBS|'$inputDBS'|g' ${crab_cfg_file}


	# sed -i 's|#config.Data.totalUnits|'config.Data.totalUnits'|g' ${crab_cfg_file}
	# sed -i 's|#totalUnits|'$maxFiles'|g' ${crab_cfg_file}

	#sed -i 's|#config.Data.lumiMask|'config.Data.lumiMask'|g' ${crab_cfg_file}
	#sed -i 's|#lumiMaskFile|'${lumiMaskFile}'|g' ${crab_cfg_file}

	python ${crab_cfg_file} | tee --append ${submit_log_file}
done
