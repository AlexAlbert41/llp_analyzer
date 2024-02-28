#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/run_printJson.sh
filesPerJob=10
ver=V1p17

listData2016=(
)
listData2017=(
)
listData2018=(
2018B
)
listData2018ABC_AOD=(
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018
)

listData2018ABC_RAW=(
Run2_displacedJetMuonNtupler_V1p17_Data2018ABC_RAW_Run2018B-v1_v1_v3
)

for year in \
Data2018ABC_RAW
do
	echo ${year}
	sampleList=list${year}[@]
        for sample in "${!sampleList}"
        do
		echo "Sample " ${sample}
		version=/${ver}/${year}/v1/
		#output=/storage/af/group/phys_exotica/delayedjets/HNL/json/${version}/v2/${sample}
		output=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonNtuple/json/${version}/sixie/EGamma/${sample}
		echo ${output}
		#inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${ver}/${year}/One_Lepton_Merged_MS_Hits/${sample}.txt

		inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/sixie/EGamma/${sample}.txt
		
		nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
		mod=`python -c "print int($nfiles.0%$filesPerJob)"`
	        if [ ${mod} -eq 0 ]
	        then
	                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	        fi
		rm -f submit/lumi_${sample}*
		rm -f log/lumi_${sample}*

		echo "job " ${maxjob}
	        jdl_file=submit/lumi_${sample}_${maxjob}.jdl
	        echo "Universe = vanilla" > ${jdl_file}
	        echo "Executable = ${job_script}" >> ${jdl_file}
	        echo "Arguments = ${inputfilelist} ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}


	        # option should always be 1, when running condor
	        echo "Log = log/lumi_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
	        echo "Output = log/lumi_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	        echo "Error = log/lumi_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

	        #echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
	        echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}

	        echo "RequestMemory = 4000" >> ${jdl_file}
	        echo "RequestCpus = 1" >> ${jdl_file}
	        echo "RequestDisk = 4" >> ${jdl_file}

	        echo "+RunAsOwner = True" >> ${jdl_file}
	        echo "+InteractiveUser = true" >> ${jdl_file}
	        echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7\"" >> ${jdl_file}
		echo '+SingularityBindCVMFS = True' >> ${jdl_file}
	        echo "run_as_owner = True" >> ${jdl_file}
	        echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
	        echo "should_transfer_files = YES" >> ${jdl_file}
	        echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	        echo "Queue ${maxjob}" >> ${jdl_file}
	        echo "condor_submit ${jdl_file}"
	        condor_submit ${jdl_file} -batch-name ${sample}


	done
done
