from pathlib import Path
from typing import List


def partition_root_files(n: int, path_input: Path) -> List[List[str]]:
    """Partition a directory of root files into n groups.

    :param n: Number of groups to split the files into.
    :type n: int
    :param path_input: Directory containing the root files to split.
    :type path_input: Path
    :return: List of n lists of str path to root file.
    :rtype: List[List[Path]]
    """
    input_files = list(path_input.glob("*.root"))
    input_files.sort()

    n_files = len(input_files)
    n_files_per_job = n_files // n
    if n_files % n != 0:
        # leftover files
        n_files_per_job += 1

    output_files = []
    for i in range(n):
        start = i * n_files_per_job
        end = (i + 1) * n_files_per_job
        input_files_group = input_files[start:end]
        # convert to string
        input_files_group = [str(file) for file in input_files_group]
        output_files.append(input_files_group)
    return output_files


if __name__ == "__main__":
    # Usage:
    # python3 submitJob_LPC_llpMerge_hadd.py --input {ntupler path} --output {output dir of hadded skim} --job-dir {dir to job script} --n {number of hadded files}
    # Example:
    # python3 submitJob_LPC_llpMerge_hadd.py --input /eos/uscms/store/user/xxxx/llp/skim_Run2023D_Muon0 --output /eos/uscms/store/user/xxxx/llp/skim_hadd_Run2023C_Muon0 --job hadd_test --n 10
    import argparse
    import os

    arg_parser = argparse.ArgumentParser(
        description="Hadd input ntupler skimmed files into ntupler files."
    )
    arg_parser.add_argument(
        "--n",
        required=True,
        type=int,
        help="Number of groups to split the input ntupler file into.",
    )
    arg_parser.add_argument(
        "--input", required=True, help="Path that contains ntupler files."
    )
    arg_parser.add_argument(
        "--output", required=True, help="Path to the output directory."
    )
    arg_parser.add_argument("--job", required=True, help="Path to the job directory.")
    args = arg_parser.parse_args()

    n = args.n
    path_output = Path(args.output)
    path_input = Path(args.input)
    path_job = Path(args.job)
    path_output.mkdir(parents=True, exist_ok=True)
    path_job.mkdir(parents=True, exist_ok=True)

    script_sh = """
#!/bin/bash
date
MAINDIR=`pwd`
ls
voms-proxy-info --all
#CMSSW from scratch (only need for root)
export CWD=${PWD}
export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export SCRAM_ARCH=slc7_amd64_gcc700
scramv1 project CMSSW CMSSW_10_6_20
cp MergeNtuples CMSSW_10_6_20/src
cd CMSSW_10_6_20/src
eval `scramv1 runtime -sh` # cmsenv
echo "Untar JEC:" 
echo "After Untar: "
echo "ls -la"
ls -la
echo "Running job..."
# jobs
"""
    files_group_list = partition_root_files(n, path_input)
    for i, files_group in enumerate(files_group_list):
        files_group = [
            f.replace("/eos/uscms/", "root://cmseos.fnal.gov//") 
            for f in files_group
        ]
        ntupler_files = " ".join(files_group)
        output_dir = path_output / f"ntupler_{i}.root"
        output_dir = str(output_dir)
        output_dir = output_dir.replace("/eos/uscms/", "root://cmseos.fnal.gov//")
        script_sh += f"hadd -f {output_dir} {ntupler_files}\n"

    # write to file
    path_script_sh = path_job / "runjob.sh"
    path_script_sh.write_text(script_sh.strip())
    os.system(f"chmod +x {path_script_sh}")

    script_jdl = """
universe = vanilla
Executable = runjob.sh
Should_Transfer_Files = YES 
WhenToTransferOutput = ON_EXIT_OR_EVICT
Output = runjob.$(Process).$(Cluster).stdout
Error  = runjob.$(Process).$(Cluster).stdout
Log    = runjob.$(Process).$(Cluster).log
Arguments = $(Process) 1
Queue 1
"""
    path_script_jdl = path_job / "runjob.jdl"
    path_script_jdl.write_text(script_jdl.strip())

    print(f"Script written to {path_job}")
