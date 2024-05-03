from pathlib import Path
from typing import List


CMSSW: str = "CMSSW_10_6_20"  # CMSSW version to use


def split_nanoAOD_txt_file(n: int, file_dir: Path) -> List[List[str]]:
    """Split a txt file that contains many nanoAOD file paths,
    one per line, into n groups.
    Note that some elements in the end might be empty.

    :param n: Number of groups to split the files into.
    :type n: int
    :param file_dir: Directory to the txt file to split.
    :type file_dir: Path
    :return: List of n lists of str path to txt file.
    :rtype: List[List[str]]
    """
    with open(file_dir, "r") as f:
        input_files = f.readlines()

    if n >= len(input_files):
        n = len(input_files)
        n_files_per_job = 1
    else:
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
        # remove newline character
        input_files_group = [file.strip() for file in input_files_group]
        output_files.append(input_files_group)
    return output_files


if __name__ == "__main__":
    # Example:
    # python3 submitJob_LPC_llpMerge_split.py --exe ./MergeNTuples --dir-ntupler /eos/uscms/store/user/xxx/llp/skim_hadd_Run2023C_Muon0 --NanoAOD lists/Run3_nanoAOD/Muon0_2023C.txt --n 50 --job-dir merge_split_test --dir-output /eos/uscms/store/user/xxx/llp/skim_hadd_Run2023C_Muon0 --copy-root
    
    # parse arguments
    import argparse

    arg_parser = argparse.ArgumentParser(
        description="Split input files for LPC llpMerge job."
    )
    arg_parser.add_argument("--exe", required=True, help="Path to the executable.")
    arg_parser.add_argument(
        "--dir-ntupler", required=True, help="Directory that contains ntupler file."
    )
    arg_parser.add_argument(
        "--NanoAOD", required=True, help="Path to the input NanoAOD txt file."
    )
    arg_parser.add_argument(
        "--copy-root",
        action="store_true",
        default=False,
        help="Whether to copy the root files (ntupler file and NanoAODs) to the CMSSW directory.",
    )
    arg_parser.add_argument(
        "--n",
        required=True,
        type=int,
        help="Number of groups to split the input NanoAOD file into.",
    )
    arg_parser.add_argument(
        "--dir-job", required=True, help="Path to the job directory."
    )
    arg_parser.add_argument(
        "--dir-output", required=True, help="Path to the resulting root file directory."
    )
    arg_parser.add_argument(
        "--memory", type=int, default=2048, help="Memory to request for the worker nodes (in MB)"
    )
    args = arg_parser.parse_args()

    # split the NanoAOD file into n groups and save them as temp lists
    path_Nano_AOD = Path(args.NanoAOD)
    dir_job = Path(args.dir_job)
    dir_job.mkdir(parents=True, exist_ok=True)
    dir_ntuplers = Path(args.dir_ntupler)
    n_ntuplers = len(list(dir_ntuplers.glob("*.root")))

    nanoAOD_group_list = split_nanoAOD_txt_file(args.n, path_Nano_AOD)
    n_tmp_list = 0  # number of temp NanoAOD lists (< n if some are empty)
    for i, nanoAOD_group in enumerate(nanoAOD_group_list):
        if len(nanoAOD_group) == 0:
            # no nanoAOD left -> break
            break
        output_file = dir_job / f"temp_nanoAOD_{i}.txt"
        with open(output_file, "w") as f:
            for j, line in enumerate(nanoAOD_group):
                if j != len(nanoAOD_group) - 1:
                    f.write(line + "\n")
                else:
                    f.write(line)
        n_tmp_list += 1

    # write job script (runjob.sh)
    exe_name = Path(args.exe).name
    path_exe = Path(args.exe).resolve()  # absolute path
    path_convertList = path_exe.parent / "convertListMerge.py"

    # set up
    script_sh = f"""
#!/bin/bash
date
MAINDIR=`pwd`
ls
voms-proxy-info --all
#CMSSW from scratch (only need for root)
export CWD=${{PWD}}
export PATH=${{PATH}}:/cvmfs/cms.cern.ch/common
export SCRAM_ARCH=slc7_amd64_gcc700
scramv1 project CMSSW {CMSSW}
echo "Inside $MAINDIR:"
ls -lah
cp {exe_name} {CMSSW}/src
cd {CMSSW}/src
eval `scramv1 runtime -sh` # cmsenv
echo "Untar JEC:" 
echo "After Untar: "
# jobs
"""

    # job to run
    dir_ntupler = (args.dir_ntupler).replace("/eos/uscms/", "root://cmseos.fnal.gov//")
    if args.copy_root:
        # copy NanoAODs and ntupler to current directory
        script_sh += f"""
python3 ${{MAINDIR}}/convertListMerge.py -i ${{MAINDIR}}/temp_nanoAOD_${{2}}.txt
xrdcp -f {dir_ntupler}/ntupler_${{1}}.root ntupler.root
echo "ls -lha"
ls -lha
echo "Running job..."
./{exe_name} ntupler.root local_list.txt ntupler_${{1}}_nanoAOD_${{2}}.root ""
""".strip()
    else:
        # run with NanoAODs and ntupler in the original directory
        script_sh += f"""
echo "ls -lha"
ls -lha
echo "Running job..."
./{exe_name} {dir_ntupler}/ntupler_${{1}}.root ${{MAINDIR}}/temp_nanoAOD_${{2}}.txt ntupler_${{1}}_nanoAOD_${{2}}.root ""
""".strip()

    # copy root files to eos and clean up
    dir_output = (args.dir_output).replace("/eos/uscms/", "root://cmseos.fnal.gov//")

    script_sh += f"""
echo "Inside $MAINDIR:"
ls -lah
echo "coping to eos: +xrdcp -f ntupler_${{1}}_nanoAOD_${{2}}.root {dir_output}"
xrdcp -f ntupler_${{1}}_nanoAOD_${{2}}.root {dir_output}
echo "DELETING..."
rm -rf {CMSSW}
rm -rf *.pdf *.C core*
cd $MAINDIR
echo "remove output local file"
rm -rf *.root
ls
date
"""
    script_sh = script_sh.strip()

    list_path_temp_nanoAODs = [f"temp_nanoAOD_{i}.txt" for i in range(n_tmp_list)]
    # write jdl
    script_jdl = f"""
universe = vanilla
request_memory = {args.memory}
Executable = runjob.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT_OR_EVICT
Transfer_Input_Files = runjob.sh,{path_exe},{path_convertList},{",".join(list_path_temp_nanoAODs)}
n1 = {n_ntuplers}
n2 = {n_tmp_list}
N = $(n1) * $(n2)
I = ($(Process) / $(n2))
J = ($(Process) % $(n2))
Output = runjob-$(Cluster)_$(Process)-$INT(I)_$INT(J).stdout
Error  = runjob-$(Cluster)_$(Process)-$INT(I)_$INT(J).stderr
Log    = runjob-$(Cluster)_$(Process)-$INT(I)_$INT(J).log
Arguments  = "$INT(I) $INT(J)"
Queue $(N)
""".strip()

    # write runjob.sh and runjob.jdl
    with open(dir_job / "runjob.sh", "w") as f:
        f.write(script_sh)
    with open(dir_job / "runjob.jdl", "w") as f:
        f.write(script_jdl)

    print(f"Job directory created at {dir_job.resolve()}")
