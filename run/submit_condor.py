#!/bin/env python

##############################################################################
#
# submit_condor
#
# a script to submit Superflow-based jobs to the Condor
# batch system from the UCI 'brick'
#
# daniel.joseph.antrim@cern.ch
# June 2017
#
##############################################################################


import os
import sys
import glob
import subprocess
import time
debug = False

ana_name = "rjigsawAna_WWBB" # name of the executable to be run
tar_name = "n0232val" # name of the directory you stored in the tarball
tar_location = "/data/uclhc/uci/user/dantrim/" # location of tarball file to take to job site
superflow_run_mode = "-c -n 10" # superflow run mode and options

###############################################################################
# output locations
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0232/f_jun5/data/Raw/" # location to dump the output ntuples
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0232/f_jun5/data/logs/" # location to dump the job logs
out_dir = "/data/uclhc/uci/user/dantrim/n0232val/SuperRest/run/testing/"
log_dir = "/data/uclhc/uci/user/dantrim/n0232val/SuperRest/run/testing/"

###############################################################################
# filelist locations
filelist_dir = "/data/uclhc/uci/user/dantrim/n0232val/filelists/" # directory of filelist directories
in_job_filelist_dir = "/n0232val/filelists/" # name of filelist dir as seen at the job site (in the tarred file)

###############################################################################
# samples to process - text files located in <filelist_dir>/<sample-name>/*.txt
samples = ["wwbb_test"]
#samples = ["ttbar"]
#samples = ["diboson_sherpa"]


###############################################################################
# sample listting - any samples in this list will have individual jobs per
# file in its filelist
# NB be sure you have provided the sumw-file with this DSID's sumw!
#samples_to_split = ["410009"] # ttbar
#samples_to_split = ["361073","361077","363356"]
samples_to_split = ["342053"]


###############################################################################
# sites to consider for processing
doBrick = True
doLocal = False
doSDSC  = False
doUC    = False

def fax_is_setup() :
    """
    check if the fax storage prefix is defined
    """

    fax_ok = True
    if os.environ.get('STORAGEPREFIX') == None :
        print "ERROR STORAGEPREFIX environment is not found, you must call 'lsetup fax' before calling this script"
        fax_ok = False
    return fax_ok

def current_dir_ok() :
    """
    we must call this script in the output directory
    """

    current_ok = True

    current_dir = str(os.environ['PWD'])
    expected_dir = str(os.path.abspath(out_dir))

    if not current_dir == expected_dir :
        print "ERROR current directory (%s) is not the output directory for the ntuples (%s)"%(current_dir, expected_dir)
        current_ok = False

    return current_ok


def environment_ready() :

    env_ok = True

    if not current_dir_ok() :
        env_ok = False

    if not fax_is_setup() :
        env_ok = False

    return env_ok

def build_condor_script(condor_script_name, executable_name, process_group, n_datasets) :

    global doBrick,doLocal,doSDSC,doUC

    brick = 'false'
    local = 'false'
    sdsc  = 'false'
    uc    = 'false'
    if doBrick : brick = 'true'
    if doLocal : local = 'true'
    if doSDSC  : sdsc  = 'true'
    if doUC    : uc    = 'true'

    f = open(condor_script_name, 'w')
    f.write('universe = vanilla\n')
    f.write('+local=%s\n'%brick)
    f.write('+site_local=%s\n'%local)
    f.write('+sdsc=%s\n'%sdsc)
    f.write('+uc=%s\n'%uc)
    f.write('executable = %s\n'%executable_name)
    f.write('arguments = $(Process) $ENV(ARGS)\n')
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    
    #f.write('output = /data/uclhc/uci/user/dantrim/n0232val/SuperRest/run/testing/test_BLAH_$(Process).txt\n')
    #f.write('transfer_output_files = test_job_output.txt\n')
    #f.write('transfer_output_remaps = "test_job_output.txt = /data/uclhc/uci/user/dantrim/n0232val/SuperRest/run/testing/test_BLAH.txt"\n') 
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.write('queue %d\n'%n_datasets)
    f.close()

def is_split_dataset(ds) :

    out_ds = ""
    is_split = False
    for dsid in samples_to_split :
        if dsid in ds :
            out_ds = str(dsid)
            is_split = True
            break
    return is_split, out_ds

def build_job_executable(executable_name, process_group, number_of_samples) :

    #print "build_jobs_executable process_group : %s, filelist_dir: %s"%(process_group, filelist_dir)
    group_list_check = "%s%s/"%(filelist_dir, process_group)

    datasets = glob.glob(group_list_check + "*.txt")

    n_additional = 0
    splits = []
    for ds in datasets :
        is_split, split_ds = is_split_dataset(ds)
        if is_split :
            lines = open(ds).readlines()
            for line in lines :
                if not line : continue
                line = line.strip()
                n_additional += 1
            splits.append(split_ds)

    split_ds_str = (",").join(splits)

    f = open(executable_name, 'w')
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " run condor %s "\n'%process_group)
    f.write('process_number=${1}\n')
    f.write('output_dir=${2}\n')
    f.write('log_dir=${3}\n')
    f.write('executable=${4}\n')
    f.write('stored_dir=${5}\n')
    f.write('group_name=${6}\n')
    f.write('split_dsids=${7}\n')
    f.write('superflow_options=${@:8}\n\n')
    f.write('echo "    process number        : ${process_number}"\n')
    f.write('echo "    output directory      : ${output_dir}"\n')
    f.write('echo "    log directory         : ${log_dir}"\n')
    f.write('echo "    superflow executable  : ${executable}"\n')
    f.write('echo "    tarred directory      : ${stored_dir}"\n')
    f.write('echo "    process group name    : ${group_name}"\n')
    f.write('echo "    split dsid list       : ${split_dsids}"\n')
    f.write('echo "    superflow options     : ${superflow_options}"\n')
    f.write('while (( "$#" )); do\n')
    f.write('    shift\n')
    f.write('done\n\n')
    f.write('work_dir=${PWD}\n')
    f.write('echo "untarring area.tgz"\n')
    f.write('tar -xvf area.tgz\n\n')
    f.write('echo "done untarring"\n')
    f.write('echo "current directory structure:"\n')
    f.write('ls -ltrh\n\n')
    f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
    f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
    f.write('echo "cd ${stored_dir}"\n')
    f.write('cd ${stored_dir}\n') 
    f.write('echo "Directory structure in untarred directory:"\n')
    f.write('ls -ltrh\n')
    f.write('lsetup fax\n')
    f.write('source susynt-read/bash/setup_root.sh\n')
    f.write('source RootCore/local_setup.sh\n')
    f.write('ls ./filelists/${group_name} > joblist_${group_name}.txt\n')
    #f.write('echo "Built in-job filelist for group ${group_name}:"\n')
    f.write('echo "python ./Superflow/run/get_filelist.py ${group_name} ${process_number} ${split_dsids} ${stored_dir} > injob_filelist_${group_name}_${process_number}.txt"\n')
    f.write('python ./Superflow/run/get_filelist.py ${group_name} ${process_number} ${split_dsids} ${stored_dir} > injob_filelist_${group_name}_${process_number}.txt\n')
    f.write('ls -ltrh\n')
    f.write('input_list_for_process_number=$(head -1 injob_filelist_${group_name}_${process_number}.txt)\n')
    f.write('echo "input list for process : ${input_list_for_process_number}"\n')
    #f.write('cat ${input_list_for_process_number}\n')

    f.write('echo "python ./Superflow/run/get_joblog_name.py ./filelists/${group_name} ${process_number} ${split_dsids} > injob_log_${group_name}_${process_number}.txt"\n')
    f.write('python ./Superflow/run/get_joblog_name.py ./filelists/${group_name} ${process_number} ${split_dsids} > injob_log_${group_name}_${process_number}.txt\n')
    #f.write('cat injob_log_${group_name}_${process_number}.txt\n')
    f.write('log_for_process=$(head -1 injob_log_${group_name}_${process_number}.txt)\n')
    f.write('split_cmd=""\n')
    f.write('n_lines=$(cat injob_filelist_${group_name}_${process_number}.txt |wc -l)\n')
    f.write('if [ ${n_lines} = "2" ]; then\n')
    f.write('    split_cmd=$(tail -1 injob_filelist_${group_name}_${process_number}.txt)\n')
    f.write('fi\n')
    #f.write('echo "SPLIT CMD = ${split_cmd}"\n')
    #f.write('echo "CATTING FILELIST FILE"\n')
    #f.write('cat injob_filelist_${group_name}_${process_number}.txt\n')
    f.write('echo "Setting log to: ${log_for_process}"\n')

    f.write('cd SuperRest/\n')
    f.write('source setRestFrames.sh\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Submitting with input filelist ${input_list_for_process_number}"\n')
    f.write('echo "Submitting with output log file ${log_for_process}"\n')
    f.write('echo "${executable} -i ${input_list_for_process_number} ${superflow_options} ${split_cmd} 2>&1 | tee ${log_for_process}"\n')
    f.write('${executable} -i ${input_list_for_process_number} ${superflow_options} ${split_cmd} 2>&1 | tee ${log_for_process}\n')

    #f.write('echo "HELLO WORLD" > test_job_output.txt\n')
    #f.write('dd if=/dev/zero of=test_job_output.txt bs=52428800 count=1\n')
    #f.write('dd if=/dev/zero of=output_file_${group_name}_${process_number}.txt bs=52428800 count=1\n')
    #f.write('touch output_file_${group_name}_${process_number}.root\n')
    #f.write('echo "HELLO WORLD" 2>&1 |tee test_job_output.txt\n')
    f.write('ls -ltrh\n')
    

def get_retry_list(sample_list) :
    return sample_list

def submit_sample(sample) :

    print 35*"- "
    print "Submitting sample : %s"%sample

    suffix = ""
    if not sample.endswith("/") : suffix = "/"
    sample_lists = glob.glob(filelist_dir + sample + suffix + "*.txt")
    if len(sample_lists) == 0 :
        print "WARNING No filelists found for sample %s in filelist dir %s!"%(sample, filelist_dir)
        sys.exit()

    #sample_lists = get_retry_list(sample_lists)

    number_of_samples = len(sample_lists)
    print "Number of datasets in sample %s : %d"%(sample, number_of_samples)

    process_group = sample

    condor_script_name = "submit_%s.condor"%process_group
    executable_name = "run_condor_%s.sh"%process_group

    total_number_of_jobs = 0
    for sl in sample_lists :
        is_split, split_ds = is_split_dataset(sl)
        if is_split :
            rfiles = open(sl).readlines()
            for rf in rfiles :
                if not rf : continue
                total_number_of_jobs += 1
        else :
            total_number_of_jobs += 1

    print "Total number of jobs to process in samples %s : %d"%(sample, total_number_of_jobs)

    build_condor_script(condor_script_name, executable_name, process_group, total_number_of_jobs)
    build_job_executable(executable_name, process_group, total_number_of_jobs)

    global superflow_run_mode

    split_list = (",").join(samples_to_split)
    if len(split_list) == 0 :
        split_list = "X"

    # build the job command
    run_cmd = "ARGS="
    run_cmd += '"'
    run_cmd += ' %s '%out_dir
    run_cmd += ' %s '%log_dir
    run_cmd += ' %s '%ana_name
    run_cmd += ' %s '%tar_name
    run_cmd += ' %s '%process_group
    run_cmd += ' %s '%split_list
    run_cmd += ' %s '%superflow_run_mode
    run_cmd += '"'
    run_cmd += ' condor_submit %s '%condor_script_name
    run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
    run_cmd += ' -append "output = %s%s" '%(log_dir, process_group + ".out")
    run_cmd += ' -append "log = %s%s" '%(log_dir, process_group + ".log")
    run_cmd += ' -append "error = %s%s" '%(log_dir, process_group + ".err")
    #run_cmd += ' -append "transfer_output_files = test_job_output.txt" '
    #run_cmd += ' -append "tansfer_output_remaps = '
    #run_cmd += '"test_job_output.txt = '
    #run_cmd += "/data/uclhc/uci/user/dantrim/n0232val/SuperRest/run/test_BLAH.txt" 
    #run_cmd += '""'
    #run_cmd += ' -append "transfer_output_files = test_job_output.txt" '
    #run_cmd += ' -append "transfer_output_remaps = '
    #run_cmd += '"test_job_output.txt = /data/uclhc/uci/user/dantrim/n0232val/SuperRest/run/test_BLAH.txt"' 

    global debug

    if debug :
        print run_cmd
    subprocess.call(run_cmd, shell=True)

    print "\n"

def main() :
    print 70*"-"
    print " submit_condor \n" 

    global debug

    debug = False
    if len(sys.argv) >= 2 :
        if sys.argv[1] == "-d" or sys.argv[1] == "--debug" :
            debug = True
    

    # check that we are situationally aware
    if not environment_ready() :
        sys.exit(1)

    for s in samples :
        submit_sample(s)


    

    

###############################################################################
if __name__ == "__main__" :
    main()


