#!/bin/env python
import os
import sys
import glob
import subprocess
import time

ana_name = "rjigsawAna_MT"
tar_location = "/data/uclhc/uci/user/dantrim/"

#n_split = sys.argv[1]

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0231/test_queue/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0231/test_queue/logs/"

filelist_dir = "/data/uclhc/uci/user/dantrim/n0231val/filelists/"
in_job_filelist_dir = "/n0231val/filelists/"
#samples = ["data16"]
#samples = ["diboson_sherpa_lvlv", "higgs"]
samples = ["higgs"]
#samples = ["diboson_sherpa_lvlv", "drellyan_sherpa", "wjets_sherpa22", "zjets_sherpa22", "ttV", "singletop"]

samples_to_split = ["410009"]

doBrick = False
doLocal = True
doSDSC  = True 
doUC    = True 


def get_retrylist() :
    #runlist = ["279867", "280464"]
    #runlist = ["363389","363373","363369","363111","363104"]
    #runlist = ["306269"]
    #runlist = ["307514"]
    runlist = []
    lines = open("retrylist.txt").readlines()
    for line in lines :
        if not line : continue
        line = line.strip()
        runlist.append(str(line))
    #runlist = ["279867","280464"]
    #retryfile = "/data/uclhc/uci/user/dantrim/n0231val/resub.txt"
    #lines = open(retryfile).readlines()
    #for line in lines :
    #    if not line : continue
    #    line = line.strip()
    #    runlist.append(line)
    return runlist

def main() :
    print "SubmitCondorSF"

    ### retry [begin]
    #retry_lines = open(retry_list).readlines()
    

#    look_for_tarball()
    #look_for_condor_script(brick_ = doBrick, local_ = doLocal, sdsc_ = doSDSC, uc_ = doUC)
    #look_for_condor_executable()

    for s in samples[:1] :
        print "Submtting sample : %s"%s
        suff = ""
        if not s.endswith("/") : suff = "/"
        sample_lists = glob.glob(filelist_dir + s + suff + "*.txt")
        if len(sample_lists) == 0 :
            print "No sample lists in filelist dir!"
            sys.exit()

        #retry_list = get_retrylist()
        #new_lists = []
        #for s_ in sample_lists :
        #    for x in retry_list[:1] :
        #        if x in s_ :
        #            new_lists.append(s_)
        #sample_lists = new_lists


        number_of_samples = len(sample_lists)
        print "number of samples: ",number_of_samples
        print s
        process_group= s


        script_name = "submitFile_%s.condor"%process_group
        executable_name = "RunCondorSF_%s.sh"%process_group
        build_condor_script(script_name, executable_name, process_group, number_of_samples, doBrick, doLocal, doSDSC, doUC)
        build_job_executable(executable_name, process_group, number_of_samples)


        # Superflow run mode
        run_mode = "-c"

        # build the job command
        run_cmd = "ARGS="
        run_cmd += '"'
        run_cmd += ' %s '%out_dir
        run_cmd += ' %s '%log_dir
        run_cmd += ' %s '%ana_name
        run_cmd += ' n0231val '
        run_cmd += ' %s '%process_group
        run_cmd += ' %s '%run_mode # put here any extra cmd line options for Superflow executable
        run_cmd += '"'
        run_cmd += ' condor_submit %s '%script_name
        run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")

        run_cmd += ' -append "output = %s%s" '%(log_dir, process_group + ".out")
        run_cmd += ' -append "log = %s%s" '%(log_dir, process_group + ".log")
        run_cmd += ' -append "error = %s%s" '%(log_dir, process_group + ".err")

        print run_cmd
        subprocess.call(run_cmd, shell=True)

        


#        for dataset in sample_lists :
#            print dataset
#            ds_injob = in_job_filelist_dir + dataset.split(in_job_filelist_dir)[1]
#            print ds_injob
#            sys.exit()

            

#            fullname = str(os.path.abspath(dataset))
#            dataset_original = dataset
#            print "    > %s"%dataset
#
#            dataset = "." + dataset[dataset.find(in_job_filelist_dir):]
#            print "    >> %s"%dataset
#
#            if not (str(os.path.abspath(out_dir)) == str(os.environ['PWD'])) :
#                print "You must call this script from the output directory where the ntuples will be stored!"
#                print " >>> Expected submission directory : %s"%os.path.abspath(out_dir)
#                sys.exit()
#
#            do_split_sample = False
#            for dsid_split in samples_to_split :
#                if dsid_split in dataset :
#                    do_split_sample = True
#
#            run_mode = "-c"

#            if not do_split_sample :
#
#                # submit the job as usual
#                run_cmd = "ARGS="
#                run_cmd += '"'
#                run_cmd += ' %s '%out_dir
#                run_cmd += ' %s '%log_dir
#                run_cmd += ' %s '%ana_name
#                #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
#                run_cmd += ' n0231val '
#                run_cmd += ' %s '%dataset
#                run_cmd += ' %s '%run_mode # any extra cmd line optino for Superflow executable
#                run_cmd += '"'
#                run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
#                lname = dataset.split("/")[-1].replace(".txt", "")
#                run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
#                run_cmd += ' -append "output = %s%s" '%(log_dir, lname + ".out")
#                run_cmd += ' -append "log = %s%s" '%(log_dir, lname + ".log")
#                run_cmd += ' -append "error = %s%s" '%(log_dir, lname + ".err")
#
#                print run_cmd
#                subprocess.call(run_cmd, shell=True)
#
#
#            elif do_split_sample :
#                split_suffix = 0
#
#                split_files = []
#                lines = open(dataset_original).readlines()
#                for line in lines :
#                    if not line : continue
#                    line = line.strip()
#                    split_files.append(line)
#
#                for split_file in split_files :
#                    print "    >>> Sub-file [%d] %s"%(split_suffix, split_file)
#
#                    run_cmd = "ARGS="
#                    run_cmd += '"'
#                    run_cmd += ' %s '%out_dir
#                    run_cmd += ' %s '%log_dir
#                    run_cmd += ' %s '%ana_name
#                    #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
#                    run_cmd += ' n0231val '
#                    run_cmd += ' %s '%split_file
#                    run_cmd += ' %s --sumw --suffix %d'%(run_mode, split_suffix) # any extra cmd line optino for Superflow executable
#                    run_cmd += '"'
#                    run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
#                    lname = dataset.split("/")[-1].replace(".txt", "")
#                    run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
#                    run_cmd += ' -append "output = %s%s" '%(log_dir, lname + "_%d.out"%split_suffix)
#                    run_cmd += ' -append "log = %s%s" '%(log_dir, lname + "_%d.log"%split_suffix)
#                    run_cmd += ' -append "error = %s%s" '%(log_dir, lname + "_%d.err"%split_suffix)
#
#                    print run_cmd
#                    subprocess.call(run_cmd, shell=True)
#
#                    split_suffix = split_suffix + 1
                    

def look_for_tarball() :
    if not os.path.isfile("area.tgz") :
        print "Tarball not found."
        sys.exit()

def build_condor_script(script_name, executable_name, proc_group_name, n_samples_in_group, doBrick, doLocal, doSDSC, doUC) :

    brick = 'false'
    local = 'false'
    sdsc  = 'false'
    uc    = 'false'
    if doBrick : brick = 'true'
    if doLocal : local = 'true'
    if doSDSC  : sdsc  = 'true'
    if doUC    : uc    = 'true'

    f = open(script_name, "w")
    f.write('universe = vanilla\n')
    f.write('+local=%s\n'%brick)
    f.write('+site_local=%s\n'%local)
    f.write('+sdsc=%s\n'%sdsc)
    f.write('+uc=%s\n'%uc)
    f.write('executable = %s\n'%executable_name)
    f.write('arguments = $(Process) $ENV(ARGS)\n')
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.write('queue %d\n'%n_samples_in_group)
    f.close()


def look_for_condor_script(brick_ = False, local_ = False, sdsc_ = False, uc_ = False) :

    brick = 'false'
    local = 'false'
    sdsc  = 'false'
    uc    = 'false'
    if brick_ : brick = 'true'
    if local_ : local = 'true'
    if sdsc_  : sdsc = 'true'
    if uc_    : uc = 'true'

    f = open('submitFile_TEMPLATE.condor', 'w')
    f.write('universe = vanilla\n')
    f.write('+local=%s\n'%brick_)
    f.write('+site_local=%s\n'%local_)
    f.write('+sdsc=%s\n'%sdsc_)
    f.write('+uc=%s\n'%uc_)
    #f.write('transfer_input_files = area.tgz.tgz\n')
    f.write('executable = RunCondorSF.sh\n')
    f.write('arguments = $(Process) $ENV(ARGS)\n')
    #f.write('arguments = $ENV(ARGS)\n')
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    #f.write('transfer_output_files = OUTFILE\n')
    #f.write('transfer_output_remaps = OUTFILE_REMAP\n')
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.write('queue 178\n')
    f.close()

def build_job_executable(executable_name, process_group, n_samples_in_group) :
    f = open(executable_name, 'w')
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " -------- RunCondorSF %s -------- "\n'%process_group)
    f.write('processs_no=${1}\n')
    f.write('output_dir=${2}\n')
    f.write('log_dir=${3}\n')
    f.write('sflow_exec=${4}\n')
    f.write('stored_dir=${5}\n')
    f.write('group_name=${6}\n')
    f.write('sflow_options=${@:7}\n\n')
    f.write('echo "    output directory    : ${output_dir}"\n')
    f.write('echo "    log directory       : ${log_dir}"\n')
    f.write('echo "    sflow executable    : ${sflow_exec}"\n')
    f.write('echo "    tarred dir          : ${stored_dir}"\n')
    f.write('echo "    process group       : ${group_name}"\n')
    f.write('echo "    sflow options       : ${sflow_options}"\n\n')
    f.write('while (( "$#" )); do\n')
    f.write('    shift\n')
    f.write('done\n\n')
    f.write('work_dir=${PWD}\n')
    f.write('echo "untarring area.tgz"\n')
    f.write('tar -xzf area.tgz\n\n')
    f.write('echo "done untarring"\n')
    f.write('echo " > current directory structure:"\n')
    f.write('ls -ltrh\n\n')
    f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
    f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
    f.write('echo "Calling : cd ${stored_dir}"\n')
    f.write('cd ${stored_dir}\n')
    f.write('echo "Directory structure:"\n')
    f.write('ls -ltrh\n')
    f.write('lsetup fax\n')
    f.write('source susynt-read/bash/setup_root.sh\n')
    f.write('echo "Calling : source RootCore/local_setup.sh"\n')
    f.write('source RootCore/local_setup.sh\n')
    f.write('ls ./filelists/${group_name} > joblist_${group_name}.txt\n')
    f.write('echo "Built in-job filelist for group ${group_name}:"\n')
    f.write('cat joblist_${group_name}.txt\n')
    f.write('input_list_for_process="./SuperRest/run/GetFileList.py ${group_name} ${process_no}"\n')
    f.write('input_list_for_process="./${stored_dir}$input_list_for_process"\n')
    f.write('echo "Found input for group: ${group_name} and process no: ${process_no}"\n')
    f.write('echo "  > ${input_list_for_process}"\n\n')
    f.write('log_for_process="./SuperRest/run/GetJobLogName.py ${input_list_for_process}"\n')
    f.write('echo "  > Run log: ${log_for_process}"\n')
    f.write('echo "Calling : cd SuperRest/"\n')
    f.write('cd SuperRest/\n')
    f.write('source setRestFrames.sh\n')
    f.write('echo "Calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Calling : ${sflow_exec} -i ${input_list_for_process} ${sflow_options} 2>&1 |tee ${log_for_process}"\n')
    f.write('${sflow_exec} -i ${input_list_for_process} ${sflow_options} 2>&1 |tee ${log_for_process}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()
    


def look_for_condor_executable() :
    f = open('RunCondorSF.sh', 'w') 
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " ------- RunCondorSF -------- "\n')
    f.write('process_no=${1}\n')
    f.write('output_dir=${2}\n')
    f.write('log_dir=${3}\n')
    f.write('sflow_exec=${4}\n')
    f.write('stored_dir=${5}\n')
    f.write('input=${6}\n')
    f.write('sflow_options=${@:7}\n\n')
    f.write('echo "    output directory   : ${output_dir}"\n')
    f.write('echo "    log directory      : ${log_dir}"\n')
    f.write('echo "    sflow executable   : ${sflow_exec}"\n')
    f.write('echo "    tarred dir         : ${stored_dir}"\n')
    f.write('echo "    sample list        : ${input}"\n')
    f.write('echo "    sflow options      : ${sflow_options}"\n\n')
    f.write('while (( "$#" )); do\n')
    f.write('    shift\n')
    f.write('done\n\n')
    f.write('work_dir=${PWD}\n')
    f.write('echo "untarring area.tgz"\n')
    f.write('tar -xzf area.tgz\n\n')
    f.write('echo "done untarring\n')
    f.write('echo "current directory structure:\n')
    f.write('ls -ltrh\n\n')
    f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
    f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
    f.write('echo "Calling : cd ${stored_dir}"\n')
    f.write('cd ${stored_dir}\n')
    f.write('echo "Directory structure:"\n')
    f.write('ls -ltrh\n')
    f.write('lsetup fax\n')
    f.write('source susynt-read/bash/setup_root.sh\n')
    f.write('echo "Calling : source RootCore/local_setup.sh"\n')
    f.write('source RootCore/local_setup.sh\n')
    f.write('cmd="./get_sample.py ${process_no}"\n')
    f.write('x="$($cmd)"\n')
    f.write('logcmd="./get_log.py ${x}"\n')
    f.write('logname="$($logcmd)"\n')
    f.write('echo "Calling : cd SuperRest/"\n')
    f.write('cd SuperRest/\n')
    f.write('source setRestFrames.sh\n')
    f.write('echo "Calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Running over sample: ${x}"\n')
    f.write('echo "Calling : ${sflow_exec} -i ./n0231val/filelists/data_toRun/${x} ${sflow_options} 2>&1 |tee ${logname}"\n')
    #f.write('echo "Calling : ${sflow_exec} -i ${input} ${sflow_options}"\n')
    f.write('${sflow_exec} -i ./n0231val/filelists/data_toRun/${x} ${sflow_options} 2>&1 |tee ${logname}\n')
    #f.write('${sflow_exec} -i ${input} ${sflow_options}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()
    
    


if __name__=="__main__" :
    main()

