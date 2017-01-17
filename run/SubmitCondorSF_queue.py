#!/bin/env python
import os
import sys
import glob
import subprocess
import time

ana_name = "rjigsawAna_TTBAR"
tar_location = "/data/uclhc/uci/user/dantrim/"

#n_split = sys.argv[1]

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0228/a_sep21/mc/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0228/a_sep21/mc/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0228/d_oct12/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0228/d_oct12/mc/logs/"

filelist_dir = "/data/uclhc/uci/user/dantrim/n0228val/filelists/"
in_job_filelist_dir = "/n0228val/filelists/"
#samples = ["diboson_sherpa_lvlv", "drellyan_sherpa", "wjets_sherpa22", "zjets_sherpa22", "ttV", "singletop"]
#samples = ["ttbar"]
#samples = ["zjets_sherpa22"]
#samples = ["data16_n0228"] #, "data16_n0228"]
#samples = ["data16_n0228"]
#samples = ["data16_remaining"]
#samples = ["bwn"]
#samples = ["missingMC"]
samples = ["ttbar"]
#samples = ["data_toRun"]

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
    #retryfile = "/data/uclhc/uci/user/dantrim/n0228val/resub.txt"
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
    look_for_condor_script(brick_ = doBrick, local_ = doLocal, sdsc_ = doSDSC, uc_ = doUC)
    look_for_condor_executable()

    for s in samples :
        print "Submtting sample : %s"%s
        suff = ""
        if not s.endswith("/") : suff = "/"
        sample_lists = glob.glob(filelist_dir + s + suff + "*.txt")
        if len(sample_lists) == 0 :
            print "No sample lists in filelist dir!"
            sys.exit()

        retry_list = get_retrylist()
        new_lists = []
        for s_ in sample_lists :
            for x in retry_list[:1] :
                if x in s_ :
                    new_lists.append(s_)
        sample_lists = new_lists


        for dataset in sample_lists :

            

            fullname = str(os.path.abspath(dataset))
            dataset_original = dataset
            print "    > %s"%dataset

            dataset = "." + dataset[dataset.find(in_job_filelist_dir):]
            print "    >> %s"%dataset

            if not (str(os.path.abspath(out_dir)) == str(os.environ['PWD'])) :
                print "You must call this script from the output directory where the ntuples will be stored!"
                print " >>> Expected submission directory : %s"%os.path.abspath(out_dir)
                sys.exit()

            do_split_sample = False
            for dsid_split in samples_to_split :
                if dsid_split in dataset :
                    do_split_sample = True

            run_mode = "-c"

            if not do_split_sample :

                # submit the job as usual
                run_cmd = "ARGS="
                run_cmd += '"'
                run_cmd += ' %s '%out_dir
                run_cmd += ' %s '%log_dir
                run_cmd += ' %s '%ana_name
                #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
                run_cmd += ' n0228val '
                run_cmd += ' %s '%dataset
                run_cmd += ' %s '%run_mode # any extra cmd line optino for Superflow executable
                run_cmd += '"'
                run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
                lname = dataset.split("/")[-1].replace(".txt", "")
                run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
                run_cmd += ' -append "output = %s%s" '%(log_dir, lname + ".out")
                run_cmd += ' -append "log = %s%s" '%(log_dir, lname + ".log")
                run_cmd += ' -append "error = %s%s" '%(log_dir, lname + ".err")

                print run_cmd
                subprocess.call(run_cmd, shell=True)


            elif do_split_sample :
                split_suffix = 0

                split_files = []
                lines = open(dataset_original).readlines()
                for line in lines :
                    if not line : continue
                    line = line.strip()
                    split_files.append(line)

                for split_file in split_files :
                    print "    >>> Sub-file [%d] %s"%(split_suffix, split_file)

                    run_cmd = "ARGS="
                    run_cmd += '"'
                    run_cmd += ' %s '%out_dir
                    run_cmd += ' %s '%log_dir
                    run_cmd += ' %s '%ana_name
                    #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
                    run_cmd += ' n0228val '
                    run_cmd += ' %s '%split_file
                    run_cmd += ' %s --sumw --suffix %d'%(run_mode, split_suffix) # any extra cmd line optino for Superflow executable
                    run_cmd += '"'
                    run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
                    lname = dataset.split("/")[-1].replace(".txt", "")
                    run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
                    run_cmd += ' -append "output = %s%s" '%(log_dir, lname + "_%d.out"%split_suffix)
                    run_cmd += ' -append "log = %s%s" '%(log_dir, lname + "_%d.log"%split_suffix)
                    run_cmd += ' -append "error = %s%s" '%(log_dir, lname + "_%d.err"%split_suffix)

                    print run_cmd
                    subprocess.call(run_cmd, shell=True)

                    split_suffix = split_suffix + 1
                    

def look_for_tarball() :
    if not os.path.isfile("area.tgz") :
        print "Tarball not found."
        sys.exit()

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
    f.write('echo "Calling : ${sflow_exec} -i ./n0228val/filelists/data_toRun/${x} ${sflow_options} 2>&1 |tee ${logname}"\n')
    #f.write('echo "Calling : ${sflow_exec} -i ${input} ${sflow_options}"\n')
    f.write('${sflow_exec} -i ./n0228val/filelists/data_toRun/${x} ${sflow_options} 2>&1 |tee ${logname}\n')
    #f.write('${sflow_exec} -i ${input} ${sflow_options}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()
    
    


if __name__=="__main__" :
    main()

