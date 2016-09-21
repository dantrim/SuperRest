#!/bin/env python
import os
import sys
import glob
import subprocess
import time

ana_name = "rjigsawAna_MT"
#ana_name = "rjigsawAna_MT"
tar_location = "/data/uclhc/uci/user/dantrim/"

#n_split = sys.argv[1]

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/forFake3/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/forFake3/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/data/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/data/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/diboson_sf/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/diboson_sf/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/ttbar/split%s/Raw/"%n_split 
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/ttbar/split%s/logs/"%n_split

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/diboson_lvlv/split%s/Raw/"%n_split
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0226/jul25/mc/diboson_lvlv/split%s/logs/"%n_split



filelist_dir = "/data/uclhc/uci/user/dantrim/n0226val/filelists/"
in_job_filelist_dir = "/n0226val/filelists/"
#samples = ['drellyan_sherpa','wjets_sherpa22','zjets_sherpa22','singletop','bwn','ttV','diboson_sherpa_noLVLV']
#samples = ['n0226_allData']
#samples = ['n0226_dataToRun']
#samples = ["diboson_sherpa_lvlv"]
samples = ["bwn"]
#samples = ["diboson_sherpa_lvlv", "bwn"]
#samples = ['drellyan_sherpa']
#samples = ['diboson_sherpa_lvlv']
#samples = ['bwn','singletop','ttV','wjets_sherpa22','zjets_sherpa22','diboson_sherpa_noLVLV']
#samples = ['ttbar_split/split%s'%n_split]
#samples = ['diboson_lvlv_split/split%s'%n_split]

doBrick = True
doLocal = False 
doSDSC  = False 
doUC    = False 


def get_retrylist() :
    runlist = ["361094", "361096"]
    retryfile = "/data/uclhc/uci/user/dantrim/n0226val/resub.txt"
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

        #new_lists = []
        ##get_these = ['280273', '284420', '284484']
        ##get_these = ['363103','363108','363113']
        ##get_these = ['410025','410026']
        ##get_these = ['363102']
        ##get_these = ['279867']#,'279932','284006']
        ##get_these = ['363366','363103','363108','363113']
        #get_these = ['301973']
        ##get_these = ['361439','361475','361482','387941','410015']
        #for s_ in sample_lists :
        #    for blah in get_these :
        #        if blah in s_ :
        #            new_lists.append(s_)
        #sample_lists = new_lists
        #retry_list = get_retrylist()
        #new_lists = []
        #for s_ in sample_lists :
        #    for x in retry_list :
        #        if x in s_ :
        #            new_lists.append(s_)
        #sample_lists = new_lists
        
            

        for dataset in sample_lists :
            fullname = str(os.path.abspath(dataset))
            print "    > %s"%dataset
           # do_this = False
           # for ds in retry_lines :
           #     ds = ds.strip()
           #     ds = ds.replace("./user", "user")
           #     ds = ds.replace(".out" , "")
           #     if ds in fullname :
           #         do_this = True
           # if not do_this :
           #     continue
            

            dataset = "." + dataset[dataset.find(in_job_filelist_dir):]
            print "    >> %s"%dataset
            

            if not (str(os.path.abspath(out_dir)) == str(os.environ['PWD'])) :
                print "You must call this script from the output directory where the ntuples will be stored!"
                print " >>> Expected submission directory : %s"%os.path.abspath(out_dir)
                sys.exit()

            run_cmd = "ARGS="
            run_cmd += '"'
            run_cmd += ' %s '%out_dir
            run_cmd += ' %s '%log_dir
            run_cmd += ' %s '%ana_name
            #run_cmd += ' %s '%(tar_location + "area3.tgz.tgz")
            run_cmd += ' n0226val '
            run_cmd += ' %s '%dataset
            run_cmd += ' -a' # any extra cmd line optino for Superflow executable
            run_cmd += '"'
            run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
            lname = dataset.split("/")[-1].replace(".txt", "")
            run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area3.tgz")
            run_cmd += ' -append "output = %s%s" '%(log_dir, lname + ".out")
            run_cmd += ' -append "log = %s%s" '%(log_dir, lname + ".log")
            run_cmd += ' -append "error = %s%s" '%(log_dir, lname + ".err")

            print run_cmd
            subprocess.call(run_cmd, shell=True)

            #time.sleep(0.5)

def look_for_tarball() :
    if not os.path.isfile("area3.tgz") :
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
    #f.write('transfer_input_files = area3.tgz.tgz\n')
    f.write('executable = RunCondorSF.sh\n')
    f.write('arguments = $ENV(ARGS)\n')
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    #f.write('transfer_output_files = OUTFILE\n')
    #f.write('transfer_output_remaps = OUTFILE_REMAP\n')
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.write('queue\n')
    f.close()

def look_for_condor_executable() :
    f = open('RunCondorSF.sh', 'w') 
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " ------- RunCondorSF -------- "\n')
    f.write('output_dir=${1}\n')
    f.write('log_dir=${2}\n')
    f.write('sflow_exec=${3}\n')
    f.write('stored_dir=${4}\n')
    f.write('input=${5}\n')
    f.write('sflow_options=${@:6}\n\n')
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
    f.write('echo "untarring area3.tgz"\n')
    f.write('tar -xzf area3.tgz\n\n')
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
    f.write('echo "Calling : cd SuperRest/"\n')
    f.write('cd SuperRest/\n')
    f.write('source setRestFrames.sh\n')
    f.write('echo "Calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Calling : ${sflow_exec} -i ${input} ${sflow_options}"\n')
    f.write('${sflow_exec} -i ${input} ${sflow_options}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()
    
    


if __name__=="__main__" :
    main()

