#!/usr/bin/env python

import sys

if __name__=="__main__" :
    input_list = sys.argv[1]

    log_name = input_list.split("/")[3]
    log_name = log_name.replace(".txt","_run.log")
    print log_name
