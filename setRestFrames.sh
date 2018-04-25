#!/bin/bash

if [ -d "${TestArea}" ]; then

    PATH=${TestArea}/../RestFrames/bin:${PATH}; export PATH
    
    LD_LIBRARY_PATH=${TestArea}/../RestFrames/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
    
    DYLD_LIBRARY_PATH=${TestArea}/../RestFrames/lib:$DYLD_LIBRARY_PATH; export DYLD_LIBRARY_PATH
    
    SHLIB_PATH=${TestArea}/../RestFrames/lib:$SHLIB_PATH; export SHLIB_PATH
    
    LIBPATH=${TestArea}/../RestFrames/lib:$LIBPATH; export LIBPATH
else
    echo "TestArea is not defined. Please setup RootCore before calling this script."
fi
