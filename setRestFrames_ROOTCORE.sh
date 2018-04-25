#!/bin/bash

if [ -d "${ROOTCOREDIR}" ]; then

    PATH=${ROOTCOREDIR}/../RestFrames/bin:$PATH; export PATH
    
    LD_LIBRARY_PATH=${ROOTCOREDIR}/../RestFrames/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
    
    DYLD_LIBRARY_PATH=${ROOTCOREDIR}/../RestFrames/lib:$DYLD_LIBRARY_PATH; export DYLD_LIBRARY_PATH
    
    SHLIB_PATH=${ROOTCOREDIR}/../RestFrames/lib:$SHLIB_PATH; export SHLIB_PATH
    
    LIBPATH=${ROOTCOREDIR}/../RestFrames/lib:$LIBPATH; export LIBPATH
else
    echo "ROOTCOREDIR is not defined. Please setup RootCore before calling this script."
fi
