// std/stl
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

// ROOT
#include "TChain.h"
#include "TVectorD.h"
#include "TRandom.h"

// SusyNtuple
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/SusyNtSys.h"
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/input_options.h"

using namespace std;
using namespace sflow;

const string analysis_name = "ntupler_jigsaw";

int main(int argc, char* argv[])
{

    SFOptions options(argc, argv);
    options.ana_name = analysis_name;
    if(!read_options(options)) {
        exit(1);
    }


    return 0;
} // main
