// rjigsawAna.cxx


// std
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
#include "SusyNtuple/TriggerTools.h"
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"

// RestFrames
#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace sflow;
using namespace RestFrames;

const string analysis_name = "rjigsawAna_MT";

////////////////////////////////////////
// Function Prototypes
////////////////////////////////////////
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, string& suffix_name_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sys_, bool& split_sumw, bool dbg);


////////////////////////////////////////
// MAIN
////////////////////////////////////////
int main(int argc, char* argv[])
{
    ///////////////////////
    // READ IN
    ///////////////////////
    int n_skip_ = 0;
    int num_events_ = -1;
    bool m_dbg = false;
    string sample_;
    string suffix_name_ = "";
    bool do_sumw_split = false;
    SuperflowRunMode run_mode_ = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, suffix_name_, run_mode_, nt_sys_, do_sumw_split, m_dbg);

    ////////////////////////////////////////////////////
    // Construct and configure the Superflow object
    ////////////////////////////////////////////////////
    Superflow* cutflow = new Superflow();
    cutflow->setAnaName(analysis_name);
    cutflow->setAnaType(AnalysisType::Ana_Stop2L);
    cutflow->setLumi(1000); // 1/fb
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode_);
    cutflow->setCountWeights(true);
    cutflow->setChain(chain);
    if(suffix_name_!="")
        cutflow->setFileSuffix(suffix_name_);
    if(do_sumw_split) {
        string sumw_file = "./n0231val/sumw_file.txt";
        //string sumw_file = "/data/uclhc/uci/user/dantrim/n0229val/sumw_file.txt"; 
        cout << analysis_name << "    Reading sumw for sample from file: " << sumw_file << endl; 
        cutflow->setUseSumwFile(sumw_file);
    }

    // initialize trigger
    cutflow->nttools().initTriggerTool(ChainHelper::firstFile(sample_, m_dbg));

    // print some useful
    cout << analysis_name << "    Total Entries    : " << chain->GetEntries() << endl;
    if(num_events_ > 0) {
    cout << analysis_name << "    Process Entries  : " << num_events_ << endl;
    } else {
    cout << analysis_name << "    Process Entries  : " << chain->GetEntries() << endl;
    }

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    //
    // Superflow methods [BEGIN]
    //
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    *cutflow << CutName("read in ") << [](Superlink* sl) -> bool { return true; };

    ////////////////////////////////////////////////////
    // Cleaning cuts
    ////////////////////////////////////////////////////
    int cutflags = 0;
    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];
        return (sl->tools->passGRL(cutflags));
    };
    *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passLarErr(cutflags));
    };
    *cutflow << CutName("Tile Error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTileErr(cutflags));
    };
    *cutflow << CutName("SCT error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passSCTErr(cutflags));
    };
    *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTTC(cutflags));
    };
    *cutflow << CutName("pass Good Vertex") << [&](Superlink * sl) -> bool {
        return (sl->tools->passGoodVtx(cutflags));
    };
    *cutflow << CutName("pass bad muon veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passBadMuon(sl->preMuons));
    };
    *cutflow << CutName("pass cosmic muon veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passCosmicMuon(sl->baseMuons));
    };
    *cutflow << CutName("pass jet cleaning") << [&](Superlink* sl) -> bool {
        return (sl->tools->passJetCleaning(sl->baseJets));
    };

//    *cutflow << CutName("pass met cleaning") << [&](Superlink* sl) -> bool {
//        return (sl->tools->passMetCleaning(sl->met));
//    };

    ///////////////////////////////////////////////////
    // Analysis Cuts
    ///////////////////////////////////////////////////
    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        return sl->leptons->size() == 2;
    };

    *cutflow << CutName("lepton pTs > 20 GeV") << [](Superlink* sl) -> bool {
        return ( (sl->leptons->at(0)->Pt()>20) && (sl->leptons->at(1)->Pt()>20) );
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
    };

    *cutflow << CutName("mll > 20 GeV") << [](Superlink* sl) -> bool {
        return ( (*sl->leptons->at(0) + *sl->leptons->at(1)).M() > 20. );
    };

    //*cutflow << CutName("veto SF Z-window (within 10 GeV)") << [](Superlink* sl) -> bool {
    //    bool pass = true;
    //    bool isSF = false;
    //    if((sl->leptons->size()==2 && (sl->electrons->size()==2 || sl->muons->size()==2))) isSF = true;
    //    if(isSF) {
    //        double mll = (*sl->leptons->at(0) + *sl->leptons->at(1)).M();
    //        if( fabs(mll-91.2) < 10. ) pass = false;
    //    }
    //    return pass;
    //};

    *cutflow << CutName("pass trigger requirement") << [&](Superlink* sl) -> bool {
        int year = sl->nt->evt()->treatAsYear;
        bool pass = false;

        bool passes_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
        bool passes_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
        bool passes_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
        bool passes_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
        bool passes_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
        bool passes_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14");

        bool trig_pass_2015 = (passes_mu18_mu8noL1 || passes_e17_lhloose_mu14 || passes_2e12_lhloose_L12EM10VH);
        bool trig_pass_2016 = (passes_2e17_lhvloose_nod0 || passes_mu22_mu8noL1 || passes_e17_lhloose_nod0_mu14);
        if( (year==2015 && trig_pass_2015==1) || (year==2016 && trig_pass_2016==1) ) {
            pass = true;
        }
        return pass;
    };
    
    //*cutflow << CutName("isSF") << [&](Superlink* sl) -> bool {
    //    bool isSF = false;
    //    bool zveto = false;
    //    if((sl->leptons->size()==2 && (sl->electrons->size()==2 || sl->muons->size()==2))) isSF = true;
    //    if(isSF) {
    //        double mll = (*sl->leptons->at(0) + *sl->leptons->at(1)).M();
    //        if(fabs(mll-91.2)>10.) zveto = true;
    //    }
    //    if(isSF && zveto) return true;
    //    else { return false; }
    //};
    //*cutflow << CutName("isDF") << [&](Superlink* sl) -> bool {
    //    bool isDF = false;
    //    if((sl->leptons->size()==2 && sl->electrons->size()==1 && sl->muons->size()==1)) isDF = true;
    //    if(isDF) return true;
    //    else { return false; }
    //};

    ///////////////////////////////////////////////////
    // Ntuple Setup
    ///////////////////////////////////////////////////


    // TRIGGERS
    bool pass_mu18_mu8noL1;
    bool pass_mu20_mu8noL1;
    bool pass_e17_lhloose_mu14;
    bool pass_2e12_lhloose_L12EM10VH;
    bool pass_2e15_lhvloose_L12EM13VH;
    // updated
    bool pass_2e17_lhvloose_nod0;
    bool pass_mu22_mu8noL1;
    bool pass_e17_lhloose_nod0_mu14;
    *cutflow << [&](Superlink* sl, var_void*) {
        pass_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
        pass_mu20_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_mu8noL1");
        pass_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
        pass_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
        pass_2e15_lhvloose_L12EM13VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e15_lhvloose_L12EM13VH");
        pass_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
        pass_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
        pass_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14"); 
    };
    *cutflow << NewVar("pass mu18_mu8noL1"); {
        *cutflow << HFTname("trig_mu18_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_mu18_mu8noL1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass mu20_mu8noL1"); {
        *cutflow << HFTname("trig_mu20_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_mu20_mu8noL1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e17_lhloose_mu14"); {
        *cutflow << HFTname("trig_e17_lhloose_mu14");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_e17_lhloose_mu14;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e12_lhloose_L12EM10VH"); {
        *cutflow << HFTname("trig_2e12_lhloose_L12EM10VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_2e12_lhloose_L12EM10VH;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e15_lhvloose_L12EM13VH"); {
        *cutflow << HFTname("trig_2e15_lhvloose_L12EM13VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_2e15_lhvloose_L12EM13VH;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass OR 2015 trigger set"); {
        *cutflow << HFTname("trig_pass2015");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return (pass_mu18_mu8noL1 || pass_e17_lhloose_mu14 || pass_2e12_lhloose_L12EM10VH);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass OR 2016 trigger set"); {
        *cutflow << HFTname("trig_pass2016");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return (pass_mu20_mu8noL1 || pass_e17_lhloose_mu14 || pass_2e15_lhvloose_L12EM13VH);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass HLT_2e17_lhvloose_nod0"); {
        *cutflow << HFTname("trig_2e17_lhvloose_nod0");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_2e17_lhvloose_nod0;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("pass mu22_mu8noL1"); {
        *cutflow << HFTname("trig_mu22_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_mu22_mu8noL1;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("pass HLT_e17_lhloose_nod0_mu14"); {
        *cutflow << HFTname("trig_e17_lhloose_nod0_mu14");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_e17_lhloose_nod0_mu14;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("pass OR UPDATED trigger set"); {
        *cutflow << HFTname("trig_pass2016update");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return (pass_2e17_lhvloose_nod0 || pass_mu22_mu8noL1 || pass_e17_lhloose_nod0_mu14); 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mcid"); {
        *cutflow << HFTname("mcid");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sl->nt->evt()->mcChannel;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("year (either 2015 or 2016 for data)"); {
        *cutflow << HFTname("year");
        *cutflow << [&](Superlink* sl, var_int*) -> int { 
            return sl->nt->evt()->treatAsYear;
        };
        *cutflow << SaveVar();
    }

    // standard variables
    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product() * sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("event weight with Sherpa V+Jets weight"); {
        *cutflow << HFTname("eventweightVJets");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            double weight = sl->weights->product() * sl->nt->evt()->wPileup;
            if(sl->nt->evt()->isSherpaVjetsSample)
                weight *= sl->nt->evt()->sherpa22VjetsWeight;
            return weight;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("event weight (no pileup)"); {
        *cutflow << HFTname("eventweightNOPUPW");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product();
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("Pile-up weight"); {
        *cutflow << HFTname("pupw");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pile-up weight (up variation)"); {
        *cutflow << HFTname("pupw_up");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup_up;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pile-up weight (down variation"); {
        *cutflow << HFTname("pupw_down");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup_dn;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Three body re-weighting -- Left polarized stop"); {
        *cutflow << HFTname("susy3BodyLeftPol");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->susy3BodyLeftPol;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Three body re-weighting -- Right polarized (70%) stop"); {
        *cutflow << HFTname("susy3BodyRightPol");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->susy3BodyRightPol;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Three body re-weighting -- Matrix element only"); {
        *cutflow << HFTname("susy3BodyOnlyMass");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->susy3BodyOnlyMass;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is MC"); {
        *cutflow << HFTname("isMC");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of primary vertices"); {
        *cutflow << HFTname("nVtx");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->nVtx; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("average interactions per b.c."); {
        *cutflow << HFTname("avgMu");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->avgMu; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("average interactions per b.c. with data scale factor applied"); {
        *cutflow << HFTname("avgMuDataSF");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->avgMuDataSF; };
        *cutflow << SaveVar();
    } 

    // lepton variables
    // lepton variables
    // lepton variables

    LeptonVector leptons;
    ElectronVector electrons;
    MuonVector muons;
    *cutflow << [&](Superlink* sl, var_void*) { leptons = *sl->leptons; };
    *cutflow << [&](Superlink* sl, var_void*) { electrons = *sl->electrons; };
    *cutflow << [&](Superlink* sl, var_void*) { muons = *sl->muons; };

    *cutflow << NewVar("number of leptons"); {
        *cutflow << HFTname("nLeptons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return leptons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of electrons"); {
        *cutflow << HFTname("nElectrons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return electrons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of muons"); {
        *cutflow << HFTname("nMuons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return muons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton charge"); {
        *cutflow << HFTname("l_q");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->q);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0"); {
        *cutflow << HFTname("l_d0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->d0);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0sig"); {
        *cutflow << HFTname("l_d0sig");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->d0sigBSCorr);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton z0sinTheta"); {
        *cutflow << HFTname("l_z0sinTheta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->z0SinTheta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton pt"); {
        *cutflow << HFTname("l_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->Pt());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("etconetopo20"); {
        *cutflow << HFTname("l_etconetopo20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->etconetopo20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("etconetopo30"); {
        *cutflow << HFTname("l_etconetopo30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->etconetopo30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptcone20"); {
        *cutflow << HFTname("l_ptcone20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptcone20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptcone30"); {
        *cutflow << HFTname("l_ptcone30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptcone30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptvarcone20"); {
        *cutflow << HFTname("l_ptvarcone20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptvarcone20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptvarcone30"); {
        *cutflow << HFTname("l_ptvarcone30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptvarcone30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton eta"); {
        *cutflow << HFTname("l_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->Eta());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton phi"); {
        *cutflow << HFTname("l_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->Phi());
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mll leptons"); {
        *cutflow << HFTname("mll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mll = -999.0;
            if(leptons.size() == 2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                mll = (*l0 + *l1).M();
            }
            return mll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dilepton pT"); {
        *cutflow << HFTname("pTll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double pTll = -999.0;
            if(leptons.size() == 2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                pTll = (*l0 + *l1).Pt();
            }
            return pTll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between to leptons"); {
        *cutflow << HFTname("dphi_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double dphi = -999.0;
            if(leptons.size() == 2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                dphi = l0.DeltaPhi(l1);
            }
            return dphi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta eta between two leptons"); {
        *cutflow << HFTname("deta_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double deta = -999.0;
            if(leptons.size() == 2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                deta = l0.Eta() - l1.Eta();
            }
            return deta;
        };
        *cutflow << SaveVar();
    }

    // jet variables
    // jet variables
    // jet variables

    JetVector jets;
    JetVector bjets;
    JetVector sjets;
    *cutflow << [&](Superlink* sl, var_void*) { jets = *sl->jets; };
    *cutflow << [&](Superlink* sl, var_void*) {
        for(int i = 0; i < jets.size(); i++) {
            Jet* j = jets[i];
            if(sl->tools->jetSelector().isB(j))  bjets.push_back(j);
            if(!sl->tools->jetSelector().isB(j)) sjets.push_back(j);
        }// i
    };
    *cutflow << NewVar("number of jets"); {
        *cutflow << HFTname("nJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return jets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets"); {
        *cutflow << HFTname("nSJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of bjets"); {
        *cutflow << HFTname("nBJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet pt"); {
        *cutflow << HFTname("j_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet pt"); {
        *cutflow << HFTname("sj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt"); {
        *cutflow << HFTname("bj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet eta"); {
        *cutflow << HFTname("j_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet eta"); {
        *cutflow << HFTname("sj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet eta"); {
        *cutflow << HFTname("bj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet phi"); {
        *cutflow << HFTname("j_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet phi"); {
        *cutflow << HFTname("sj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet phi"); {
        *cutflow << HFTname("bj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between dilepton system and leading jet"); {
        *cutflow << HFTname("dphi_j0_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(jets.size()>0) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = jets.at(0)->DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading sjet"); {
        *cutflow << HFTname("dphi_j0_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(jets.size()>0) {
                out = jets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading sjet"); {
        *cutflow << HFTname("dphi_sj0_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(sjets.size()>0) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = sjets.at(0)->DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading sjet"); {
        *cutflow << HFTname("dphi_sj0_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(sjets.size()>0) {
                out = sjets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading bjet"); {
        *cutflow << HFTname("dphi_bj0_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(bjets.size()>0) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = bjets.at(0)->DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading bjet"); {
        *cutflow << HFTname("dphi_bj0_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(bjets.size()>0) {
                out = bjets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    // met variables
    // met variables
    // met variables
    Met met;
    *cutflow << [&](Superlink* sl, var_void*) { met = *sl->met; };
    *cutflow << NewVar("transverse missing energy (Etmiss)"); {
        *cutflow << HFTname("met");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Pt(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("phi coord. of Etmiss"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Phi(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met TST"); {
        *cutflow << HFTname("metTST");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.softTerm_et; };
        *cutflow << SaveVar();
    }
//    *cutflow << NewVar("delta phi TST and MET"); {
//        *cutflow << HFTname("dphi_met_tst");
//        *cutflow << [&](Superlink* sl, var_float*) -> double { return fabs(TVector2::Phi_mpi_pi(met.softTerm_phi - met.lv().Phi())); };
//        *cutflow << SaveVar();
//    }

    *cutflow << NewVar("delta phi between dilepton system and met"); {
        *cutflow << HFTname("dphi_met_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector l0, l1, ll;
            l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
            ll = l0 + l1;
            return met.lv().DeltaPhi(ll);
        };
        *cutflow << SaveVar();
    }
//    *cutflow << NewVar("delta phi between lead lepton and met"); {
//        *cutflow << HFTname("dphi_met_l0");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            TLorentzVector l0;
//            l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
//            return met.lv().DeltaPhi(l0);
//        };
//        *cutflow << SaveVar();
//    }

    *cutflow << NewVar("mt2"); {
        *cutflow << HFTname("mt2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mt2 = -999.0;
            if(leptons.size() == 2) {
                mt2 = kin::getMT2(*sl->leptons, *sl->met);
            }
            return mt2;
        };
        *cutflow << SaveVar();
    }
    double meff;
    *cutflow << NewVar("meff : scalar sum pt of all jets, leptons, and met"); {
        *cutflow << HFTname("meff");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            meff = 0.0;
            // met
            meff += met.lv().Pt();
            // jets
            for(unsigned int ij = 0; ij < jets.size(); ij++){
                meff += jets.at(ij)->Pt();
            }
            // leptons
            for(unsigned int il=0; il < leptons.size(); il++){
                meff += leptons.at(il)->Pt();
            }
            return meff;
        };
        *cutflow << SaveVar();
    }
    double meff_S2L;
    *cutflow << NewVar("meff S2L : scalar sum pt of leptons, met, and up to two jets"); {
        *cutflow << HFTname("meff_S2L");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            meff_S2L = 0.0;
            // met
            meff_S2L += met.lv().Pt();
            // leptons
            for(int il=0; il < (int)leptons.size(); il++){
                meff_S2L += leptons.at(il)->Pt();
            }
            // jets
            int n_j = 0;
            for(int ij = 0; ij < (int)jets.size(); ij++){
                if(n_j < 2) {
                    meff_S2L += jets.at(ij)->Pt();
                    n_j++;
                }
            }
            return meff_S2L;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R1 : met / meff"); {
        *cutflow << HFTname("R1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R1 = -999.0;
            if(meff>0.0) {
                R1 = met.lv().Pt() / meff * 1.0;
            }
            return R1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R1 S2L : met / meff_S2L"); {
        *cutflow << HFTname("R1_S2L");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R1_S2L = -999.0;
            if(meff_S2L>0.0) {
                R1_S2L = met.lv().Pt() / (meff_S2L * 1.0);
            }
            return R1_S2L;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R2 : met / (met + l0pt + l1pt)"); {
        *cutflow << HFTname("R2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R2 = -999.0;
            if(leptons.size() == 2) {
                double denom = met.lv().Pt() + leptons.at(0)->Pt() + leptons.at(1)->Pt();
                R2 = met.lv().Pt() / denom * 1.0;
            }
            return R2;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("cosThetaB (WW-like)"); {
        *cutflow << HFTname("cosThetaB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double cosThetaB = -999;
            if(leptons.size()==2) {
                TLorentzVector lp, lm;
                for(int il = 0; il < (int)leptons.size(); il++) {
                    if(leptons.at(il)->q < 0) lm = *leptons.at(il);
                    else if(leptons.at(il)->q > 0) lp = *leptons.at(il);
                } // il
                TLorentzVector ll = lp+lm;
                TVector3 boost = ll.BoostVector();
                lp.Boost(-boost);
                lm.Boost(-boost);
                cosThetaB = tanh((lp.Eta()-lm.Eta())/2.);
            }
            return cosThetaB;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton centrality (OS-VBF) - sjets"); {
        *cutflow << HFTname("lep_centrality_sjet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(leptons.size()==2 && sjets.size()>=2) {
                float l_eta1 = leptons.at(0)->Eta();
                float l_eta2 = leptons.at(1)->Eta();
                float j_eta1 = sjets.at(0)->Eta();
                float j_eta2 = sjets.at(1)->Eta();

                float min_l_etas = std::min(l_eta1, l_eta2);
                float min_j_etas = std::min(j_eta1, j_eta2);
                float max_l_etas = std::max(l_eta1, l_eta2);
                float max_j_etas = std::max(j_eta1, j_eta2);

                out = std::min( (min_l_etas - min_j_etas), (max_j_etas - max_l_etas) );
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton centrality (OS-VBF) - jets"); {
        *cutflow << HFTname("lep_centrality_jet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(leptons.size()==2 && jets.size()>=2) {
                float l_eta1 = leptons.at(0)->Eta();
                float l_eta2 = leptons.at(1)->Eta();
                float j_eta1 = jets.at(0)->Eta();
                float j_eta2 = jets.at(1)->Eta();

                float min_l_etas = std::min(l_eta1, l_eta2);
                float min_j_etas = std::min(j_eta1, j_eta2);
                float max_l_etas = std::max(l_eta1, l_eta2);
                float max_j_etas = std::max(j_eta1, j_eta2);

                out = std::min( (min_l_etas - min_j_etas), (max_j_etas - max_l_etas) );
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("f_recoil_sjets"); {
        *cutflow << HFTname("f_recoil_sjet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            float numerator = 0.0;
            TVector2 f_num;
            for(auto& jet : sjets) {
                if(jet->Pt()>10. && std::abs(jet->Eta())<4.5 && jet->jvt!=-1) {
                    float jvt = std::abs(jet->jvt); 
                    TVector2 jet_pt(jvt * jet->Px(), jvt * jet->Py());
                    f_num += jet_pt;
                }
            } // jet
            numerator = sqrt(f_num.Px()*f_num.Px() + f_num.Py()*f_num.Py());
            float denom = (*leptons.at(0) + *leptons.at(1)).Pt();
            return (numerator/denom*1.0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("f_recoil_jets"); {
        *cutflow << HFTname("f_recoil_jet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            float numerator = 0.0;
            TVector2 f_num;
            for(auto& jet : jets) {
                if(jet->Pt()>10. && std::abs(jet->Eta())<4.5 && jet->jvt!=-1) {
                    float jvt = std::abs(jet->jvt); 
                    TVector2 jet_pt(jvt * jet->Px(), jvt * jet->Py());
                    f_num += jet_pt;
                }
            } // jet
            numerator = sqrt(f_num.Px()*f_num.Px() + f_num.Py()*f_num.Py());
            float denom = (*leptons.at(0) + *leptons.at(1)).Pt();
            return (numerator/denom*1.0);
        };
        *cutflow << SaveVar();
    }



    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //
    // RESTFRAMES - 0 bjets, 2 leptons
    // RESTFRAMES - 0 bjets, 2 leptons
    // RESTFRAMES - 0 bjets, 2 leptons
    //
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    double H_11_SS;
    double H_21_SS;
    double H_12_SS;
    double H_22_SS;
    double H_11_S1;
    double H_11_SS_T;
    double H_21_SS_T;
    double H_22_SS_T;
    double H_11_S1_T;
    double shat;
    double pTT_T;
    double pTT_Z;
    double RPT;
    double RPT_H_11_SS;
    double RPT_H_21_SS;
    double RPT_H_22_SS;
    double RPZ_H_11_SS;
    double RPZ_H_21_SS;
    double RPZ_H_22_SS;
    double RPT_H_11_SS_T;
    double RPT_H_21_SS_T;
    double RPT_H_22_SS_T;
    double RPZ;
    double RPZ_H_11_SS_T;
    double RPZ_H_21_SS_T;
    double RPZ_H_22_SS_T;
    double gamInvRp1;
    double MDR;
    double costheta_SS;
    double dphi_v_SS;
    double DPB_vSS;
    double cosB_1;
    double cosB_2;
    double cosB_3;
    double cosB_4;
    double dphi_v1_i1_ss;
    double dphi_s1_s2_ss;
    double dphiS_I_ss;
    double dphiS_I_s1;

    *cutflow << [&](Superlink* sl, var_void*) {
    ///////////////////////////////////////////

        // declare the frames
        LabRecoFrame lab("lab", "lab");
        DecayRecoFrame ss("ss", "ss");
        DecayRecoFrame s1("s1", "s1");
        DecayRecoFrame s2("s2", "s2");
        VisibleRecoFrame v1("v1", "v1");
        VisibleRecoFrame v2("v2", "v2");
        InvisibleRecoFrame i1("i1", "i1");
        InvisibleRecoFrame i2("i2", "i2");

        // connect the frames
        lab.SetChildFrame(ss);
        ss.AddChildFrame(s1);
        ss.AddChildFrame(s2);
        s1.AddChildFrame(v1);
        s1.AddChildFrame(i1);
        s2.AddChildFrame(v2);
        s2.AddChildFrame(i2);

        // check that the decay tree is connected properly
        if(!lab.InitializeTree()) {
            cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ << ")    Unable to initialize tree from lab frame. Exitting." << endl;
            exit(1);
        }

        // define groups

        InvisibleGroup inv("inv", "invsible group jigsaws");
        inv.AddFrame(i1);
        inv.AddFrame(i2);

        CombinatoricGroup vis("vis", "visible object jigsaws");
        vis.AddFrame(v1);
        vis.SetNElementsForFrame(v1, 1, false);
        vis.AddFrame(v2);
        vis.SetNElementsForFrame(v2, 1, false);

        SetMassInvJigsaw MinMassJigsaw("MinMass", "Invisible system mass jigsaw");
        inv.AddJigsaw(MinMassJigsaw);

        SetRapidityInvJigsaw RapidityJigsaw("RapidityJigsaw", "invisible system rapidity jigsaw");
        inv.AddJigsaw(RapidityJigsaw);
        RapidityJigsaw.AddVisibleFrames(lab.GetListVisibleFrames());

        ContraBoostInvJigsaw ContraBoostJigsaw("ContraBoostJigsaw", "ContraBoost Invariant Jigsaw");
        inv.AddJigsaw(ContraBoostJigsaw);
        ContraBoostJigsaw.AddVisibleFrames((s1.GetListVisibleFrames()), 0);
        ContraBoostJigsaw.AddVisibleFrames((s2.GetListVisibleFrames()), 1);
        ContraBoostJigsaw.AddInvisibleFrame(i1, 0);
        ContraBoostJigsaw.AddInvisibleFrame(i2, 1);

        MinMassesCombJigsaw HemiJigsaw("hemi_jigsaw", "Minimize m_{v_{1,2}} jigsaw");
        vis.AddJigsaw(HemiJigsaw);
        HemiJigsaw.AddFrame(v1, 0);
        HemiJigsaw.AddFrame(v2, 1);

        // check that the jigsaws are in place
        if(!lab.InitializeAnalysis()) {
            cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exitting." << endl;
            exit(1);
        }

        // clear the event for sho
        lab.ClearEvent();

        // set the met
        TVector3 met3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
        inv.SetLabFrameThreeVector(met3vector);

        // add leptons to the visible group
        vis.AddLabFrameFourVector(*leptons.at(0));
        vis.AddLabFrameFourVector(*leptons.at(1));

        // analayze that
        lab.AnalyzeEvent();

        //////////////////////////////
        // HT variables -- SS frame
        //////////////////////////////

        // H_1_1^SS
        TLorentzVector tlv_v1_ss = v1.GetFourVector(ss);
        TLorentzVector tlv_v2_ss = v2.GetFourVector(ss);
        TLorentzVector tlv_i1_ss = i1.GetFourVector(ss);
        TLorentzVector tlv_i2_ss = i2.GetFourVector(ss);

        TVector3 p_v1_ss = tlv_v1_ss.Vect();
        TVector3 p_v2_ss = tlv_v2_ss.Vect();
        TVector3 p_i1_ss = tlv_i1_ss.Vect();
        TVector3 p_i2_ss = tlv_i2_ss.Vect();

        TVector3 p_v_ss = p_v1_ss + p_v2_ss;
        TVector3 p_i_ss = p_i1_ss + p_i2_ss;

        H_11_SS = p_v_ss.Mag() + p_i_ss.Mag();

        // H_2_1^SS
        H_21_SS = p_v1_ss.Mag() + p_v2_ss.Mag() + p_i_ss.Mag();

        // H_1_2^SS
        H_12_SS = p_v_ss.Mag() + p_i1_ss.Mag() + p_i2_ss.Mag();

        // H_2_2^SS
        H_22_SS = p_v1_ss.Mag() + p_v2_ss.Mag() + p_i1_ss.Mag() + p_i2_ss.Mag();

        //////////////////////////////
        // HT variables -- S1 frame
        //////////////////////////////
        TLorentzVector tlv_v1_s1 = v1.GetFourVector(s1);
        TLorentzVector tlv_i1_s1 = i1.GetFourVector(s1);

        TVector3 p_v1_s1 = tlv_v1_s1.Vect();
        TVector3 p_i1_s1 = tlv_i1_s1.Vect();

        H_11_S1 = p_v1_s1.Mag() + p_i1_s1.Mag();

        ///////////////////////////////
        // transverse scale variables
        ///////////////////////////////
        TVector3 tp_v1_ss = tlv_v1_ss.Vect(); tp_v1_ss.SetZ(0.);
        TVector3 tp_v2_ss = tlv_v2_ss.Vect(); tp_v2_ss.SetZ(0.);
        TVector3 tp_i1_ss = tlv_i1_ss.Vect(); tp_i1_ss.SetZ(0.);
        TVector3 tp_i2_ss = tlv_i2_ss.Vect(); tp_i2_ss.SetZ(0.);
        TVector3 tp_v1_s1 = tlv_v1_s1.Vect(); tp_v1_s1.SetZ(0.);
        TVector3 tp_i1_s1 = tlv_i1_s1.Vect(); tp_i1_s1.SetZ(0.);

        H_11_SS_T = (tp_v1_ss + tp_v2_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();
        H_21_SS_T = tp_v1_ss.Mag() + tp_v2_ss.Mag() + (tp_i1_ss + tp_i2_ss).Mag();
        H_22_SS_T = tp_v1_ss.Mag() + tp_v2_ss.Mag() + tp_i1_ss.Mag() + tp_i2_ss.Mag();
        H_11_S1_T = tp_v1_s1.Mag() + tp_i1_s1.Mag();

        /// system mass
        shat = ss.GetMass();

        //////////////////////
        // RATIO OF CM pT
        TVector3 vPTT = ss.GetFourVector(lab).Vect();
        pTT_T = vPTT.Pt();
        pTT_Z = vPTT.Pz();
        RPT = vPTT.Pt() / (vPTT.Pt() + shat / 4.);
        RPZ = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + shat / 4.);

        RPT_H_11_SS = vPTT.Pt() / (vPTT.Pt() + H_11_SS/4.);
        RPT_H_21_SS = vPTT.Pt() / (vPTT.Pt() + H_21_SS/4.);
        RPT_H_22_SS = vPTT.Pt() / (vPTT.Pt() + H_22_SS/4.);
        RPZ_H_11_SS = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + H_11_SS/4.);
        RPZ_H_21_SS = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + H_21_SS/4.);
        RPZ_H_22_SS = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + H_22_SS/4.);

        RPT_H_11_SS_T = vPTT.Pt() / (vPTT.Pt() + H_11_SS_T/4.);
        RPT_H_21_SS_T = vPTT.Pt() / (vPTT.Pt() + H_21_SS_T/4.);
        RPT_H_22_SS_T = vPTT.Pt() / (vPTT.Pt() + H_22_SS_T/4.);
        RPZ_H_11_SS_T = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + H_11_SS_T/4.);
        RPZ_H_21_SS_T = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + H_21_SS_T/4.);
        RPZ_H_22_SS_T = fabs(vPTT.Pz()) / (fabs(vPTT.Pz()) + H_22_SS_T/4.);

        //////////////////////
        // shapes
        gamInvRp1 = ss.GetVisibleShape();

        //////////////////////
        // MDR
        MDR = 2.0 * v1.GetEnergy(s1);

        /////////////////////
        // ANGLES
        costheta_SS = ss.GetCosDecayAngle();
        dphi_v_SS = ss.GetDeltaPhiVisible();

        // costhetaB emulatur
        TVector3 v_s = s1.GetFourVector(ss).Vect().Unit();
        TVector3 v_v = v1.GetFourVector(s1).Vect().Unit();
        cosB_1 = v_s.Dot(v_v);

        cosB_2 = v1.GetCosDecayAngle(s1);

        cosB_3 = v1.GetCosDecayAngle(ss);

        TVector3 v_v2 = v1.GetFourVector(ss).Vect().Unit();
        cosB_4 = v_s.Dot(v_v2);

        // angle between invisible
        dphi_v1_i1_ss = -1.;//v1.GetFourVector(ss).DeltaPhi(i1.GetFourVector(ss));
        dphi_s1_s2_ss = -1.;//s1.GetFourVector(ss).DeltaPhi(s2.GetFourVector(ss)); 


        dphiS_I_ss = -1.;//s1.GetFourVector(ss).DeltaPhi(i1.GetFourVector(ss));
        dphiS_I_s1 = -1.;//s1.GetFourVector(ss).DeltaPhi(i1.GetFourVector(s1));
        
        

        ////////////////////
        // BOOST ANGLES
        DPB_vSS = ss.GetDeltaPhiBoostVisible();



    ///////////////////////////////////////////
    }; //end void


    *cutflow << NewVar("HT : H_11_SS"); {
        *cutflow << HFTname("H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("HT : H_21_SS"); {
        *cutflow << HFTname("H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("HT : H_12_SS"); {
        *cutflow << HFTname("H_12_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_12_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("HT : H_22_SS"); {
        *cutflow << HFTname("H_22_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_22_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("HT : H_11_S1"); {
        *cutflow << HFTname("H_11_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_11_S1;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("H_11_SS_T"); {
        *cutflow << HFTname("H_11_SS_T");
        *cutflow <<[&](Superlink* sl, var_float*) -> double {
            return H_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("H_21_SS_T"); {
        *cutflow << HFTname("H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("H_22_SS_T"); {
        *cutflow << HFTname("H_22_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_22_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("H_11_S1_T"); {
        *cutflow << HFTname("H_11_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return H_11_S1_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("shat"); {
        *cutflow << HFTname("shat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return shat;
        };
        *cutflow << SaveVar();  
    }
    *cutflow << NewVar("pTT_T"); {
        *cutflow << HFTname("pTT_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return pTT_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pTT_Z"); {
        *cutflow << HFTname("pTT_Z");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return pTT_Z;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPT"); {
        *cutflow << HFTname("RPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ"); {
        *cutflow << HFTname("RPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPT_H_11_SS"); {
        *cutflow << HFTname("RPT_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_H_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPT_H_21_SS"); {
        *cutflow << HFTname("RPT_H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_H_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPT_H_22_SS"); {
        *cutflow << HFTname("RPT_H_22_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_H_22_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ_H_11_SS"); {
        *cutflow << HFTname("RPZ_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_H_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ_H_21_SS"); {
        *cutflow << HFTname("RPZ_H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_H_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ_H_22_SS"); {
        *cutflow << HFTname("RPZ_H_22_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_H_22_SS;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("RPT_H_11_SS_T"); {
        *cutflow << HFTname("RPT_H_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_H_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPT_H_21_SS_T"); {
        *cutflow << HFTname("RPT_H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_H_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPT_H_22_SS_T"); {
        *cutflow << HFTname("RPT_H_22_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_H_22_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ_H_11_SS_T"); {
        *cutflow << HFTname("RPZ_H_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_H_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ_H_21_SS_T"); {
        *cutflow << HFTname("RPZ_H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_H_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RPZ_H_22_SS_T"); {
        *cutflow << HFTname("RPZ_H_22_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_H_22_SS_T;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("gamInvRp1"); {
        *cutflow << HFTname("gamInvRp1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return gamInvRp1;
        };
        *cutflow << SaveVar();  
    }
    *cutflow << NewVar("MDR"); {
        *cutflow << HFTname("MDR");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("costheta_SS"); {
        *cutflow << HFTname("costheta_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return costheta_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dphi_v_SS"); {
        *cutflow << HFTname("dphi_v_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_v_SS;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dphiS_I_SS"); {
        *cutflow << HFTname("dphiS_I_ss");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphiS_I_ss;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dphiS_I_s1"); {
        *cutflow << HFTname("dphiS_I_s1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphiS_I_s1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("deltaX (WW-like)"); {
        *cutflow << HFTname("deltaX");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double deltaX = -999.0;
            if(leptons.size()==2) {
                double sqrtS = 13000.0;
                double num = leptons.at(0)->Pz() + leptons.at(1)->Pz();
                deltaX = num / sqrtS * 1.0;
            }
            return deltaX;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("cosThetaB emulator 1"); {
        *cutflow << HFTname("cosB_1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return cosB_1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("cosThetaB emulator 2"); {
        *cutflow << HFTname("cosB_2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return cosB_2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("cosThetaB emulator 3"); {
        *cutflow << HFTname("cosB_3");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return cosB_3;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("cosThetaB emulator 4"); {
        *cutflow << HFTname("cosB_4");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return cosB_4;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between visible & invisible in SS frame"); {
        *cutflow << HFTname("dphi_v1_i1_ss");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_v1_i1_ss;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between s1 and s2 in SS frame"); {
        *cutflow << HFTname("dphi_s1_s2_ss");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_s1_s2_ss;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("DPB_vSS"); {
        *cutflow << HFTname("DPB_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPB_vSS;
        };
        *cutflow << SaveVar();
    }
//
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    //
//    // RESTFRAMES - 2 jets, 2 leptons
//    // RESTFRAMES - 2 jets, 2 leptons
//    // RESTFRAMES - 2 jets, 2 leptons
//    // RESTFRAMES - 2 jets, 2 leptons
//    //
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
/////*
//    double xshat;
//    double xgaminv;
//    double xRPT;
//    double xRPZ;
//    double xcosSS;
//    double xdphiLSS;
//    double xMS;
//    double xPS;
//    double xMSS;
//    double xgaminvSS;
//    double xDeltaBetaSS;
//    double xDPD_vSS;
//    double xDPB_vSS;
//    int xNV[2]; // number of visible objects in hemisphere
//    double xcosS[2]; // cosine stop decay angle
//    double xcosC[2]; // cosine intermediate child decay angle
//    double xdphiSC[2]; // cosine between stop and child decay planes
//    double xRCS[2]; // ratio of child and stop masses (w/ WIMP masses subtracted);
//    double xjet1PT[2]; // first leading jet pT associated with this hemisphere
//    double xjet2PT[2]; // second leading jet pT associated with this hemisphere
//    double xPinv[2]; // Pinv / HS
//    double xH_11_SS;
//    double xH_21_SS;
//    double xH_41_SS;
//    double xH_42_SS;
//    double xH_11_S1;
//    double xH_21_S1;
//    double xH_11_SS_T;
//    double xH_21_SS_T;
//    double xH_41_SS_T;
//    double xH_42_SS_T;
//    double xH_11_S1_T;
//    double xH_21_S1_T;
//
//    double xRPT_H_11_SS; 
//    double xRPT_H_21_SS; 
//    double xRPT_H_41_SS; 
//    double xRPT_H_42_SS; 
//    double xRPZ_H_11_SS; 
//    double xRPZ_H_21_SS; 
//    double xRPZ_H_41_SS; 
//    double xRPZ_H_42_SS; 
//    double xRPT_H_11_SS_T; 
//    double xRPT_H_21_SS_T; 
//    double xRPT_H_41_SS_T; 
//    double xRPT_H_42_SS_T; 
//    double xRPZ_H_11_SS_T; 
//    double xRPZ_H_21_SS_T; 
//    double xRPZ_H_41_SS_T; 
//    double xRPZ_H_42_SS_T; 
//
//    double xdphiVS_I[2];
//
//    // feb 1 vars
//    //double xCosP1;
//    //double xCosP2;
//
//
//    *cutflow << [&](Superlink* sl, var_void*) {
//        if(sjets.size()>=2 && bjets.size()>=1) {
//
//            // setup the analysis tree
//            LabRecoFrame xlab("xlab", "xlab");
//            DecayRecoFrame xss("xss", "xss");
//            DecayRecoFrame xs1("xs1", "xs1");
//            DecayRecoFrame xs2("xs2", "xs2");
//            DecayRecoFrame xc1("xc1", "xc1");
//            DecayRecoFrame xc2("xc2", "xc2");
//            VisibleRecoFrame xv1s("xv1s", "xv1s");
//            VisibleRecoFrame xv2s("xv2s", "xv2s");
//            InvisibleRecoFrame xi1("xi1", "xi1");
//            VisibleRecoFrame xv1c("xv1c", "xv1c");
//            VisibleRecoFrame xv2c("xv2c", "xv2c");
//            InvisibleRecoFrame xi2("xi2", "xi2");
//
//            xlab.SetChildFrame(xss);
//            xss.AddChildFrame(xs1);
//            xss.AddChildFrame(xs2);
//            xs1.AddChildFrame(xc1);
//            xs1.AddChildFrame(xv1s);
//            xc1.AddChildFrame(xi1);
//            xc1.AddChildFrame(xv1c);
//            xs2.AddChildFrame(xc2);
//            xs2.AddChildFrame(xv2s);
//            xc2.AddChildFrame(xi2);
//            xc2.AddChildFrame(xv2c);
//
//            // check that the decay tree is connected properly
//            if(!xlab.InitializeTree()) {
//                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ <<")    Unable to initialize tree from lab frame. Exiting." << endl;
//                exit(1);
//            }
//
//
//            // define groupes
//
//            InvisibleGroup xinv("xinv", "x-invisible gruop jigsaws");
//            xinv.AddFrame(xi1);
//            xinv.AddFrame(xi2);
//
//            CombinatoricGroup xvis("xvis", "x-visible object jigsaws");
//
//            // visible frames in first decay step must always have at least one element
//            xvis.AddFrame(xv1s);
//            xvis.AddFrame(xv2s);
//            xvis.SetNElementsForFrame(xv1s, 1, false);
//            xvis.SetNElementsForFrame(xv2s, 1, false);
//            // visible frames in second decay step can have zero elements
//            xvis.AddFrame(xv1c);
//            xvis.AddFrame(xv2c);
//            xvis.SetNElementsForFrame(xv1c, 0, false);
//            xvis.SetNElementsForFrame(xv2c, 0, false);
//
//            // define jigsaws
//            SetMassInvJigsaw xMinMassJigsaw("xminmass", "x-Invisible system mass jigsaw");
//            xinv.AddJigsaw(xMinMassJigsaw);
//            SetRapidityInvJigsaw xRapidityJigsaw("xRapidity", "x-Invisible system rapidity jigsaw");
//            xinv.AddJigsaw(xRapidityJigsaw);
//            xRapidityJigsaw.AddVisibleFrames((xlab.GetListVisibleFrames()));
//            ContraBoostInvJigsaw xContraBoostJigsaw("xContra", "x-Contraboost invariant jigsaw");
//            xinv.AddJigsaw(xContraBoostJigsaw);
//            xContraBoostJigsaw.AddVisibleFrames((xs1.GetListVisibleFrames()), 0);
//            xContraBoostJigsaw.AddVisibleFrames((xs2.GetListVisibleFrames()), 1);
//            xContraBoostJigsaw.AddInvisibleFrames((xs1.GetListInvisibleFrames()), 0);
//            xContraBoostJigsaw.AddInvisibleFrames((xs2.GetListInvisibleFrames()), 1);
//
//            MinMassesCombJigsaw xHemiJigsaw("xHemi", "x-Minimize m_{V_{1,2}} jigsaw");
//            xvis.AddJigsaw(xHemiJigsaw);
//            xHemiJigsaw.AddFrame(xv1s, 0);
//            xHemiJigsaw.AddFrame(xv2s, 1);
//            xHemiJigsaw.AddFrame(xv1c, 0);
//            xHemiJigsaw.AddFrame(xv2c, 1);
//
//            MinMassesCombJigsaw x1HemiJigsaw("x-1 Hemi", "x-1 Minimize m_{C_{a}} jigsaw");
//            xvis.AddJigsaw(x1HemiJigsaw);
//            x1HemiJigsaw.AddFrame(xv1s, 0);
//            x1HemiJigsaw.AddFrame(xv1c, 1);
//            x1HemiJigsaw.AddFrame(xi1, 1);
//
//            MinMassesCombJigsaw x2HemiJigsaw("x-2 Hemi", "x-2 Minimize m_{C_{b}} jigsaw");
//            xvis.AddJigsaw(x2HemiJigsaw);
//            x2HemiJigsaw.AddFrame(xv2s, 0);
//            x2HemiJigsaw.AddFrame(xv2c, 1);
//            x2HemiJigsaw.AddFrame(xi2, 1);
//
//            // check that the jigsaws are in place
//            if(!xlab.InitializeAnalysis()) {
//                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
//                exit(1);
//            }
//
//            // clear the event for sho
//            xlab.ClearEvent();
//
//            // set the met
//            TVector3 xmet3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
//            xinv.SetLabFrameThreeVector(xmet3vector);
//
//            // set up the jets and leptons
//            //vector<TLorentzVector> xjets;
//            //TLorentzVector l1, l2;
//
//            //l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
//            //l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
//            //for(int ij = 0; ij < (int)jets.size(); ij++) {
//            //    TLorentzVector xj;
//            //    xj.SetPtEtaPhiM(jets.at(ij)->Pt(), jets.at(ij)->Eta(), jets.at(ij)->Phi(), jets.at(ij)->M());
//            //    xjets.push_back(xj);
//            //}
//
//            //vector<RFKey> jetID;
//            //jetID.push_back(xvis.AddLabFrameFourVector(l1));
//            //jetID.push_back(xvis.AddLabFrameFourVector(l2));
//            //for(int ij = 0; ij < (int)xjets.size(); ij++) {
//            //    jetID.push_back(xvis.AddLabFrameFourVector(xjets.at(ij)));
//            //}
//
//
//            TLorentzVector j1, j2, l1, l2;
//            j1.SetPtEtaPhiM(jets.at(0)->Pt(), jets.at(0)->Eta(), jets.at(0)->Phi(), jets.at(0)->M());
//            j2.SetPtEtaPhiM(jets.at(1)->Pt(), jets.at(1)->Eta(), jets.at(1)->Phi(), jets.at(1)->M());
//            l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
//            l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Eta(), leptons.at(1)->M());
//
//            vector<RFKey> jetID;
//            jetID.push_back(xvis.AddLabFrameFourVector(j1));
//            jetID.push_back(xvis.AddLabFrameFourVector(j2));
//
//            vector<RFKey> lepID;
//            lepID.push_back(xvis.AddLabFrameFourVector(l1));
//            lepID.push_back(xvis.AddLabFrameFourVector(l2));
//
//
//            // analyze that
//            xlab.AnalyzeEvent();
//
//            DecayRecoFrame* S[2];
//            DecayRecoFrame* C[2];
//            VisibleRecoFrame* VS[2];
//            VisibleRecoFrame* VC[2];
//            InvisibleRecoFrame* I[2];
//            //randomize
//            int flip = (gRandom->Rndm() > 0.5);
//            S[flip] = &xs1;
//            S[(flip+1)%2] = &xs2;
//            C[flip] = &xc1;
//            C[(flip+1)%2] = &xc2;
//            VS[flip] = &xv1s;
//            VS[(flip+1)%2] = &xv2s;
//            VC[flip] = &xv1c;
//            VC[(flip+1)%2] = &xv2c;
//            I[flip] = &xi1;
//            I[(flip+1)%2] = &xi2;
//
//
//            //////////////////////////////////////
//            // Observables
//            //////////////////////////////////////
//
//
//            // total CM mass
//            xshat = xss.GetMass();
//            // 'mass-less' stop assumption gamma in CM frame
//            xgaminv = xss.GetVisibleShape();
//
//            TVector3 xvPSS = xss.GetFourVector(xlab).Vect();
//
//            // ratio of CM pT to CM mass
//            xRPT = xvPSS.Pt() / (xvPSS.Pt() + xshat/4.);
//            // ratio of CM pz to CM mass
//            xRPZ = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xshat/4.);
//            // cos decay angle of ss system
//            xcosSS = xss.GetCosDecayAngle();
//            // delta phi between lab and SS decay planes
//            xdphiLSS = xlab.GetDeltaPhiDecayPlanes(xss);
//
//
//            TLorentzVector vVS1 = S[0]->GetVisibleFourVector(*S[0]);
//            TLorentzVector vVS2 = S[1]->GetVisibleFourVector(*S[1]);
//
//            // stop mass
//            xMS = (vVS1.M2() - vVS2.M2()) / (2.*(vVS1.E() - vVS2.E()));
//
//            xPS = S[0]->GetMomentum(xss);
//            xMSS = 2.*sqrt(xPS*xPS + xMS*xMS);
//            xgaminvSS = 2.*xMS/xMSS;
//            double beta = sqrt(1.-xgaminv*xgaminv);
//            double betaSS = sqrt(1.-xgaminvSS*xgaminvSS);
//
//            // velocity difference between 'massive' and 'mass-less'
//            xDeltaBetaSS = -(betaSS-beta)/(1.-betaSS*beta);
//
//            // dleta phi between SS visible decay products and SS decay axis
//            xDPD_vSS = xss.GetDeltaPhiDecayVisible();
//            // delta phi between SS visible decay products and SS momentum
//            xDPB_vSS = xss.GetDeltaPhiBoostVisible();
//
//
//            // number of visible objects in hemisphere
//            for(int i = 0; i < 2; i++) {
//                xNV[i] = xvis.GetNElementsInFrame(*VS[i]); 
//                xNV[i] += xvis.GetNElementsInFrame(*VC[i]);
//
//                TVector3 xvP1 = VS[i]->GetFourVector(*S[i]).Vect();
//                TVector3 xvP2 = VC[i]->GetFourVector(*S[i]).Vect();
//                xPinv[i] = 2.*(xvP1+xvP2).Mag()/(xvP1.Mag()+xvP2.Mag() + (xvP1+xvP2).Mag());
//
//                xcosS[i] = S[i]->GetCosDecayAngle();
//
//                int N = jetID.size();
//                double pTmax[2]; pTmax[0] = -1.; pTmax[1] = -1.;
//                for(int j = 0; j < N; j++) {
//                    const RestFrame& frame = xvis.GetFrame(jetID[j]);
//                    if(VS[i]->IsSame(frame) || VC[i]->IsSame(frame)) { // jet is in hemisphere 'i'
//                        double pT_ = xvis.GetLabFrameFourVector(jetID[j]).Pt();
//                        if(pT_ > pTmax[0]) {
//                            pTmax[1] = pTmax[0];
//                            pTmax[0] = pT_;
//                        } else {
//                            if(pT_ > pTmax[1]) pTmax[1] = pT_;
//                        }
//                    }
//                } // j
//
//                xjet1PT[i] = pTmax[0]; // lead visible object pT in hemisphere i
//                xjet2PT[i] = pTmax[1]; // sub lead visible object pT in hemisphere i
//
//                xdphiVS_I[i] = VS[i]->GetFourVector(*S[i]).DeltaPhi(I[i]->GetFourVector(*S[i]));
//
//                if(xNV[i] > 1) {
//                    xcosS[i] = C[i]->GetCosDecayAngle();
//                    xdphiSC[i] = S[i]->GetDeltaPhiDecayPlanes(*C[i]);
//                    xRCS[i] = (C[i]->GetMass() - I[i]->GetMass())/(S[i]->GetMass()-I[i]->GetMass());
//
//                } else {
//                    xcosS[i] = -2;
//                    xdphiSC[i] = -5;
//                    xRCS[i] = -5;
//                }
//                
//            } // i
//
//            ////////////////////////////////////////////////////////////
//            // scale variables
//            ////////////////////////////////////////////////////////////
//
//            TLorentzVector v_v1s_ss = VS[0]->GetFourVector(xss);
//            TLorentzVector v_v2s_ss = VS[1]->GetFourVector(xss);
//            TLorentzVector v_v1c_ss = VC[0]->GetFourVector(xss);
//            TLorentzVector v_v2c_ss = VC[1]->GetFourVector(xss);
//            TLorentzVector v_i1_ss  = I[0]->GetFourVector(xss);
//            TLorentzVector v_i2_ss  = I[1]->GetFourVector(xss);
//            if(v_i1_ss.Vect().Mag() > 1.0e4) {
//                v_i1_ss.SetPtEtaPhiM(0, 0, 0, 0);
//            }
//            if(v_i2_ss.Vect().Mag() > 1.0e4) {
//                v_i2_ss.SetPtEtaPhiM(0, 0, 0, 0);
//            }
//
//            TLorentzVector v_v1s_s1 = VS[0]->GetFourVector(xs1);
//            TLorentzVector v_v1c_s1 = VC[0]->GetFourVector(xs1);
//            TLorentzVector v_i1_s1  = I[0]->GetFourVector(xs1);
//
//            // H_11_SS
//            TVector3 p_vis_H11SS = (v_v1s_ss + v_v2s_ss + v_v1c_ss + v_v2c_ss).Vect();
//            TVector3 p_invis_H11SS = (v_i1_ss + v_i2_ss).Vect();
//            xH_11_SS = p_vis_H11SS.Mag() + p_invis_H11SS.Mag();
//
//            // H_21_SS
//            TVector3 p_vis_1_H21SS = (v_v1s_ss + v_v1c_ss).Vect();
//            TVector3 p_vis_2_H21SS = (v_v2s_ss + v_v2c_ss).Vect();
//            TVector3 p_invis_H21SS = (v_i1_ss + v_i2_ss).Vect();
//            xH_21_SS = p_vis_1_H21SS.Mag() + p_vis_2_H21SS.Mag() + p_invis_H21SS.Mag();
//
//            // H_41_SS
//            xH_41_SS = v_v1s_ss.Vect().Mag() + v_v2s_ss.Vect().Mag() + v_v1c_ss.Vect().Mag() + v_v2c_ss.Vect().Mag() + (v_i1_ss + v_i2_ss).Vect().Mag();
//
//            // H_42_SS
//            xH_42_SS = v_v1s_ss.Vect().Mag() + v_v2s_ss.Vect().Mag() + v_v1c_ss.Vect().Mag() + v_v2c_ss.Vect().Mag() + v_i1_ss.Vect().Mag() + v_i2_ss.Vect().Mag();
//
//            // H_11_S1
//            xH_11_S1 = (v_v1s_s1+v_v1c_s1).Vect().Mag() + v_i1_s1.Vect().Mag();
//
//            // H_21_S1
//            xH_21_S1 = v_v1s_s1.Vect().Mag() + v_v1c_s1.Vect().Mag() + v_i1_s1.Vect().Mag();
//
//            ////////////////
//            // transverse scale variables
//            ////////////////
//            TVector3 tp_v1s_ss = v_v1s_ss.Vect(); tp_v1s_ss.SetZ(0.);
//            TVector3 tp_v2s_ss = v_v2s_ss.Vect(); tp_v2s_ss.SetZ(0.);
//            TVector3 tp_v1c_ss = v_v1c_ss.Vect(); tp_v1c_ss.SetZ(0.);
//            TVector3 tp_v2c_ss = v_v2c_ss.Vect(); tp_v2c_ss.SetZ(0.);
//            TVector3 tp_i1_ss  = v_i1_ss.Vect();  tp_i1_ss.SetZ(0.);
//            TVector3 tp_i2_ss  = v_i2_ss.Vect();  tp_i2_ss.SetZ(0.);
//
//            TVector3 tp_v1s_s1 = v_v1s_s1.Vect(); tp_v1s_s1.SetZ(0.);
//            TVector3 tp_v1c_s1 = v_v1c_s1.Vect(); tp_v1c_s1.SetZ(0.);
//            TVector3 tp_i1_s1  = v_i1_s1.Vect();  tp_i1_s1.SetZ(0.);
//
//            // H_11_SS_T
//            xH_11_SS_T = (tp_v1s_ss + tp_v2s_ss + tp_v1c_ss + tp_v2c_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();
//
//            // H_21_SS_T
//            xH_21_SS_T = (tp_v1s_ss + tp_v1c_ss).Mag() + (tp_v2s_ss + tp_v2c_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();
//
//            // H_41_SS_T
//            xH_41_SS_T = tp_v1s_ss.Mag() + tp_v2s_ss.Mag() + tp_v1c_ss.Mag() + tp_v2c_ss.Mag() + (tp_i1_ss+tp_i2_ss).Mag();
//
//            // H_42_SS_T
//            xH_42_SS_T = tp_v1s_ss.Mag() + tp_v2s_ss.Mag() + tp_v1c_ss.Mag() + tp_v2c_ss.Mag() + tp_i1_ss.Mag() + tp_i2_ss.Mag();
//
//            // H_11_S1_T
//            xH_11_S1_T = (tp_v1s_s1 + tp_v1c_s1).Mag() + tp_i1_s1.Mag();
//
//            // H_21_S1_T
//            xH_21_S1_T = tp_v1s_s1.Mag() + tp_v1c_s1.Mag() + tp_i1_s1.Mag();
//        
//
//            xRPT_H_11_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_11_SS/4.);
//            xRPT_H_21_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_21_SS/4.);
//            xRPT_H_41_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_41_SS/4.);
//            xRPT_H_42_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_42_SS/4.);
//            xRPZ_H_11_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_11_SS/4.);
//            xRPZ_H_21_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_21_SS/4.);
//            xRPZ_H_41_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_41_SS/4.);
//            xRPZ_H_42_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_42_SS/4.);
//
//            xRPT_H_11_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_11_SS_T/4.);
//            xRPT_H_21_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_21_SS_T/4.);
//            xRPT_H_41_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_41_SS_T/4.);
//            xRPT_H_42_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_42_SS_T/4.);
//            xRPZ_H_11_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_11_SS_T/4.);
//            xRPZ_H_21_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_21_SS_T/4.);
//            xRPZ_H_41_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_41_SS_T/4.);
//            xRPZ_H_42_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_42_SS_T/4.);
//
//
//            ///////////////////////////
//            // feb 1 variables
//            //xCosP1 = S[0]->GetCosDecayAngle(I[0]);
//            //xCosP2 = S[1]->GetCosDecayAngle(I[1]);
//
//        } // njets == 2
//    };
//
//    *cutflow << NewVar("xshat"); {
//        *cutflow << HFTname("xshat");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xshat;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xgaminv"); {
//        *cutflow << HFTname("xgaminv");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xgaminv;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT"); {
//        *cutflow << HFTname("xRPT");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ"); {
//        *cutflow << HFTname("xRPZ");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xcosSS"); {
//        *cutflow << HFTname("xcosSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xcosSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xdphiLSS"); {
//        *cutflow << HFTname("xdphiLSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xdphiLSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xMS"); {
//        *cutflow << HFTname("xMS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xMS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xPS"); {
//        *cutflow << HFTname("xPS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xPS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xMSS"); {
//        *cutflow << HFTname("xMSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xMSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xgaminvSS"); {
//        *cutflow << HFTname("xgaminvSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xgaminvSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xDeltaBetaSS"); {
//        *cutflow << HFTname("xDeltaBetaSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xDeltaBetaSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xDPD_vSS"); {
//        *cutflow << HFTname("xDPD_vSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xDPD_vSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xDPB_vSS"); {
//        *cutflow << HFTname("xDPB_vSS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xDPB_vSS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xNV"); {
//        *cutflow << HFTname("xNV");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xNV[0]);
//            out.push_back(xNV[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xNV_0"); {
//        *cutflow << HFTname("xNV_0");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            return xNV[0];
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xNV_1"); {
//        *cutflow << HFTname("xNV_1");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            return xNV[1];
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xcosS"); {
//        *cutflow << HFTname("xcosS");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xcosS[0]);
//            out.push_back(xcosS[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xcosC"); {
//        *cutflow << HFTname("xcosC");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xcosC[0]);
//            out.push_back(xcosC[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xdphiSC"); {
//        *cutflow << HFTname("xdphiSC");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xdphiSC[0]);
//            out.push_back(xdphiSC[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRCS"); {
//        *cutflow << HFTname("xRCS");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xRCS[0]);
//            out.push_back(xRCS[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xjet1PT"); {
//        *cutflow << HFTname("xjet1PT");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xjet1PT[0]);
//            out.push_back(xjet1PT[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xjet2PT"); {
//        *cutflow << HFTname("xjet2PT");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xjet2PT[0]);
//            out.push_back(xjet2PT[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xPinv"); {
//        *cutflow << HFTname("xPinv");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xPinv[0]);
//            out.push_back(xPinv[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_11_SS"); {
//        *cutflow << HFTname("xH_11_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_11_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_21_SS"); {
//        *cutflow << HFTname("xH_21_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_21_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_41_SS"); {
//        *cutflow << HFTname("xH_41_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_41_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_42_SS"); {
//        *cutflow << HFTname("xH_42_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_42_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_11_S1"); {
//        *cutflow << HFTname("xH_11_S1");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_11_S1;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_21_S1"); {
//        *cutflow << HFTname("xH_21_S1");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_21_S1;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_11_SS_T"); {
//        *cutflow << HFTname("xH_11_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_11_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_21_SS_T"); {
//        *cutflow << HFTname("xH_21_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_21_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_41_SS_T"); {
//        *cutflow << HFTname("xH_41_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_41_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_42_SS_T"); {
//        *cutflow << HFTname("xH_42_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_42_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_11_S1_T"); {
//        *cutflow << HFTname("xH_11_S1_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_11_S1_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xH_21_S1_T"); {
//        *cutflow << HFTname("xH_21_S1_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xH_21_S1_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_11_SS"); {
//        *cutflow << HFTname("xRPT_H_11_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_11_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_21_SS"); {
//        *cutflow << HFTname("xRPT_H_21_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_21_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_41_SS"); {
//        *cutflow << HFTname("xRPT_H_41_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_41_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_42_SS"); {
//        *cutflow << HFTname("xRPT_H_42_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_42_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    *cutflow << NewVar("xRPZ_H_11_SS"); {
//        *cutflow << HFTname("xRPZ_H_11_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_11_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_21_SS"); {
//        *cutflow << HFTname("xRPZ_H_21_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_21_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_41_SS"); {
//        *cutflow << HFTname("xRPZ_H_41_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_41_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_42_SS"); {
//        *cutflow << HFTname("xRPZ_H_42_SS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_42_SS;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_11_SS_T"); {
//        *cutflow << HFTname("xRPT_H_11_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_11_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_21_SS_T"); {
//        *cutflow << HFTname("xRPT_H_21_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_21_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_41_SS_T"); {
//        *cutflow << HFTname("xRPT_H_41_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_41_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPT_H_42_SS_T"); {
//        *cutflow << HFTname("xRPT_H_42_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPT_H_42_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_11_SS_T"); {
//        *cutflow << HFTname("xRPZ_H_11_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_11_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_21_SS_T"); {
//        *cutflow << HFTname("xRPZ_H_21_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_21_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_41_SS_T"); {
//        *cutflow << HFTname("xRPZ_H_41_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_41_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xRPZ_H_42_SS_T"); {
//        *cutflow << HFTname("xRPZ_H_42_SS_T");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            return xRPZ_H_42_SS_T;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("xdphiVS_I"); {
//        *cutflow << HFTname("xdphiVS_I");
//        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
//            vector<double> out;
//            out.clear();
//            out.push_back(xdphiVS_I[0]);
//            out.push_back(xdphiVS_I[1]);
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//*/
/*
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // RESTFRAMES - mW
    // RESTFRAMES - mW
    // RESTFRAMES - mW
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    double w_shat;
    double w_gaminv;
    double w_pTCM;
    double w_pZCM;
    double w_RPT;
    double w_RPZ;
    double w_cosSS;
    double w_dphiLSS;
    double w_dphiSS_SA;
    double w_dphiSS_SB;
    double w_dphiSA_SB;
    double w_MS;
    double w_DPD_vSS;
    double w_DPB_vSS;
    int w_N_VA1;
    int w_N_VA2;
    int w_N_VB1;
    int w_N_VB2;
    int w_N_VA;
    int w_N_VB;
    int w_N_V1;
    int w_N_V2;
    int w_Nsj_VA1;
    int w_Nsj_VA2;
    int w_Nsj_VB1;
    int w_Nsj_VB2;
    int w_Nbj_VA1;
    int w_Nbj_VA2;
    int w_Nbj_VB1;
    int w_Nbj_VB2;
    int w_Nsj_A ;
    int w_Nsj_B ;
    int w_Nsj_V1;
    int w_Nsj_V2;
    int w_Nbj_A ;
    int w_Nbj_B ;
    int w_Nbj_V1;
    int w_Nbj_V2;
    double w_H_11_SS;
    double w_H_21_SS;
    double w_H_41_SS;
    double w_H_22_SS;
    double w_H_42_SS;
    double w_H_11_SS_T;
    double w_H_21_SS_T;
    double w_H_41_SS_T;
    double w_H_22_SS_T;
    double w_H_42_SS_T;
    double w_H_11_SA;
    double w_H_21_SA;
    double w_H_11_SA_T;
    double w_H_21_SA_T;
    double w_MDR_VA1_SA;

    double w_gaminv_SA;
    double w_H_11_SS_lep;
    double w_H_21_SS_lep;
    double w_H_41_SS_lep;
    double w_H_22_SS_lep;
    double w_H_42_SS_lep;
    double w_H_61_SS_lep;
    double w_H_62_SS_lep;
    double w_H_11_SS_lep_T;
    double w_H_21_SS_lep_T;
    double w_H_41_SS_lep_T;
    double w_H_22_SS_lep_T;
    double w_H_42_SS_lep_T;
    double w_H_61_SS_lep_T;
    double w_H_62_SS_lep_T;
    double w_H_11_SA_lep;
    double w_H_21_SA_lep;
    double w_H_31_SA_lep;
    double w_H_11_SA_lep_T;
    double w_H_21_SA_lep_T;
    double w_H_31_SA_lep_T;
    double w_MDR_CLA_CA;


    *cutflow << [&](Superlink* sl, var_void*) {
        if(jets.size()>=2 && leptons.size()==2) {

            // setup the tree
            LabRecoFrame wLAB("LAB", "LAB");
            DecayRecoFrame wSS("SS", "SS");
            DecayRecoFrame wSA("SA", "SA");
            DecayRecoFrame wSB("SB", "SB");
            DecayRecoFrame wCA("CA", "CA");
            DecayRecoFrame wCB("CB", "CB");
            VisibleRecoFrame wVA1("VA1", "VA1");
            VisibleRecoFrame wVB1("VB1", "VB1");
            InvisibleRecoFrame wIA("IA", "IA");
            InvisibleRecoFrame wIB("IB", "IB");
            DecayRecoFrame wCLA("CLA", "CLA");
            DecayRecoFrame wCLB("CLB", "CLB");
            VisibleRecoFrame wVA2("VA2", "VA2");
            VisibleRecoFrame wLA("LA", "LA");
            VisibleRecoFrame wVB2("VB2", "VB2");
            VisibleRecoFrame wLB("LB", "LB");

            wLAB.SetChildFrame(wSS);

            wSS.AddChildFrame(wSA);
            wSS.AddChildFrame(wSB);

            wSA.AddChildFrame(wCA);
            wSA.AddChildFrame(wVA1);
            wSB.AddChildFrame(wCB);
            wSB.AddChildFrame(wVB1);

            wCA.AddChildFrame(wCLA);
            wCA.AddChildFrame(wIA);
            wCB.AddChildFrame(wCLB);
            wCB.AddChildFrame(wIB);

            wCLA.AddChildFrame(wVA2);
            wCLA.AddChildFrame(wLA);
            wCLB.AddChildFrame(wVB2);
            wCLB.AddChildFrame(wLB);

            if(!wLAB.InitializeTree()) {
                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ << ")    Unable to initialize tree from lab frame. Exiting." << endl;
                exit(1);
            }

            // define groups
            InvisibleGroup wINV("wINV", "w-Invisible group jigsaws");
            wINV.AddFrame(wIA);
            wINV.AddFrame(wIB);

            CombinatoricGroup wVIS("wVIS", "w-Visible object jigsaws");

            wVIS.AddFrame(wVA1);
            wVIS.AddFrame(wVB1);
            wVIS.AddFrame(wVA2);
            wVIS.AddFrame(wVB2);
            wVIS.SetNElementsForFrame(wVA1, 1, false);
            wVIS.SetNElementsForFrame(wVB1, 1, false);
            wVIS.SetNElementsForFrame(wVA2, 0, false); // don't require any jets in 2nd layer 
            wVIS.SetNElementsForFrame(wVB2, 0, false); // don't require any jets in 2nd layer  

            // define jigsaws
            SetMassInvJigsaw wMinMassJigsaw("wminmass", "w-Invisible system mass jigsaw");
            wINV.AddJigsaw(wMinMassJigsaw);

            SetRapidityInvJigsaw wRapidityJigsaw("wRapidity", "w-Invisible system rapidity jigsaw");
            wINV.AddJigsaw(wRapidityJigsaw);
            wRapidityJigsaw.AddVisibleFrames((wLAB.GetListVisibleFrames()));

            ContraBoostInvJigsaw wContraBoostJigsaw("wContra", "w-ContraBoost Invariant Jigsaw");
            wINV.AddJigsaw(wContraBoostJigsaw);
            wContraBoostJigsaw.AddVisibleFrames((wSA.GetListVisibleFrames()), 0);
            wContraBoostJigsaw.AddVisibleFrames((wSB.GetListVisibleFrames()), 1);
            wContraBoostJigsaw.AddInvisibleFrames((wSA.GetListInvisibleFrames()), 0);
            wContraBoostJigsaw.AddInvisibleFrames((wSB.GetListInvisibleFrames()), 1);

            MinMassesCombJigsaw wHemiJigsaw("wHemi", "w-Minimize m_{V_{a,b}} jigsaw");
            wVIS.AddJigsaw(wHemiJigsaw);
            wHemiJigsaw.AddFrame(wVA1, 0);
            wHemiJigsaw.AddFrame(wVA2, 0);
            wHemiJigsaw.AddFrame(wVB1, 1);
            wHemiJigsaw.AddFrame(wVB2, 1);

            MinMassesCombJigsaw waHemiJigsaw("wHemi-a", "w-Minimize m_{C_{a}} jigsaw");
            wVIS.AddJigsaw(waHemiJigsaw);
            waHemiJigsaw.AddFrame(wVA1, 0);
            waHemiJigsaw.AddFrame(wVA2, 1);
            waHemiJigsaw.AddFrame(wIA,  1);

            MinMassesCombJigsaw wbHemiJigsaw("wHemi-b", "w-Minimize m_{C_{b}} jigsaw");
            wVIS.AddJigsaw(wbHemiJigsaw);
            wbHemiJigsaw.AddFrame(wVB1, 0);
            wbHemiJigsaw.AddFrame(wVB2, 1);
            wbHemiJigsaw.AddFrame(wIB,  1);

            // lepton jigsaws
            CombinatoricGroup wLEPS("wLEPS", "w-Lepton jigsaws");
            wLEPS.AddFrame(wLA);
            wLEPS.AddFrame(wLB);
            wLEPS.SetNElementsForFrame(wLA, 1, true); // in our decay tree implementation we have ==1 lepton per leg
            wLEPS.SetNElementsForFrame(wLB, 1, true); // in our decay tree implementation we have ==1 lepton per leg

            MinMassesCombJigsaw wLepHemiJigsaw("wHemi-lep", "w-Minimize m_{l,l} jigsaw");
            wLEPS.AddJigsaw(wLepHemiJigsaw);
            wLepHemiJigsaw.AddFrame(wLA, 0);
            wLepHemiJigsaw.AddFrame(wLB, 1);      

            // check that jigsaws are connected properly
            if(!wLAB.InitializeAnalysis()) {
                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
                exit(1);
            }

            // clear the event
            wLAB.ClearEvent();

            // set the met
            TVector3 tmet3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
            wINV.SetLabFrameThreeVector(tmet3vector);

            // build up our ish
            vector<TLorentzVector> wJETS, wSJETS, wBJETS, wLEPTONS;
            for(int i = 0; i < (int)jets.size(); i++) {
                Jet* j = jets[i];
                if(sl->tools->jetSelector().isB(j)) {
                    TLorentzVector b;
                    b.SetPtEtaPhiM(j->Pt(), j->Eta(), j->Phi(), j->M());
                    wBJETS.push_back(b);
                }
                else if(!(sl->tools->jetSelector().isB(j))) {
                    TLorentzVector sj;
                    sj.SetPtEtaPhiM(j->Pt(), j->Eta(), j->Phi(), j->M());
                    wSJETS.push_back(sj);
                }
            }
            for(int i = 0; i < (int)leptons.size(); i++) {
                TLorentzVector l;
                l.SetPtEtaPhiM(leptons[i]->Pt(), leptons[i]->Eta(), leptons[i]->Phi(), leptons[i]->M());
                wLEPTONS.push_back(l);
            }

            vector<RFKey> sjetID, bjetID;
            for(int i = 0; i < (int)wBJETS.size(); i++) {
                bjetID.push_back(wVIS.AddLabFrameFourVector(wBJETS[i]));
            }
            for(int i = 0; i < (int)wSJETS.size(); i++) {
                sjetID.push_back(wVIS.AddLabFrameFourVector(wSJETS[i]));
            }

            // add leptons
            wLEPS.AddLabFrameFourVector(wLEPTONS[0]);
            wLEPS.AddLabFrameFourVector(wLEPTONS[1]);

            // analyze
            wLAB.AnalyzeEvent();

            ////////////////////////////////////////////////////
            // Observables
            ////////////////////////////////////////////////////

            w_shat = wSS.GetMass();
            w_gaminv = wSS.GetVisibleShape();
            w_gaminv_SA = wSA.GetVisibleShape();

            TVector3 wvSS = wSS.GetFourVector(wLAB).Vect();
            w_pTCM = wvSS.Pt();
            w_pZCM = fabs(wvSS.Pz());
            w_RPT = w_pTCM / (w_pTCM + w_shat /4.);
            w_RPZ = fabs(w_pZCM) / (fabs(w_pZCM) + w_shat /4.);

            w_cosSS = wSS.GetCosDecayAngle();
            w_dphiLSS = wLAB.GetDeltaPhiDecayPlanes(wSS);
            if(w_dphiLSS > acos(-1)) w_dphiLSS = w_dphiLSS - 2*acos(-1);
            w_dphiSS_SA = wSS.GetDeltaPhiDecayPlanes(wSA);
            if(w_dphiSS_SA > acos(-1)) w_dphiSS_SA = w_dphiSS_SA - 2*acos(-1);
            w_dphiSS_SB = wSS.GetDeltaPhiDecayPlanes(wSB);
            if(w_dphiSS_SB > acos(-1)) w_dphiSS_SB = w_dphiSS_SB - 2*acos(-1);
            w_dphiSA_SB = wSA.GetDeltaPhiDecayPlanes(wSB);
            if(w_dphiSA_SB > acos(-1)) w_dphiSA_SB = w_dphiSA_SB - 2*acos(-1);

            TLorentzVector pSA = wSA.GetVisibleFourVector(wSA);
            TLorentzVector pSB = wSB.GetVisibleFourVector(wSB);
            w_MS = (pSA.M2() - pSB.M2()) / (2.*(pSA.E() - pSB.E()));

            // delta phi between SS visible decay products and SS axis
            w_DPD_vSS = wSS.GetDeltaPhiDecayVisible();
            // delta phi between SS visible decay products and SS boost
            w_DPB_vSS = wSS.GetDeltaPhiBoostVisible();

            // where does everything end up?
            w_N_VA1 = wVIS.GetNElementsInFrame(wVA1);
            w_N_VA2 = wVIS.GetNElementsInFrame(wVA2);
            w_N_VB1 = wVIS.GetNElementsInFrame(wVB1);
            w_N_VB2 = wVIS.GetNElementsInFrame(wVB2);

            w_N_VA = w_N_VA1 + w_N_VA2;
            w_N_VB = w_N_VB1 + w_N_VB2;
            w_N_V1 = w_N_VA1 + w_N_VB1;
            w_N_V2 = w_N_VA2 + w_N_VB2;

            w_Nsj_VA1 = 0;
            w_Nsj_VA2 = 0;
            w_Nsj_VB1 = 0;
            w_Nsj_VB2 = 0;


            int nS = sjetID.size();
            for(int i = 0; i < nS; i++) {
                if(wVIS.GetFrame(sjetID[i]) == wVA1) {
                    w_Nsj_VA1++;
                }
                else if(wVIS.GetFrame(sjetID[i]) == wVA2) {
                    w_Nsj_VA2++;
                }
                else if(wVIS.GetFrame(sjetID[i]) == wVB1) {
                    w_Nsj_VB1++;
                }
                else if(wVIS.GetFrame(sjetID[i]) == wVB2) {
                    w_Nsj_VB2++;
                }
            }

            w_Nsj_A  = w_Nsj_VA1 + w_Nsj_VA2;
            w_Nsj_B  = w_Nsj_VB1 + w_Nsj_VB2;
            w_Nsj_V1 = w_Nsj_VA1 + w_Nsj_VB1;
            w_Nsj_V2 = w_Nsj_VA2 + w_Nsj_VB2;

            w_Nbj_VA1 = 0;
            w_Nbj_VA2 = 0;
            w_Nbj_VB1 = 0;
            w_Nbj_VB2 = 0;

            int nB = bjetID.size();
            for(int i = 0; i < nB; i++) {
                if(wVIS.GetFrame(bjetID[i]) == wVA1) {
                    w_Nbj_VA1++;
                }
                else if(wVIS.GetFrame(bjetID[i]) == wVA2) {
                    w_Nbj_VA2++;
                }
                else if(wVIS.GetFrame(bjetID[i]) == wVB1) {
                    w_Nbj_VB1++;
                }
                else if(wVIS.GetFrame(bjetID[i]) == wVB2) {
                    w_Nbj_VB2++;
                }
            }

            w_Nbj_A  = w_Nbj_VA1 + w_Nbj_VA2;
            w_Nbj_B  = w_Nbj_VB1 + w_Nbj_VB2;
            w_Nbj_V1 = w_Nbj_VA1 + w_Nbj_VB1;
            w_Nbj_V2 = w_Nbj_VA2 + w_Nbj_VB2;


            ///////////////////////////////
            // Scale Variables (jets only)
            ///////////////////////////////

            // SS frame
            TLorentzVector vP_VA1_SS = wVA1.GetFourVector(wSS);
            TLorentzVector vP_VA2_SS = wVA2.GetFourVector(wSS);
            TLorentzVector vP_VB1_SS = wVB1.GetFourVector(wSS);
            TLorentzVector vP_VB2_SS = wVB2.GetFourVector(wSS);
            TLorentzVector vP_IA_SS  = wIA.GetFourVector(wSS);
            TLorentzVector vP_IB_SS  = wIB.GetFourVector(wSS);

            w_H_11_SS = (vP_VA1_SS + vP_VA2_SS + vP_VB1_SS + vP_VB2_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            w_H_21_SS = (vP_VA1_SS + vP_VA2_SS).P() + (vP_VB1_SS + vP_VB2_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            w_H_41_SS = vP_VA1_SS.P() + vP_VA2_SS.P() + vP_VB1_SS.P() + vP_VB2_SS.P() + (vP_IA_SS + vP_IB_SS).P();
            w_H_22_SS = (vP_VA1_SS + vP_VA2_SS).P() + (vP_VB1_SS + vP_VB2_SS).P() + vP_IA_SS.P() + vP_IB_SS.P();
            w_H_42_SS = vP_VA1_SS.P() + vP_VA2_SS.P() + vP_VB1_SS.P() + vP_VB2_SS.P() + vP_IA_SS.P() + vP_IB_SS.P();

            //   > transverse
            w_H_11_SS_T = (vP_VA1_SS + vP_VA2_SS + vP_VB1_SS + vP_VB2_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            w_H_21_SS_T = (vP_VA1_SS + vP_VA2_SS).Pt() + (vP_VB1_SS + vP_VB2_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            w_H_41_SS_T = vP_VA1_SS.Pt() + vP_VA2_SS.Pt() + vP_VB1_SS.Pt() + vP_VB2_SS.Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            w_H_22_SS_T = (vP_VA1_SS + vP_VA2_SS).Pt() + (vP_VB1_SS + vP_VB2_SS).Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt();
            w_H_42_SS_T = vP_VA1_SS.Pt() + vP_VA2_SS.Pt() + vP_VB1_SS.Pt() + vP_VB2_SS.Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt();

            // S1 frame
            TLorentzVector vP_VA1_SA = wVA1.GetFourVector(wSA);
            TLorentzVector vP_VA2_SA = wVA2.GetFourVector(wSA);
            TLorentzVector vP_IA_SA  = wIA.GetFourVector(wSA);

            w_H_11_SA = (vP_VA1_SA + vP_VA2_SA).P() + vP_IA_SA.P();
            w_H_21_SA = vP_VA1_SA.P() + vP_VA2_SA.P() + vP_IA_SA.P();

            //   > transverse
            w_H_11_SA_T = (vP_VA1_SA + vP_VA2_SA).Pt() + vP_IA_SA.Pt();
            w_H_21_SA_T = vP_VA1_SA.Pt() + vP_VA2_SA.Pt() + vP_IA_SA.Pt();


            ///////////////////////////////
            // Scale Variables (with leptons)
            ///////////////////////////////
            TLorentzVector vP_LA_SS = wLA.GetFourVector(wSS);
            TLorentzVector vP_LB_SS = wLB.GetFourVector(wSS);

            w_H_11_SS_lep = (vP_VA1_SS + vP_VA2_SS + vP_VB1_SS + vP_VB2_SS + vP_LA_SS + vP_LB_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            w_H_21_SS_lep = (vP_VA1_SS + vP_VA2_SS + vP_LA_SS).P() + (vP_VB1_SS + vP_VB2_SS + vP_LB_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            w_H_41_SS_lep = vP_VA1_SS.P() + (vP_VA2_SS + vP_LA_SS).P() + vP_VB1_SS.P() + (vP_VB2_SS + vP_LB_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            w_H_22_SS_lep = (vP_VA1_SS + vP_VA2_SS + vP_LA_SS).P() + (vP_VB1_SS + vP_VB2_SS + vP_LB_SS).P() + vP_IA_SS.P() + vP_IB_SS.P(); 
            w_H_42_SS_lep = vP_VA1_SS.P() + (vP_VA2_SS + vP_LA_SS).P() + vP_VB1_SS.P() + (vP_VB2_SS + vP_LB_SS).P() + vP_IA_SS.P() + vP_IB_SS.P();
            w_H_61_SS_lep = vP_VA1_SS.P() + vP_VA2_SS.P() + vP_LA_SS.P() + vP_VB1_SS.P() + vP_VB2_SS.P() + vP_LB_SS.P() + (vP_IA_SS + vP_IB_SS).P(); 
            w_H_62_SS_lep = vP_VA1_SS.P() + vP_VA2_SS.P() + vP_LA_SS.P() + vP_VB1_SS.P() + vP_VB2_SS.P() + vP_LB_SS.P() + vP_IA_SS.P() + vP_IB_SS.P(); 

            //  > transverse
            w_H_11_SS_lep_T = (vP_VA1_SS + vP_VA2_SS + vP_VB1_SS + vP_VB2_SS + vP_LA_SS + vP_LB_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            w_H_21_SS_lep_T = (vP_VA1_SS + vP_VA2_SS + vP_LA_SS).Pt() + (vP_VB1_SS + vP_VB2_SS + vP_LB_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            w_H_41_SS_lep_T = vP_VA1_SS.Pt() + (vP_VA2_SS + vP_LA_SS).Pt() + vP_VB1_SS.Pt() + (vP_VB2_SS + vP_LB_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            w_H_22_SS_lep_T = (vP_VA1_SS + vP_VA2_SS + vP_LA_SS).Pt() + (vP_VB1_SS + vP_VB2_SS + vP_LB_SS).Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt(); 
            w_H_42_SS_lep_T = vP_VA1_SS.Pt() + (vP_VA2_SS + vP_LA_SS).Pt() + vP_VB1_SS.Pt() + (vP_VB2_SS + vP_LB_SS).Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt();
            w_H_61_SS_lep_T = vP_VA1_SS.Pt() + vP_VA2_SS.Pt() + vP_LA_SS.Pt() + vP_VB1_SS.Pt() + vP_VB2_SS.Pt() + vP_LB_SS.Pt() + (vP_IA_SS + vP_IB_SS).Pt(); 
            w_H_62_SS_lep_T = vP_VA1_SS.Pt() + vP_VA2_SS.Pt() + vP_LA_SS.Pt() + vP_VB1_SS.Pt() + vP_VB2_SS.Pt() + vP_LB_SS.Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt(); 

            // S1 frame
            TLorentzVector vP_LA_SA = wLA.GetFourVector(wSA);

            w_H_11_SA_lep = (vP_VA1_SA + vP_VA2_SA + vP_LA_SA).P() + vP_IA_SA.P();
            w_H_21_SA_lep = (vP_VA2_SA + vP_LA_SA).P() + vP_VA1_SA.P();
            w_H_31_SA_lep = vP_VA1_SA.P() + vP_VA2_SA.P() + vP_LA_SA.P() + vP_IA_SA.P();

            //  > transverse
            w_H_11_SA_lep_T = (vP_VA1_SA + vP_VA2_SA + vP_LA_SA).Pt() + vP_IA_SA.Pt();
            w_H_21_SA_lep_T = (vP_VA2_SA + vP_LA_SA).Pt() + vP_VA1_SA.Pt();
            w_H_31_SA_lep_T = vP_VA1_SA.Pt() + vP_VA2_SA.Pt() + vP_LA_SA.Pt() + vP_IA_SA.Pt();

            //////////////////////////////////////
            // Energy variables
            //////////////////////////////////////

            w_MDR_VA1_SA = 2.0 * wVA1.GetEnergy(wSA);
            w_MDR_CLA_CA = 2.0 * wCLA.GetEnergy(wCA);
            


        } // 2 signal leptons and >=2 jets
    }; // end RESTFRAMES mT setup
    
    *cutflow << NewVar("w_shat"); {
        *cutflow << HFTname("w_shat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_shat;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_gaminv"); {
        *cutflow << HFTname("w_gaminv");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_gaminv;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_pTCM"); {
        *cutflow << HFTname("w_pTCM");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_pTCM;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_pZCM"); {
        *cutflow << HFTname("w_pZCM");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_pZCM;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_RPT"); {
        *cutflow << HFTname("w_RPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_RPT;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_RPZ"); {
        *cutflow << HFTname("w_RPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_RPZ;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_cosSS"); {
        *cutflow << HFTname("w_cosSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_cosSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_dphiLSS"); {
        *cutflow << HFTname("w_dphiLSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_dphiLSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_dphiSS_SA"); {
        *cutflow << HFTname("w_dphiSS_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_dphiSS_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_dphiSS_SB"); {
        *cutflow << HFTname("w_dphiSS_SB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_dphiSS_SB;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_dphiSA_SB"); {
        *cutflow << HFTname("w_dphiSA_SB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_dphiSA_SB;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_MS"); {
        *cutflow << HFTname("w_MS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_MS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_DPD_vSS"); {
        *cutflow << HFTname("w_DPD_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_DPD_vSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_DPB_vSS"); {
        *cutflow << HFTname("w_DPB_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_DPB_vSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_VA1"); {
        *cutflow << HFTname("w_N_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_VA2"); {
        *cutflow << HFTname("w_N_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_VB1"); {
        *cutflow << HFTname("w_N_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_VB2"); {
        *cutflow << HFTname("w_N_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_VA"); {
        *cutflow << HFTname("w_N_VA");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_VA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_VB"); {
        *cutflow << HFTname("w_N_VB");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_VB;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_V1"); {
        *cutflow << HFTname("w_N_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_N_V2"); {
        *cutflow << HFTname("w_N_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_N_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_VA1"); {
        *cutflow << HFTname("w_Nsj_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_VA2"); {
        *cutflow << HFTname("w_Nsj_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_VB1"); {
        *cutflow << HFTname("w_Nsj_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_VB2"); {
        *cutflow << HFTname("w_Nsj_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_VA1"); {
        *cutflow << HFTname("w_Nbj_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_VA2"); {
        *cutflow << HFTname("w_Nbj_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_VB1"); {
        *cutflow << HFTname("w_Nbj_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_VB2"); {
        *cutflow << HFTname("w_Nbj_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_A"); {
        *cutflow << HFTname("w_Nsj_A");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_A;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_B"); {
        *cutflow << HFTname("w_Nsj_B");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_B;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_V1"); {
        *cutflow << HFTname("w_Nsj_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nsj_V2"); {
        *cutflow << HFTname("w_Nsj_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nsj_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_A"); {
        *cutflow << HFTname("w_Nbj_A");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_A;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_B"); {
        *cutflow << HFTname("w_Nbj_B");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_B;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_V1"); {
        *cutflow << HFTname("w_Nbj_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_Nbj_V2"); {
        *cutflow << HFTname("w_Nbj_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_Nbj_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SS"); {
        *cutflow << HFTname("w_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SS"); {
        *cutflow << HFTname("w_H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_41_SS"); {
        *cutflow << HFTname("w_H_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_41_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_22_SS"); {
        *cutflow << HFTname("w_H_22_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_22_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_42_SS"); {
        *cutflow << HFTname("w_H_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_42_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SS_T"); {
        *cutflow << HFTname("w_H_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SS_T"); {
        *cutflow << HFTname("w_H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_41_SS_T"); {
        *cutflow << HFTname("w_H_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_41_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_22_SS_T"); {
        *cutflow << HFTname("w_H_22_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_22_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_42_SS_T"); {
        *cutflow << HFTname("w_H_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_42_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SA"); {
        *cutflow << HFTname("w_H_11_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SA"); {
        *cutflow << HFTname("w_H_21_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SA_T"); {
        *cutflow << HFTname("w_H_11_SA_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SA_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SA_T"); {
        *cutflow << HFTname("w_H_21_SA_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SA_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_MDR_VA1_SA"); {
        *cutflow << HFTname("w_MDR_VA1_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_MDR_VA1_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_gaminv_SA"); {
        *cutflow << HFTname("w_gaminv_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_gaminv_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SS_lep"); {
        *cutflow << HFTname("w_H_11_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SS_lep"); {
        *cutflow << HFTname("w_H_21_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_41_SS_lep"); {
        *cutflow << HFTname("w_H_41_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_41_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_22_SS_lep"); {
        *cutflow << HFTname("w_H_22_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_22_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_42_SS_lep"); {
        *cutflow << HFTname("w_H_42_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_42_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_61_SS_lep"); {
        *cutflow << HFTname("w_H_61_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_61_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_62_SS_lep"); {
        *cutflow << HFTname("w_H_62_SS_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_62_SS_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SS_lep_T"); {
        *cutflow << HFTname("w_H_11_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SS_lep_T"); {
        *cutflow << HFTname("w_H_21_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_41_SS_lep_T"); {
        *cutflow << HFTname("w_H_41_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_41_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_22_SS_lep_T"); {
        *cutflow << HFTname("w_H_22_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_22_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_42_SS_lep_T"); {
        *cutflow << HFTname("w_H_42_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_42_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_61_SS_lep_T"); {
        *cutflow << HFTname("w_H_61_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_61_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_62_SS_lep_T"); {
        *cutflow << HFTname("w_H_62_SS_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_62_SS_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SA_lep"); {
        *cutflow << HFTname("w_H_11_SA_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SA_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SA_lep"); {
        *cutflow << HFTname("w_H_21_SA_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SA_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_31_SA_lep"); {
        *cutflow << HFTname("w_H_31_SA_lep");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_31_SA_lep;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_11_SA_lep_T"); {
        *cutflow << HFTname("w_H_11_SA_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_11_SA_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_21_SA_lep_T"); {
        *cutflow << HFTname("w_H_21_SA_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_21_SA_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_H_31_SA_lep_T"); {
        *cutflow << HFTname("w_H_31_SA_lep_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_H_31_SA_lep_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("w_MDR_CLA_CA"); {
        *cutflow << HFTname("w_MDR_CLA_CA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = w_MDR_CLA_CA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 


*/

//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    //
//    // RESTFRAMES -- COMPRESSED
//    // RESTFRAMES -- COMPRESSED
//    // RESTFRAMES -- COMPRESSED
//    //
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////////////////////////
//
//    int cN_ISR;
//    int cN_V;
//    int cNj_ISR;
//    int cNj_V;
//    int cNb_ISR;
//    int cNb_V;
//    double cpT_ISR;
//    double cpT_I;
//    double cpT_S;
//    double cR_ISR;
//    double cR_ISR_S;
//    double cCosS;
//    double cCosS_I;
//    double cMS;
//    double cMV;
//    double cMISR;
//
//
//    *cutflow << [&](Superlink* sl, var_void*) {
//    ///////////////////////////////////////////
//        if(sjets.size()>=2 && bjets.size()>=1) {
//
//            // setup the analysis tree
//            LabRecoFrame        cLAB("cLAB", "cLAB");
//            DecayRecoFrame      cCM("cCM", "cCM");
//            DecayRecoFrame      cS("cS", "cS");
//            VisibleRecoFrame    cISR("cISR", "cISR");
//            VisibleRecoFrame    cV("cV", "cV");
//            DecayRecoFrame      cINVLL("cINVLL", "cINVLL");
//            DecayRecoFrame      cLL("cLL", "cLL");
//            InvisibleRecoFrame  cI("cI", "cI");
//            VisibleRecoFrame    cLEP1("cLEP1", "cLEP1");
//            VisibleRecoFrame    cLEP2("cLEP2", "cLEP2");
//
//            // connect the tree according to our interpretation
//            cLAB.SetChildFrame(cCM);
//            cCM.AddChildFrame(cISR);
//            cCM.AddChildFrame(cS);
//            cS.AddChildFrame(cV);
//            cS.AddChildFrame(cINVLL);
//            cINVLL.AddChildFrame(cLL);
//            cINVLL.AddChildFrame(cI);
//            cLL.AddChildFrame(cLEP1);
//            cLL.AddChildFrame(cLEP2);
//
//            // check that the decay tree is connected properly
//            if(!cLAB.InitializeTree()) {
//                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ << ")    Unable to initialize tree from lab frame. Exiting." << endl;
//                exit(1);
//            }
//
//            // define groups
//            InvisibleGroup cINV("cINV", "c-Invisible system");
//            cINV.AddFrame(cI);
//
//            // combinatoric for all of our jets (sjets & bjets)
//            // put loosest restrictions on NElements (can always cut on these)
//            CombinatoricGroup cVIS("cVIS", "c-Visible jets");
//            cVIS.AddFrame(cISR);
//            cVIS.SetNElementsForFrame(cISR,1,false);
//            cVIS.AddFrame(cV);
//            cVIS.SetNElementsForFrame(cV,0,false);
//
//            // set the invisible system mass to zero
//            SetMassInvJigsaw cInvMass("cInvMass", "c-Invisible system mass to zero jigsaw");
//            cINV.AddJigsaw(cInvMass);
//           
//            // define the rule to partition objects between "ISR" and "V"
//            MinMassesCombJigsaw cSplitVis("cSplitVis", "Minimize M_{ISR} and M_{V} Jigsaw"); 
//            cVIS.AddJigsaw(cSplitVis);
//            // "0" group (ISR)
//            cSplitVis.AddFrame(cISR, 0);
//            // "1" group (V + I + LL)
//            cSplitVis.AddFrame(cV,1);
//            cSplitVis.AddFrame(cI,1);
//            cSplitVis.AddFrame(cLL,1);
//
//            // check that the jigsaws are in place
//            if(!cLAB.InitializeAnalysis()) {
//                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
//                exit(1);
//            }
//
//            // clear the event for sho
//            cLAB.ClearEvent();
//
//            // set the met
//            TVector3 cMET(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
//
//            // vector of sjets
//            vector<TLorentzVector> sJETS;
//            for(int ij = 0; ij < (int)sjets.size(); ij++) {
//                TLorentzVector j;
//                j.SetPtEtaPhiM(sjets.at(ij)->Pt(), sjets.at(ij)->Eta(), sjets.at(ij)->Phi(), sjets.at(ij)->M());
//                sJETS.push_back(j);
//            }
//            // vector of bjets
//            vector<TLorentzVector> bJETS;
//            for(int ib = 0; ib < (int)bjets.size(); ib++) {
//                TLorentzVector b;
//                b.SetPtEtaPhiM(bjets.at(ib)->Pt(), bjets.at(ib)->Eta(), bjets.at(ib)->Phi(), bjets.at(ib)->M());
//                bJETS.push_back(b);
//            }
//
//            // add the jets (in transverse view) to the combinatoric group and index each type
//            vector<RFKey> sjetID;
//            vector<RFKey> bjetID;
//            for(int ij = 0; ij < (int)sJETS.size(); ij++) {
//                TLorentzVector jet_ = sJETS[ij];
//                jet_.SetPtEtaPhiM(jet_.Pt(), 0.0, jet_.Phi(), jet_.M());
//                sjetID.push_back(cVIS.AddLabFrameFourVector(jet_));
//            }
//            for(int ib = 0; ib < (int)bJETS.size(); ib++) {
//                TLorentzVector jet_ = bJETS[ib];
//                jet_.SetPtEtaPhiM(jet_.Pt(), 0.0, jet_.Phi(), jet_.M());
//                bjetID.push_back(cVIS.AddLabFrameFourVector(jet_));
//            }
//
//            // add the met vector to the invisible system
//            cINV.SetLabFrameThreeVector(cMET);
//
//            // add the leptons (in transverse view) to their frames
//            TLorentzVector lep1, lep2;
//            lep1.SetPtEtaPhiM(leptons.at(0)->Pt(), 0.0, leptons.at(0)->Phi(), leptons.at(0)->M());
//            lep2.SetPtEtaPhiM(leptons.at(1)->Pt(), 0.0, leptons.at(1)->Phi(), leptons.at(1)->M());
//            cLEP1.SetLabFrameFourVector(lep1);
//            cLEP2.SetLabFrameFourVector(lep2);
//            
//
//            // analyze the event
//            if(!cLAB.AnalyzeEvent()) {
//                cout << analysis_name << "    RestFrames::AnalyzeEvent ERROR (" << __LINE__ << ")    Unable to analyze event from lab frame. Exiting." << endl;
//                exit(1);
//            }
//
//            cN_ISR = cVIS.GetNElementsInFrame(cISR);
//            cN_V   = cVIS.GetNElementsInFrame(cV);
//            cNj_ISR = 0;
//            cNj_V = 0;
//            cNb_ISR = 0;
//            cNb_V = 0;
//
//
//            int nS = sjetID.size();
//            for(int i = 0; i < nS; i++) {
//                if(cVIS.GetFrame(sjetID[i]) == cISR) {
//                    cNj_ISR++;
//                }
//                else if(cVIS.GetFrame(sjetID[i]) == cV) {
//                    cNj_V++;
//                }
//            }
//            int nB = bjetID.size();
//            for(int i = 0; i < nB; i++) {
//                if(cVIS.GetFrame(bjetID[i]) == cISR) {
//                    cNb_ISR++;
//                }
//                else if(cVIS.GetFrame(bjetID[i]) == cV) {
//                    cNb_V++;
//                }
//            }
//
//            TVector3 vP_ISR = cISR.GetFourVector(cCM).Vect();
//            TVector3 vP_S = cS.GetFourVector(cCM).Vect();
//            TVector3 vP_I   = cI.GetFourVector(cCM).Vect();
//
//            cpT_ISR = vP_ISR.Mag();
//            cpT_S = vP_S.Mag();
//            cpT_I = vP_I.Mag();
//            cR_ISR = fabs(vP_I.Dot(vP_ISR.Unit())) / cpT_ISR;
//            cR_ISR_S = fabs(vP_I.Dot(vP_S.Unit())) / cpT_S;
//
//            cCosS = cS.GetCosDecayAngle(); // decay angle relative to cINVLL decay
//            cCosS_I = cS.GetCosDecayAngle(cI); // decay angle relative to the invisible system (a la two-body decay to cINV and cLL)
//
//            cMS = cS.GetMass();
//            cMV = cV.GetMass();
//            cMISR = cISR.GetMass();
//        } // sjets >= 2 && bjest >= 1
//
//    }; // RESTFRAMES -- COMPRESSED [END]
//
//    *cutflow << NewVar("cN_ISR -- number of visible objects in the cISR frame"); {
//        *cutflow << HFTname("cN_ISR");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            int out = -1;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cN_ISR;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cN_V -- number of visible objects in the cV frame"); {
//        *cutflow << HFTname("cN_V");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            int out = -1;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cN_V;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cNj_ISR -- number of sjets in the cISR frame"); {
//        *cutflow << HFTname("cNj_ISR");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            int out = -1;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cNj_ISR;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cNj_V -- number of sjets in the cV frame"); {
//        *cutflow << HFTname("cNj_V");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            int out = -1;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cNj_V;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cNb_ISR -- number of bjets in the cISR frame"); {
//        *cutflow << HFTname("cNb_ISR");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            int out = -1;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cNb_ISR;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cNb_V -- number of bjets in the cV frame"); {
//        *cutflow << HFTname("cNb_V");
//        *cutflow << [&](Superlink* sl, var_int*) -> int {
//            int out = -1;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cNb_V;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cpT_ISR -- pT of cISR frame in the cCM restframe"); {
//        *cutflow << HFTname("cpT_ISR");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -10.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cpT_ISR;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cpT_I -- pT of cI frame in the cCM restframe"); {
//        *cutflow << HFTname("cpT_I");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -10.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cpT_I;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cpT_S -- pT of cS frame in the cCM restframe"); {
//        *cutflow << HFTname("cpT_S");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -10.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cpT_S;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cR_ISR -- ratio of inv. system momentum to ISR momentum"); {
//        *cutflow << HFTname("cR_ISR");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -5.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cR_ISR;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cR_ISR_S -- ratio of inv. system momtnum to S momentum"); {
//        *cutflow << HFTname("cR_ISR_S");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -5.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cR_ISR_S;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cCosS -- cosine of decay of cS relative to the cINVLL decay frame"); {
//        *cutflow << HFTname("cCosS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -5.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cCosS;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cCosS_I -- cosine of decay of cS relative to the cI frame"); {
//        *cutflow << HFTname("cCosS_I");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -5.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cCosS_I;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cMS -- mass of the sparticle decay system cS"); {
//        *cutflow << HFTname("cMS");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -10.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cMS;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cMV -- mass of the cV system"); {
//        *cutflow << HFTname("cMV");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -10.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cMV;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }
//    *cutflow << NewVar("cMISR -- mass of the cISR system"); {
//        *cutflow << HFTname("cMISR");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            double out = -10.;
//            if(sjets.size()>=2 && bjets.size()>=1) {
//                out = cMISR;
//            }
//            return out;
//        };
//        *cutflow << SaveVar();
//    }






    // clear the wectors
    *cutflow << [&](Superlink* sl, var_void*) { leptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { electrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { muons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { jets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { met.clear(); };


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Sysystematics [BEGIN]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////
    // weight systematics
    ////////////////////////////////////

    // electron eff
    *cutflow << NewSystematic("shift in electron ID efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ID_UP, SupersysWeight::EL_EFF_ID_DN);
        *cutflow << TreeName("EL_EFF_ID");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron ISO efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ISO_UP, SupersysWeight::EL_EFF_ISO_DN);
        *cutflow << TreeName("EL_EFF_Iso");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron RECO efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_RECO_UP, SupersysWeight::EL_EFF_RECO_DN);
        *cutflow << TreeName("EL_EFF_Reco");
        *cutflow << SaveSystematic();
    }

    // muon eff
    *cutflow << NewSystematic("muon eff stat uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_UP, SupersysWeight::MUON_EFF_STAT_DN);
        *cutflow << TreeName("MUON_EFF_STAT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff stat uncertainty (low pt)"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_LOWPT_UP, SupersysWeight::MUON_EFF_STAT_LOWPT_DN);
        *cutflow << TreeName("MUON_EFF_STAT_LOWPT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff syst uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYS_UP, SupersysWeight::MUON_EFF_SYS_DN);
        *cutflow << TreeName("MUON_EFF_SYS");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff syst uncertainty (low pt"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYS_LOWPT_UP, SupersysWeight::MUON_EFF_SYS_LOWPT_DN);
        *cutflow << TreeName("MUON_EFF_SYS_LOWPT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff iso stat uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_STAT_UP, SupersysWeight::MUON_EFF_ISO_STAT_DN);
        *cutflow << TreeName("MUON_ISO_STAT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff iso syst uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_SYS_UP, SupersysWeight::MUON_EFF_ISO_SYS_DN);
        *cutflow << TreeName("MUON_ISO_SYS");
        *cutflow << SaveSystematic();
    }

    // flavor tagging eff
    *cutflow << NewSystematic("shift in b-tag efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_B_UP, SupersysWeight::FT_EFF_B_DN);
        *cutflow << TreeName("FT_EFF_B");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in c-tag efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_C_UP, SupersysWeight::FT_EFF_C_DN);
        *cutflow << TreeName("FT_EFF_C");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in light tag (i.e. mis-tag) efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_LT_UP, SupersysWeight::FT_EFF_LT_DN);
        *cutflow << TreeName("FT_EFF_Light");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift flavor tagging extrapolation ?"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAP_UP, SupersysWeight::FT_EFF_EXTRAP_DN);
        *cutflow << TreeName("FT_EFF_extrapolation");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift flavor tagging extrapolation - charm ?"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAPC_UP, SupersysWeight::FT_EFF_EXTRAPC_DN);
        *cutflow << TreeName("FT_EFF_extrapolation_charm");
        *cutflow << SaveSystematic();
    }

    // jvt eff
    *cutflow << NewSystematic("shift in JVT efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::JVT_EFF_UP, SupersysWeight::JVT_EFF_DN);
        *cutflow << TreeName("JET_JVTEff");
        *cutflow << SaveSystematic();
    }

    // pileup
    *cutflow << NewSystematic("shift in data mu (pile-up)"); {
        *cutflow << WeightSystematic(SupersysWeight::PILEUP_UP, SupersysWeight::PILEUP_DN);
        *cutflow << TreeName("PILEUP");
        *cutflow << SaveSystematic();
    }


    ////////////////////////////////////
    // shape systematics
    ////////////////////////////////////

    // egamma
    *cutflow << NewSystematic("shift in e-gamma resolution (UP)"); {
        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_UP);
        *cutflow << TreeName("EG_RESOLUTION_ALL_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma resolution (DOWN)"); {
        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_DN);
        *cutflow << TreeName("EG_RESOLUTION_ALL_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma scale (UP)"); {
        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_UP);
        *cutflow << TreeName("EG_SCALE_ALL_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma scale (DOWN)"); {
        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_DN);
        *cutflow << TreeName("EG_SCALE_ALL_DN");
        *cutflow << SaveSystematic();
    }
    // muon
    *cutflow << NewSystematic("muon ID (UP)"); {
        *cutflow << EventSystematic(NtSys::MUON_ID_UP);
        *cutflow << TreeName("MUON_ID_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon ID (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUON_ID_DN);
        *cutflow << TreeName("MUON_ID_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (UP)"); {
        *cutflow << EventSystematic(NtSys::MUON_MS_UP);
        *cutflow << TreeName("MUON_MS_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUON_MS_DN);
        *cutflow << TreeName("MUON_MS_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (UP)"); {
        *cutflow << EventSystematic(NtSys::MUON_SCALE_UP);
        *cutflow << TreeName("MUON_SCALE_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (DN)"); {
        *cutflow << EventSystematic(NtSys::MUON_SCALE_DN);
        *cutflow << TreeName("MUON_SCALE_DN");
        *cutflow << SaveSystematic();
    }

    // jet
    *cutflow << NewSystematic("JER"); {
        *cutflow << EventSystematic(NtSys::JER);
        *cutflow << TreeName("JER");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 1 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_UP);
        *cutflow << TreeName("JET_GroupedNP_1_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 1 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_DN);
        *cutflow << TreeName("JET_GroupedNP_1_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 2 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_2_UP);
        *cutflow << TreeName("JET_GroupedNP_2_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 2 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_2_DN);
        *cutflow << TreeName("JET_GroupedNP_2_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 3 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_3_UP);
        *cutflow << TreeName("JET_GroupedNP_3_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 3 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_3_DN);
        *cutflow << TreeName("JET_GroupedNP_3_DN");
        *cutflow << SaveSystematic();
    }

    // met
 //   *cutflow << NewSystematic("MET TST Soft-Term resolution (parallel)"); {
 //       *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPara);
 //       *cutflow << TreeName("MET_SoftTrk_ResoPara");
 //       *cutflow << SaveSystematic();
 //   }
 //   *cutflow << NewSystematic("MET TST Soft-Term resolution (perpendicular)"); {
 //       *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPerp);
 //       *cutflow << TreeName("MET_SoftTrk_ResoPerp");
 //       *cutflow << SaveSystematic();
 //   }
 //   *cutflow << NewSystematic("MET TST Soft-Term shift in scale (UP)"); {
 //       *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleUp);
 //       *cutflow << TreeName("MET_SoftTrk_ScaleUp");
 //       *cutflow << SaveSystematic();
 //   }
 //   *cutflow << NewSystematic("MET TST Soft-Term shift in scale (DOWN)"); {
 //       *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleDown);
 //       *cutflow << TreeName("MET_SoftTrk_ScaleDown");
 //       *cutflow << SaveSystematic();
 //   }


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Superflow methods [END]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // initialize the cutflow and start the event loop
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);
    delete cutflow;
    delete chain;
    cout << "La Fin." << endl;
    exit(0);

}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
    string& suffix_name_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sys, bool& do_sumw_split, bool dbg)
{
    bool nominal_ = false;
    bool nominal_and_weight_sys_ = false;
    bool all_sys_ = false;
    string input;

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-n") == 0)
            num_events_ = atoi(argv[++i]);
        else if(strcmp(argv[i], "-i") == 0) {
            input = argv[++i];
        }
        else if(strcmp(argv[i], "-c") == 0) {
            nominal_ = true;
        }
        else if(strcmp(argv[i], "-w") == 0) {
            nominal_and_weight_sys_ = true;
        }
        else if(strcmp(argv[i], "-a") == 0) {
            all_sys_ = true;
        }
        else if(strcmp(argv[i], "-d") == 0) {
            dbg = true;
        }
        else if(strcmp(argv[i], "--suffix") == 0) {
            suffix_name_ = argv[++i];
        }
        else if(strcmp(argv[i], "--sumw") == 0) {
            do_sumw_split = true;
        }
        else {
            cout << analysis_name << "    Error (fatal) : Bad arguments." << endl;
            exit(1);
        }
    } // i

    sample_ = input;

    cout << analysis_name << "    input: " << input << endl;
    cout << analysis_name << "    input: " << input << endl;
    cout << analysis_name << "    input: " << input << endl;
    ChainHelper::addInput(chain, input, dbg); 

    Long64_t tot_num_events = chain->GetEntries();
    num_events_ = (num_events_ < 0 ? tot_num_events : num_events_);
    if(nominal_) {
        run_mode_ = SuperflowRunMode::nominal;
        cout << analysis_name << "    run mode: SuperflowRunMode::nominal" << endl;
    }
    if(nominal_and_weight_sys_) {
        run_mode_ = SuperflowRunMode::nominal_and_weight_syst;
        cout << analysis_name << "    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if(all_sys_) {
        run_mode_ = SuperflowRunMode::all_syst;
        cout << analysis_name << "    run mode: SuperflowRunMode::all_syst" << endl;
    }
    if(suffix_name_!="") {
        cout << analysis_name << "    Setting output file suffix to: " << suffix_name_ << endl;
    }
}

