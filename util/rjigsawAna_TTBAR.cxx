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
    cutflow->setLumi(3209); // 3.34/fb
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode_);
    cutflow->setCountWeights(true);
    cutflow->setChain(chain);
    if(suffix_name_!="")
        cutflow->setFileSuffix(suffix_name_);
    if(do_sumw_split) {
        //string sumw_file = "/data/uclhc/uci/user/dantrim/n0228val/sumw_file.txt"; 
        string sumw_file = "./n0228val/sumw_file.txt";
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

    
    *cutflow << NewVar("jet jvt"); {
        *cutflow << HFTname("j_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->jvt);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet jvt"); {
        *cutflow << HFTname("sj_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->jvt);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet jvt"); {
        *cutflow << HFTname("bj_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->jvt);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet nTrk"); {
        *cutflow << HFTname("j_nTrk");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->nTracks);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet nTrk"); {
        *cutflow << HFTname("sj_nTrk");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->nTracks);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet nTrk"); {
        *cutflow << HFTname("bj_nTrk");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->nTracks);
            }
            return out;
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
/*
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
*/
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    //  RESTFRAMES TOP
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    double Mtt;
    double pT_TT;
    double Mt;
    double Eb_Ta;
    double Eb_Tb;
    double M_top;
    double M_bl;
    double E_bl;

    *cutflow << [&](Superlink* sl, var_void*) {
        if(jets.size()>=2 && bjets.size()>0) {

            LabRecoFrame LAB("LAB", "LAB");
            DecayRecoFrame TT("TT", "TT");
            DecayRecoFrame Ta("Ta", "T_{a}");
            DecayRecoFrame Tb("Tb", "T_{b}");
            DecayRecoFrame Wa("Wa", "W_{a}");
            DecayRecoFrame Wb("Wb", "W_{b}");
            VisibleRecoFrame Ba("Ba", "b_{a}");
            VisibleRecoFrame La("La", "l_{a}");
            InvisibleRecoFrame Na("Na", "n_{a}");
            VisibleRecoFrame Bb("Bb", "b_{b}");
            VisibleRecoFrame Lb("Lb", "l_{b}");
            InvisibleRecoFrame Nb("Nb", "n_{b}");

            LAB.SetChildFrame(TT);
            TT.AddChildFrame(Ta);
            TT.AddChildFrame(Tb);
            Ta.AddChildFrame(Ba);
            Ta.AddChildFrame(Wa);
            Tb.AddChildFrame(Bb);
            Tb.AddChildFrame(Wb);
            Wa.AddChildFrame(La);
            Wa.AddChildFrame(Na);
            Wb.AddChildFrame(Lb);
            Wb.AddChildFrame(Nb);

            if(!(LAB.InitializeTree())) {
                cout << "failed to initailize tree" << endl;
                exit(1);
            }

            InvisibleGroup INV("INV", "INV");
            INV.AddFrame(Na);
            INV.AddFrame(Nb);
            // combinatoric group for b's
            CombinatoricGroup BVIS("BVIS", "BVIS");
            BVIS.AddFrame(Ba);
            BVIS.AddFrame(Bb);
            // bjet frames must have at least 1 element
            BVIS.SetNElementsForFrame(Ba,1);
            BVIS.SetNElementsForFrame(Bb,1);

            // minimize equal top masses neutrino jigsaws
            SetMassInvJigsaw InvMassJigsaw("INVMass", "INVMass");
            INV.AddJigsaw(InvMassJigsaw);
            SetRapidityInvJigsaw InvRapidityJigsaw("INVRapidity", "INVRapidity");
            INV.AddJigsaw(InvRapidityJigsaw);
            InvRapidityJigsaw.AddVisibleFrames(La+Ba+Lb+Bb);

            ContraBoostInvJigsaw MinMt("MinMt", "MinMt");
            INV.AddJigsaw(MinMt);
            MinMt.AddVisibleFrames(La+Ba,0);
            MinMt.AddVisibleFrames(Lb+Bb,1);
            MinMt.AddInvisibleFrame(Na, 0);
            MinMt.AddInvisibleFrame(Nb,1);


            // bjet combinatoric jigsaws
            MinMassesCombJigsaw MinBL("MinBL", "MinBL");
            BVIS.AddJigsaw(MinBL);
            MinBL.AddFrames(La+Ba,0);
            MinBL.AddFrames(Lb+Bb,1);

            if(!(LAB.InitializeAnalysis())) {
                cout << "failed to initialize analysis" << endl;
                exit(1);
            }


            LAB.ClearEvent();

            TVector3 met3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
            INV.SetLabFrameThreeVector(met3vector);

            La.SetLabFrameFourVector(*sl->leptons->at(0));
            Lb.SetLabFrameFourVector(*sl->leptons->at(1));
            for(int ijet=0; ijet < (int)sl->jets->size(); ijet++) {
                BVIS.AddLabFrameFourVector(*sl->jets->at(ijet));
            }

            LAB.AnalyzeEvent();

            Mtt = TT.GetMass();
            pT_TT = TT.GetFourVector(LAB).Pt();
            Mt = Ta.GetMass()/2;
            Eb_Ta = Ba.GetFourVector(Ta).E();
            Eb_Tb = Bb.GetFourVector(Tb).E();
            M_top = 0.5*(Ba.GetFourVector(LAB) + Na.GetFourVector(LAB) + La.GetFourVector(LAB)).M();
            M_bl = 0.5*(Ba.GetFourVector(LAB) + La.GetFourVector(LAB)).M();
            E_bl = (Ba.GetFourVector(Ta)+La.GetFourVector(Ta)).E();



/*            // setup the analysis tree
            LabRecoFrame LAB("lab", "lab"); 
            DecayRecoFrame TT("tt", "tt");
            DecayRecoFrame T1("t1", "t1");
            DecayRecoFrame T2("t2", "t2");

            // T1 side
            DecayRecoFrame W1("w1", "w1");
            VisibleRecoFrame L1("l1", "l1");
            InvisibleRecoFrame N1("n1", "n1");
            VisibleRecoFrame B1("b1", "b1");

            // T2 side
            DecayRecoFrame W2("w2", "w2");
            VisibleRecoFrame L2("l2", "l2");
            InvisibleRecoFrame N2("n2", "n2");
            VisibleRecoFrame B2("b2", "b2");

            // connect the tree parts
            LAB.SetChildFrame(TT);
            TT.AddChildFrame(T1);
            T1.AddChildFrame(B1);
            T1.AddChildFrame(W1);
            W1.AddChildFrame(L1);
            W1.AddChildFrame(N1);

            TT.AddChildFrame(T2);
            T2.AddChildFrame(B2);
            T2.AddChildFrame(W2);
            W2.AddChildFrame(L2);
            W2.AddChildFrame(N2);

            // check that things are OK
            if(!LAB.InitializeTree()) {
                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ << ")    Unable to initialize tree from lab frame. Exiting." << endl;
                exit(1);
            }

            // define groups
            InvisibleGroup INV("INV", "neutrino invisible jigsaws");
            INV.AddFrame(N1);
            INV.AddFrame(N2);

            CombinatoricGroup VIS("VIS", "visible object jigsaws");
            VIS.AddFrame(B1);
            VIS.AddFrame(B2);
            VIS.AddFrame(L1);
            VIS.AddFrame(L2);

            VIS.SetNElementsForFrame(L1, 1, true);
            VIS.SetNElementsForFrame(L2, 1, true);
            VIS.SetNElementsForFrame(B1, 1, false);
            VIS.SetNElementsForFrame(B2, 1, false);

            // define jigsaws
            SetMassInvJigsaw MinMassJigsaw("minmass", "invisible system mass jigsaw");
            INV.AddJigsaw(MinMassJigsaw);
            SetRapidityInvJigsaw RapidityJigsaw("rapidity", "invisible system rapidity jigsaw");
            INV.AddJigsaw(RapidityJigsaw);
            RapidityJigsaw.AddVisibleFrames(LAB.GetListVisibleFrames());
            ContraBoostInvJigsaw ContraBoostJigsaw("contra", "contraboost invariant jigsaw");
            INV.AddJigsaw(ContraBoostJigsaw);
            ContraBoostJigsaw.AddVisibleFrames(T1.GetListVisibleFrames(), 0);
            ContraBoostJigsaw.AddVisibleFrames(T2.GetListVisibleFrames(), 1);
            ContraBoostJigsaw.AddInvisibleFrames(T1.GetListInvisibleFrames(), 0);
            ContraBoostJigsaw.AddInvisibleFrames(T2.GetListInvisibleFrames(), 1);

            // group the leptons and (b-)jets to appropriate hemisphere
            MinMassesCombJigsaw HemiJigsaw("hemi", "minimize m_{V_{1,2}} jigsaw");
            VIS.AddJigsaw(HemiJigsaw);
            HemiJigsaw.AddFrame(B1, 0);
            HemiJigsaw.AddFrame(B2, 1);
            HemiJigsaw.AddFrame(L1, 0);
            HemiJigsaw.AddFrame(L2, 1);

            // define rule to partition objects between B and L
            MinMassesCombJigsaw T1SplitVis("T1splitvis", "Minimize M_{B1,L1} jigsaw");
            VIS.AddJigsaw(T1SplitVis);
            T1SplitVis.AddFrame(B1,0);
            T1SplitVis.AddFrame(N1,1);
            T1SplitVis.AddFrame(L1,1);

            MinMassesCombJigsaw T2SplitVis("T2splitvis", "Minimize M_{B2,L2} jigsaw");
            VIS.AddJigsaw(T2SplitVis);
            T2SplitVis.AddFrame(B2,0);
            T2SplitVis.AddFrame(N2,1);
            T2SplitVis.AddFrame(L2,1);


            // check that jigsaws are correctly in place
            if(!LAB.InitializeAnalysis()) {
                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
                exit(1);
            }

            // clear the envet
            LAB.ClearEvent();

            // set the met
            TVector3 met3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
            INV.SetLabFrameThreeVector(met3vector);

*/
        }
    };

    *cutflow << NewVar("Mtt"); {
        *cutflow << HFTname("Mtt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return Mtt;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pT_TT"); {
        *cutflow << HFTname("pT_TT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return pT_TT;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Mt"); {
        *cutflow << HFTname("Mt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return Mt;
        };
        *cutflow << SaveVar();  
    }
    *cutflow << NewVar("Eb_Ta"); {
        *cutflow << HFTname("Eb_Ta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return Eb_Ta;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Eb_Tb"); {
        *cutflow << HFTname("Eb_Tb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return Eb_Tb;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("M_top"); {
        *cutflow << HFTname("M_top");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return M_top;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("M_bl"); {
        *cutflow << HFTname("M_bl");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return M_bl;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("E_bl"); {
        *cutflow << HFTname("E_bl");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return E_bl;
        };
        *cutflow << SaveVar();
    }




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
        *cutflow << EventSystematic(NtSys::MUONS_ID_UP);
        *cutflow << TreeName("MUONS_ID_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon ID (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUONS_ID_DN);
        *cutflow << TreeName("MUONS_ID_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (UP)"); {
        *cutflow << EventSystematic(NtSys::MUONS_MS_UP);
        *cutflow << TreeName("MUONS_MS_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUONS_MS_DN);
        *cutflow << TreeName("MUONS_MS_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (UP)"); {
        *cutflow << EventSystematic(NtSys::MUONS_SCALE_UP);
        *cutflow << TreeName("MUONS_SCALE_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (DN)"); {
        *cutflow << EventSystematic(NtSys::MUONS_SCALE_DN);
        *cutflow << TreeName("MUONS_SCALE_DN");
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

