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

const string analysis_name = "rjigsawAna";

////////////////////////////////////////
// Function Prototypes
////////////////////////////////////////
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sys_, bool dbg);


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
    SuperflowRunMode run_mode_ = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode_, nt_sys_, m_dbg);

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
        cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
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

    ///////////////////////////////////////////////////
    // Analysis Cuts
    ///////////////////////////////////////////////////
    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        return sl->leptons->size() == 2;
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
    };
    

    ///////////////////////////////////////////////////
    // Ntuple Setup
    ///////////////////////////////////////////////////

    // standard variables
    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
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
    *cutflow << NewVar("delta phi between lead lepton and met"); {
        *cutflow << HFTname("dphi_met_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector l0;
            l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            return met.lv().DeltaPhi(l0);
        };
        *cutflow << SaveVar();
    }

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
        dphi_v1_i1_ss = v1.GetFourVector(ss).DeltaPhi(i1.GetFourVector(ss));
        dphi_s1_s2_ss = s1.GetFourVector(ss).DeltaPhi(s2.GetFourVector(ss)); 


        dphiS_I_ss = s1.GetFourVector(ss).DeltaPhi(i1.GetFourVector(ss));
        dphiS_I_s1 = s1.GetFourVector(ss).DeltaPhi(i1.GetFourVector(s1));
        
        

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
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //
    // RESTFRAMES - 2 jets, 2 leptons
    // RESTFRAMES - 2 jets, 2 leptons
    // RESTFRAMES - 2 jets, 2 leptons
    // RESTFRAMES - 2 jets, 2 leptons
    //
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    double xshat;
    double xgaminv;
    double xRPT;
    double xRPZ;
    double xcosSS;
    double xdphiLSS;
    double xMS;
    double xPS;
    double xMSS;
    double xgaminvSS;
    double xDeltaBetaSS;
    double xDPD_vSS;
    double xDPB_vSS;
    int xNV[2]; // number of visible objects in hemisphere
    double xcosS[2]; // cosine stop decay angle
    double xcosC[2]; // cosine intermediate child decay angle
    double xdphiSC[2]; // cosine between stop and child decay planes
    double xRCS[2]; // ratio of child and stop masses (w/ WIMP masses subtracted);
    double xjet1PT[2]; // first leading jet pT associated with this hemisphere
    double xjet2PT[2]; // second leading jet pT associated with this hemisphere
    double xPinv[2]; // Pinv / HS
    double xH_11_SS;
    double xH_21_SS;
    double xH_41_SS;
    double xH_42_SS;
    double xH_11_S1;
    double xH_21_S1;
    double xH_11_SS_T;
    double xH_21_SS_T;
    double xH_41_SS_T;
    double xH_42_SS_T;
    double xH_11_S1_T;
    double xH_21_S1_T;

    double xRPT_H_11_SS; 
    double xRPT_H_21_SS; 
    double xRPT_H_41_SS; 
    double xRPT_H_42_SS; 
    double xRPZ_H_11_SS; 
    double xRPZ_H_21_SS; 
    double xRPZ_H_41_SS; 
    double xRPZ_H_42_SS; 
    double xRPT_H_11_SS_T; 
    double xRPT_H_21_SS_T; 
    double xRPT_H_41_SS_T; 
    double xRPT_H_42_SS_T; 
    double xRPZ_H_11_SS_T; 
    double xRPZ_H_21_SS_T; 
    double xRPZ_H_41_SS_T; 
    double xRPZ_H_42_SS_T; 

    double xdphiVS_I[2];

    // feb 1 vars
    //double xCosP1;
    //double xCosP2;


    *cutflow << [&](Superlink* sl, var_void*) {
        if(sjets.size()>=2 && bjets.size()>=1) {

            // setup the analysis tree
            LabRecoFrame xlab("xlab", "xlab");
            DecayRecoFrame xss("xss", "xss");
            DecayRecoFrame xs1("xs1", "xs1");
            DecayRecoFrame xs2("xs2", "xs2");
            DecayRecoFrame xc1("xc1", "xc1");
            DecayRecoFrame xc2("xc2", "xc2");
            VisibleRecoFrame xv1s("xv1s", "xv1s");
            VisibleRecoFrame xv2s("xv2s", "xv2s");
            InvisibleRecoFrame xi1("xi1", "xi1");
            VisibleRecoFrame xv1c("xv1c", "xv1c");
            VisibleRecoFrame xv2c("xv2c", "xv2c");
            InvisibleRecoFrame xi2("xi2", "xi2");

            xlab.SetChildFrame(xss);
            xss.AddChildFrame(xs1);
            xss.AddChildFrame(xs2);
            xs1.AddChildFrame(xc1);
            xs1.AddChildFrame(xv1s);
            xc1.AddChildFrame(xi1);
            xc1.AddChildFrame(xv1c);
            xs2.AddChildFrame(xc2);
            xs2.AddChildFrame(xv2s);
            xc2.AddChildFrame(xi2);
            xc2.AddChildFrame(xv2c);

            // check that the decay tree is connected properly
            if(!xlab.InitializeTree()) {
                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ <<")    Unable to initialize tree from lab frame. Exiting." << endl;
                exit(1);
            }


            // define groupes

            InvisibleGroup xinv("xinv", "x-invisible gruop jigsaws");
            xinv.AddFrame(xi1);
            xinv.AddFrame(xi2);

            CombinatoricGroup xvis("xvis", "x-visible object jigsaws");

            // visible frames in first decay step must always have at least one element
            xvis.AddFrame(xv1s);
            xvis.AddFrame(xv2s);
            xvis.SetNElementsForFrame(xv1s, 1, false);
            xvis.SetNElementsForFrame(xv2s, 1, false);
            // visible frames in second decay step can have zero elements
            xvis.AddFrame(xv1c);
            xvis.AddFrame(xv2c);
            xvis.SetNElementsForFrame(xv1c, 0, false);
            xvis.SetNElementsForFrame(xv2c, 0, false);

            // define jigsaws
            SetMassInvJigsaw xMinMassJigsaw("xminmass", "x-Invisible system mass jigsaw");
            xinv.AddJigsaw(xMinMassJigsaw);
            SetRapidityInvJigsaw xRapidityJigsaw("xRapidity", "x-Invisible system rapidity jigsaw");
            xinv.AddJigsaw(xRapidityJigsaw);
            xRapidityJigsaw.AddVisibleFrames((xlab.GetListVisibleFrames()));
            ContraBoostInvJigsaw xContraBoostJigsaw("xContra", "x-Contraboost invariant jigsaw");
            xinv.AddJigsaw(xContraBoostJigsaw);
            xContraBoostJigsaw.AddVisibleFrames((xs1.GetListVisibleFrames()), 0);
            xContraBoostJigsaw.AddVisibleFrames((xs2.GetListVisibleFrames()), 1);
            xContraBoostJigsaw.AddInvisibleFrames((xs1.GetListInvisibleFrames()), 0);
            xContraBoostJigsaw.AddInvisibleFrames((xs2.GetListInvisibleFrames()), 1);

            MinMassesCombJigsaw xHemiJigsaw("xHemi", "x-Minimize m_{V_{1,2}} jigsaw");
            xvis.AddJigsaw(xHemiJigsaw);
            xHemiJigsaw.AddFrame(xv1s, 0);
            xHemiJigsaw.AddFrame(xv2s, 1);
            xHemiJigsaw.AddFrame(xv1c, 0);
            xHemiJigsaw.AddFrame(xv2c, 1);

            MinMassesCombJigsaw x1HemiJigsaw("x-1 Hemi", "x-1 Minimize m_{C_{a}} jigsaw");
            xvis.AddJigsaw(x1HemiJigsaw);
            x1HemiJigsaw.AddFrame(xv1s, 0);
            x1HemiJigsaw.AddFrame(xv1c, 1);
            x1HemiJigsaw.AddFrame(xi1, 1);

            MinMassesCombJigsaw x2HemiJigsaw("x-2 Hemi", "x-2 Minimize m_{C_{b}} jigsaw");
            xvis.AddJigsaw(x2HemiJigsaw);
            x2HemiJigsaw.AddFrame(xv2s, 0);
            x2HemiJigsaw.AddFrame(xv2c, 1);
            x2HemiJigsaw.AddFrame(xi2, 1);

            // check that the jigsaws are in place
            if(!xlab.InitializeAnalysis()) {
                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
                exit(1);
            }

            // clear the event for sho
            xlab.ClearEvent();

            // set the met
            TVector3 xmet3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
            xinv.SetLabFrameThreeVector(xmet3vector);

            // set up the jets and leptons
            //vector<TLorentzVector> xjets;
            //TLorentzVector l1, l2;

            //l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            //l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
            //for(int ij = 0; ij < (int)jets.size(); ij++) {
            //    TLorentzVector xj;
            //    xj.SetPtEtaPhiM(jets.at(ij)->Pt(), jets.at(ij)->Eta(), jets.at(ij)->Phi(), jets.at(ij)->M());
            //    xjets.push_back(xj);
            //}

            //vector<RFKey> jetID;
            //jetID.push_back(xvis.AddLabFrameFourVector(l1));
            //jetID.push_back(xvis.AddLabFrameFourVector(l2));
            //for(int ij = 0; ij < (int)xjets.size(); ij++) {
            //    jetID.push_back(xvis.AddLabFrameFourVector(xjets.at(ij)));
            //}


            TLorentzVector j1, j2, l1, l2;
            j1.SetPtEtaPhiM(jets.at(0)->Pt(), jets.at(0)->Eta(), jets.at(0)->Phi(), jets.at(0)->M());
            j2.SetPtEtaPhiM(jets.at(1)->Pt(), jets.at(1)->Eta(), jets.at(1)->Phi(), jets.at(1)->M());
            l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Eta(), leptons.at(1)->M());

            vector<RFKey> jetID;
            jetID.push_back(xvis.AddLabFrameFourVector(j1));
            jetID.push_back(xvis.AddLabFrameFourVector(j2));

            vector<RFKey> lepID;
            lepID.push_back(xvis.AddLabFrameFourVector(l1));
            lepID.push_back(xvis.AddLabFrameFourVector(l2));


            // analyze that
            xlab.AnalyzeEvent();

            DecayRecoFrame* S[2];
            DecayRecoFrame* C[2];
            VisibleRecoFrame* VS[2];
            VisibleRecoFrame* VC[2];
            InvisibleRecoFrame* I[2];
            //randomize
            int flip = (gRandom->Rndm() > 0.5);
            S[flip] = &xs1;
            S[(flip+1)%2] = &xs2;
            C[flip] = &xc1;
            C[(flip+1)%2] = &xc2;
            VS[flip] = &xv1s;
            VS[(flip+1)%2] = &xv2s;
            VC[flip] = &xv1c;
            VC[(flip+1)%2] = &xv2c;
            I[flip] = &xi1;
            I[(flip+1)%2] = &xi2;


            //////////////////////////////////////
            // Observables
            //////////////////////////////////////


            // total CM mass
            xshat = xss.GetMass();
            // 'mass-less' stop assumption gamma in CM frame
            xgaminv = xss.GetVisibleShape();

            TVector3 xvPSS = xss.GetFourVector(xlab).Vect();

            // ratio of CM pT to CM mass
            xRPT = xvPSS.Pt() / (xvPSS.Pt() + xshat/4.);
            // ratio of CM pz to CM mass
            xRPZ = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xshat/4.);
            // cos decay angle of ss system
            xcosSS = xss.GetCosDecayAngle();
            // delta phi between lab and SS decay planes
            xdphiLSS = xlab.GetDeltaPhiDecayPlanes(xss);


            TLorentzVector vVS1 = S[0]->GetVisibleFourVector(*S[0]);
            TLorentzVector vVS2 = S[1]->GetVisibleFourVector(*S[1]);

            // stop mass
            xMS = (vVS1.M2() - vVS2.M2()) / (2.*(vVS1.E() - vVS2.E()));

            xPS = S[0]->GetMomentum(xss);
            xMSS = 2.*sqrt(xPS*xPS + xMS*xMS);
            xgaminvSS = 2.*xMS/xMSS;
            double beta = sqrt(1.-xgaminv*xgaminv);
            double betaSS = sqrt(1.-xgaminvSS*xgaminvSS);

            // velocity difference between 'massive' and 'mass-less'
            xDeltaBetaSS = -(betaSS-beta)/(1.-betaSS*beta);

            // dleta phi between SS visible decay products and SS decay axis
            xDPD_vSS = xss.GetDeltaPhiDecayVisible();
            // delta phi between SS visible decay products and SS momentum
            xDPB_vSS = xss.GetDeltaPhiBoostVisible();


            // number of visible objects in hemisphere
            for(int i = 0; i < 2; i++) {
                xNV[i] = xvis.GetNElementsInFrame(*VS[i]); 
                xNV[i] += xvis.GetNElementsInFrame(*VC[i]);

                TVector3 xvP1 = VS[i]->GetFourVector(*S[i]).Vect();
                TVector3 xvP2 = VC[i]->GetFourVector(*S[i]).Vect();
                xPinv[i] = 2.*(xvP1+xvP2).Mag()/(xvP1.Mag()+xvP2.Mag() + (xvP1+xvP2).Mag());

                xcosS[i] = S[i]->GetCosDecayAngle();

                int N = jetID.size();
                double pTmax[2]; pTmax[0] = -1.; pTmax[1] = -1.;
                for(int j = 0; j < N; j++) {
                    const RestFrame& frame = xvis.GetFrame(jetID[j]);
                    if(VS[i]->IsSame(frame) || VC[i]->IsSame(frame)) { // jet is in hemisphere 'i'
                        double pT_ = xvis.GetLabFrameFourVector(jetID[j]).Pt();
                        if(pT_ > pTmax[0]) {
                            pTmax[1] = pTmax[0];
                            pTmax[0] = pT_;
                        } else {
                            if(pT_ > pTmax[1]) pTmax[1] = pT_;
                        }
                    }
                } // j

                xjet1PT[i] = pTmax[0]; // lead visible object pT in hemisphere i
                xjet2PT[i] = pTmax[1]; // sub lead visible object pT in hemisphere i

                xdphiVS_I[i] = VS[i]->GetFourVector(*S[i]).DeltaPhi(I[i]->GetFourVector(*S[i]));

                if(xNV[i] > 1) {
                    xcosS[i] = C[i]->GetCosDecayAngle();
                    xdphiSC[i] = S[i]->GetDeltaPhiDecayPlanes(*C[i]);
                    xRCS[i] = (C[i]->GetMass() - I[i]->GetMass())/(S[i]->GetMass()-I[i]->GetMass());

                } else {
                    xcosS[i] = -999;
                    xdphiSC[i] = -999;
                    xRCS[i] = -999;
                }
                
            } // i

            ////////////////////////////////////////////////////////////
            // scale variables
            ////////////////////////////////////////////////////////////

            TLorentzVector v_v1s_ss = VS[0]->GetFourVector(xss);
            TLorentzVector v_v2s_ss = VS[1]->GetFourVector(xss);
            TLorentzVector v_v1c_ss = VC[0]->GetFourVector(xss);
            TLorentzVector v_v2c_ss = VC[1]->GetFourVector(xss);
            TLorentzVector v_i1_ss  = I[0]->GetFourVector(xss);
            TLorentzVector v_i2_ss  = I[1]->GetFourVector(xss);
            if(v_i1_ss.Vect().Mag() > 1.0e4) {
                v_i1_ss.SetPtEtaPhiM(0, 0, 0, 0);
            }
            if(v_i2_ss.Vect().Mag() > 1.0e4) {
                v_i2_ss.SetPtEtaPhiM(0, 0, 0, 0);
            }

            TLorentzVector v_v1s_s1 = VS[0]->GetFourVector(xs1);
            TLorentzVector v_v1c_s1 = VC[0]->GetFourVector(xs1);
            TLorentzVector v_i1_s1  = I[0]->GetFourVector(xs1);

            // H_11_SS
            TVector3 p_vis_H11SS = (v_v1s_ss + v_v2s_ss + v_v1c_ss + v_v2c_ss).Vect();
            TVector3 p_invis_H11SS = (v_i1_ss + v_i2_ss).Vect();
            xH_11_SS = p_vis_H11SS.Mag() + p_invis_H11SS.Mag();

            // H_21_SS
            TVector3 p_vis_1_H21SS = (v_v1s_ss + v_v1c_ss).Vect();
            TVector3 p_vis_2_H21SS = (v_v2s_ss + v_v2c_ss).Vect();
            TVector3 p_invis_H21SS = (v_i1_ss + v_i2_ss).Vect();
            xH_21_SS = p_vis_1_H21SS.Mag() + p_vis_2_H21SS.Mag() + p_invis_H21SS.Mag();

            // H_41_SS
            xH_41_SS = v_v1s_ss.Vect().Mag() + v_v2s_ss.Vect().Mag() + v_v1c_ss.Vect().Mag() + v_v2c_ss.Vect().Mag() + (v_i1_ss + v_i2_ss).Vect().Mag();

            // H_42_SS
            xH_42_SS = v_v1s_ss.Vect().Mag() + v_v2s_ss.Vect().Mag() + v_v1c_ss.Vect().Mag() + v_v2c_ss.Vect().Mag() + v_i1_ss.Vect().Mag() + v_i2_ss.Vect().Mag();

            // H_11_S1
            xH_11_S1 = (v_v1s_s1+v_v1c_s1).Vect().Mag() + v_i1_s1.Vect().Mag();

            // H_21_S1
            xH_21_S1 = v_v1s_s1.Vect().Mag() + v_v1c_s1.Vect().Mag() + v_i1_s1.Vect().Mag();

            ////////////////
            // transverse scale variables
            ////////////////
            TVector3 tp_v1s_ss = v_v1s_ss.Vect(); tp_v1s_ss.SetZ(0.);
            TVector3 tp_v2s_ss = v_v2s_ss.Vect(); tp_v2s_ss.SetZ(0.);
            TVector3 tp_v1c_ss = v_v1c_ss.Vect(); tp_v1c_ss.SetZ(0.);
            TVector3 tp_v2c_ss = v_v2c_ss.Vect(); tp_v2c_ss.SetZ(0.);
            TVector3 tp_i1_ss  = v_i1_ss.Vect();  tp_i1_ss.SetZ(0.);
            TVector3 tp_i2_ss  = v_i2_ss.Vect();  tp_i2_ss.SetZ(0.);

            TVector3 tp_v1s_s1 = v_v1s_s1.Vect(); tp_v1s_s1.SetZ(0.);
            TVector3 tp_v1c_s1 = v_v1c_s1.Vect(); tp_v1c_s1.SetZ(0.);
            TVector3 tp_i1_s1  = v_i1_s1.Vect();  tp_i1_s1.SetZ(0.);

            // H_11_SS_T
            xH_11_SS_T = (tp_v1s_ss + tp_v2s_ss + tp_v1c_ss + tp_v2c_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();

            // H_21_SS_T
            xH_21_SS_T = (tp_v1s_ss + tp_v1c_ss).Mag() + (tp_v2s_ss + tp_v2c_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();

            // H_41_SS_T
            xH_41_SS_T = tp_v1s_ss.Mag() + tp_v2s_ss.Mag() + tp_v1c_ss.Mag() + tp_v2c_ss.Mag() + (tp_i1_ss+tp_i2_ss).Mag();

            // H_42_SS_T
            xH_42_SS_T = tp_v1s_ss.Mag() + tp_v2s_ss.Mag() + tp_v1c_ss.Mag() + tp_v2c_ss.Mag() + tp_i1_ss.Mag() + tp_i2_ss.Mag();

            // H_11_S1_T
            xH_11_S1_T = (tp_v1s_s1 + tp_v1c_s1).Mag() + tp_i1_s1.Mag();

            // H_21_S1_T
            xH_21_S1_T = tp_v1s_s1.Mag() + tp_v1c_s1.Mag() + tp_i1_s1.Mag();
        

            xRPT_H_11_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_11_SS/4.);
            xRPT_H_21_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_21_SS/4.);
            xRPT_H_41_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_41_SS/4.);
            xRPT_H_42_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_42_SS/4.);
            xRPZ_H_11_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_11_SS/4.);
            xRPZ_H_21_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_21_SS/4.);
            xRPZ_H_41_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_41_SS/4.);
            xRPZ_H_42_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_42_SS/4.);

            xRPT_H_11_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_11_SS_T/4.);
            xRPT_H_21_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_21_SS_T/4.);
            xRPT_H_41_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_41_SS_T/4.);
            xRPT_H_42_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_42_SS_T/4.);
            xRPZ_H_11_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_11_SS_T/4.);
            xRPZ_H_21_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_21_SS_T/4.);
            xRPZ_H_41_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_41_SS_T/4.);
            xRPZ_H_42_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_42_SS_T/4.);


            ///////////////////////////
            // feb 1 variables
            //xCosP1 = S[0]->GetCosDecayAngle(I[0]);
            //xCosP2 = S[1]->GetCosDecayAngle(I[1]);

        } // njets == 2
    };

    *cutflow << NewVar("xshat"); {
        *cutflow << HFTname("xshat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xshat;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xgaminv"); {
        *cutflow << HFTname("xgaminv");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xgaminv;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT"); {
        *cutflow << HFTname("xRPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ"); {
        *cutflow << HFTname("xRPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosSS"); {
        *cutflow << HFTname("xcosSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xcosSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiLSS"); {
        *cutflow << HFTname("xdphiLSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xdphiLSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xMS"); {
        *cutflow << HFTname("xMS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xMS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xPS"); {
        *cutflow << HFTname("xPS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xPS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xMSS"); {
        *cutflow << HFTname("xMSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xMSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xgaminvSS"); {
        *cutflow << HFTname("xgaminvSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xgaminvSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDeltaBetaSS"); {
        *cutflow << HFTname("xDeltaBetaSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xDeltaBetaSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDPD_vSS"); {
        *cutflow << HFTname("xDPD_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xDPD_vSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDPB_vSS"); {
        *cutflow << HFTname("xDPB_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xDPB_vSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV"); {
        *cutflow << HFTname("xNV");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xNV[0]);
            out.push_back(xNV[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV_0"); {
        *cutflow << HFTname("xNV_0");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return xNV[0];
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV_1"); {
        *cutflow << HFTname("xNV_1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return xNV[1];
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosS"); {
        *cutflow << HFTname("xcosS");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xcosS[0]);
            out.push_back(xcosS[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosC"); {
        *cutflow << HFTname("xcosC");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xcosC[0]);
            out.push_back(xcosC[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiSC"); {
        *cutflow << HFTname("xdphiSC");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xdphiSC[0]);
            out.push_back(xdphiSC[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRCS"); {
        *cutflow << HFTname("xRCS");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xRCS[0]);
            out.push_back(xRCS[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xjet1PT"); {
        *cutflow << HFTname("xjet1PT");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xjet1PT[0]);
            out.push_back(xjet1PT[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xjet2PT"); {
        *cutflow << HFTname("xjet2PT");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xjet2PT[0]);
            out.push_back(xjet2PT[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xPinv"); {
        *cutflow << HFTname("xPinv");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xPinv[0]);
            out.push_back(xPinv[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_SS"); {
        *cutflow << HFTname("xH_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_SS"); {
        *cutflow << HFTname("xH_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_41_SS"); {
        *cutflow << HFTname("xH_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_41_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_42_SS"); {
        *cutflow << HFTname("xH_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_42_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_S1"); {
        *cutflow << HFTname("xH_11_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_S1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_S1"); {
        *cutflow << HFTname("xH_21_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_S1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_SS_T"); {
        *cutflow << HFTname("xH_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_SS_T"); {
        *cutflow << HFTname("xH_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_41_SS_T"); {
        *cutflow << HFTname("xH_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_41_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_42_SS_T"); {
        *cutflow << HFTname("xH_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_42_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_S1_T"); {
        *cutflow << HFTname("xH_11_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_S1_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_S1_T"); {
        *cutflow << HFTname("xH_21_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_S1_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_11_SS"); {
        *cutflow << HFTname("xRPT_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_11_SS;
        };
        *cutflow << SaveVar();
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //
    // RESTFRAMES - 2 jets, 2 leptons
    // RESTFRAMES - 2 jets, 2 leptons
    // RESTFRAMES - 2 jets, 2 leptons
    // RESTFRAMES - 2 jets, 2 leptons
    //
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
/*
    double xshat;
    double xgaminv;
    double xRPT;
    double xRPZ;
    double xcosSS;
    double xdphiLSS;
    double xMS;
    double xPS;
    double xMSS;
    double xgaminvSS;
    double xDeltaBetaSS;
    double xDPD_vSS;
    double xDPB_vSS;
    int xNV[2]; // number of visible objects in hemisphere
    double xcosS[2]; // cosine stop decay angle
    double xcosC[2]; // cosine intermediate child decay angle
    double xdphiSC[2]; // cosine between stop and child decay planes
    double xRCS[2]; // ratio of child and stop masses (w/ WIMP masses subtracted);
    double xjet1PT[2]; // first leading jet pT associated with this hemisphere
    double xjet2PT[2]; // second leading jet pT associated with this hemisphere
    double xPinv[2]; // Pinv / HS
    double xH_11_SS;
    double xH_21_SS;
    double xH_41_SS;
    double xH_42_SS;
    double xH_11_S1;
    double xH_21_S1;
    double xH_11_SS_T;
    double xH_21_SS_T;
    double xH_41_SS_T;
    double xH_42_SS_T;
    double xH_11_S1_T;
    double xH_21_S1_T;

    double xRPT_H_11_SS; 
    double xRPT_H_21_SS; 
    double xRPT_H_41_SS; 
    double xRPT_H_42_SS; 
    double xRPZ_H_11_SS; 
    double xRPZ_H_21_SS; 
    double xRPZ_H_41_SS; 
    double xRPZ_H_42_SS; 
    double xRPT_H_11_SS_T; 
    double xRPT_H_21_SS_T; 
    double xRPT_H_41_SS_T; 
    double xRPT_H_42_SS_T; 
    double xRPZ_H_11_SS_T; 
    double xRPZ_H_21_SS_T; 
    double xRPZ_H_41_SS_T; 
    double xRPZ_H_42_SS_T; 

    double xdphiVS_I[2];

    // feb 1 vars
    //double xCosP1;
    //double xCosP2;


    *cutflow << [&](Superlink* sl, var_void*) {
        if(sjets.size()>=2 && bjets.size()>=1) {

            // setup the analysis tree
            LabRecoFrame xlab("xlab", "xlab");
            DecayRecoFrame xss("xss", "xss");
            DecayRecoFrame xs1("xs1", "xs1");
            DecayRecoFrame xs2("xs2", "xs2");
            DecayRecoFrame xc1("xc1", "xc1");
            DecayRecoFrame xc2("xc2", "xc2");
            VisibleRecoFrame xv1s("xv1s", "xv1s");
            VisibleRecoFrame xv2s("xv2s", "xv2s");
            InvisibleRecoFrame xi1("xi1", "xi1");
            VisibleRecoFrame xv1c("xv1c", "xv1c");
            VisibleRecoFrame xv2c("xv2c", "xv2c");
            InvisibleRecoFrame xi2("xi2", "xi2");

            xlab.SetChildFrame(xss);
            xss.AddChildFrame(xs1);
            xss.AddChildFrame(xs2);
            xs1.AddChildFrame(xc1);
            xs1.AddChildFrame(xv1s);
            xc1.AddChildFrame(xi1);
            xc1.AddChildFrame(xv1c);
            xs2.AddChildFrame(xc2);
            xs2.AddChildFrame(xv2s);
            xc2.AddChildFrame(xi2);
            xc2.AddChildFrame(xv2c);

            // check that the decay tree is connected properly
            if(!xlab.InitializeTree()) {
                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ <<")    Unable to initialize tree from lab frame. Exiting." << endl;
                exit(1);
            }


            // define groupes

            InvisibleGroup xinv("xinv", "x-invisible gruop jigsaws");
            xinv.AddFrame(xi1);
            xinv.AddFrame(xi2);

            CombinatoricGroup xvis("xvis", "x-visible object jigsaws");

            // visible frames in first decay step must always have at least one element
            xvis.AddFrame(xv1s);
            xvis.AddFrame(xv2s);
            xvis.SetNElementsForFrame(xv1s, 1, false);
            xvis.SetNElementsForFrame(xv2s, 1, false);
            // visible frames in second decay step can have zero elements
            xvis.AddFrame(xv1c);
            xvis.AddFrame(xv2c);
            xvis.SetNElementsForFrame(xv1c, 0, false);
            xvis.SetNElementsForFrame(xv2c, 0, false);

            // define jigsaws
            SetMassInvJigsaw xMinMassJigsaw("xminmass", "x-Invisible system mass jigsaw");
            xinv.AddJigsaw(xMinMassJigsaw);
            SetRapidityInvJigsaw xRapidityJigsaw("xRapidity", "x-Invisible system rapidity jigsaw");
            xinv.AddJigsaw(xRapidityJigsaw);
            xRapidityJigsaw.AddVisibleFrames((xlab.GetListVisibleFrames()));
            ContraBoostInvJigsaw xContraBoostJigsaw("xContra", "x-Contraboost invariant jigsaw");
            xinv.AddJigsaw(xContraBoostJigsaw);
            xContraBoostJigsaw.AddVisibleFrames((xs1.GetListVisibleFrames()), 0);
            xContraBoostJigsaw.AddVisibleFrames((xs2.GetListVisibleFrames()), 1);
            xContraBoostJigsaw.AddInvisibleFrames((xs1.GetListInvisibleFrames()), 0);
            xContraBoostJigsaw.AddInvisibleFrames((xs2.GetListInvisibleFrames()), 1);

            MinMassesCombJigsaw xHemiJigsaw("xHemi", "x-Minimize m_{V_{1,2}} jigsaw");
            xvis.AddJigsaw(xHemiJigsaw);
            xHemiJigsaw.AddFrame(xv1s, 0);
            xHemiJigsaw.AddFrame(xv2s, 1);
            xHemiJigsaw.AddFrame(xv1c, 0);
            xHemiJigsaw.AddFrame(xv2c, 1);

            MinMassesCombJigsaw x1HemiJigsaw("x-1 Hemi", "x-1 Minimize m_{C_{a}} jigsaw");
            xvis.AddJigsaw(x1HemiJigsaw);
            x1HemiJigsaw.AddFrame(xv1s, 0);
            x1HemiJigsaw.AddFrame(xv1c, 1);
            x1HemiJigsaw.AddFrame(xi1, 1);

            MinMassesCombJigsaw x2HemiJigsaw("x-2 Hemi", "x-2 Minimize m_{C_{b}} jigsaw");
            xvis.AddJigsaw(x2HemiJigsaw);
            x2HemiJigsaw.AddFrame(xv2s, 0);
            x2HemiJigsaw.AddFrame(xv2c, 1);
            x2HemiJigsaw.AddFrame(xi2, 1);

            // check that the jigsaws are in place
            if(!xlab.InitializeAnalysis()) {
                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
                exit(1);
            }

            // clear the event for sho
            xlab.ClearEvent();

            // set the met
            TVector3 xmet3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
            xinv.SetLabFrameThreeVector(xmet3vector);

            // set up the jets and leptons
            //vector<TLorentzVector> xjets;
            //TLorentzVector l1, l2;

            //l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            //l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
            //for(int ij = 0; ij < (int)jets.size(); ij++) {
            //    TLorentzVector xj;
            //    xj.SetPtEtaPhiM(jets.at(ij)->Pt(), jets.at(ij)->Eta(), jets.at(ij)->Phi(), jets.at(ij)->M());
            //    xjets.push_back(xj);
            //}

            //vector<RFKey> jetID;
            //jetID.push_back(xvis.AddLabFrameFourVector(l1));
            //jetID.push_back(xvis.AddLabFrameFourVector(l2));
            //for(int ij = 0; ij < (int)xjets.size(); ij++) {
            //    jetID.push_back(xvis.AddLabFrameFourVector(xjets.at(ij)));
            //}


            TLorentzVector j1, j2, l1, l2;
            j1.SetPtEtaPhiM(jets.at(0)->Pt(), jets.at(0)->Eta(), jets.at(0)->Phi(), jets.at(0)->M());
            j2.SetPtEtaPhiM(jets.at(1)->Pt(), jets.at(1)->Eta(), jets.at(1)->Phi(), jets.at(1)->M());
            l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Eta(), leptons.at(1)->M());

            vector<RFKey> jetID;
            jetID.push_back(xvis.AddLabFrameFourVector(j1));
            jetID.push_back(xvis.AddLabFrameFourVector(j2));

            vector<RFKey> lepID;
            lepID.push_back(xvis.AddLabFrameFourVector(l1));
            lepID.push_back(xvis.AddLabFrameFourVector(l2));


            // analyze that
            xlab.AnalyzeEvent();

            DecayRecoFrame* S[2];
            DecayRecoFrame* C[2];
            VisibleRecoFrame* VS[2];
            VisibleRecoFrame* VC[2];
            InvisibleRecoFrame* I[2];
            //randomize
            int flip = (gRandom->Rndm() > 0.5);
            S[flip] = &xs1;
            S[(flip+1)%2] = &xs2;
            C[flip] = &xc1;
            C[(flip+1)%2] = &xc2;
            VS[flip] = &xv1s;
            VS[(flip+1)%2] = &xv2s;
            VC[flip] = &xv1c;
            VC[(flip+1)%2] = &xv2c;
            I[flip] = &xi1;
            I[(flip+1)%2] = &xi2;


            //////////////////////////////////////
            // Observables
            //////////////////////////////////////


            // total CM mass
            xshat = xss.GetMass();
            // 'mass-less' stop assumption gamma in CM frame
            xgaminv = xss.GetVisibleShape();

            TVector3 xvPSS = xss.GetFourVector(xlab).Vect();

            // ratio of CM pT to CM mass
            xRPT = xvPSS.Pt() / (xvPSS.Pt() + xshat/4.);
            // ratio of CM pz to CM mass
            xRPZ = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xshat/4.);
            // cos decay angle of ss system
            xcosSS = xss.GetCosDecayAngle();
            // delta phi between lab and SS decay planes
            xdphiLSS = xlab.GetDeltaPhiDecayPlanes(xss);


            TLorentzVector vVS1 = S[0]->GetVisibleFourVector(*S[0]);
            TLorentzVector vVS2 = S[1]->GetVisibleFourVector(*S[1]);

            // stop mass
            xMS = (vVS1.M2() - vVS2.M2()) / (2.*(vVS1.E() - vVS2.E()));

            xPS = S[0]->GetMomentum(xss);
            xMSS = 2.*sqrt(xPS*xPS + xMS*xMS);
            xgaminvSS = 2.*xMS/xMSS;
            double beta = sqrt(1.-xgaminv*xgaminv);
            double betaSS = sqrt(1.-xgaminvSS*xgaminvSS);

            // velocity difference between 'massive' and 'mass-less'
            xDeltaBetaSS = -(betaSS-beta)/(1.-betaSS*beta);

            // dleta phi between SS visible decay products and SS decay axis
            xDPD_vSS = xss.GetDeltaPhiDecayVisible();
            // delta phi between SS visible decay products and SS momentum
            xDPB_vSS = xss.GetDeltaPhiBoostVisible();


            // number of visible objects in hemisphere
            for(int i = 0; i < 2; i++) {
                xNV[i] = xvis.GetNElementsInFrame(*VS[i]); 
                xNV[i] += xvis.GetNElementsInFrame(*VC[i]);

                TVector3 xvP1 = VS[i]->GetFourVector(*S[i]).Vect();
                TVector3 xvP2 = VC[i]->GetFourVector(*S[i]).Vect();
                xPinv[i] = 2.*(xvP1+xvP2).Mag()/(xvP1.Mag()+xvP2.Mag() + (xvP1+xvP2).Mag());

                xcosS[i] = S[i]->GetCosDecayAngle();

                int N = jetID.size();
                double pTmax[2]; pTmax[0] = -1.; pTmax[1] = -1.;
                for(int j = 0; j < N; j++) {
                    const RestFrame& frame = xvis.GetFrame(jetID[j]);
                    if(VS[i]->IsSame(frame) || VC[i]->IsSame(frame)) { // jet is in hemisphere 'i'
                        double pT_ = xvis.GetLabFrameFourVector(jetID[j]).Pt();
                        if(pT_ > pTmax[0]) {
                            pTmax[1] = pTmax[0];
                            pTmax[0] = pT_;
                        } else {
                            if(pT_ > pTmax[1]) pTmax[1] = pT_;
                        }
                    }
                } // j

                xjet1PT[i] = pTmax[0]; // lead visible object pT in hemisphere i
                xjet2PT[i] = pTmax[1]; // sub lead visible object pT in hemisphere i

                xdphiVS_I[i] = VS[i]->GetFourVector(*S[i]).DeltaPhi(I[i]->GetFourVector(*S[i]));

                if(xNV[i] > 1) {
                    xcosS[i] = C[i]->GetCosDecayAngle();
                    xdphiSC[i] = S[i]->GetDeltaPhiDecayPlanes(*C[i]);
                    xRCS[i] = (C[i]->GetMass() - I[i]->GetMass())/(S[i]->GetMass()-I[i]->GetMass());

                } else {
                    xcosS[i] = -999;
                    xdphiSC[i] = -999;
                    xRCS[i] = -999;
                }
                
            } // i

            ////////////////////////////////////////////////////////////
            // scale variables
            ////////////////////////////////////////////////////////////

            TLorentzVector v_v1s_ss = VS[0]->GetFourVector(xss);
            TLorentzVector v_v2s_ss = VS[1]->GetFourVector(xss);
            TLorentzVector v_v1c_ss = VC[0]->GetFourVector(xss);
            TLorentzVector v_v2c_ss = VC[1]->GetFourVector(xss);
            TLorentzVector v_i1_ss  = I[0]->GetFourVector(xss);
            TLorentzVector v_i2_ss  = I[1]->GetFourVector(xss);
            if(v_i1_ss.Vect().Mag() > 1.0e4) {
                v_i1_ss.SetPtEtaPhiM(0, 0, 0, 0);
            }
            if(v_i2_ss.Vect().Mag() > 1.0e4) {
                v_i2_ss.SetPtEtaPhiM(0, 0, 0, 0);
            }

            TLorentzVector v_v1s_s1 = VS[0]->GetFourVector(xs1);
            TLorentzVector v_v1c_s1 = VC[0]->GetFourVector(xs1);
            TLorentzVector v_i1_s1  = I[0]->GetFourVector(xs1);

            // H_11_SS
            TVector3 p_vis_H11SS = (v_v1s_ss + v_v2s_ss + v_v1c_ss + v_v2c_ss).Vect();
            TVector3 p_invis_H11SS = (v_i1_ss + v_i2_ss).Vect();
            xH_11_SS = p_vis_H11SS.Mag() + p_invis_H11SS.Mag();

            // H_21_SS
            TVector3 p_vis_1_H21SS = (v_v1s_ss + v_v1c_ss).Vect();
            TVector3 p_vis_2_H21SS = (v_v2s_ss + v_v2c_ss).Vect();
            TVector3 p_invis_H21SS = (v_i1_ss + v_i2_ss).Vect();
            xH_21_SS = p_vis_1_H21SS.Mag() + p_vis_2_H21SS.Mag() + p_invis_H21SS.Mag();

            // H_41_SS
            xH_41_SS = v_v1s_ss.Vect().Mag() + v_v2s_ss.Vect().Mag() + v_v1c_ss.Vect().Mag() + v_v2c_ss.Vect().Mag() + (v_i1_ss + v_i2_ss).Vect().Mag();

            // H_42_SS
            xH_42_SS = v_v1s_ss.Vect().Mag() + v_v2s_ss.Vect().Mag() + v_v1c_ss.Vect().Mag() + v_v2c_ss.Vect().Mag() + v_i1_ss.Vect().Mag() + v_i2_ss.Vect().Mag();

            // H_11_S1
            xH_11_S1 = (v_v1s_s1+v_v1c_s1).Vect().Mag() + v_i1_s1.Vect().Mag();

            // H_21_S1
            xH_21_S1 = v_v1s_s1.Vect().Mag() + v_v1c_s1.Vect().Mag() + v_i1_s1.Vect().Mag();

            ////////////////
            // transverse scale variables
            ////////////////
            TVector3 tp_v1s_ss = v_v1s_ss.Vect(); tp_v1s_ss.SetZ(0.);
            TVector3 tp_v2s_ss = v_v2s_ss.Vect(); tp_v2s_ss.SetZ(0.);
            TVector3 tp_v1c_ss = v_v1c_ss.Vect(); tp_v1c_ss.SetZ(0.);
            TVector3 tp_v2c_ss = v_v2c_ss.Vect(); tp_v2c_ss.SetZ(0.);
            TVector3 tp_i1_ss  = v_i1_ss.Vect();  tp_i1_ss.SetZ(0.);
            TVector3 tp_i2_ss  = v_i2_ss.Vect();  tp_i2_ss.SetZ(0.);

            TVector3 tp_v1s_s1 = v_v1s_s1.Vect(); tp_v1s_s1.SetZ(0.);
            TVector3 tp_v1c_s1 = v_v1c_s1.Vect(); tp_v1c_s1.SetZ(0.);
            TVector3 tp_i1_s1  = v_i1_s1.Vect();  tp_i1_s1.SetZ(0.);

            // H_11_SS_T
            xH_11_SS_T = (tp_v1s_ss + tp_v2s_ss + tp_v1c_ss + tp_v2c_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();

            // H_21_SS_T
            xH_21_SS_T = (tp_v1s_ss + tp_v1c_ss).Mag() + (tp_v2s_ss + tp_v2c_ss).Mag() + (tp_i1_ss + tp_i2_ss).Mag();

            // H_41_SS_T
            xH_41_SS_T = tp_v1s_ss.Mag() + tp_v2s_ss.Mag() + tp_v1c_ss.Mag() + tp_v2c_ss.Mag() + (tp_i1_ss+tp_i2_ss).Mag();

            // H_42_SS_T
            xH_42_SS_T = tp_v1s_ss.Mag() + tp_v2s_ss.Mag() + tp_v1c_ss.Mag() + tp_v2c_ss.Mag() + tp_i1_ss.Mag() + tp_i2_ss.Mag();

            // H_11_S1_T
            xH_11_S1_T = (tp_v1s_s1 + tp_v1c_s1).Mag() + tp_i1_s1.Mag();

            // H_21_S1_T
            xH_21_S1_T = tp_v1s_s1.Mag() + tp_v1c_s1.Mag() + tp_i1_s1.Mag();
        

            xRPT_H_11_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_11_SS/4.);
            xRPT_H_21_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_21_SS/4.);
            xRPT_H_41_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_41_SS/4.);
            xRPT_H_42_SS = xvPSS.Pt() / (xvPSS.Pt() + xH_42_SS/4.);
            xRPZ_H_11_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_11_SS/4.);
            xRPZ_H_21_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_21_SS/4.);
            xRPZ_H_41_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_41_SS/4.);
            xRPZ_H_42_SS = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_42_SS/4.);

            xRPT_H_11_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_11_SS_T/4.);
            xRPT_H_21_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_21_SS_T/4.);
            xRPT_H_41_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_41_SS_T/4.);
            xRPT_H_42_SS_T = xvPSS.Pt() / (xvPSS.Pt() + xH_42_SS_T/4.);
            xRPZ_H_11_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_11_SS_T/4.);
            xRPZ_H_21_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_21_SS_T/4.);
            xRPZ_H_41_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_41_SS_T/4.);
            xRPZ_H_42_SS_T = fabs(xvPSS.Pz()) / (fabs(xvPSS.Pz()) + xH_42_SS_T/4.);


            ///////////////////////////
            // feb 1 variables
            //xCosP1 = S[0]->GetCosDecayAngle(I[0]);
            //xCosP2 = S[1]->GetCosDecayAngle(I[1]);

        } // njets == 2
    };

    *cutflow << NewVar("xshat"); {
        *cutflow << HFTname("xshat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xshat;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xgaminv"); {
        *cutflow << HFTname("xgaminv");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xgaminv;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT"); {
        *cutflow << HFTname("xRPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ"); {
        *cutflow << HFTname("xRPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosSS"); {
        *cutflow << HFTname("xcosSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xcosSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiLSS"); {
        *cutflow << HFTname("xdphiLSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xdphiLSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xMS"); {
        *cutflow << HFTname("xMS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xMS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xPS"); {
        *cutflow << HFTname("xPS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xPS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xMSS"); {
        *cutflow << HFTname("xMSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xMSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xgaminvSS"); {
        *cutflow << HFTname("xgaminvSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xgaminvSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDeltaBetaSS"); {
        *cutflow << HFTname("xDeltaBetaSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xDeltaBetaSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDPD_vSS"); {
        *cutflow << HFTname("xDPD_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xDPD_vSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDPB_vSS"); {
        *cutflow << HFTname("xDPB_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xDPB_vSS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV"); {
        *cutflow << HFTname("xNV");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xNV[0]);
            out.push_back(xNV[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV_0"); {
        *cutflow << HFTname("xNV_0");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return xNV[0];
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV_1"); {
        *cutflow << HFTname("xNV_1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return xNV[1];
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosS"); {
        *cutflow << HFTname("xcosS");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xcosS[0]);
            out.push_back(xcosS[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosC"); {
        *cutflow << HFTname("xcosC");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xcosC[0]);
            out.push_back(xcosC[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiSC"); {
        *cutflow << HFTname("xdphiSC");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xdphiSC[0]);
            out.push_back(xdphiSC[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRCS"); {
        *cutflow << HFTname("xRCS");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xRCS[0]);
            out.push_back(xRCS[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xjet1PT"); {
        *cutflow << HFTname("xjet1PT");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xjet1PT[0]);
            out.push_back(xjet1PT[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xjet2PT"); {
        *cutflow << HFTname("xjet2PT");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xjet2PT[0]);
            out.push_back(xjet2PT[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xPinv"); {
        *cutflow << HFTname("xPinv");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xPinv[0]);
            out.push_back(xPinv[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_SS"); {
        *cutflow << HFTname("xH_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_SS"); {
        *cutflow << HFTname("xH_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_41_SS"); {
        *cutflow << HFTname("xH_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_41_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_42_SS"); {
        *cutflow << HFTname("xH_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_42_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_S1"); {
        *cutflow << HFTname("xH_11_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_S1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_S1"); {
        *cutflow << HFTname("xH_21_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_S1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_SS_T"); {
        *cutflow << HFTname("xH_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_SS_T"); {
        *cutflow << HFTname("xH_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_41_SS_T"); {
        *cutflow << HFTname("xH_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_41_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_42_SS_T"); {
        *cutflow << HFTname("xH_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_42_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_S1_T"); {
        *cutflow << HFTname("xH_11_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_11_S1_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_S1_T"); {
        *cutflow << HFTname("xH_21_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xH_21_S1_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_11_SS"); {
        *cutflow << HFTname("xRPT_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_21_SS"); {
        *cutflow << HFTname("xRPT_H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_41_SS"); {
        *cutflow << HFTname("xRPT_H_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_41_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_42_SS"); {
        *cutflow << HFTname("xRPT_H_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_42_SS;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("xRPZ_H_11_SS"); {
        *cutflow << HFTname("xRPZ_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_11_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_21_SS"); {
        *cutflow << HFTname("xRPZ_H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_21_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_41_SS"); {
        *cutflow << HFTname("xRPZ_H_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_41_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_42_SS"); {
        *cutflow << HFTname("xRPZ_H_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_42_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_11_SS_T"); {
        *cutflow << HFTname("xRPT_H_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_21_SS_T"); {
        *cutflow << HFTname("xRPT_H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_41_SS_T"); {
        *cutflow << HFTname("xRPT_H_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_41_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT_H_42_SS_T"); {
        *cutflow << HFTname("xRPT_H_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPT_H_42_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_11_SS_T"); {
        *cutflow << HFTname("xRPZ_H_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_11_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_21_SS_T"); {
        *cutflow << HFTname("xRPZ_H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_21_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_41_SS_T"); {
        *cutflow << HFTname("xRPZ_H_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_41_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ_H_42_SS_T"); {
        *cutflow << HFTname("xRPZ_H_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return xRPZ_H_42_SS_T;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiVS_I"); {
        *cutflow << HFTname("xdphiVS_I");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.clear();
            out.push_back(xdphiVS_I[0]);
            out.push_back(xdphiVS_I[1]);
            return out;
        };
        *cutflow << SaveVar();
    }
*/
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // RESTFRAMES - mT
    // RESTFRAMES - mT
    // RESTFRAMES - mT
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    double t_shat;
    double t_gaminv;
    double t_pTCM;
    double t_pZCM;
    double t_RPT;
    double t_RPZ;
    double t_cosSS;
    double t_dphiLSS;
    double t_dphiSS_SA;
    double t_dphiSS_SB;
    double t_dphiSA_SB;
    double t_MS;
    double t_DPD_vSS;
    double t_DPB_vSS;
    int t_N_VA1;
    int t_N_VA2;
    int t_N_VB1;
    int t_N_VB2;
    int t_N_VA;
    int t_N_VB;
    int t_N_V1;
    int t_N_V2;
    int t_Nsj_VA1;
    int t_Nsj_VA2;
    int t_Nsj_VB1;
    int t_Nsj_VB2;
    int t_Nbj_VA1;
    int t_Nbj_VA2;
    int t_Nbj_VB1;
    int t_Nbj_VB2;
    int t_Nlep_VA1;
    int t_Nlep_VA2;
    int t_Nlep_VB1;
    int t_Nlep_VB2;
    int t_Nsj_A ;
    int t_Nsj_B ;
    int t_Nsj_V1;
    int t_Nsj_V2;
    int t_Nbj_A ;
    int t_Nbj_B ;
    int t_Nbj_V1;
    int t_Nbj_V2;
    int t_Nlep_A ;
    int t_Nlep_B ;
    int t_Nlep_V1;
    int t_Nlep_V2;
    double t_H_11_SS;
    double t_H_21_SS;
    double t_H_41_SS;
    double t_H_22_SS;
    double t_H_42_SS;
    double t_H_11_SS_T;
    double t_H_21_SS_T;
    double t_H_41_SS_T;
    double t_H_22_SS_T;
    double t_H_42_SS_T;
    double t_H_11_SA;
    double t_H_21_SA;
    double t_H_11_SA_T;
    double t_H_21_SA_T;
    double t_MDR_VA1_SA;
    double t_MDR_VA2_CA;

    *cutflow << [&](Superlink* sl, var_void*) {
        if(jets.size()>=2 && leptons.size()==2) {

            // setup the tree
            LabRecoFrame tLAB("LAB", "LAB");
            DecayRecoFrame tSS("SS", "SS");
            DecayRecoFrame tSA("SA", "SA");
            DecayRecoFrame tSB("SB", "SB");
            DecayRecoFrame tCA("CA", "CA");
            DecayRecoFrame tCB("CB", "CB");
            VisibleRecoFrame tVA1("VA1", "VA1");
            VisibleRecoFrame tVB1("VB1", "VB1");
            InvisibleRecoFrame tIA("IA", "IA");
            InvisibleRecoFrame tIB("IB", "IB");
            VisibleRecoFrame tVA2("VA2", "VA2");
            VisibleRecoFrame tVB2("VB2", "VB2");

            tLAB.SetChildFrame(tSS);

            tSS.AddChildFrame(tSA);
            tSS.AddChildFrame(tSB);

            tSA.AddChildFrame(tCA);
            tSA.AddChildFrame(tVA1);
            tSB.AddChildFrame(tCB);
            tSB.AddChildFrame(tVB1);

            tCA.AddChildFrame(tIA);
            tCA.AddChildFrame(tVA2);
            tCB.AddChildFrame(tIB);
            tCB.AddChildFrame(tVB2);

            if(!tLAB.InitializeTree()) {
                cout << analysis_name << "    RestFrames::InitializeTree ERROR (" << __LINE__ << ")    Unable to initialize tree from lab frame. Exiting." << endl;
                exit(1);
            }

            // define groups
            InvisibleGroup tINV("tINV", "t-Invisible group jigsaws");
            tINV.AddFrame(tIA);
            tINV.AddFrame(tIB);

            CombinatoricGroup tVIS("tVIS", "t-Visible object jigsaws");

            tVIS.AddFrame(tVA1);
            tVIS.AddFrame(tVB1);
            tVIS.AddFrame(tVA2);
            tVIS.AddFrame(tVB2);
            tVIS.SetNElementsForFrame(tVA1, 1, false);
            tVIS.SetNElementsForFrame(tVB1, 1, false);
            tVIS.SetNElementsForFrame(tVA2, 1, false);
            tVIS.SetNElementsForFrame(tVB2, 1, false);

            // define jigsaws
            SetMassInvJigsaw tMinMassJigsaw("tminmass", "t-Invisible system mass jigsaw");
            tINV.AddJigsaw(tMinMassJigsaw);

            SetRapidityInvJigsaw tRapidityJigsaw("tRapidity", "t-Invisible system rapidity jigsaw");
            tINV.AddJigsaw(tRapidityJigsaw);
            tRapidityJigsaw.AddVisibleFrames((tLAB.GetListVisibleFrames()));

            ContraBoostInvJigsaw tContraBoostJigsaw("tContra", "t-ContraBoost Invariant Jigsaw");
            tINV.AddJigsaw(tContraBoostJigsaw);
            tContraBoostJigsaw.AddVisibleFrames((tSA.GetListVisibleFrames()), 0);
            tContraBoostJigsaw.AddVisibleFrames((tSB.GetListVisibleFrames()), 1);
            tContraBoostJigsaw.AddInvisibleFrames((tSA.GetListInvisibleFrames()), 0);
            tContraBoostJigsaw.AddInvisibleFrames((tSB.GetListInvisibleFrames()), 1);

            MinMassesCombJigsaw tHemiJigsaw("tHemi", "t-Minimize m_{V_{a,b}} jigsaw");
            tVIS.AddJigsaw(tHemiJigsaw);
            tHemiJigsaw.AddFrame(tVA1, 0);
            tHemiJigsaw.AddFrame(tVA2, 0);
            tHemiJigsaw.AddFrame(tVB1, 1);
            tHemiJigsaw.AddFrame(tVB2, 1);

            MinMassesCombJigsaw taHemiJigsaw("tHemi-a", "t-Minimize m_{C_{a}} jigsaw");
            tVIS.AddJigsaw(taHemiJigsaw);
            taHemiJigsaw.AddFrame(tVA1, 0);
            taHemiJigsaw.AddFrame(tVA2, 1);
            taHemiJigsaw.AddFrame(tIA,  1);

            MinMassesCombJigsaw tbHemiJigsaw("tHemi-b", "t-Minimize m_{C_{b}} jigsaw");
            tVIS.AddJigsaw(tbHemiJigsaw);
            tbHemiJigsaw.AddFrame(tVB1, 0);
            tbHemiJigsaw.AddFrame(tVB2, 1);
            tbHemiJigsaw.AddFrame(tIB,  1);

            // check that jigsaws are connected properly
            if(!tLAB.InitializeAnalysis()) {
                cout << analysis_name << "    RestFrames::InitializeAnalysis ERROR (" << __LINE__ << ")    Unable to initialize analysis from lab frame. Exiting." << endl;
                exit(1);
            }

            // clear the event
            tLAB.ClearEvent();

            // set the met
            TVector3 tmet3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
            tINV.SetLabFrameThreeVector(tmet3vector);

            // build up our ish
            vector<TLorentzVector> tJETS, tSJETS, tBJETS, tLEPTONS;
            for(int i = 0; i < (int)jets.size(); i++) {
                Jet* j = jets[i];
                if(sl->tools->jetSelector().isB(j)) {
                    TLorentzVector b;
                    b.SetPtEtaPhiM(j->Pt(), j->Eta(), j->Phi(), j->M());
                    tBJETS.push_back(b);
                }
                else if(!(sl->tools->jetSelector().isB(j))) {
                    TLorentzVector sj;
                    sj.SetPtEtaPhiM(j->Pt(), j->Eta(), j->Phi(), j->M());
                    tSJETS.push_back(sj);
                }
            }
            for(int i = 0; i < (int)leptons.size(); i++) {
                TLorentzVector l;
                l.SetPtEtaPhiM(leptons[i]->Pt(), leptons[i]->Eta(), leptons[i]->Phi(), leptons[i]->M());
                tLEPTONS.push_back(l);
            }

            vector<RFKey> sjetID, bjetID, lepID;
            for(int i = 0; i < (int)tBJETS.size(); i++) {
                bjetID.push_back(tVIS.AddLabFrameFourVector(tBJETS[i]));
            }
            for(int i = 0; i < (int)tSJETS.size(); i++) {
                sjetID.push_back(tVIS.AddLabFrameFourVector(tSJETS[i]));
            }
            for(int i = 0; i < (int)tLEPTONS.size(); i++) {
                lepID.push_back(tVIS.AddLabFrameFourVector(tLEPTONS[i]));
            }

            // analyze
            tLAB.AnalyzeEvent();

            ////////////////////////////////////////////////////
            // Observables
            ////////////////////////////////////////////////////

            t_shat = tSS.GetMass();
            t_gaminv = tSS.GetVisibleShape();

            TVector3 tvSS = tSS.GetFourVector(tLAB).Vect();
            t_pTCM = tvSS.Pt();
            t_pZCM = fabs(tvSS.Pz());
            t_RPT = t_pTCM / (t_pTCM + t_shat /4.);
            t_RPZ = fabs(t_pZCM) / (fabs(t_pZCM) + t_shat /4.);

            t_cosSS = tSS.GetCosDecayAngle();
            t_dphiLSS = tLAB.GetDeltaPhiDecayPlanes(tSS);
            if(t_dphiLSS > acos(-1)) t_dphiLSS = t_dphiLSS - 2*acos(-1);
            t_dphiSS_SA = tSS.GetDeltaPhiDecayPlanes(tSA);
            if(t_dphiSS_SA > acos(-1)) t_dphiSS_SA = t_dphiSS_SA - 2*acos(-1);
            t_dphiSS_SB = tSS.GetDeltaPhiDecayPlanes(tSB);
            if(t_dphiSS_SB > acos(-1)) t_dphiSS_SB = t_dphiSS_SB - 2*acos(-1);
            t_dphiSA_SB = tSA.GetDeltaPhiDecayPlanes(tSB);
            if(t_dphiSA_SB > acos(-1)) t_dphiSA_SB = t_dphiSA_SB - 2*acos(-1);

            TLorentzVector pSA = tSA.GetVisibleFourVector(tSA);
            TLorentzVector pSB = tSB.GetVisibleFourVector(tSB);
            t_MS = (pSA.M2() - pSB.M2()) / (2.*(pSA.E() - pSB.E()));

            // delta phi between SS visible decay products and SS axis
            t_DPD_vSS = tSS.GetDeltaPhiDecayVisible();
            // delta phi between SS visible decay products and SS boost
            t_DPB_vSS = tSS.GetDeltaPhiBoostVisible();

            // where does everything end up?
            t_N_VA1 = tVIS.GetNElementsInFrame(tVA1);
            t_N_VA2 = tVIS.GetNElementsInFrame(tVA2);
            t_N_VB1 = tVIS.GetNElementsInFrame(tVB1);
            t_N_VB2 = tVIS.GetNElementsInFrame(tVB2);

            t_N_VA = t_N_VA1 + t_N_VA2;
            t_N_VB = t_N_VB1 + t_N_VB2;
            t_N_V1 = t_N_VA1 + t_N_VB1;
            t_N_V2 = t_N_VA2 + t_N_VB2;

            t_Nsj_VA1 = 0;
            t_Nsj_VA2 = 0;
            t_Nsj_VB1 = 0;
            t_Nsj_VB2 = 0;


            int nS = sjetID.size();
            for(int i = 0; i < nS; i++) {
                if(tVIS.GetFrame(sjetID[i]) == tVA1) {
                    t_Nsj_VA1++;
                }
                else if(tVIS.GetFrame(sjetID[i]) == tVA2) {
                    t_Nsj_VA2++;
                }
                else if(tVIS.GetFrame(sjetID[i]) == tVB1) {
                    t_Nsj_VB1++;
                }
                else if(tVIS.GetFrame(sjetID[i]) == tVB2) {
                    t_Nsj_VB2++;
                }
            }

            t_Nsj_A  = t_Nsj_VA1 + t_Nsj_VA2;
            t_Nsj_B  = t_Nsj_VB1 + t_Nsj_VB2;
            t_Nsj_V1 = t_Nsj_VA1 + t_Nsj_VB1;
            t_Nsj_V2 = t_Nsj_VA2 + t_Nsj_VB2;

            t_Nbj_VA1 = 0;
            t_Nbj_VA2 = 0;
            t_Nbj_VB1 = 0;
            t_Nbj_VB2 = 0;

            int nB = bjetID.size();
            for(int i = 0; i < nB; i++) {
                if(tVIS.GetFrame(bjetID[i]) == tVA1) {
                    t_Nbj_VA1++;
                }
                else if(tVIS.GetFrame(bjetID[i]) == tVA2) {
                    t_Nbj_VA2++;
                }
                else if(tVIS.GetFrame(bjetID[i]) == tVB1) {
                    t_Nbj_VB1++;
                }
                else if(tVIS.GetFrame(bjetID[i]) == tVB2) {
                    t_Nbj_VB2++;
                }
            }

            t_Nbj_A  = t_Nbj_VA1 + t_Nbj_VA2;
            t_Nbj_B  = t_Nbj_VB1 + t_Nbj_VB2;
            t_Nbj_V1 = t_Nbj_VA1 + t_Nbj_VB1;
            t_Nbj_V2 = t_Nbj_VA2 + t_Nbj_VB2;

            t_Nlep_VA1 = 0;
            t_Nlep_VA2 = 0;
            t_Nlep_VB1 = 0;
            t_Nlep_VB2 = 0;

            int nL = lepID.size();
            for(int i = 0; i < nL; i++) {
                if(tVIS.GetFrame(lepID[i]) == tVA1) {
                    t_Nlep_VA1++;
                }
                else if(tVIS.GetFrame(lepID[i]) == tVA2) {
                    t_Nlep_VA2++;
                }
                else if(tVIS.GetFrame(lepID[i]) == tVB1) {
                    t_Nlep_VB1++;
                }
                else if(tVIS.GetFrame(lepID[i]) == tVB2) {
                    t_Nlep_VB2++;
                }
            }
                 
            t_Nlep_A  = t_Nlep_VA1 + t_Nlep_VA2;
            t_Nlep_B  = t_Nlep_VB1 + t_Nlep_VB2;
            t_Nlep_V1 = t_Nlep_VA1 + t_Nlep_VB1;
            t_Nlep_V2 = t_Nlep_VA2 + t_Nlep_VB2;

            ///////////////////////////////
            // Scale Variables
            ///////////////////////////////

            // SS frame
            TLorentzVector vP_VA1_SS = tVA1.GetFourVector(tSS);
            TLorentzVector vP_VA2_SS = tVA2.GetFourVector(tSS);
            TLorentzVector vP_VB1_SS = tVB1.GetFourVector(tSS);
            TLorentzVector vP_VB2_SS = tVB2.GetFourVector(tSS);
            TLorentzVector vP_IA_SS  = tIA.GetFourVector(tSS);
            TLorentzVector vP_IB_SS  = tIB.GetFourVector(tSS);

            t_H_11_SS = (vP_VA1_SS + vP_VA2_SS + vP_VB1_SS + vP_VB2_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            t_H_21_SS = (vP_VA1_SS + vP_VA2_SS).P() + (vP_VB1_SS + vP_VB2_SS).P() + (vP_IA_SS + vP_IB_SS).P();
            t_H_41_SS = vP_VA1_SS.P() + vP_VA2_SS.P() + vP_VB1_SS.P() + vP_VB2_SS.P() + (vP_IA_SS + vP_IB_SS).P();
            t_H_22_SS = (vP_VA1_SS + vP_VA2_SS).P() + (vP_VB1_SS + vP_VB2_SS).P() + vP_IA_SS.P() + vP_IB_SS.P();
            t_H_42_SS = vP_VA1_SS.P() + vP_VA2_SS.P() + vP_VB1_SS.P() + vP_VB2_SS.P() + vP_IA_SS.P() + vP_IB_SS.P();

            //   > transverse
            t_H_11_SS_T = (vP_VA1_SS + vP_VA2_SS + vP_VB1_SS + vP_VB2_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            t_H_21_SS_T = (vP_VA1_SS + vP_VA2_SS).Pt() + (vP_VB1_SS + vP_VB2_SS).Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            t_H_41_SS_T = vP_VA1_SS.Pt() + vP_VA2_SS.Pt() + vP_VB1_SS.Pt() + vP_VB2_SS.Pt() + (vP_IA_SS + vP_IB_SS).Pt();
            t_H_22_SS_T = (vP_VA1_SS + vP_VA2_SS).Pt() + (vP_VB1_SS + vP_VB2_SS).Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt();
            t_H_42_SS_T = vP_VA1_SS.Pt() + vP_VA2_SS.Pt() + vP_VB1_SS.Pt() + vP_VB2_SS.Pt() + vP_IA_SS.Pt() + vP_IB_SS.Pt();

            // S1 frame
            TLorentzVector vP_VA1_SA = tVA1.GetFourVector(tSA);
            TLorentzVector vP_VA2_SA = tVA2.GetFourVector(tSA);
            TLorentzVector vP_IA_SA  = tIA.GetFourVector(tSA);

            t_H_11_SA = (vP_VA1_SA + vP_VA2_SA).P() + vP_IA_SA.P();
            t_H_21_SA = vP_VA1_SA.P() + vP_VA2_SA.P() + vP_IA_SA.P();

            //   > transverse
            t_H_11_SA_T = (vP_VA1_SA + vP_VA2_SA).Pt() + vP_IA_SA.Pt();
            t_H_21_SA_T = vP_VA1_SA.Pt() + vP_VA2_SA.Pt() + vP_IA_SA.Pt();

            //////////////////////////////////////
            // Energy variables
            //////////////////////////////////////

            t_MDR_VA1_SA = 2.0 * tVA1.GetEnergy(tSA);
            t_MDR_VA2_CA = 2.0 * tVA2.GetEnergy(tCA);


        } // 2 signal leptons and >=2 jets
    }; // end RESTFRAMES mT setup
    
    *cutflow << NewVar("t_shat"); {
        *cutflow << HFTname("t_shat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_shat;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_gaminv"); {
        *cutflow << HFTname("t_gaminv");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_gaminv;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_pTCM"); {
        *cutflow << HFTname("t_pTCM");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_pTCM;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_pZCM"); {
        *cutflow << HFTname("t_pZCM");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_pZCM;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_RPT"); {
        *cutflow << HFTname("t_RPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_RPT;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_RPZ"); {
        *cutflow << HFTname("t_RPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_RPZ;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_cosSS"); {
        *cutflow << HFTname("t_cosSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_cosSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_dphiLSS"); {
        *cutflow << HFTname("t_dphiLSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_dphiLSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_dphiSS_SA"); {
        *cutflow << HFTname("t_dphiSS_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_dphiSS_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_dphiSS_SB"); {
        *cutflow << HFTname("t_dphiSS_SB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_dphiSS_SB;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_dphiSA_SB"); {
        *cutflow << HFTname("t_dphiSA_SB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_dphiSA_SB;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_MS"); {
        *cutflow << HFTname("t_MS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_MS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_DPD_vSS"); {
        *cutflow << HFTname("t_DPD_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_DPD_vSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_DPB_vSS"); {
        *cutflow << HFTname("t_DPB_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_DPB_vSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_VA1"); {
        *cutflow << HFTname("t_N_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_VA2"); {
        *cutflow << HFTname("t_N_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_VB1"); {
        *cutflow << HFTname("t_N_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_VB2"); {
        *cutflow << HFTname("t_N_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_VA"); {
        *cutflow << HFTname("t_N_VA");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_VA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_VB"); {
        *cutflow << HFTname("t_N_VB");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_VB;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_V1"); {
        *cutflow << HFTname("t_N_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_N_V2"); {
        *cutflow << HFTname("t_N_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_N_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_VA1"); {
        *cutflow << HFTname("t_Nsj_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_VA2"); {
        *cutflow << HFTname("t_Nsj_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_VB1"); {
        *cutflow << HFTname("t_Nsj_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_VB2"); {
        *cutflow << HFTname("t_Nsj_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_VA1"); {
        *cutflow << HFTname("t_Nbj_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_VA2"); {
        *cutflow << HFTname("t_Nbj_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_VB1"); {
        *cutflow << HFTname("t_Nbj_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_VB2"); {
        *cutflow << HFTname("t_Nbj_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_VA1"); {
        *cutflow << HFTname("t_Nlep_VA1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_VA1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_VA2"); {
        *cutflow << HFTname("t_Nlep_VA2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_VA2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_VB1"); {
        *cutflow << HFTname("t_Nlep_VB1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_VB1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_VB2"); {
        *cutflow << HFTname("t_Nlep_VB2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_VB2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_A"); {
        *cutflow << HFTname("t_Nsj_A");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_A;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_B"); {
        *cutflow << HFTname("t_Nsj_B");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_B;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_V1"); {
        *cutflow << HFTname("t_Nsj_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nsj_V2"); {
        *cutflow << HFTname("t_Nsj_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nsj_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_A"); {
        *cutflow << HFTname("t_Nbj_A");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_A;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_B"); {
        *cutflow << HFTname("t_Nbj_B");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_B;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_V1"); {
        *cutflow << HFTname("t_Nbj_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nbj_V2"); {
        *cutflow << HFTname("t_Nbj_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nbj_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_A"); {
        *cutflow << HFTname("t_Nlep_A");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_A;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_B"); {
        *cutflow << HFTname("t_Nlep_B");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_B;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_V1"); {
        *cutflow << HFTname("t_Nlep_V1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_V1;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_Nlep_V2"); {
        *cutflow << HFTname("t_Nlep_V2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -1;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_Nlep_V2;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_11_SS"); {
        *cutflow << HFTname("t_H_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_11_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_21_SS"); {
        *cutflow << HFTname("t_H_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_21_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_41_SS"); {
        *cutflow << HFTname("t_H_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_41_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_22_SS"); {
        *cutflow << HFTname("t_H_22_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_22_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_42_SS"); {
        *cutflow << HFTname("t_H_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_42_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_11_SS_T"); {
        *cutflow << HFTname("t_H_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_11_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_21_SS_T"); {
        *cutflow << HFTname("t_H_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_21_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_41_SS_T"); {
        *cutflow << HFTname("t_H_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_41_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_22_SS_T"); {
        *cutflow << HFTname("t_H_22_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_22_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_42_SS_T"); {
        *cutflow << HFTname("t_H_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_42_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_11_SA"); {
        *cutflow << HFTname("t_H_11_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_11_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_21_SA"); {
        *cutflow << HFTname("t_H_21_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_21_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_11_SA_T"); {
        *cutflow << HFTname("t_H_11_SA_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_11_SA_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_H_21_SA_T"); {
        *cutflow << HFTname("t_H_21_SA_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_H_21_SA_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_MDR_VA1_SA"); {
        *cutflow << HFTname("t_MDR_VA1_SA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_MDR_VA1_SA;
            }
            return out;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("t_MDR_VA2_CA"); {
        *cutflow << HFTname("t_MDR_VA2_CA");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            int out = -99.0;
            if(jets.size()>=2 && leptons.size()==2) {
                out = t_MDR_VA2_CA;
            }
            return out;
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
    SuperflowRunMode& run_mode_, SusyNtSys& nt_sys, bool dbg)
{
    bool nominal_ = true;
    string input;

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-n") == 0)
            num_events_ = atoi(argv[++i]);
        else if(strcmp(argv[i], "-i") == 0) {
            input = argv[++i];
        }
        else if(strcmp(argv[i], "-d") == 0) {
            dbg = true;
        }
        else {
            cout << analysis_name << "    Error (fatal) : Bad arguments." << endl;
            exit(1);
        }
    } // i

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
}
