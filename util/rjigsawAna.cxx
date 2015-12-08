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
    cutflow->setLumi(3340); // 3.34/fb
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

    *cutflow << NewVar("is MC"); {
        *cutflow << HFTname("isMC");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
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
    double shat;
    double pTT_T;
    double pTT_Z;
    double RPT;
    double RPZ;
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



        /// system mass
        shat = ss.GetMass();

        //////////////////////
        // RATIO OF CM pT
        TVector3 vPTT = ss.GetFourVector(lab).Vect();
        pTT_T = vPTT.Pt();
        pTT_Z = vPTT.Pz();
        RPT = vPTT.Pt() / (vPTT.Pt() + shat / 4.);
        RPZ = vPTT.Pz() / (vPTT.Pz() + shat / 4.);


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
    

    *cutflow << [&](Superlink* sl, var_void*) {
        if(jets.size()>=2) {

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
            TLorentzVector j1, j2, l1, l2;
            j1.SetPtEtaPhiM(jets.at(0)->Pt(), jets.at(0)->Eta(), jets.at(0)->Phi(), jets.at(0)->M());
            j2.SetPtEtaPhiM(jets.at(1)->Pt(), jets.at(1)->Eta(), jets.at(1)->Phi(), jets.at(1)->M());
            l1.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            l2.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Eta(), leptons.at(1)->M());

            vector<RFKey> jetID;
            jetID.push_back(xvis.AddLabFrameFourVector(j1));
            jetID.push_back(xvis.AddLabFrameFourVector(j2));
            jetID.push_back(xvis.AddLabFrameFourVector(l1));
            jetID.push_back(xvis.AddLabFrameFourVector(l2)); 


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
            xRPZ = xvPSS.Pz() / (xvPSS.Pz() + xshat/4.);
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
        

        } // njets == 2
    };

    *cutflow << NewVar("xshat"); {
        *cutflow << HFTname("xshat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xshat;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xgaminv"); {
        *cutflow << HFTname("xgaminv");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xgaminv;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPT"); {
        *cutflow << HFTname("xRPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xRPT;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRPZ"); {
        *cutflow << HFTname("xRPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xRPZ;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosSS"); {
        *cutflow << HFTname("xcosSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xcosSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiLSS"); {
        *cutflow << HFTname("xdphiLSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xdphiLSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xMS"); {
        *cutflow << HFTname("xMS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xMS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xPS"); {
        *cutflow << HFTname("xPS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xPS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xMSS"); {
        *cutflow << HFTname("xMSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xMSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xgaminvSS"); {
        *cutflow << HFTname("xgaminvSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xgaminvSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDeltaBetaSS"); {
        *cutflow << HFTname("xDeltaBetaSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xDeltaBetaSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDPD_vSS"); {
        *cutflow << HFTname("xDPD_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xDPD_vSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xDPB_vSS"); {
        *cutflow << HFTname("xDPB_vSS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xDPB_vSS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV"); {
        *cutflow << HFTname("xNV");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999);
            out.push_back(-999);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xNV[0]);
                out.push_back(xNV[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV_0"); {
        *cutflow << HFTname("xNV_0");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -999;
            if(jets.size()>=2) {
                out = xNV[0];
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xNV_1"); {
        *cutflow << HFTname("xNV_1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            int out = -999;
            if(jets.size()>=2) {
                out = xNV[1];
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosS"); {
        *cutflow << HFTname("xcosS");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xcosS[0]);
                out.push_back(xcosS[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xcosC"); {
        *cutflow << HFTname("xcosC");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xcosC[0]);
                out.push_back(xcosC[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xdphiSC"); {
        *cutflow << HFTname("xdphiSC");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xdphiSC[0]);
                out.push_back(xdphiSC[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xRCS"); {
        *cutflow << HFTname("xRCS");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xRCS[0]);
                out.push_back(xRCS[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xjet1PT"); {
        *cutflow << HFTname("xjet1PT");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xjet1PT[0]);
                out.push_back(xjet1PT[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xjet2PT"); {
        *cutflow << HFTname("xjet2PT");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xjet2PT[0]);
                out.push_back(xjet2PT[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xPinv"); {
        *cutflow << HFTname("xPinv");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            out.push_back(-999.);
            out.push_back(-999.);
            if(jets.size()>=2) {
                out.clear();
                out.push_back(xPinv[0]);
                out.push_back(xPinv[1]);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_SS"); {
        *cutflow << HFTname("xH_11_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_11_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_SS"); {
        *cutflow << HFTname("xH_21_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_21_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_41_SS"); {
        *cutflow << HFTname("xH_41_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_41_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_42_SS"); {
        *cutflow << HFTname("xH_42_SS");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_42_SS;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_S1"); {
        *cutflow << HFTname("xH_11_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_11_S1;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_S1"); {
        *cutflow << HFTname("xH_21_S1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_21_S1;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_SS_T"); {
        *cutflow << HFTname("xH_11_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_11_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_SS_T"); {
        *cutflow << HFTname("xH_21_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_21_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_41_SS_T"); {
        *cutflow << HFTname("xH_41_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_41_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_42_SS_T"); {
        *cutflow << HFTname("xH_42_SS_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_42_SS_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_11_S1_T"); {
        *cutflow << HFTname("xH_11_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_11_S1_T;
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xH_21_S1_T"); {
        *cutflow << HFTname("xH_21_S1_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999.;
            if(jets.size()>=2) {
                out = xH_21_S1_T;
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
