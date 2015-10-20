#pragma once

#include <vector>
#include <Math/Vector4D.h>
#include <cp3_llbb/HHAnalysis/interface/Enums.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> LorentzVector;

namespace HH {
    struct Lepton {
        LorentzVector p4;
        int8_t charge;
        unsigned int idx;
        int8_t hlt_idx = -1; // Index to the matched HLT object. -1 if no match. 
                             // Example : t->Draw("hh_leptons.p4.Pt() - hlt_object_p4[hh_leptons.hlt_idx].Pt()","hh_leptons.hlt_idx != -1","")
        bool hlt_already_tried_matching = false; // do the matching only once, even when the lepton is in several Dilepton
        bool isMu;
        bool isEl;
        bool id_L; // Loose
        bool id_T; // Tight
        bool iso_L; // Loose
        bool iso_T; // Tight
    };
    struct Dilepton {
        LorentzVector p4;
        std::pair<unsigned int, unsigned int> idxs; // indices in the collection of HH::Lepton
        unsigned int ilep1; // index in the corresponding framework collection
        unsigned int ilep2; // index in the corresponding framework collection
        std::pair<int8_t, int8_t> hlt_idxs = std::make_pair(-1,-1); // Stores indices of matched online objects. (-1,-1) if no match
        bool isOS; // Opposite Sign
        bool isMuMu;
        bool isElEl;
        bool isElMu;
        bool isMuEl;
        bool isSF; // Same Flavour
        bool id_LL;
        bool id_LT;
        bool id_TL;
        bool id_TT;
        bool iso_LL;
        bool iso_LT;
        bool iso_TL;
        bool iso_TT;
        float DR_l_l;
        float DPhi_l_l;
    };
    struct Met {
        LorentzVector p4;
        bool isNoHF;
    };
    struct DileptonMet : public Dilepton, public Met {
        LorentzVector p4;
        unsigned int ill; // index in the HH::Dilepton collection
        unsigned int imet; // index in the HH::Met collection
        float DPhi_ll_met;
        float minDPhi_l_met;
        float maxDPhi_l_met;
        float MT;
        float MT_formula;
        float projectedMet;
    };
    struct Jet {
        LorentzVector p4;
        unsigned idx;
        bool btagL;
        bool btagM;
        bool btagT;
    };
    struct Dijet {
        LorentzVector p4;
        std::pair<unsigned int, unsigned int> idxs; // indices in the collection of HH::Jet
        unsigned int ijet1; // indices in the framework collection
        unsigned int ijet2;
        bool btag_LL;
        bool btag_LM;
        bool btag_LT;
        bool btag_ML;
        bool btag_MM;
        bool btag_MT;
        bool btag_TL;
        bool btag_TM;
        bool btag_TT;
        float DR_j_j;
        float DPhi_j_j;
    };
    struct DileptonMetDijet : public DileptonMet, public Dijet {
        LorentzVector p4;
        unsigned int illmet; // index in the HH::DileptonMet collection
        unsigned int ijj; // index in the HH::Dijet collection
        unsigned int imet; // index in the HH::Met collection
        float DPhi_jj_met;
        float minDPhi_j_met;
        float maxDPhi_j_met;
        float maxDR_l_j;
        float minDR_l_j;
        float DR_ll_jj;
        float DPhi_ll_jj;
        float DR_llmet_jj;
        float DPhi_llmet_jj;
        float cosThetaStar_CS;
    };
}


