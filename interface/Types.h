#pragma once

#include <vector>
#include <Math/Vector4D.h>
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> LorentzVector;

namespace HH {
    typedef std::vector< std::vector<unsigned int> > MapType;
    struct Lepton {
        LorentzVector p4;
        int8_t charge;
        unsigned int idx;
        bool isMu;
        bool isEl;
        bool isID_L; // Loose
        bool isID_T; // Tight
    };
    struct Dilepton {
        LorentzVector p4;
        std::pair<unsigned int, unsigned int> idxs; // indices in the collection of HH::Lepton
        unsigned int ilep1; // index in the corresponding framework collection
        unsigned int ilep2; // index in the corresponding framework collection
        bool isOS; // Opposite Sign
        bool isMuMu;
        bool isElEl;
        bool isElMu;
        bool isMuEl;
        bool isSF; // Same Flavour
        bool isID_LL;
        bool isID_LT;
        bool isID_TL;
        bool isID_TT;
        float DR;
        float DPhi;
    };
    struct Met {
        LorentzVector p4;
        bool isNoHF;
    };
    struct DileptonMet : Dilepton {
        LorentzVector p4;
        unsigned int ill; // index in the HH::Dilepton collection
        unsigned int imet; // index in the HH::Met collection
        bool isNoHF;
        float DPhi_ll_met;
        float minDPhi_l_met;
        float maxDPhi_l_met;
        float MT;
        float MT_formula;
        float projectedMet;
    };
}


