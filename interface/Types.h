#pragma once

#include <Math/Vector4D.h>
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> LorentzVector;

struct Lepton { 
    LorentzVector p4; 
    int8_t charge;
    unsigned int idx; 
    bool isMu; 
    bool isEl; 
};  

struct Dilepton { 
    LorentzVector p4; 
    std::pair<unsigned int, unsigned int> idxs; 
    int chargeProduct;
    bool isMuMu; 
    bool isElEl; 
    bool isElMu; 
    bool isMuEl; 
};  


