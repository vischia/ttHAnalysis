#include <string>

#include <cp3_llbb/HHAnalysis/interface/Indices.h>

namespace HHAnalysis {
  
  // Combination of Lepton ID + Lepton Isolation for a single lepton
  uint16_t LepIDIso(const LepID::LepID& id, const LepIso::LepIso& iso){
    return LepIso::Count * id + iso;
  }
  std::string LepIDIsoStr(const LepID::LepID& id, const LepIso::LepIso& iso){
    return "ID" + LepID::map.at(id) + "_Iso" + LepIso::map.at(iso);
  }

  // Combination of Lepton ID for a DiLepton object
  uint16_t LepLepID(const LepID::LepID& id1, const LepID::LepID& id2){
    return LepID::Count * id1 +  id2;
  }
  std::string LepLepIDStr(const LepID::LepID& id1, const LepID::LepID& id2){
    return LepID::map.at(id1) + LepID::map.at(id2);
  }

  // Combination of Lepton Isolation for a DiLepton object
  uint16_t LepLepIso(const LepIso::LepIso& iso1, const LepIso::LepIso& iso2){
    return LepIso::Count * iso1 + iso2;
  }
  std::string LepLepIsoStr(const LepIso::LepIso& iso1, const LepIso::LepIso& iso2){
    return LepIso::map.at(iso1) + LepIso::map.at(iso2);
  }

  // Combination of Lepton ID + Lepton Isolation for a DiLepton object
  uint16_t LepLepIDIso(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2){
    return LepIso::Count*LepID::Count*LepIso::Count * id1 + LepID::Count*LepIso::Count * iso1 + LepIso::Count * id2 + iso2;
  }
  std::string LepLepIDIsoStr(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2){
    return "ID" + LepID::map.at(id1) + LepID::map.at(id2) + "_Iso" + LepIso::map.at(iso1) + LepIso::map.at(iso2);
  }

  // Combination of Jet IDs for two jets (NOTE: NOT USED FOR NOW)
  uint16_t JetJetID(const JetID::JetID& id1, const JetID::JetID& id2){
    return JetID::Count * id1 + id2;
  }
  std::string JetJetIDStr(const JetID::JetID& id1, const JetID::JetID& id2){
    return JetID::map.at(id1) + JetID::map.at(id2);
  }
  
  // Combination of Jet ID and B-tagging working point (NOTE: NOT USED FOR NOW)
  uint16_t JetIDBWP(const JetID::JetID& id, const BWP::BWP& wp){
    return BWP::Count * id + wp;
  }
  std::string JetIDBWPStr(const JetID::JetID& id, const BWP::BWP& wp){
    return "ID" + JetID::map.at(id) + "_B" + BWP::map.at(wp);
  }
  
  // Combination of Lepton ID + Lepton Isolation (one lepton) and B-tagging working point for one jet
  uint16_t LepIDIsoJetBWP(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp){
    return LepIso::Count*BWP::Count * id + BWP::Count * iso + wp;
  }
  std::string LepIDIsoJetBWPStr(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp){
    return "ID" + LepID::map.at(id) + "_Iso" + LepIso::map.at(iso) + "_B" + BWP::map.at(wp);
  }

  // Combination of B-tagging working points for two jets
  uint16_t JetJetBWP(const BWP::BWP& wp1, const BWP::BWP& wp2){
    return BWP::Count * wp1 + wp2;
  }
  std::string JetJetBWPStr(const BWP::BWP& wp1, const BWP::BWP& wp2){
    return BWP::map.at(wp1) + BWP::map.at(wp2);
  }

  // Combination of Jet ID and B-tagging working points for two jets 
  uint16_t JetJetIDBWP(const JetID::JetID& id1, const BWP::BWP& wp1, const JetID::JetID& id2, const BWP::BWP& wp2){
    return 
      BWP::Count*JetID::Count*BWP::Count * id1 + 
                 JetID::Count*BWP::Count * wp1 + 
                              BWP::Count * id2 + 
                                           wp2 ;
  }
  std::string JetJetIDBWPStr(const JetID::JetID& id1, const JetID::JetID& id2, const BWP::BWP wp1, const BWP::BWP wp2){
    return "ID" + JetID::map.at(id1) + JetID::map.at(id2) + "_B" + BWP::map.at(wp1) + BWP::map.at(wp2);
  }
  
  // Combination of Lepton ID + Lepton Isolation (one lepton) and B-tagging working points for two jets
  uint16_t LepIDIsoJetJetBWP(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp1, const BWP::BWP& wp2){
    return LepIso::Count*BWP::Count*BWP::Count * id + BWP::Count*BWP::Count * iso + BWP::Count * wp1 + wp2;
  }
  std::string LepIDIsoJetJetBWPStr(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp1, const BWP::BWP& wp2){
    return "ID" + LepID::map.at(id) + "_Iso" + LepIso::map.at(iso) + "_B" + BWP::map.at(wp1) + BWP::map.at(wp2);
  }
  
  // Combination of Lepton ID, Lepton Isolation, Jet ID, B-tagging working points and jetPair ordering for a two-lepton-two-b-jets object
  uint16_t LepLepIDIsoJetJetIDBWPPair(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2, const JetID::JetID& jetid1, const JetID::JetID& jetid2, const BWP::BWP& wp1, const BWP::BWP& wp2, const JetPair::JetPair& jetpair){
    return 
      LepIso::Count*LepID::Count*LepIso::Count*JetID::Count*JetID::Count*BWP::Count*BWP::Count*JetPair::Count * id1     + 
                    LepID::Count*LepIso::Count*JetID::Count*JetID::Count*BWP::Count*BWP::Count*JetPair::Count * iso1    + 
                                 LepIso::Count*JetID::Count*JetID::Count*BWP::Count*BWP::Count*JetPair::Count * id2     + 
                                               JetID::Count*JetID::Count*BWP::Count*BWP::Count*JetPair::Count * iso2    + 
                                                            JetID::Count*BWP::Count*BWP::Count*JetPair::Count * jetid1  + 
                                                                         BWP::Count*BWP::Count*JetPair::Count * jetid2  + 
                                                                                    BWP::Count*JetPair::Count * wp1     + 
                                                                                               JetPair::Count * wp2     + 
                                                                                                                jetpair ; 
  }
  std::string LepLepIDIsoJetJetIDBWPPairStr(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2, const JetID::JetID& jetid1, const JetID::JetID& jetid2, const BWP::BWP& wp1, const BWP::BWP& wp2, const JetPair::JetPair& jetpair){
    return "Lep_ID" + LepID::map.at(id1) + LepID::map.at(id2) + "_Iso" + LepIso::map.at(iso1) + LepIso::map.at(iso2) + "_Jet_ID" + JetID::map.at(jetid2) + JetID::map.at(jetid2) + "_B" + BWP::map.at(wp1) + BWP::map.at(wp2) + "_Ordered" + JetPair::map.at(jetpair);
  }

}
