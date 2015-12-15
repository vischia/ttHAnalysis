#pragma once

#include <map>
#include <array>
#include <string>

namespace HHAnalysis {
  
  // Lepton IDs
  namespace LepID {
    enum LepID{ L, T, Count };
    // Ugly way to allow iterating over all items in the enumeration ( for(const LepID::LepID& id: LepID::it) )
    const std::array<LepID, Count> it = {{ L, T }};
    // Is useful in categories to construct cut strings out of each working point
    const std::map<LepID, std::string> map = { {L, "L"}, {T, "T"} };
  }
  
  // Lepton Isolation
  namespace LepIso {
    enum LepIso{ no, L, T, Count };
    const std::array<LepIso, Count> it = {{ no, L, T }};
    const std::map<LepIso, std::string> map = { {no, "no"}, {L, "L"}, {T, "T"} };
  }

  // Combination of Lepton ID + Lepton Isolation for a single lepton
  uint16_t LepIDIso(const LepID::LepID& id, const LepIso::LepIso& iso);
  std::string LepIDIsoStr(const LepID::LepID& id, const LepIso::LepIso& iso);

  // Combination of Lepton ID for a DiLepton object
  uint16_t LepLepID(const LepID::LepID& id1, const LepID::LepID& id2);
  std::string LepLepIDStr(const LepID::LepID& id1, const LepID::LepID& id2);

  // Combination of Lepton Isolation for a DiLepton object
  uint16_t LepLepIso(const LepIso::LepIso& iso1, const LepIso::LepIso& iso2);
  std::string LepLepIsoStr(const LepIso::LepIso& iso1, const LepIso::LepIso& iso2);

  // Combination of Lepton ID + Lepton Isolation for a DiLepton object
  uint16_t LepLepIDIso(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2);
  std::string LepLepIDIsoStr(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2);

  // Jet ID
  namespace JetID {
    enum JetID{ L, T, TLV, no, Count };
    const std::array<JetID, Count> it = {{ L, T, TLV, no }};
    const std::map<JetID, std::string> map = { {L, "L"}, {T, "T"}, {TLV, "TLV"}, {no, "no"} };
  }

  // Combination of Jet IDs for two jets (NOTE: NOT USED FOR NOW)
  uint16_t JetJetID(const JetID::JetID& id1, const JetID::JetID& id2);
  std::string JetJetIDStr(const JetID::JetID& id1, const JetID::JetID& id2);
  
  // B-tagging working points
  namespace BWP {
    enum BWP{ no, L, M, T, Count };
    const std::array<BWP, Count> it = {{ no, L, M, T }};
    const std::map<BWP, std::string> map = { {no, "no"}, {L, "L"}, {M, "M"}, {T, "T"} };
  }

  // Jet combinatoric
  namespace JetPair {
    enum JetPair { ht, mh, pt, csv, jp, ptOverM, Count };
    const std::array<JetPair, Count> it = {{ ht, mh, pt, csv, jp, ptOverM }};
    const std::map<JetPair, std::string> map = { {ht, "ht"}, {mh, "mh"}, {pt, "pt"}, {csv, "csv"}, {jp, "jp"}, {ptOverM, "ptOverM"} }; 
  }

  // Combination of Jet ID and B-tagging working point
  uint16_t JetIDBWP(const JetID::JetID& id, const BWP::BWP& wp);
  std::string JetIDBWPStr(const JetID::JetID& id, const BWP::BWP& wp);
  
  // Combination of Lepton ID + Lepton Isolation (one lepton) and B-tagging working point for one jet
  uint16_t LepIDIsoJetBWP(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp);
  std::string LepIDIsoJetBWPStr(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp);

  // Combination of B-tagging working points for two jets
  uint16_t JetJetBWP(const BWP::BWP& wp1, const BWP::BWP& wp2);
  std::string JetJetBWPStr(const BWP::BWP& wp1, const BWP::BWP& wp2);

  // Combination of Jet ID and B-tagging working points for two jets 
  uint16_t JetJetIDBWP(const JetID::JetID& id1, const JetID::JetID& id2, const BWP::BWP& wp1, const BWP::BWP& wp2);
  std::string JetJetIDBWPStr(const JetID::JetID& id1, const JetID::JetID& id2, const BWP::BWP& wp1, const BWP::BWP& wp2);

  // Combination of Lepton ID, Lepton Isolation, Jet ID, B-tagging working points and jetPair ordering for a two-lepton-two-b-jets object
  uint16_t LepLepIDIsoJetJetIDBWPPair(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2, const JetID::JetID& jetid1, const JetID::JetID& jetid2, const BWP::BWP& wp1, const BWP::BWP& wp2, const JetPair::JetPair& jetpair);
  std::string LepLepIDIsoJetJetIDBWPPairStr(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2, const JetID::JetID& jetid1, const JetID::JetID& jetid2, const BWP::BWP& wp1, const BWP::BWP& wp2, const JetPair::JetPair& jetpair);

}
