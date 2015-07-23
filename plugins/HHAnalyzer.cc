#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
// from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
//    enum StatusBits {
//0      kIsPrompt = 0,
//1      kIsDecayedLeptonHadron,
//2      kIsTauDecayProduct,
//3      kIsPromptTauDecayProduct,
//4      kIsDirectTauDecayProduct,
//5      kIsDirectPromptTauDecayProduct,
//6      kIsDirectHadronDecayProduct,
//7      kIsHardProcess,
//8      kFromHardProcess,
//9      kIsHardProcessTauDecayProduct,
//10      kIsDirectHardProcessTauDecayProduct,
//11      kFromHardProcessBeforeFSR,
//12      kIsFirstCopy,
//13      kIsLastCopy,
//14      kIsLastCopyBeforeFSR
//    };

void HHAnalyzer::analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager& producers) {

    const GenParticlesProducer& gp = dynamic_cast<const GenParticlesProducer&>(producers.get("gen_particles"));
    const JetsProducer& jets = dynamic_cast<const JetsProducer&>(producers.get("jets"));

/*    LorentzVector gen_ll(0.,0.,0.,0.);
    LorentzVector gen_llFSR(0.,0.,0.,0.);
    LorentzVector gen_met(0.,0.,0.,0.);
    LorentzVector gen_llmet(0.,0.,0.,0.);
    LorentzVector gen_llFSRmet(0.,0.,0.,0.);
*/
    int iHiggs1, iHiggs2;
    iHiggs1 = iHiggs2 = -1;
/*
    int iB1, iB2;
    int iJet1, iJet2;
    int iL1, iL2;
    int iNu1, iNu2;
    iB1 = iB2 = -1;
    iJet1 = iJet2 = -1;
    iL1 = iL2 = -1;
    iNu1 = iNu2 = -1;
*/
    // Construct the signal gen-level info
    // Get the Higgses !
    for (unsigned int ip = 0 ; ip < gp.pruned_p4.size() ; ip++) {
        if( gp.pruned_pdg_id[ip] != 25 ) continue; // just take care of the Higgses for now
        std::bitset<15> flags (gp.pruned_status_flags[ip]);
        if( flags.test(13) ) continue; // must be the last copy
        if( iHiggs1 == -1 )
            iHiggs1 = ip;
        else if( iHiggs2 == -1 )
            iHiggs2 = ip;
        else
            std::cout << "ERROR, more than 3 Higgses in the event!" << std::endl;
    }
    std::cout << "\tiHiggs1= " << iHiggs1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iHiggs1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iHiggs1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iHiggs1].Pt() << " , " << gp.pruned_p4[iHiggs1].Eta() << " , " << gp.pruned_p4[iHiggs1].Phi() << " , " << gp.pruned_p4[iHiggs1].E() << " , " << gp.pruned_p4[iHiggs1].M() << ")" << std::endl;
    std::cout << "\tiHiggs2= " << iHiggs2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iHiggs2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iHiggs2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iHiggs2].Pt() << " , " << gp.pruned_p4[iHiggs2].Eta() << " , " << gp.pruned_p4[iHiggs2].Phi() << " , " << gp.pruned_p4[iHiggs2].E() << " , " << gp.pruned_p4[iHiggs2].M() << ")" << std::endl;

    // Get the Higgses decay products
    for (unsigned int ip = 0 ; ip < gp.pruned_p4.size() ; ip++) {
//        if( gp.pruned_pdg_id[ip] == 25 ) continue; // skip the Higgs this time
        std::bitset<15> flags (gp.pruned_status_flags[ip]);
        if( !(flags.test(13) && flags.test(8)) ) continue; // take the last copies coming from the hard process
        std::cout << "\t\tip= " << ip << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[ip] << "\tflags= " << flags << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[ip].Pt() << " , " << gp.pruned_p4[ip].Eta() << " , " << gp.pruned_p4[ip].Phi() << " , " << gp.pruned_p4[ip].E() << " , " << gp.pruned_p4[ip].M() << ")" << std::endl;
        for(unsigned int imother = 0 ; imother < gp.pruned_mothers_index[ip].size() ; imother++)
            std::cout << "\t\t\timother= " << imother << "\tindex= " << gp.pruned_mothers_index[ip][imother] << std::endl;
    }
/*
    for (unsigned int ip = 0 ; ip < gp.pruned_p4.size() ; ip++) {
        std::bitset<15> flags (gp.pruned_status_flags[ip]);
        if( !(gp.pruned_pdg_id[ip] == 25 || abs(gp.pruned_pdg_id[ip]) == 24) ) continue;
//        if( !((flags.test(0) == 1 || flags.test(7) == 1 || flags.test(8) == 1) && (flags.test(13) == 1)) ) continue; // only care of prompt from hard process, only care of the last copy (after FSR and all that)
//        if( abs(gp.pruned_pdg_id[ip]) < 11 || (abs(gp.pruned_pdg_id[ip]) > 16 && gp.pruned_pdg_id[ip]!= 22) ) continue; // leptons and neutrinos
        std::cout << "ip= " << ip << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[ip] << "\tflags= " << flags << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[ip].Pt() << " , " << gp.pruned_p4[ip].Eta() << " , " << gp.pruned_p4[ip].Phi() << " , " << gp.pruned_p4[ip].E() << " , " << gp.pruned_p4[ip].M() << ")" << std::endl;
        if( abs(gp.pruned_pdg_id[ip]) == 11 || abs(gp.pruned_pdg_id[ip]) == 13 )
        {
            gen_ll += gp.pruned_p4[ip];
            gen_llFSR += gp.pruned_p4[ip];
        }
        else if( gp.pruned_pdg_id[ip] == 22 )
        {
            gen_llFSR += gp.pruned_p4[ip];
        }
        else
        {
            gen_met += gp.pruned_p4[ip];
        }
        gen_llmet = gen_ll + gen_met;
        gen_llFSRmet = gen_llFSR + gen_met;
    }
*/
/*    std::cout << "\tGEN LL(pt, eta, phi, e, mass)= (" << gen_ll.Pt() << " , " << gen_ll.Eta() << " , " << gen_ll.Phi() << " , " << gen_ll.E() << " , " << gen_ll.M() << ")" << std::endl;
    std::cout << "\tGEN LL + FSR(pt, eta, phi, e, mass)= (" << gen_llFSR.Pt() << " , " << gen_llFSR.Eta() << " , " << gen_llFSR.Phi() << " , " << gen_llFSR.E() << " , " << gen_llFSR.M() << ")" << std::endl;
    std::cout << "\tGEN MET(pt, eta, phi, e, mass)= (" << gen_met.Pt() << " , " << gen_met.Eta() << " , " << gen_met.Phi() << " , " << gen_met.E() << " , " << gen_met.M() << ")" << std::endl;
    std::cout << "\tGEN LLMET(pt, eta, phi, e, mass)= (" << gen_llmet.Pt() << " , " << gen_llmet.Eta() << " , " << gen_llmet.Phi() << " , " << gen_llmet.E() << " , " << gen_llmet.M() << ")" << std::endl;
    std::cout << "\tGEN LLMET + FSR(pt, eta, phi, e, mass)= (" << gen_llFSRmet.Pt() << " , " << gen_llFSRmet.Eta() << " , " << gen_llFSRmet.Phi() << " , " << gen_llFSRmet.E() << " , " << gen_llFSRmet.M() << ")" << std::endl;
*/
/*
    for (auto p4: jets.p4) {
        std::cout << "Jet pt: " << p4.P() << std::endl;
    }
*/
}
