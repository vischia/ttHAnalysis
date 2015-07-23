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
    int iH1, iH2;
    int iB1, iB2;
    int iL1, iL2;
    int iNu1, iNu2;
    iH1 = iH2 = -1;
    iB1 = iB2 = -1;
    iL1 = iL2 = -1;
    iNu1 = iNu2 = -1;

    // Construct the signal gen-level info
    for (unsigned int ip = 0 ; ip < gp.pruned_p4.size() ; ip++) {
        std::bitset<15> flags (gp.pruned_status_flags[ip]);
        if( !(flags.test(13) && flags.test(8)) ) continue; // take the last copies coming from the hard process
        if( abs(gp.pruned_pdg_id[ip]) == 25 )
        {
            if( iH1 == -1 )
                iH1 = ip;
            else if( iH2 == -1 )
                iH2 = ip;
            else
                std::cout << "ERROR, more than 2 Higgs in the hard process!" << std::endl;
        }
        else if( abs(gp.pruned_pdg_id[ip]) == 5 )
        {
            if( iB1 == -1 )
                iB1 = ip;
            else if( iB2 == -1 )
                iB2 = ip;
            else
                std::cout << "ERROR, more than 2 bjets in the hard process!" << std::endl;
        }
        else if( abs(gp.pruned_pdg_id[ip]) == 11 || abs(gp.pruned_pdg_id[ip]) == 13 )
        {
            if( iL1 == -1 )
                iL1 = ip;
            else if( iL2 == -1 )
                iL2 = ip;
            else
                std::cout << "ERROR, more than 2 leptons in the hard process!" << std::endl;
        }
        else if( abs(gp.pruned_pdg_id[ip]) == 12 || abs(gp.pruned_pdg_id[ip]) == 14 || abs(gp.pruned_pdg_id[ip]) == 16 )
        {
            if( iNu1 == -1 )
                iNu1 = ip;
            else if( iNu2 == -1 )
                iNu2 = ip;
            else
                std::cout << "ERROR, more than 2 neutrinos in the hard process!" << std::endl;
        }
//FIXME if the lepton is not in the hard process itself, then look for FSR photons
//FIXME keep radion / graviton / W / Z when present
    } // end of loop over pruned gen particles
    std::cout << "\tiH1= " << iH1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iH1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iH1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iH1].Pt() << " , " << gp.pruned_p4[iH1].Eta() << " , " << gp.pruned_p4[iH1].Phi() << " , " << gp.pruned_p4[iH1].E() << " , " << gp.pruned_p4[iH1].M() << ")" << std::endl;
    std::cout << "\tiH2= " << iH2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iH2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iH2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iH2].Pt() << " , " << gp.pruned_p4[iH2].Eta() << " , " << gp.pruned_p4[iH2].Phi() << " , " << gp.pruned_p4[iH2].E() << " , " << gp.pruned_p4[iH2].M() << ")" << std::endl;
    std::cout << "\tiB1= " << iB1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iB1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iB1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iB1].Pt() << " , " << gp.pruned_p4[iB1].Eta() << " , " << gp.pruned_p4[iB1].Phi() << " , " << gp.pruned_p4[iB1].E() << " , " << gp.pruned_p4[iB1].M() << ")" << std::endl;
    std::cout << "\tiB2= " << iB2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iB2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iB2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iB2].Pt() << " , " << gp.pruned_p4[iB2].Eta() << " , " << gp.pruned_p4[iB2].Phi() << " , " << gp.pruned_p4[iB2].E() << " , " << gp.pruned_p4[iB2].M() << ")" << std::endl;
    std::cout << "\tiL1= " << iL1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iL1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iL1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iL1].Pt() << " , " << gp.pruned_p4[iL1].Eta() << " , " << gp.pruned_p4[iL1].Phi() << " , " << gp.pruned_p4[iL1].E() << " , " << gp.pruned_p4[iL1].M() << ")" << std::endl;
    std::cout << "\tiL2= " << iL2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iL2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iL2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iL2].Pt() << " , " << gp.pruned_p4[iL2].Eta() << " , " << gp.pruned_p4[iL2].Phi() << " , " << gp.pruned_p4[iL2].E() << " , " << gp.pruned_p4[iL2].M() << ")" << std::endl;
    std::cout << "\tiNu1= " << iNu1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iNu1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iNu1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iNu1].Pt() << " , " << gp.pruned_p4[iNu1].Eta() << " , " << gp.pruned_p4[iNu1].Phi() << " , " << gp.pruned_p4[iNu1].E() << " , " << gp.pruned_p4[iNu1].M() << ")" << std::endl;
    std::cout << "\tiNu2= " << iNu2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iNu2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iNu2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iNu2].Pt() << " , " << gp.pruned_p4[iNu2].Eta() << " , " << gp.pruned_p4[iNu2].Phi() << " , " << gp.pruned_p4[iNu2].E() << " , " << gp.pruned_p4[iNu2].M() << ")" << std::endl;


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
