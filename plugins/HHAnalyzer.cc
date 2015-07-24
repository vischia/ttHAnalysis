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

/*
    LorentzVector gen_B1(0.,0.,0.,0.);
    LorentzVector gen_B2(0.,0.,0.,0.);
    LorentzVector gen_Nu1(0.,0.,0.,0.);
    LorentzVector gen_Nu2(0.,0.,0.,0.);
    LorentzVector gen_L1(0.,0.,0.,0.);
    LorentzVector gen_L2(0.,0.,0.,0.);
    LorentzVector gen_LL(0.,0.,0.,0.);
    LorentzVector gen_BB(0.,0.,0.,0.);
    LorentzVector gen_L1FSRNu(0.,0.,0.,0.);
    LorentzVector gen_L2FSRNu(0.,0.,0.,0.);
    LorentzVector gen_L1FSR(0.,0.,0.,0.);
    LorentzVector gen_L2FSR(0.,0.,0.,0.);
    LorentzVector gen_LLFSR(0.,0.,0.,0.);
    LorentzVector gen_NuNu(0.,0.,0.,0.);
    LorentzVector gen_LLNuNu(0.,0.,0.,0.);
    LorentzVector gen_LLFSRNuNu(0.,0.,0.,0.);
    LorentzVector gen_LLFSRNuNuBB(0.,0.,0.,0.);
*/
    int iX; 
    int iH1, iH2;
    int iV1, iV2;
    int iB1, iB2;
    int iL1, iL2;
    int iNu1, iNu2;
    std::vector<int> iG1, iG2;
    iX = -1;
    iH1 = iH2 = -1;
    iV1 = iV2 = -1;
    iB1 = iB2 = -1;
    iL1 = iL2 = -1;
    iNu1 = iNu2 = -1;
    iG1.clear();
    iG2.clear();
    bool isThereFSRforL1 = false;
    bool isThereFSRforL2 = false;

    // Construct the signal gen-level info
    for (unsigned int ip = 0 ; ip < gp.pruned_p4.size() ; ip++) {
        std::bitset<15> flags (gp.pruned_status_flags[ip]);
        if( !flags.test(13) ) continue; // take the last copies
        if( flags.test(8) ) // first look at the hard process
        {
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
            // if the lepton is not in the hard process itself, then look for FSR photons later on
                if( iL1 == -1 )
                {
                    if( flags.test(8) && !flags.test(7) )
                        isThereFSRforL1 = true;
                    iL1 = ip;
                }
                else if( iL2 == -1 )
                {
                    if( flags.test(8) && !flags.test(7) )
                        isThereFSRforL2 = true;
                    iL2 = ip;
                }
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
            else if( abs(gp.pruned_pdg_id[ip]) == 35 || abs(gp.pruned_pdg_id[ip]) == 39 )
            {
                if( iX == -1 )
                    iX = ip;
                else
                    std::cout << "ERROR, more than 1 radion/graviton in the hard process!" << std::endl;
            }
            else if( abs(gp.pruned_pdg_id[ip]) == 23 || abs(gp.pruned_pdg_id[ip]) == 24 )
            {
                if( iV1 == -1 )
                    iV1 = ip;
                else if( iV2 == -1 )
                    iV2 = ip;
                else
                    std::cout << "ERROR, more than 2 vector bosons in the hard process!" << std::endl;
            }
        } // end if coming from hard process
    } // end of loop over pruned gen particles
    // new loop to find FSR photons
    if( isThereFSRforL1 || isThereFSRforL2 )
    {
        for (unsigned int ip = 0 ; ip < gp.pruned_p4.size() ; ip++) {
            if( gp.pruned_pdg_id[ip] != 22 ) continue;
            std::bitset<15> flags (gp.pruned_status_flags[ip]);
            if( !flags.test(13) ) continue; // take the last copies
            for(unsigned int imother = 0; imother < gp.pruned_mothers_index[ip].size() ; imother++)
            {
                if( isThereFSRforL1 )
                {
                    for(unsigned int jmother = 0; jmother < gp.pruned_mothers_index[iL1].size() ; jmother++)
                    {
                        if( gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[iL1][jmother] )
                            iG1.push_back(ip);
                    }
                }
                if( isThereFSRforL2 ) 
                {
                    for(unsigned int jmother = 0; jmother < gp.pruned_mothers_index[iL2].size() ; jmother++)
                    {
                        if( gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[iL2][jmother] )
                            iG2.push_back(ip);
                    }
                }
            }


/*            
                std::cout << "\tip= " << ip << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[ip] << "\tflags= " << flags << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[ip].Pt() << " , " << gp.pruned_p4[ip].Eta() << " , " << gp.pruned_p4[ip].Phi() << " , " << gp.pruned_p4[ip].E() << " , " << gp.pruned_p4[ip].M() << ")" << std::endl;
                for(unsigned int imother = 0; imother < gp.pruned_mothers_index[ip].size() ; imother++)
                    std::cout << "\timother= " << imother << "\tgp.pruned_mothers_index[" << ip << "][" << imother << "]= " << gp.pruned_mothers_index[ip][imother] << std::endl;
            
*/
        } // end of loop over pruned gen particles
    } // end of FSR loop
    // Sanity checks
    if(isThereFSRforL1) assert(iG1.size() > 0);
    if(isThereFSRforL2) assert(iG2.size() > 0);

    if( iX != -1 )
        std::cout << "\tiX= " << iX << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iX] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iX]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iX].Pt() << " , " << gp.pruned_p4[iX].Eta() << " , " << gp.pruned_p4[iX].Phi() << " , " << gp.pruned_p4[iX].E() << " , " << gp.pruned_p4[iX].M() << ")" << std::endl;
    std::cout << "\tiH1= " << iH1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iH1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iH1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iH1].Pt() << " , " << gp.pruned_p4[iH1].Eta() << " , " << gp.pruned_p4[iH1].Phi() << " , " << gp.pruned_p4[iH1].E() << " , " << gp.pruned_p4[iH1].M() << ")" << std::endl;
    std::cout << "\tiH2= " << iH2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iH2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iH2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iH2].Pt() << " , " << gp.pruned_p4[iH2].Eta() << " , " << gp.pruned_p4[iH2].Phi() << " , " << gp.pruned_p4[iH2].E() << " , " << gp.pruned_p4[iH2].M() << ")" << std::endl;
    std::cout << "\tiB1= " << iB1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iB1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iB1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iB1].Pt() << " , " << gp.pruned_p4[iB1].Eta() << " , " << gp.pruned_p4[iB1].Phi() << " , " << gp.pruned_p4[iB1].E() << " , " << gp.pruned_p4[iB1].M() << ")" << std::endl;
    std::cout << "\tiB2= " << iB2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iB2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iB2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iB2].Pt() << " , " << gp.pruned_p4[iB2].Eta() << " , " << gp.pruned_p4[iB2].Phi() << " , " << gp.pruned_p4[iB2].E() << " , " << gp.pruned_p4[iB2].M() << ")" << std::endl;
    if( iV1 != -1 )
        std::cout << "\tiV1= " << iV1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iV1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iV1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iV1].Pt() << " , " << gp.pruned_p4[iV1].Eta() << " , " << gp.pruned_p4[iV1].Phi() << " , " << gp.pruned_p4[iV1].E() << " , " << gp.pruned_p4[iV1].M() << ")" << std::endl;
    if( iV2 != -1 )
        std::cout << "\tiV2= " << iV2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iV2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iV2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iV2].Pt() << " , " << gp.pruned_p4[iV2].Eta() << " , " << gp.pruned_p4[iV2].Phi() << " , " << gp.pruned_p4[iV2].E() << " , " << gp.pruned_p4[iV2].M() << ")" << std::endl;
    std::cout << "\tiL1= " << iL1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iL1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iL1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iL1].Pt() << " , " << gp.pruned_p4[iL1].Eta() << " , " << gp.pruned_p4[iL1].Phi() << " , " << gp.pruned_p4[iL1].E() << " , " << gp.pruned_p4[iL1].M() << ")" << std::endl;
    std::cout << "\tiL2= " << iL2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iL2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iL2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iL2].Pt() << " , " << gp.pruned_p4[iL2].Eta() << " , " << gp.pruned_p4[iL2].Phi() << " , " << gp.pruned_p4[iL2].E() << " , " << gp.pruned_p4[iL2].M() << ")" << std::endl;
    std::cout << "\tiNu1= " << iNu1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iNu1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iNu1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iNu1].Pt() << " , " << gp.pruned_p4[iNu1].Eta() << " , " << gp.pruned_p4[iNu1].Phi() << " , " << gp.pruned_p4[iNu1].E() << " , " << gp.pruned_p4[iNu1].M() << ")" << std::endl;
    std::cout << "\tiNu2= " << iNu2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iNu2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iNu2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iNu2].Pt() << " , " << gp.pruned_p4[iNu2].Eta() << " , " << gp.pruned_p4[iNu2].Phi() << " , " << gp.pruned_p4[iNu2].E() << " , " << gp.pruned_p4[iNu2].M() << ")" << std::endl;

    for(unsigned int iG = 0; iG < iG1.size() ; iG++)
    {
        std::cout << "\tiG1[" << iG << "]= " << iG1[iG] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iG1[iG]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iG1[iG]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iG1[iG]].Pt() << " , " << gp.pruned_p4[iG1[iG]].Eta() << " , " << gp.pruned_p4[iG1[iG]].Phi() << " , " << gp.pruned_p4[iG1[iG]].E() << " , " << gp.pruned_p4[iG1[iG]].M() << ")" << std::endl;
    }
    for(unsigned int iG = 0; iG < iG2.size() ; iG++)
    {
        std::cout << "\tiG2[" << iG << "]= " << iG2[iG] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[iG2[iG]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[iG2[iG]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[iG2[iG]].Pt() << " , " << gp.pruned_p4[iG2[iG]].Eta() << " , " << gp.pruned_p4[iG2[iG]].Phi() << " , " << gp.pruned_p4[iG2[iG]].E() << " , " << gp.pruned_p4[iG2[iG]].M() << ")" << std::endl;
    }




/*
    for (auto p4: jets.p4) {
        std::cout << "Jet pt: " << p4.P() << std::endl;
    }
*/
}
