#pragma once

// Code from https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h

#include <bitset>

struct GenStatusFlags {

    public:

        enum StatusBits {
            kIsPrompt = 0,
            kIsDecayedLeptonHadron,
            kIsTauDecayProduct,
            kIsPromptTauDecayProduct,
            kIsDirectTauDecayProduct,
            kIsDirectPromptTauDecayProduct,
            kIsDirectHadronDecayProduct,
            kIsHardProcess,
            kFromHardProcess,
            kIsHardProcessTauDecayProduct,
            kIsDirectHardProcessTauDecayProduct,
            kFromHardProcessBeforeFSR,
            kIsFirstCopy,
            kIsLastCopy,
            kIsLastCopyBeforeFSR
        };

        GenStatusFlags(int16_t flags):
            flags_(flags) {
            // Empty
        }

        /////////////////////////////////////////////////////////////////////////////
        //these are robust, generator-independent functions for categorizing
        //mainly final state particles, but also intermediate hadrons/taus    

        //is particle prompt (not from hadron, muon, or tau decay)
        bool isPrompt() const { return flags_[kIsPrompt]; }

        //is particle a decayed hadron, muon, or tau (does not include resonance decays like W,Z,Higgs,top,etc)
        //This flag is equivalent to status 2 in the current HepMC standard
        //but older generators (pythia6, herwig6) predate this and use status 2 also for other intermediate
        //particles/states    
        bool isDecayedLeptonHadron() const { return flags_[kIsDecayedLeptonHadron]; }

        //this particle is a direct or indirect tau decay product
        bool isTauDecayProduct() const { return flags_[kIsTauDecayProduct]; }

        //this particle is a direct or indirect decay product of a prompt tau
        bool isPromptTauDecayProduct() const { return flags_[kIsPromptTauDecayProduct]; }

        //this particle is a direct tau decay product
        bool isDirectTauDecayProduct() const { return flags_[kIsDirectTauDecayProduct]; }

        //this particle is a direct decay product from a prompt tau 
        bool isDirectPromptTauDecayProduct() const { return flags_[kIsDirectPromptTauDecayProduct]; }

        //this particle is a direct decay product from a hadron
        bool isDirectHadronDecayProduct() const { return flags_[kIsDirectHadronDecayProduct]; }

        /////////////////////////////////////////////////////////////////////////////
        //these are generator history-dependent functions for tagging particles
        //associated with the hard process
        //Currently implemented for Pythia 6 and Pythia 8 status codes and history   
        //and may not have 100% consistent meaning across all types of processes
        //Users are strongly encouraged to stick to the more robust flags above    

        //this particle is part of the hard process
        bool isHardProcess() const { return flags_[kIsHardProcess]; }

        //this particle is the direct descendant of a hard process particle of the same pdg id
        bool fromHardProcess() const { return flags_[kFromHardProcess]; }

        //this particle is a direct or indirect decay product of a tau
        //from the hard process
        bool isHardProcessTauDecayProduct() const { return flags_[kIsHardProcessTauDecayProduct]; }

        //this particle is a direct decay product of a tau
        //from the hard process
        bool isDirectHardProcessTauDecayProduct() const { return flags_[kIsDirectHardProcessTauDecayProduct]; }

        //this particle is the direct descendant of a hard process particle of the same pdg id
        //For outgoing particles the kinematics are those before QCD or QED FSR
        //This corresponds roughly to status code 3 in pythia 6    
        bool fromHardProcessBeforeFSR() const { return flags_[kFromHardProcessBeforeFSR]; }

        //this particle is the first copy of the particle in the chain with the same pdg id 
        bool isFirstCopy() const { return flags_[kIsFirstCopy]; }

        //this particle is the last copy of the particle in the chain with the same pdg id
        //(and therefore is more likely, but not guaranteed, to carry the final physical momentum)    
        bool isLastCopy() const { return flags_[kIsLastCopy]; }

        //this particle is the last copy of the particle in the chain with the same pdg id
        //before QED or QCD FSR
        //(and therefore is more likely, but not guaranteed, to carry the momentum after ISR)  
        bool isLastCopyBeforeFSR() const { return flags_[kIsLastCopyBeforeFSR]; }

        void dump() const {
            std::cout << "Generator status flags:" << std::endl;
            std::cout << "\tisPrompt: " << isPrompt() << std::endl;
            std::cout << "\tisDecayedLeptonHadron: " << isDecayedLeptonHadron() << std::endl;
            std::cout << "\tisTauDecayProduct: " << isTauDecayProduct() << std::endl;
            std::cout << "\tisPromptTauDecayProduct: " << isPromptTauDecayProduct() << std::endl;
            std::cout << "\tisDirectTauDecayProduct: " << isDirectTauDecayProduct() << std::endl;
            std::cout << "\tisDirectPromptTauDecayProduct: " << isDirectPromptTauDecayProduct() << std::endl;
            std::cout << "\tisDirectHadronDecayProduct: " << isDirectHadronDecayProduct() << std::endl;
            std::cout << "\tisHardProcess: " << isHardProcess() << std::endl;
            std::cout << "\tfromHardProcess: " << fromHardProcess() << std::endl;
            std::cout << "\tisHardProcessTauDecayProduct: " << isHardProcessTauDecayProduct() << std::endl;
            std::cout << "\tisDirectHardProcessTauDecayProduct: " << isDirectHardProcessTauDecayProduct() << std::endl;
            std::cout << "\tfromHardProcessBeforeFSR: " << fromHardProcessBeforeFSR() << std::endl;
            std::cout << "\tisFirstCopy: " << isFirstCopy() << std::endl;
            std::cout << "\tisLastCopy: " << isLastCopy() << std::endl;
            std::cout << "\tisLastCopyBeforeFSR: " << isLastCopyBeforeFSR() << std::endl;
        }


    private:
        std::bitset<15> flags_;
};
