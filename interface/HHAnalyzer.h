#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/Category.h>

struct lepton { 
    LorentzVector p4; 
    unsigned int idx; 
    bool isMu; 
    bool isEl; 
};  

struct dilepton { 
    LorentzVector p4; 
    std::pair<unsigned int, unsigned int> idxs; 
    bool isMuMu; 
    bool isElEl; 
    bool isElMu; 
    bool isMuEl; 
};  

class HHAnalyzer: public Framework::Analyzer {
    public:
        HHAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),
            m_electronIsoCut( config.getUntrackedParameter<double>("electronIsoCut") ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut") ),
            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut") ),
            m_muonIsoCut( config.getUntrackedParameter<double>("muonIsoCut") ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut") ),
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut") )
        {
            m_electron_loose_wp_name = config.getUntrackedParameter<std::string>("electrons_loose_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-loose");
            m_electron_tight_wp_name = config.getUntrackedParameter<std::string>("electrons_tight_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-tight");
        }

        std::vector<lepton> Leptons;
        std::vector<dilepton> diLeptons;

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet& config) override;
        
        BRANCH(selectedjets, std::vector<LorentzVector>);
        BRANCH(selectedbjets, std::vector<LorentzVector>);
        BRANCH(dijets, LorentzVector);
        BRANCH(dibjets, LorentzVector);

        BRANCH(Leptons_p4, std::vector<LorentzVector>); // make a list of leptons sorted by pt 
        //std::pair<LorentzVector, LorentzVector>& Lepton_p4 = tree["Lepton_p4"].write<std::pair<LorentzVector, LorentzVector>>();
        BRANCH(Leptons_isMu, std::vector<bool>);
        BRANCH(Leptons_isEl, std::vector<bool>);
        BRANCH(Leptons_idx, std::vector<unsigned int>); // if one need more informations then the p4, still possible to get it back with the idx of the lepton (we also know which collection it belongs to) 
        BRANCH(diLeptons_p4, std::vector<LorentzVector>);
        BRANCH(diLeptons_isMuMu, std::vector<bool>);
        BRANCH(diLeptons_isElEl, std::vector<bool>);
        BRANCH(diLeptons_isElMu, std::vector<bool>);
        BRANCH(diLeptons_isMuEl, std::vector<bool>);
        BRANCH(diLeptons_dRll, std::vector<float>);
        BRANCH(diLeptons_dPhill, std::vector<float>);
        BRANCH(diLeptons_Mll, std::vector<float>);
            
        BRANCH(isolatedElectrons, std::vector<unsigned int>); 
        BRANCH(isolatedMuons, std::vector<unsigned int>);
        BRANCH(tightElectrons, std::vector<unsigned int>);
        BRANCH(tightMuons, std::vector<unsigned int>);
        BRANCH(looseElectrons, std::vector<unsigned int>);
        BRANCH(looseMuons, std::vector<unsigned int>);
        BRANCH(selectedElectrons, std::vector<unsigned int>);
        BRANCH(selectedMuons, std::vector<unsigned int>);
        
    private:
        const float m_electronIsoCut, m_electronEtaCut, m_electronPtCut;
        const float m_muonIsoCut, m_muonEtaCut, m_muonPtCut;

        std::string m_electron_loose_wp_name;
        std::string m_electron_tight_wp_name;

};

#endif
