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
            Analyzer(name, tree_, config)
        {
            m_muonIsoCut = config.getUntrackedParameter<double>("muonIsoCut", 0.12 );
            m_muonEtaCut = config.getUntrackedParameter<double>("muonEtaCut", 2.4 );
            m_muonPtCut = config.getUntrackedParameter<double>("muonPtCut", 20 );

            m_electronIsoCut = config.getUntrackedParameter<double>("electronIsoCut", 0.11 );
            m_electronEtaCut = config.getUntrackedParameter<double>("electronEtaCut", 2.5 );
            m_electronPtCut = config.getUntrackedParameter<double>("electronPtCut", 20 );
            m_electron_loose_wp_name = config.getUntrackedParameter<std::string>("electrons_loose_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-loose");
            m_electron_tight_wp_name = config.getUntrackedParameter<std::string>("electrons_tight_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-tight");

            m_jetEtaCut = config.getUntrackedParameter<double>("jetEtaCut", 2.4);
            m_jetPtCut = config.getUntrackedParameter<double>("jetPtCut", 20);
            m_jet_bDiscrName = config.getUntrackedParameter<std::string>("discr_name", "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            m_jet_bDiscrCut = config.getUntrackedParameter<double>("discr_cut", 0.89);
        }

        std::vector<lepton> Leptons;
        std::vector<dilepton> diLeptons;

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet& config) override;
        
        BRANCH(selectedjets_p4, std::vector<LorentzVector>);
        BRANCH(selectedjets_idx, std::vector<unsigned int>);
        
        BRANCH(dijets_p4, std::vector<LorentzVector>);
        BRANCH(dijets_idx, std::vector<std::pair<unsigned int, unsigned int>>);// NB : this index refers, so far, to the entries in selectedjets_p4
        BRANCH(dijets_Ptjj, std::vector<float>);
        BRANCH(dijets_DRjj, std::vector<float>);
        BRANCH(dijets_DPhijj, std::vector<float>);
        BRANCH(dijets_Mjj, std::vector<float>);

        BRANCH(h_dijet_idx, unsigned int);

        BRANCH(selectedbjets_p4, std::vector<LorentzVector>);
        BRANCH(selectedbjets_idx, std::vector<unsigned int>);

        BRANCH(dibjets_p4, std::vector<LorentzVector>);
        BRANCH(dibjets_idx, std::vector<std::pair<unsigned int, unsigned int>>); // NB : this index refers, so far, to the entries in selectedbjets_p4
        BRANCH(dibjets_Ptbb, std::vector<float>);
        BRANCH(dibjets_DRbb, std::vector<float>);
        BRANCH(dibjets_DPhibb, std::vector<float>);
        BRANCH(dibjets_Mbb, std::vector<float>);

        BRANCH(h_dibjet_idx, unsigned int);

        BRANCH(selectedElectrons, std::vector<unsigned int>);
        BRANCH(selectedMuons, std::vector<unsigned int>);

        BRANCH(Leptons_p4, std::vector<LorentzVector>); // list of leptons p4 sorted by pt 
        BRANCH(Leptons_isMu, std::vector<bool>);
        BRANCH(Leptons_isEl, std::vector<bool>);
        BRANCH(Leptons_idx, std::vector<unsigned int>);  

        BRANCH(diLeptons_p4, std::vector<LorentzVector>);
        BRANCH(diLeptons_idx, std::vector<std::pair<unsigned int, unsigned int>>);  // refers to Lepton indices
        BRANCH(diLeptons_isMuMu, std::vector<bool>);
        BRANCH(diLeptons_isElEl, std::vector<bool>);
        BRANCH(diLeptons_isElMu, std::vector<bool>);
        BRANCH(diLeptons_isMuEl, std::vector<bool>);
        BRANCH(diLeptons_Ptll, std::vector<float>);
        BRANCH(diLeptons_DRll, std::vector<float>);
        BRANCH(diLeptons_DPhill, std::vector<float>);
        BRANCH(diLeptons_Mll, std::vector<float>);

        BRANCH(DR_ll_jj, float);
        BRANCH(DPhi_ll_jj, float);
        BRANCH(Pt_lljj, float);
        BRANCH(M_lljj, float);

        BRANCH(DR_ll_bb, float);
        BRANCH(DPhi_ll_bb, float);
        BRANCH(Pt_llbb, float);
        BRANCH(M_llbb, float);

        BRANCH(minDR_jl, float);
        BRANCH(minDR_bl, float);

        BRANCH(DPhi_ll_met, float);
        BRANCH(minDPhi_l_met, float);
        BRANCH(MT, float);
        BRANCH(MT_formula, float);
        BRANCH(projectedMet, float);

        BRANCH(nJet, unsigned int);
        BRANCH(nbJet, unsigned int);
        BRANCH(nMu, unsigned int);
        BRANCH(nEl, unsigned int);
        BRANCH(nLep, unsigned int);

        float m_electronIsoCut, m_electronEtaCut, m_electronPtCut;
        float m_muonIsoCut, m_muonEtaCut, m_muonPtCut;
        float m_jetEtaCut, m_jetPtCut, m_jet_bDiscrCut;
        std::string m_jet_bDiscrName;
        std::string m_electron_loose_wp_name;
        std::string m_electron_tight_wp_name;

};

#endif
