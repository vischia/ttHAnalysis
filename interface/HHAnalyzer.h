#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/HHAnalysis/interface/Types.h>

#include <Math/VectorUtil.h>

using namespace HH;

class HHAnalyzer: public Framework::Analyzer {
    public:
        HHAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config)
        {
            m_muonIsoCut = config.getUntrackedParameter<double>("muonIsoCut", 0.12);
            m_muonEtaCut = config.getUntrackedParameter<double>("muonEtaCut", 2.4);
            m_muonPtCut = config.getUntrackedParameter<double>("muonPtCut", 20);

            m_electronIsoCut = config.getUntrackedParameter<double>("electronIsoCut", 0.11);
            m_electronEtaCut = config.getUntrackedParameter<double>("electronEtaCut", 2.5);
            m_electronPtCut = config.getUntrackedParameter<double>("electronPtCut", 20);
            m_electron_loose_wp_name = config.getUntrackedParameter<std::string>("electrons_loose_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-loose");
            m_electron_tight_wp_name = config.getUntrackedParameter<std::string>("electrons_tight_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-tight");

            m_jetEtaCut = config.getUntrackedParameter<double>("jetEtaCut", 2.4);
            m_jetPtCut = config.getUntrackedParameter<double>("jetPtCut", 20);
            m_jet_bDiscrName = config.getUntrackedParameter<std::string>("discr_name", "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            m_jet_bDiscrCut = config.getUntrackedParameter<double>("discr_cut", 0.89);
        }

        BRANCH(leptons, std::vector<Lepton>);
        BRANCH(ll, std::vector<Dilepton>);

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet& config) override;

        // jets and dijets stuff
        BRANCH(jets_p4, std::vector<LorentzVector>);
        BRANCH(jets_idx, std::vector<unsigned int>);

        BRANCH(jj_p4, std::vector<LorentzVector>);
        BRANCH(jj_idx, std::vector<std::pair<unsigned int, unsigned int>>);// NB : this index refers, so far, to the entries in jets_p4
        BRANCH(jj_DR, std::vector<float>);
        BRANCH(jj_DPhi, std::vector<float>);
        BRANCH(jj_DPhi_met, std::vector<float>);
        BRANCH(jj_minDPhi_jmet, std::vector<float>);
        BRANCH(jj_maxDPhi_jmet, std::vector<float>);

        BRANCH(h_dijet_idx, unsigned int);

        BRANCH(bjets_p4, std::vector<LorentzVector>);
        BRANCH(bjets_idx, std::vector<unsigned int>);

        BRANCH(bb_p4, std::vector<LorentzVector>);
        BRANCH(bb_idx, std::vector<std::pair<unsigned int, unsigned int>>); // NB : this index refers, so far, to the entries in bjets_p4
        BRANCH(bb_DR, std::vector<float>);
        BRANCH(bb_DPhi, std::vector<float>);
        BRANCH(bb_DPhi_met, std::vector<float>);
        BRANCH(bb_minDPhi_jmet, std::vector<float>);
        BRANCH(bb_maxDPhi_jmet, std::vector<float>);

        BRANCH(h_dibjet_idx, unsigned int);

        // leptons and dileptons stuff
        BRANCH(electrons, std::vector<unsigned int>);
        BRANCH(muons, std::vector<unsigned int>);

        BRANCH(leptons_p4, std::vector<LorentzVector>); // list of leptons p4 sorted by pt
        BRANCH(leptons_isMu, std::vector<bool>);
        BRANCH(leptons_isEl, std::vector<bool>);
        BRANCH(leptons_idx, std::vector<unsigned int>);

        BRANCH(ll_p4, std::vector<LorentzVector>);
        BRANCH(llmet_p4, std::vector<LorentzVector>);
        BRANCH(ll_idx, std::vector<std::pair<unsigned int, unsigned int>>);  // refers to leptons indices
        BRANCH(ll_isMuMu, std::vector<bool>);
        BRANCH(ll_isElEl, std::vector<bool>);
        BRANCH(ll_isElMu, std::vector<bool>);
        BRANCH(ll_isMuEl, std::vector<bool>);
        BRANCH(ll_DR, std::vector<float>);
        BRANCH(ll_DPhi, std::vector<float>);
        BRANCH(ll_DPhi_met, std::vector<float>);
        BRANCH(ll_minDPhi_lmet, std::vector<float>);
        BRANCH(ll_maxDPhi_lmet, std::vector<float>);
        BRANCH(ll_MT, std::vector<float>);
        BRANCH(ll_MT_formula, std::vector<float>);
        BRANCH(ll_projectedMet, std::vector<float>);

        // lljj and llbb stuff
        BRANCH(lljj_p4, std::vector<LorentzVector>);
        BRANCH(lljj_idx, std::vector<std::pair<unsigned int, unsigned int>>);  // refers to ll and jj indices
        BRANCH(lljj_DR, std::vector<float>);
        BRANCH(lljj_DPhi, std::vector<float>);
        BRANCH(lljj_minDR_lj, std::vector<float>);
        BRANCH(lljj_maxDR_lj, std::vector<float>);

        BRANCH(llbb_p4, std::vector<LorentzVector>);
        BRANCH(llbb_idx, std::vector<std::pair<unsigned int, unsigned int>>);  // refers to ll and bb indices
        BRANCH(llbb_DR, std::vector<float>);
        BRANCH(llbb_DPhi, std::vector<float>);
        BRANCH(llbb_minDR_lb, std::vector<float>);
        BRANCH(llbb_maxDR_lb, std::vector<float>);

        // lljjmet and llbbmet stuff
        // as there is only one met, all the following vectors are in sync with lljj vectors
        // i.e. no need to store ll and jj indices
        BRANCH(lljjmet_p4, std::vector<LorentzVector>);
        BRANCH(lljjmet_DR, std::vector<float>);
        BRANCH(lljjmet_DPhi, std::vector<float>);
        BRANCH(lljjmet_cosThetaStar_CS, std::vector<float>);

        BRANCH(llbbmet_p4, std::vector<LorentzVector>);
        BRANCH(llbbmet_DR, std::vector<float>);
        BRANCH(llbbmet_DPhi, std::vector<float>);
        BRANCH(llbbmet_cosThetaStar_CS, std::vector<float>);

        // global event stuff (selected objects multiplicity)
        BRANCH(nJets, unsigned int);
        BRANCH(nBJets, unsigned int);
        BRANCH(nMuons, unsigned int);
        BRANCH(nElectrons, unsigned int);
        BRANCH(nLeptons, unsigned int);

        float m_electronIsoCut, m_electronEtaCut, m_electronPtCut;
        float m_muonIsoCut, m_muonEtaCut, m_muonPtCut;
        float m_jetEtaCut, m_jetPtCut, m_jet_bDiscrCut;
        std::string m_jet_bDiscrName;
        std::string m_electron_loose_wp_name;
        std::string m_electron_tight_wp_name;

        // utilities
        float getCosThetaStar_CS(const LorentzVector & h1, const LorentzVector & h2, float ebeam = 6500)
        {// cos theta star angle in the Collins Soper frame
            LorentzVector p1, p2;
            p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
            p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

            LorentzVector hh = h1 + h2;
            ROOT::Math::Boost boost(-hh.X() / hh.T(), -hh.Y() / hh.T(), -hh.Z() / hh.T());
            p1 = boost(p1);
            p2 = boost(p2);
            LorentzVector newh1 = boost(h1);
            ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> CSaxis(p1.Vect().Unit() - p2.Vect().Unit());

            return cos(ROOT::Math::VectorUtil::Angle(CSaxis.Unit(), newh1.Vect().Unit()));
        }
};

#endif
