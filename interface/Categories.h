#ifndef HH_CATEGORIES_H
#define HH_CATEGORIES_H

#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

class DileptonCategory: public Category {
    public:
        virtual void configure(const edm::ParameterSet& conf) override {
            m_mll_cut = conf.getUntrackedParameter<double>("mll_cut", 20);
            m_mll_lowerZcut = conf.getUntrackedParameter<double>("mll_lowerZcut", 85);
        }
        const std::vector<Lepton>& getLeptons(const AnalyzersManager& analyzers) const ;
        const std::vector<Dilepton>& getDileptons(const AnalyzersManager& analyzers) const ;
        const unsigned int getNJets(const AnalyzersManager& analyzers) const ;
        const unsigned int getNBJets(const AnalyzersManager& analyzers) const ;

    protected:
        float m_mll_cut;
        float m_mll_lowerZcut;
};

class MuMuCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class ElElCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class ElMuCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class MuElCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

#endif
