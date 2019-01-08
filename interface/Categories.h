#ifndef ttH_CATEGORIES_H
#define ttH_CATEGORIES_H

#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/ttHAnalysis/interface/ttHAnalyzer.h>

class DileptonCategory: public Category {
    public:
        const std::vector<ttH::Lepton>& getLeptons(const AnalyzersManager& analyzers) const ;
        const std::vector<ttH::Dilepton>& getDileptons(const AnalyzersManager& analyzers) const ;
        const std::vector<ttH::DileptonMetDijet>& getDileptonMetDijets(const AnalyzersManager& analyzers) const ;
        virtual void configure(const edm::ParameterSet& conf) override {
            m_analyzer_name = conf.getUntrackedParameter<std::string>("m_analyzer_name", "hh_analyzer");
        }
    private:
        std::string m_analyzer_name;

    protected:
        float m_leadingLeptonPtCut;
        float m_subleadingLeptonPtCut;
};

class MuMuCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void configure(const edm::ParameterSet& conf) override;
};

class ElElCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void configure(const edm::ParameterSet& conf) override;
};

class ElMuCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void configure(const edm::ParameterSet& conf) override;
};

class MuElCategory: public DileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void configure(const edm::ParameterSet& conf) override;
};

#endif
