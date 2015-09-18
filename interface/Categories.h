#ifndef CATEGORIES_H
#define CATEGORIES_H

#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

class DiJetCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
};

class HHDileptonCategory: public Category {
    public:
        //virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
        //virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
        //virtual void register_cuts(CutManager& manager) override;
        //virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
        virtual void configure(const edm::ParameterSet& conf) override {
            m_mll_cut = conf.getUntrackedParameter<double>("mll_cut", 20);
        }
        const std::vector<lepton>& getLeptons(const AnalyzersManager& analyzers) const ;
        const std::vector<dilepton>& getdiLeptons(const AnalyzersManager& analyzers) const ;

    protected:
        float m_mll_cut;
};

class HHMuMuCategory: public HHDileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    const std::vector<lepton>& getLeptons(const AnalyzersManager& analyzers) const;
};

class HHMuElCategory: public HHDileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
};

class HHElMuCategory: public HHDileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
};


class HHElElCategory: public HHDileptonCategory {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
};

#endif
