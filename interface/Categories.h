#ifndef CATEGORIES_H
#define CATEGORIES_H

#include <cp3_llbb/Framework/interface/Category.h>

class MuMuCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class MuElCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class ElMuCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class ElElCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class DiJetCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const override;
};

#endif
