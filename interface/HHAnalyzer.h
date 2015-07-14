#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>

class TwoMuonsCategory: public Category {
    virtual bool event_in_category() const override {
        return rand() / (float) RAND_MAX > 0.5;
    };

    virtual void register_cuts(CutManager& manager) override {
        manager.new_cut("cut_1", "cut test");
        manager.new_cut("cut_2", "cut test");
    };

    virtual void evaluate_cuts(CutManager& manager) const override {
        if (rand() / (float) RAND_MAX > 0.5)
            manager.pass_cut("cut_1");

        if (rand() / (float) RAND_MAX > 0.5)
            manager.pass_cut("cut_2");
    }
};

class HHAnalyzer: public Framework::Analyzer {
    public:
        HHAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config) {

        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&) override;

        virtual void registerCategories(CategoryManager& manager) {
            manager.new_category("two_muons", "At least two muons category", &two_muons_category);
        }

    private:
        TwoMuonsCategory two_muons_category;

};

#endif
