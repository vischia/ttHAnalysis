#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

class MuMuCategory: public Category {
    virtual bool event_in_category(const ProducersManager& producers) const override {
        const MuonsProducer& muons = dynamic_cast<const MuonsProducer&>(producers.get("muons"));
        const ElectronsProducer& electrons = dynamic_cast<const ElectronsProducer&>(producers.get("electrons"));
        if( muons.p4.size() >= 2 )
        {
            if( electrons.p4.size() >= 1 ) // if there is electrons at all, check the muons are the leading leptons
            {
                if( muons.p4[1].Pt() > electrons.p4[0].Pt() )
                    return true;
            } else {
                return true;
            }
        }
        return false;
    };
    virtual void register_cuts(CutManager& manager) override {
        manager.new_cut("muon_1_pt", "pt > 10");
    };
    virtual void evaluate_cuts(CutManager& manager, const ProducersManager& producers) const override {
        const MuonsProducer& muons = dynamic_cast<const MuonsProducer&>(producers.get("muons"));
        if (muons.p4[0].Pt() > 10)
            manager.pass_cut("muon_1_pt");
    }
};

class MuElCategory: public Category {
    virtual bool event_in_category(const ProducersManager& producers) const override {
        const MuonsProducer& muons = dynamic_cast<const MuonsProducer&>(producers.get("muons"));
        const ElectronsProducer& electrons = dynamic_cast<const ElectronsProducer&>(producers.get("electrons"));
        if( (muons.p4.size() == 1) && (electrons.p4.size() >= 1) )
        {
            if( muons.p4[0].Pt() > electrons.p4[0].Pt() )
                return true;
        }
        else if( (muons.p4.size() >= 1) && (electrons.p4.size() >= 1) )
        {
            if( (muons.p4[0].Pt() > electrons.p4[0].Pt()) && (muons.p4[1].Pt() < electrons.p4[0].Pt()) )
              return true;
        }
        return false;
    };
    virtual void register_cuts(CutManager& manager) override {
        manager.new_cut("muon_1_pt", "pt > 10");
    };
    virtual void evaluate_cuts(CutManager& manager, const ProducersManager& producers) const override {
        const MuonsProducer& muons = dynamic_cast<const MuonsProducer&>(producers.get("muons"));
        if (muons.p4[0].Pt() > 10)
            manager.pass_cut("muon_1_pt");
    }
};

class ElMuCategory: public Category {
    virtual bool event_in_category(const ProducersManager& producers) const override {
        const MuonsProducer& muons = dynamic_cast<const MuonsProducer&>(producers.get("muons"));
        const ElectronsProducer& electrons = dynamic_cast<const ElectronsProducer&>(producers.get("electrons"));
        if( (electrons.p4.size() == 1) && (muons.p4.size() >= 1) )
        {
            if( electrons.p4[0].Pt() > muons.p4[0].Pt() )
                return true;
        }
        else if( (electrons.p4.size() >= 1) && (muons.p4.size() >= 1) )
        {
            if( (electrons.p4[0].Pt() > muons.p4[0].Pt()) && (electrons.p4[1].Pt() < muons.p4[0].Pt()) )
              return true;
        }
        return false;
    };
    virtual void register_cuts(CutManager& manager) override {
        manager.new_cut("electron_1_pt", "pt > 10");
    };
    virtual void evaluate_cuts(CutManager& manager, const ProducersManager& producers) const override {
        const ElectronsProducer& electrons = dynamic_cast<const ElectronsProducer&>(producers.get("electrons"));
        if (electrons.p4[0].Pt() > 10)
            manager.pass_cut("electron_1_pt");
    }
};

class ElElCategory: public Category {
    virtual bool event_in_category(const ProducersManager& producers) const override {
        const MuonsProducer& muons = dynamic_cast<const MuonsProducer&>(producers.get("muons"));
        const ElectronsProducer& electrons = dynamic_cast<const ElectronsProducer&>(producers.get("electrons"));
        if( electrons.p4.size() >= 2 )
        {
            if( muons.p4.size() >= 1 ) // if there is muons at all, check the electrons are the leading leptons
            {
                if( electrons.p4[1].Pt() > muons.p4[0].Pt() )
                    return true;
            } else {
                return true;
            }
        }
        return false;
    };
    virtual void register_cuts(CutManager& manager) override {
        manager.new_cut("electron_1_pt", "pt > 10");
    };
    virtual void evaluate_cuts(CutManager& manager, const ProducersManager& producers) const override {
        const ElectronsProducer& electrons = dynamic_cast<const ElectronsProducer&>(producers.get("electrons"));
        if (electrons.p4[0].Pt() > 10)
            manager.pass_cut("electron_1_pt");
    }
};

class DiJetCategory: public Category {
    virtual bool event_in_category(const ProducersManager& producers) const override {
        const JetsProducer& jets = dynamic_cast<const JetsProducer&>(producers.get("jets"));
        if( jets.p4.size() >= 2 )
            return true;
        else
            return false;
    };
    virtual void register_cuts(CutManager& manager) override {
        manager.new_cut("jet_1_pt", "pt > 20");
        manager.new_cut("jet_2_pt", "pt > 20");
    };
    virtual void evaluate_cuts(CutManager& manager, const ProducersManager& producers) const override {
        const JetsProducer& jets = dynamic_cast<const JetsProducer&>(producers.get("jets"));
        if (jets.p4[0].Pt() > 20)
            manager.pass_cut("jet_1_pt");
        if (jets.p4[1].Pt() > 20)
            manager.pass_cut("jet_2_pt");
    }
};


class HHAnalyzer: public Framework::Analyzer {
    public:
        HHAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config) {

        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&) override;

        virtual void registerCategories(CategoryManager& manager) {
            manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons");
            manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons");
            manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, electron");
            manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, muon");
            manager.new_category<DiJetCategory>("dijet", "Category with at least two jets");
        }

    private:

};

#endif
