#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>
#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return muons.p4.size() >= 2 ;
};

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.dimuons.size() > 0;
};

void MuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void MuMuCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    if( (muons.p4[0] + muons.p4[1]).M() > 20.)
        manager.pass_cut("ll_mass");
}

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    if( hh_analyzer.dimuons[0].M() > 20. )
        manager.pass_cut("ll_mass");
}

// ***** ***** *****
// Dilepton Mu-E category
// ***** ***** *****
bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return (muons.p4.size() >= 1) && (electrons.p4.size() >= 1);
};

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.dileptons_mue.size() > 0;
};

void MuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void MuElCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    if( (muons.p4[0] + electrons.p4[0]).M() > 20.)
        manager.pass_cut("ll_mass");
}

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    if( hh_analyzer.dileptons_mue[0].M() > 20. )
        manager.pass_cut("ll_mass");
}

// ***** ***** *****
// Dilepton E-Mu category
// ***** ***** *****
bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return (muons.p4.size() >= 1) && (electrons.p4.size() >= 1);
};

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.dileptons_emu.size() > 0;
};

void ElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void ElMuCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    if( (muons.p4[0] + electrons.p4[0]).M() > 20.)
        manager.pass_cut("ll_mass");
}

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    if( hh_analyzer.dileptons_emu[0].M() > 20. )
        manager.pass_cut("ll_mass");
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return electrons.p4.size() >= 2;
};

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.dielectrons.size() > 0;
};

void ElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void ElElCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    if( (electrons.p4[0] + electrons.p4[1]).M() > 20.)
        manager.pass_cut("ll_mass");
}

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    if( hh_analyzer.dielectrons[0].M() > 20. )
        manager.pass_cut("ll_mass");
}

// ***** ***** *****
// Asking a minimum of two jets in the event
// ***** ***** *****
bool DiJetCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return jets.p4.size() >= 2;
};
bool DiJetCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    return true;
};
void DiJetCategory::register_cuts(CutManager& manager) {
    // empty
};
void DiJetCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    // empty
}
