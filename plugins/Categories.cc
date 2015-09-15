#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>
#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

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
