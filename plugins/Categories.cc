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

const std::vector<lepton>& HHDileptonCategory::getLeptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hhanalyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hhanalyzer.Leptons;
}

const std::vector<dilepton>& HHDileptonCategory::getdiLeptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hhanalyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hhanalyzer.diLeptons;
}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool HHMuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return muons.p4.size() >= 2 ;
};

bool HHMuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isMuMu) return true;
    }
    return false;
};

void HHMuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void HHMuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].p4.M() > m_mll_cut ) manager.pass_cut("ll_mass");
    }
}

//// ***** ***** *****
//// Dilepton Mu-E category
//// ***** ***** *****
//bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
//    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
//    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
//    return (muons.p4.size() >= 1) && (electrons.p4.size() >= 1);
//};
//
//bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//    const DileptonAnalyzer& dilepton_analyzer = analyzers.get<DileptonAnalyzer>("dilepton");
//    return dilepton_analyzer.muel.size() > 0;
//};
//
//void MuElCategory::register_cuts(CutManager& manager) {
//    manager.new_cut("ll_mass", "mll > 20");
//};
//
//void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//    const DileptonAnalyzer& dilepton_analyzer = analyzers.get<DileptonAnalyzer>("dilepton");
//    for(unsigned int idilepton = 0; idilepton < dilepton_analyzer.muel.size(); idilepton++)
//        if( dilepton_analyzer.muel[idilepton].M() > m_mll_cut)
//        {
//            manager.pass_cut("ll_mass");
//            break;
//        }
//}
//
//// ***** ***** *****
//// Dilepton E-Mu category
//// ***** ***** *****
//bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
//    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
//    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
//    return (muons.p4.size() >= 1) && (electrons.p4.size() >= 1);
//};
//
//bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//    const DileptonAnalyzer& dilepton_analyzer = analyzers.get<DileptonAnalyzer>("dilepton");
//    return dilepton_analyzer.elmu.size() > 0;
//};
//
//void ElMuCategory::register_cuts(CutManager& manager) {
//    manager.new_cut("ll_mass", "mll > 20");
//};
//
//void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//    const DileptonAnalyzer& dilepton_analyzer = analyzers.get<DileptonAnalyzer>("dilepton");
//    for(unsigned int idilepton = 0; idilepton < dilepton_analyzer.elmu.size(); idilepton++)
//        if( dilepton_analyzer.elmu[idilepton].M() > m_mll_cut)
//        {
//            manager.pass_cut("ll_mass");
//            break;
//        }
//}
//
//// ***** ***** *****
//// Dilepton El-El category
//// ***** ***** *****
//bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
//    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
//    return electrons.p4.size() >= 2;
//};
//
//bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//    const DileptonAnalyzer& dilepton_analyzer = analyzers.get<DileptonAnalyzer>("dilepton");
//    return dilepton_analyzer.elel.size() > 0;
//};
//
//void ElElCategory::register_cuts(CutManager& manager) {
//    manager.new_cut("ll_mass", "mll > 20");
//};
//
//void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//    const DileptonAnalyzer& dilepton_analyzer = analyzers.get<DileptonAnalyzer>("dilepton");
//    for(unsigned int idilepton = 0; idilepton < dilepton_analyzer.elel.size(); idilepton++)
//        if( dilepton_analyzer.elel[idilepton].M() > m_mll_cut)
//        {
//            manager.pass_cut("ll_mass");
//            break;
//        }
//}
//
//
//
//
//
