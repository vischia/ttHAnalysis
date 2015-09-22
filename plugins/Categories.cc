#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>
#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

// ***** ***** *****
// Asking a minimum of two jets in the event
// ***** ***** *****
bool HHDiJetCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return jets.p4.size() >= 2;
};
bool HHDiJetCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hhanalyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hhanalyzer.selectedjets_p4.size() >= 2;
};
void HHDiJetCategory::register_cuts(CutManager& manager) {
    // empty
};
void HHDiJetCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    // empty
}

// ***** ***** *****
// Asking a minimum of two b-jets in the event
// ***** ***** *****
bool HHDibJetCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return jets.p4.size() >= 2;
};
bool HHDibJetCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hhanalyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hhanalyzer.selectedbjets_p4.size() >= 2;
};
void HHDibJetCategory::register_cuts(CutManager& manager) {
    // empty
};
void HHDibJetCategory::evaluate_cuts_pre_analyzers(CutManager& manager, const ProducersManager& producers) const {
    // empty
}

// ***** ***** *****
// Dilepton categories
// ***** ***** *****

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
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
};

void HHMuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if ( tempdiLeptons[idilep].isMuMu ) {
            if (  tempdiLeptons[idilep].p4.M() > m_mll_cut  ) manager.pass_cut("ll_mass");
            if (  tempdiLeptons[idilep].p4.M() > m_mll_lowerZcut  ) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool HHElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return electrons.p4.size() >= 2 ;
};

bool HHElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isElEl) return true;
    }
    return false;
};

void HHElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
};

void HHElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if ( tempdiLeptons[idilep].isElEl ) {
            if (  tempdiLeptons[idilep].p4.M() > m_mll_cut  ) manager.pass_cut("ll_mass");
            if (  tempdiLeptons[idilep].p4.M() > m_mll_lowerZcut  ) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool HHElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return electrons.p4.size() + muons.p4.size() >= 2 ;
};

bool HHElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isElMu) return true;
    }
    return false;
};

void HHElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
};

void HHElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if ( tempdiLeptons[idilep].isElMu ) {
            if (  tempdiLeptons[idilep].p4.M() > m_mll_cut  ) manager.pass_cut("ll_mass");
            if (  tempdiLeptons[idilep].p4.M() > m_mll_lowerZcut  ) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool HHMuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return electrons.p4.size() + muons.p4.size() >= 2 ;
};

bool HHMuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isMuEl) return true;
    }
    return false;
};

void HHMuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
};

void HHMuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if ( tempdiLeptons[idilep].isMuEl ) {
            if (  tempdiLeptons[idilep].p4.M() > m_mll_cut  ) manager.pass_cut("ll_mass");
            if (  tempdiLeptons[idilep].p4.M() > m_mll_lowerZcut  ) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
}

