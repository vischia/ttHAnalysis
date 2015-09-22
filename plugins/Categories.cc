#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>
#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

// NB : All di-lepton categories require also two jets, as we do not want to store events without at least two leptons and two jets.
//      The criteria "having two b-jets" is passed as a cut in the category in question.

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

const unsigned int HHDileptonCategory::getNumberOfJet(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hhanalyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hhanalyzer.selectedjets_p4.size();
}

const unsigned int HHDileptonCategory::getNumberOfBJet(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hhanalyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hhanalyzer.selectedbjets_p4.size();
}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool HHMuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return (muons.p4.size() >= 2 && jets.p4.size() >= 2);
};

bool HHMuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    bool isMuMu = false;
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isMuMu) isMuMu = true;
    }
    const unsigned int numberOfJet = getNumberOfJet(analyzers);
    return (isMuMu && numberOfJet >= 2);
};

void HHMuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
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
    const unsigned int numberOfBJet = getNumberOfBJet(analyzers);
    if ( numberOfBJet >=2 ) manager.pass_cut("has_two_bJets"); 
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool HHElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return (electrons.p4.size() >= 2  && jets.p4.size() >= 2);
};

bool HHElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    bool isElEl = false;
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isElEl) isElEl = true;
    }
    const unsigned int numberOfJet = getNumberOfJet(analyzers);
    return (isElEl && numberOfJet >= 2);
};

void HHElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
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
    const unsigned int numberOfBJet = getNumberOfBJet(analyzers);
    if ( numberOfBJet >=2 ) manager.pass_cut("has_two_bJets"); 
}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool HHElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return ((electrons.p4.size() + muons.p4.size() >= 2)  && jets.p4.size() >= 2);
};

bool HHElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    bool isElMu = false;
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isElMu) isElMu = true;
    }
    const unsigned int numberOfJet = getNumberOfJet(analyzers);
    return (isElMu && numberOfJet >= 2);
};

void HHElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
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
    const unsigned int numberOfBJet = getNumberOfBJet(analyzers);
    if ( numberOfBJet >=2 ) manager.pass_cut("has_two_bJets"); 
}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool HHMuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    return ((electrons.p4.size() + muons.p4.size() >= 2)  && jets.p4.size() >= 2);
};

bool HHMuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<dilepton>& tempdiLeptons = getdiLeptons(analyzers);
    bool isMuEl = false;
    for (unsigned int idilep = 0; idilep < tempdiLeptons.size(); idilep++) 
    {
        if (tempdiLeptons[idilep].isMuEl) isMuEl = true;
    }
    const unsigned int numberOfJet = getNumberOfJet(analyzers);
    return (isMuEl && numberOfJet >= 2);
};

void HHMuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
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
    const unsigned int numberOfBJet = getNumberOfBJet(analyzers);
    if ( numberOfBJet >=2 ) manager.pass_cut("has_two_bJets"); 
}

