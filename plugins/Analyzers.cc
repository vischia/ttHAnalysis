#include <FWCore/PluginManager/interface/PluginFactory.h>

#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>

DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, HHAnalyzer, "hh_analyzer");
