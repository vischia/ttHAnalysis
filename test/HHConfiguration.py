
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer

runOnData = False

if runOnData :
    globalTag = '74X_dataRun2_v2'
else : 
    globalTag = '74X_mcRun2_asymptotic_v2'

process = Framework.create(runOnData, eras.Run2_25ns, globalTag, cms.PSet(

    hh_analyzer = cms.PSet(
        type = cms.string('hh_analyzer'),
        prefix = cms.string('hh_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            mll_cut = cms.untracked.double(10)
            ),
        parameters = cms.PSet(
            # Here are the default value (just to show what is configurable)
            electronIsoCut_EB_Loose = cms.untracked.double(0.0893), # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            electronIsoCut_EE_Loose = cms.untracked.double(0.121),
            electronIsoCut_EB_Tight = cms.untracked.double(0.0354),
            electronIsoCut_EE_Tight = cms.untracked.double(0.0646),
            electronPtCut = cms.untracked.double(20),
            electronEtaCut = cms.untracked.double(2.5),
            muonIsoCut = cms.untracked.double(.20), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonPtCut = cms.untracked.double(20),
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),
            discr_name =  cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            discr_cut =  cms.untracked.double(0.89)
            ),
        )
    ), 
    
    redoJEC=True,

    )

# Add PUPPI MET
process.framework.producers.puppimet = cms.PSet(METProducer.default_configuration.clone())
process.framework.producers.puppimet.prefix = cms.string('puppimet_')
process.framework.producers.puppimet.parameters.met = cms.untracked.InputTag('slimmedMETsPuppi')

# Custom the cuts
process.framework.producers.jets.parameters.cut = cms.untracked.string("pt > 20")
#process.framework.producers.muons.parameters.cut = cms.untracked.string("pt > 20")
#process.framework.producers.electrons.parameters.cut = cms.untracked.string("pt > 20")

Framework.schedule(process, ['hh_analyzer'])

if runOnData : 
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2015B/DoubleMuon/MINIAOD/17Jul2015-v1/30000/D8ED75E7-C12E-E511-8CBF-0025905A608C.root'
        )
else : 
    process.source.fileNames = cms.untracked.vstring(
        'file:///home/fynu/sbrochet/storage/MINIAODSIM/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_reduced.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
