
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer

process = Framework.create(False, eras.Run2_25ns, '74X_mcRun2_asymptotic_v2', cms.PSet(
#    dilepton = cms.PSet(
#        type = cms.string('dilepton_analyzer'),
#        prefix = cms.string('dilepton_'),
#        enable = cms.bool(True),
#        categories_parameters = cms.PSet(
#            mll_cut = cms.untracked.double(10)
#            ),
#        parameters = cms.PSet(
#            standalone = cms.untracked.bool(True),
#            muons_wp = cms.untracked.string('loose'),
#            electrons_wp = cms.untracked.string('loose')
#            )
#        ),

    bTagsLoose = cms.PSet(
        type = cms.string('btags_analyzer'),
        prefix = cms.string('btags_CSVv2_loose'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            discr_name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            discr_cut = cms.untracked.double(0.605),
            eta_cut = cms.untracked.double(2.4),
            pt_cut = cms.untracked.double(30)
            )
        ),

    bTagsMedium = cms.PSet(
        type = cms.string('btags_analyzer'),
        prefix = cms.string('btags_CSVv2_medium'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            discr_name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            discr_cut = cms.untracked.double(0.89),
            eta_cut = cms.untracked.double(2.4),
            pt_cut = cms.untracked.double(30)
            )
        ),

    bTagsTight = cms.PSet(
        type = cms.string('btags_analyzer'),
        prefix = cms.string('btags_CSVv2_tight'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            discr_name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            discr_cut = cms.untracked.double(0.97),
            eta_cut = cms.untracked.double(2.4),
            pt_cut = cms.untracked.double(30)
            )
        ),

    hh_analyzer = cms.PSet(
        type = cms.string('hh_analyzer'),
        prefix = cms.string('hh_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            mll_cut = cms.untracked.double(10)
            ),
        parameters = cms.PSet(
            electronIsoCut = cms.untracked.double(.11),
            electronPtCut = cms.untracked.double(20),
            electronEtaCut = cms.untracked.double(2.5),
            muonIsoCut = cms.untracked.double(.12),
            muonPtCut = cms.untracked.double(20),
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-tight")
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

Framework.schedule(process, ['bTagsLoose', 'bTagsMedium', 'bTagsTight', 'hh_analyzer'])

process.source.fileNames = cms.untracked.vstring(
        'file:///home/fynu/sbrochet/storage/MINIAODSIM/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_reduced.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
