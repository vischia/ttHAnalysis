
import FWCore.ParameterSet.Config as cms

process = cms.Process("HHAna")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#process.GlobalTag.globaltag = "MCRUN2_74_V9" # for Spring15 MC with 25ns asymptotic scenario
process.GlobalTag.globaltag = "74X_dataRun2_Prompt_v1" # for Run2015B/C (50/25ns) data taking

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))
process.source = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring(
#        '/store/mc/RunIISpring15DR74/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/50C1DCD3-7013-E511-8A98-008CFA010718.root'
#        '/store/mc/RunIISpring15DR74/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/7846F40B-D912-E511-832C-0CC47A4D99D6.root'
        'file:///nfs/user/obondu/MINIAODSIM/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph.root'
#        'file:///home/fynu/sbrochet/storage/MINIAODSIM/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9.root'
#        'file:///home/fynu/sbrochet/storage/MINIAODSIM/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A.root'
        )
##### #####
# Electron ID recipe
##### #####
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V1_cff','RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

##### #####
# Services
##### #####
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Producers
from cp3_llbb.Framework import EventProducer
#from cp3_llbb.Framework import GenParticlesProducer
from cp3_llbb.Framework import HLTProducer
from cp3_llbb.Framework import JetsProducer
from cp3_llbb.Framework import METProducer
from cp3_llbb.Framework import MuonsProducer
from cp3_llbb.Framework import ElectronsProducer
from cp3_llbb.Framework import VerticesProducer

process.framework = cms.EDProducer("ExTreeMaker",
        output = cms.string('output_mc.root'),

        producers = cms.PSet(

            event = EventProducer.default_configuration,

            hlt = HLTProducer.default_configuration,

#            gen_particles = GenParticlesProducer.default_configuration,

            jets = JetsProducer.default_configuration,

            met = METProducer.default_configuration,

            muons = MuonsProducer.default_configuration,

            electrons = ElectronsProducer.default_configuration,

            vertices = VerticesProducer.default_configuration,
            ),

        analyzers = cms.PSet(
            hh_analyzer = cms.PSet(
                type = cms.string('hh_analyzer'),
                prefix = cms.string('hh_'),
                enable = cms.bool(True)
                )
            )
        )

process.framework.producers.jets.parameters.cut = cms.untracked.string("pt > 30")
process.framework.producers.muons.parameters.cut = cms.untracked.string("pt > 20")
process.framework.producers.electrons.parameters.cut = cms.untracked.string("pt > 25")
process.framework.producers.met.parameters.cut = cms.untracked.string("pt > 20")

process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.framework
        )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.options.allowUnscheduled = cms.untracked.bool(True)
