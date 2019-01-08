import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer
from cp3_llbb.Framework.CmdLine import CmdLine

options = CmdLine()
runOnData = options.runOnData == 1

print('Now starting')
# Not sure it's the correct way to add a new era
options.changeDefaults(era='2016')
print('Changed era')

# 2017 tags from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis
if runOnData :
    options.changeDefaults(globalTag='94X_dataRun2_ReReco_EOY17_v6', process='RECO')
else:
    options.changeDefaults(globalTag='94X_mc2017_realistic_v14', process='PAT')

print('Changed defaults')

framework = Framework.Framework(options)

print('Created framework object')
framework.addAnalyzer('tth_analyzer', cms.PSet(
        type = cms.string('tth_analyzer'),
        prefix = cms.string('tth_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            # Per-category lepton pt cuts
            mumu_leadingLeptonPtCut = cms.untracked.double(20), # muon
            mumu_subleadingLeptonPtCut = cms.untracked.double(10), # muon
            elel_leadingLeptonPtCut = cms.untracked.double(25), # electron
            elel_subleadingLeptonPtCut = cms.untracked.double(15), # electron
            muel_leadingLeptonPtCut = cms.untracked.double(25), # muon
            muel_subleadingLeptonPtCut = cms.untracked.double(15), # electron
            elmu_leadingLeptonPtCut = cms.untracked.double(25), # electron
            elmu_subleadingLeptonPtCut = cms.untracked.double(10), # muon
        ),
        parameters = cms.PSet(
            # Producers
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),
            nohfMETProducer = cms.string('nohf_met'),

            # Pre-selection pt cut, applied to all leptons
            leadingElectronPtCut = cms.untracked.double(25),
            subleadingElectronPtCut = cms.untracked.double(15),
            leadingMuonPtCut = cms.untracked.double(20),
            subleadingMuonPtCut = cms.untracked.double(10),

            electronEtaCut = cms.untracked.double(2.5),
            muonLooseIsoCut = cms.untracked.double(.25), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonTightIsoCut = cms.untracked.double(.15), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-loose"),
            electrons_medium_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-medium"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-tight"),
            electrons_hlt_safe_wp_name = cms.untracked.string("cutBasedElectronHLTPreselection-Summer16-V1"),
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),

            # BTAG INFO
            # Working points from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
            discr_name =  cms.untracked.string("pfCombinedMVAV2BJetTags"),
            discr_cut_loose =  cms.untracked.double(-0.5884),
            discr_cut_medium =  cms.untracked.double(0.4432),
            discr_cut_tight =  cms.untracked.double(0.9432),

            minDR_l_j_Cut = cms.untracked.double(0.3),
            hltDRCut = cms.untracked.double(0.1),
            hltDPtCut = cms.untracked.double(0.5),  # cut will be DPt/Pt < hltDPtCut
            applyBJetRegression = cms.untracked.bool(False), # BE SURE TO ACTIVATE computeRegression FLAG BELOW

            hlt_efficiencies = cms.untracked.PSet(

                    IsoMu17leg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Muon_DoubleIsoMu17Mu8_IsoMu17leg.json'),
                    IsoMu8orIsoTkMu8leg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg.json'),

                    DoubleEleHighPtleg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Electron_IsoEle23Leg.json'),
                    DoubleEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Electron_IsoEle12Leg.json'),

                    EleMuHighPtleg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Electron_IsoEle23Leg.json'),
                    MuEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Electron_IsoEle12Leg.json'),

                    IsoMu8leg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Muon_XPathIsoMu8leg.json'),
                    IsoMu23leg = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/Efficiencies/Muon_XPathIsoMu23leg.json'),
            )
        )
    )
)

# Remove fat jets
framework.removeProducer('fat_jets')

framework.getProducer('hlt').parameters.triggers = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/triggers.xml')
# framework.getProducer('jets').parameters.cut = cms.untracked.string("pt > 20")
#framework.getProducer('jets').parameters.computeRegression = cms.untracked.bool(True)

framework.getProducer('electrons').parameters.scale_factors.id_mediumplushltsafe_hh = cms.untracked.FileInPath('cp3_llbb/ttHAnalysis/data/ScaleFactors/Electron_MediumPlusHLTSafeID_moriond17.json')

if runOnData:
    framework.redoJEC()

# Must activate Rochester correction (input file missing)
# framework.applyMuonCorrection('rochester', input='')

# Missing regressionWeights_cfi
# framework.applyElectronRegression()
# framework.applyElectronSmearing()

if not runOnData:
    framework.smearJets(resolutionFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt', scaleFactorFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_SF_AK4PFchs.txt')
    framework.doSystematics(['jec', 'jer'], jec={'uncertaintiesFile': 'cp3_llbb/ttHAnalysis/data/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt', 'splitBySources': True})

process = framework.create()

if runOnData: 
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2017C/DoubleMuon/MINIAOD/31Mar2018-v1/00000/1E30A24E-8E39-E811-8AD3-8CDCD4A99E08.root'
        )
else: 
    process.framework.treeFlushSize = cms.untracked.uint64(5 * 1024 * 1024)

    process.source.fileNames = cms.untracked.vstring(
        # Signal
        '/store/mc/RunIIFall17MiniAODv2/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/00000/10AE80AC-C94C-E811-9E46-FA163E6C5A08.root'
        
        # TT
        # '/store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_widthx1p15_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/110000/42928F17-F6B2-E811-B466-FA163E6E1BE8.root'

        # DY
        # '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/110000/F2DF6219-3FBA-E811-BEF5-02163E019FB4.root'
        )

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.source.skipEvents = cms.untracked.uint32(10)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
