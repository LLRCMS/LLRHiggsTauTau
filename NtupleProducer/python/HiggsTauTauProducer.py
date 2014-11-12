import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")

#set this cut in the cfg file
try: IsMC
except NameError:
    IsMC=True
try: LEPTON_SETUP
except NameError:
    LEPTON_SETUP="2012"
try: ELEREGRESSION
except NameError:
    ELEREGRESSION="Paper"
try: ELECORRTYPE
except NameError:
    ELECORRTYPE="Paper"
try:
    GOODLEPTON
except NameError:
    GOODLEPTON = "mass>0"
print GOODLEPTON


### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'PLS170_V6AN1::All'#'GR_70_V2_AN1::All'   # data in 70X, cf https://twiki.cern.ch/twiki/bin/view/CMS/MiniAOD
print process.GlobalTag.globaltag

### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

### ----------------------------------------------------------------------
### Trigger bit Requests 
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi 
# !!!!!!!!!!!!
print "Trigger part"
## LUCA: added trigger path for tet purposes
process.hltFilterDiMu = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiMu.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiMu.HLTPaths = ["HLT_*"] #["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"] # to run on DATA/MC 2012
process.triggerDiMu = cms.Path(process.hltFilterDiMu)
# !!!!!!!!!!!!

#MC stuff

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.drawTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticlesPruned"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False) )


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("genParticlesPruned")
                                   )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       ),
                                                   )


process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  #src = cms.InputTag("offlinePrimaryVertices"),
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

### Mu Ghost cleaning
process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("slimmedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))


process.softMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("cleanedMu"),
    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) &&" +
                     "pt>5 && abs(eta)<2.4")
#    Lowering pT cuts
#    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) &&" +
#                     "pt>3 && p>3.5 && abs(eta)<2.4")
)


# MC matching. As the genParticles are no more available in cmg, we re-match with genParticlesPruned.
process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                   src     = cms.InputTag("softMuons"), # RECO objects to match  
                                   matched = cms.InputTag("genParticlesPruned"),   # mc-truth particle collection
                                   mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
                                   checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
                                   mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
                                   maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
                                   resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(False), # False = just match input in order; True = pick lowest deltaR pair first
                                   )

process.muons =  cms.Sequence(process.cleanedMu + process.softMuons)
    

###
### Electrons
###
##--- Electron regression+calibrarion must be applied after BDT is recomputed
## NOTE patElectronsWithRegression->eleRegressionEnergy;  calibratedElectrons-> calibratedPatElectrons
## Default: NEW ECAL regression + NEW calibration + NEW combination
#process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
#process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('patElectronsWithTrigger')
#process.eleRegressionEnergy.energyRegressionType = 2 ## 1: ECAL regression w/o subclusters 2 (default): ECAL regression w/ subclusters)
##process.eleRegressionEnergy.vertexCollection = cms.InputTag('goodPrimaryVertices')
#process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
#process.calibratedPatElectrons.correctionsType = 2 # 1 = old regression, 2 = new regression, 3 = no regression, 0 = nothing
#process.calibratedPatElectrons.combinationType = 3
#process.calibratedPatElectrons.lumiRatio = cms.double(1.0)
#process.calibratedPatElectrons.isMC    = IsMC
#process.calibratedPatElectrons.synchronization = cms.bool(False)
#
##if (LEPTON_SETUP == 2011):
##   process.eleRegressionEnergy.rhoCollection = cms.InputTag('kt6PFJetsForIso:rho')
##   if (IsMC):
##       process.calibratedPatElectrons.inputDataset = "Fall11"
##   else :
##       process.calibratedPatElectrons.inputDataset = "Jan16ReReco"
##else :
##if (IsMC):
#process.calibratedPatElectrons.inputDataset = "Summer12_LegacyPaper"
##   else :
##process.calibratedPatElectrons.inputDataset = "22Jan2013ReReco"

    

process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("slimmedElectrons"),#"calibratedPatElectrons"),
   cut = cms.string("pt>7 && abs(eta)<2.5 &&" +
                    "gsfTrack.trackerExpectedHitsInner.numberOfHits<=1"
                    )
   )

#process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons) #+ process.softElectrons)
process.electrons = cms.Sequence(process.bareSoftElectrons) #+ process.softElectrons)

# Handle special cases
if ELEREGRESSION == "None" and ELECORRTYPE == "None" :   # No correction at all. Skip correction modules.
    process.bareSoftElectrons.src = cms.InputTag('slimmedElectrons')#patElectronsWithTrigger')#RH
    process.electrons = cms.Sequence(process.bareSoftElectrons + process.softElectrons)

elif ELEREGRESSION == "Moriond" and ELECORRTYPE == "Moriond" : # Moriond corrections: OLD ECAL regression + OLD calibration + OLD combination 
    if (LEPTON_SETUP == 2011):
        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2011Weights_V1.root")
    else :
        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")
    process.eleRegressionEnergy.energyRegressionType = 1
    process.calibratedPatElectrons.correctionsType   = 1
    process.calibratedPatElectrons.combinationType   = 1

elif ELEREGRESSION == "Paper" and ELECORRTYPE == "None" : # NEW ECAL regression + NO calibration + NO combination
    process.eleRegressionEnergy.energyRegressionType = 2
    process.calibratedPatElectrons.correctionsType   = 0
    process.calibratedPatElectrons.combinationType   = 0

elif ELEREGRESSION == "Paper" and ELECORRTYPE == "PaperNoComb" : # NEW ECAL regression + NEW calibration + NO combination
    process.eleRegressionEnergy.energyRegressionType = 2
    process.calibratedPatElectrons.correctionsType   = 1
    process.calibratedPatElectrons.combinationType   = 0
    

process.electronMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                       src         = cms.InputTag("bareSoftElectrons"), # RECO objects to match
                                       matched     = cms.InputTag("genParticlesPruned"), # mc-truth particle collection
                                       mcPdgId     = cms.vint32(11),               # one or more PDG ID (11 = electron); absolute values (see below)
                                       checkCharge = cms.bool(True),               # True = require RECO and MC objects to have the same charge
                                       mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                       maxDeltaR   = cms.double(0.5),              # Minimum deltaR for the match
                                       maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
                                       resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                       resolveByMatchQuality = cms.bool(False),    # False = just match input in order; True = pick lowest deltaR pair first
                                       )


### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("(isGlobalMuon || userFloat('isPFMuon'))"), #
           deltaR              = cms.double(0.05),  
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)

##
## Taus
##
process.taus = cms.EDFilter("PATTauRefSelector",
   src = cms.InputTag("slimmedTaus"),
   cut = cms.string("pt>7 && abs(eta)<2.5")
   )

#Leptons
process.softLeptons = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("slimmedMuons"), cms.InputTag("slimmedElectrons"),cms.InputTag("slimmedTaus"))
#    src = cms.VInputTag(cms.InputTag("appendPhotons:muons"), cms.InputTag("appendPhotons:electrons"),cms.InputTag("taus"))
)

#
#Jets
#
process.jets = cms.EDFilter("PATJetRefSelector",
                           src = cms.InputTag("slimmedJets"),
                           cut = cms.string("pt>10")
)
##
## Build ll candidates (here OS)
##
process.barellCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("softLeptons@+ softLeptons@-"),
                                    cut = cms.string("mass > 0"),
                                    checkCharge = cms.bool(True)
)


##
## SV fit
##
#process.SVllCand = cms.EDProducer("SVfitter",
#                                  src = cms.InputTag("barellCand"))

##
## Paths
##
process.PVfilter = cms.Path(process.goodPrimaryVertices)

# Prepare lepton collections
process.Candidates = cms.Path(
       process.muons             +
       process.electrons         +  #process.cleanSoftElectrons +
       #process.fsrPhotonSequence + process.appendPhotons     +#RH
       process.taus              +
       process.softLeptons       + process.barellCand +
       process.jets              #+
# Build dilepton candidates
       #process.SVllCand
  )

SkimPaths = cms.vstring('PVfilter') #Do not apply skim 

