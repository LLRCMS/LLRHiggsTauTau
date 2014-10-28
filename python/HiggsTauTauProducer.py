import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")

#flags (move to cfg?)
try:
    IsMC
except NameError:
    IsMC = True
try:
    LEPTON_SETUP
except NameError:
    LEPTON_SETUP = 2012 # define the set of effective areas, rho corrections, etc.
try:
    APPLYMUCORR
except NameError:
    APPLYMUCORR = True

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_70_V2_AN1::All'   # data in 70X, cf https://twiki.cern.ch/twiki/bin/view/CMS/MiniAOD
print process.GlobalTag.globaltag

### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### ----------------------------------------------------------------------
### Source
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/patTuple_1.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

### ----------------------------------------------------------------------
### Trigger bit Requests 
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi 

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

##
## DEFINE CUTS HERE
##

#Muon stuff from HZZ4L, in case revisit it
# Mu e-scale corrections (MuScleFit)
process.calibratedMuons = cms.EDProducer("MuScleFitPATMuonCorrector", 
                         src = cms.InputTag("slimmedMuons"),#patMuonsWithTrigger"), 
                         debug = cms.bool(False), 
                         identifier = cms.string("Summer12_DR53X_smearReReco"), 
                         applySmearing = cms.bool(IsMC), 
                         fakeSmearing = cms.bool(False)
                         )

# Set correct identifier for MuScleFit muon corrections
if LEPTON_SETUP == 2011:
    if IsMC:
        process.calibratedMuons.identifier = cms.string("Fall11_START42")
    else:
        process.calibratedMuons.identifier = cms.string("Data2011_42X")
else:
    if IsMC:
        process.calibratedMuons.identifier = cms.string("Summer12_DR53X_smearReReco")
    else:
        process.calibratedMuons.identifier = cms.string("Data2012_53X_ReReco")
        

### Mu Ghost cleaning
process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("calibratedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
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

process.softMuons = cms.EDProducer("MuFiller",
    src = cms.InputTag("bareSoftMuons"),
    sampleType = cms.int32(SAMPLE_TYPE),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
#    cut = cms.string("userFloat('SIP')<100"),
    cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1."),
    flags = cms.PSet(
        ID = cms.string("userFloat('isPFMuon')" ), # PF ID
        isGood = cms.string(GOODLEPTON)
    )
)


if APPLYMUCORR :
    process.muons =  cms.Sequence(process.calibratedMuons + process.cleanedMu + process.bareSoftMuons + process.softMuons)
else:
    process.cleanedMu.src = src = cms.InputTag("slimmedMuons")#patMuonsWithTrigger")
    process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons + process.softMuons)
    

##
## Electrons
##
#--- Electron regression+calibrarion must be applied after BDT is recomputed
# NOTE patElectronsWithRegression->eleRegressionEnergy;  calibratedElectrons-> calibratedPatElectrons
# Default: NEW ECAL regression + NEW calibration + NEW combination
process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('patElectronsWithTrigger')
process.eleRegressionEnergy.energyRegressionType = 2 ## 1: ECAL regression w/o subclusters 2 (default): ECAL regression w/ subclusters)
#process.eleRegressionEnergy.vertexCollection = cms.InputTag('goodPrimaryVertices')
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
process.calibratedPatElectrons.correctionsType = 2 # 1 = old regression, 2 = new regression, 3 = no regression, 0 = nothing
process.calibratedPatElectrons.combinationType = 3
process.calibratedPatElectrons.lumiRatio = cms.double(1.0)
process.calibratedPatElectrons.isMC    = IsMC
process.calibratedPatElectrons.synchronization = cms.bool(False)

if (LEPTON_SETUP == 2011):
   process.eleRegressionEnergy.rhoCollection = cms.InputTag('kt6PFJetsForIso:rho')
   if (IsMC):
       process.calibratedPatElectrons.inputDataset = "Fall11"
   else :
       process.calibratedPatElectrons.inputDataset = "Jan16ReReco"
else :
   if (IsMC):
       process.calibratedPatElectrons.inputDataset = "Summer12_LegacyPaper"
   else :
       process.calibratedPatElectrons.inputDataset = "22Jan2013ReReco"

    

process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("calibratedPatElectrons"),
   cut = cms.string("pt>7 && abs(eta)<2.5 &&" +
                    "gsfTrack.trackerExpectedHitsInner.numberOfHits<=1"
                    )
   )

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
   sampleType = cms.int32(SAMPLE_TYPE),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
#    cut = cms.string("userFloat('SIP')<100"),
   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"), # BDT MVA ID
        isGood = cms.string(GOODLEPTON)
        )
   )


process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)

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
   cut = cms.string("pt>7 && abs(eta)<2.5 &&")
   )

#Leptons
process.softLeptons = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag("softMuons"), cms.InputTag("cleanSoftElectrons"))
    src = cms.VInputTag(cms.InputTag("appendPhotons:muons"), cms.InputTag("appendPhotons:electrons"),cms.InputTag("taus"))
)


##
## Build ll candidates (here OS)
##

process.barellCand = cms.EDProducer("CandViewShallowCombiner",
                                    decay = cms.string("softLeptons@+ softLeptons@-"),
                                    cut = cms.string("mass > 0"),
                                    checkCharge = cms.bool(True)
)


##
## SV fit
##
process.SVllCand = cms.EDProducer("SVfitter",
                                  src = cms.VInputTag("barellCand"))

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
       process.softLeptons       +
# Build dilepton candidates
       process.SVllCand
       #process.bareEECand        + process.EECand   +  
       #process.bareMMCand        + process.MMCand   + 
       #process.bareMMMMCand      + process.MMMMCand  +
       #process.bareEEEECand      + process.EEEECand  +
       #process.bareEEMMCand      + process.EEMMCand  +
       #process.ZZCand
  )

SkimPaths = cms.vstring('PVfilter') #Do not apply skim 
