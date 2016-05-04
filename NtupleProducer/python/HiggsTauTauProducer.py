import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

#set this cut in the cfg file
try: APPLYELECORR
except NameError:
    APPLYELECORR="None"
ELECORRTYPE=APPLYELECORR
try: IsMC
except NameError:
    IsMC=True
try: LEPTON_SETUP
except NameError:
    LEPTON_SETUP=2012
try: APPLYFSR
except NameError:
    APPLYFSR=False
try: BUILDONLYOS
except NameError:
    BUILDONLYOS=False
try: Is25ns
except NameError:
    Is25ns=True

try: USE_NOHFMET
except NameError:
    USE_NOHFMET=False

PFMetName = "slimmedMETs"
if USE_NOHFMET: PFMetName = "slimmedMETsNoHF"

try: APPLYMETCORR
except NameError:
    APPLYMETCORR=TRUE

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")    
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
if IsMC:
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1' #MC 25 ns miniAODv2    
else :
    process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'
print process.GlobalTag.globaltag

nanosec="25"
if not Is25ns: nanosec="50"

METfiltersProcess = "PAT" if IsMC else "RECO" # NB! this is not guaranteed to be true! the following is valid on 2015 Run C + Run D data. Check:
# NB: for MET filters, use PAT or RECO depending if the miniAOD was generated simultaneously with RECO or in a separated step
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters

### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

### ----------------------------------------------------------------------
### Counters 
### ----------------------------------------------------------------------
process.nEventsTotal = cms.EDProducer("EventCountProducer")       # don't change producer name
process.nEventsPassTrigger = cms.EDProducer("EventCountProducer") # these names are then "hard-coded" inside the ntuplizer plugin

### ----------------------------------------------------------------------
### Trigger bit Requests - filter 
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = TRIGGERLIST,
    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False) #if True: throws exception if a trigger path is invalid  
)

#MC stuff

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.drawTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False) )


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       ),
                                                   )


process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string(PVERTEXCUT),
  filter = cms.bool(True),
)

### Mu Ghost cleaning
process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("slimmedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("cleanedMu"),
    cut = cms.string(MUCUT)
#    Lowering pT cuts
#    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) &&" +
#                     "pt>3 && p>3.5 && abs(eta)<2.4")
)


# MC matching. As the genParticles are no more available in cmg, we re-match with prunedGenParticles.
process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                   src     = cms.InputTag("softMuons"), # RECO objects to match  
                                   matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
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
    genCollection = cms.InputTag("prunedGenParticles"),
    rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
    vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    sampleType = cms.int32(LEPTON_SETUP),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
#    cut = cms.string("userFloat('SIP')<100"),
#    cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1."),
    cut = cms.string(""),
    flags = cms.PSet(
        ID = cms.string("userFloat('isPFMuon')" ), # PF ID
        isGood = cms.string(MUCUT)
    )
)

process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons+ process.softMuons)
    

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


# START ELECTRON CUT BASED ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
#process.load("RecoEgamma.ElectronIdentification.ElectronRegressionValueMapProducer_cfi")

#**********************
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
#**********************

process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# Define which IDs we want to produce
# Each of these two example IDs contains all four standard 
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V1_miniAOD_cff']
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']
my_id_modules =[
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',    # both 25 and 50 ns cutbased ids produced
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',                 # recommended for both 50 and 25 ns
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff', # will not be produced for 50 ns, triggering still to come
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',    # 25 ns trig
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_50ns_Trig_V1_cff',    # 50 ns trig
] 
#['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_'+nanosec+'ns_V1_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
#Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    

#process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
#   src = cms.InputTag("slimmedElectrons"),#"calibratedPatElectrons"),
#   cut = cms.string(ELECUT)
#   )


process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("slimmedElectrons"),
   rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
   vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   genCollection = cms.InputTag("prunedGenParticles"),
   sampleType = cms.int32(LEPTON_SETUP),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.

   #CUT BASED ELE ID
   #electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-veto"),
   #electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-tight"),
   #electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-medium"),
   #electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-loose"),
   #electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
   #electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
   #electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
   #electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
   electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-veto"),
   electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-tight"),
   electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-medium"),
   electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-loose"),

   #MVA ELE ID (only for 25ns right now)
   #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80"),
   #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90"),
   #mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
   #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories"),
   eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
   eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
   mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
   mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),

#    cut = cms.string("userFloat('SIP')<100"),
#   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),
   cut = cms.string(ELECUT),
   flags = cms.PSet(
        ID = cms.string("userInt('isBDT')"), # BDT MVA ID
        isGood = cms.string("")
        )
   )

#process.electrons = cms.Sequence(process.egmGsfElectronIDSequence * process.softElectrons)#process.bareSoftElectrons

egmMod = 'egmGsfElectronIDs'
mvaMod = 'electronMVAValueMapProducer'
regMod = 'electronRegressionValueMapProducer'
egmSeq = 'egmGsfElectronIDSequence'
setattr(process,egmMod,process.egmGsfElectronIDs.clone())
setattr(process,mvaMod,process.electronMVAValueMapProducer.clone())
setattr(process,regMod,process.electronRegressionValueMapProducer.clone())
setattr(process,egmSeq,cms.Sequence(getattr(process,mvaMod)*getattr(process,egmMod)*getattr(process,regMod)))
process.electrons = cms.Sequence(getattr(process,mvaMod)*getattr(process,egmMod)*getattr(process,regMod) * process.softElectrons)#process.bareSoftElectrons

# Handle special cases
#if ELECORRTYPE == "None" :   # No correction at all. Skip correction modules.
#    process.bareSoftElectrons.src = cms.InputTag('slimmedElectrons')#patElectronsWithTrigger')#RH
#    process.electrons = cms.Sequence(process.bareSoftElectrons + process.softElectrons)
#
#elif ELECORRTYPE == "Moriond" : # Moriond corrections: OLD ECAL regression + OLD calibration + OLD combination 
#    if (LEPTON_SETUP == 2011):
#        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2011Weights_V1.root")
#    else :
#        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")
#    process.eleRegressionEnergy.energyRegressionType = 1
#    process.calibratedPatElectrons.correctionsType   = 1
#    process.calibratedPatElectrons.combinationType   = 1
#
#elif ELECORRTYPE == "Paper" : # NEW ECAL regression + NO calibration + NO combination
#    process.eleRegressionEnergy.energyRegressionType = 2
#    process.calibratedPatElectrons.correctionsType   = 0
#    process.calibratedPatElectrons.combinationType   = 0
#    

process.electronMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                       src         = cms.InputTag("bareSoftElectrons"), # RECO objects to match
                                       matched     = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
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
process.bareTaus = cms.EDFilter("PATTauRefSelector",
   src = cms.InputTag("slimmedTaus"),
   cut = cms.string(TAUCUT)
   )

##NOT USED FOR NOW, TBD Later
process.cleanTaus = cms.EDProducer("PATTauCleaner",
    src = cms.InputTag("bareTaus"),
    # preselection (any string-based cut on pat::Tau)
    preselection = cms.string(
            'tauID("decayModeFinding") > 0.5 &'
            ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
            ' tauID("againstMuonTight") > 0.5 &'
            ' tauID("againstElectronMedium") > 0.5'
        ),
    
   # overlap checking configurables
   checkOverlaps = cms.PSet(
      muons = cms.PSet(
          src       = cms.InputTag("cleanPatMuons"),
          algorithm = cms.string("byDeltaR"),
          preselection        = cms.string(""),
          deltaR              = cms.double(0.3),
          checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
          pairCut             = cms.string(""),
          requireNoOverlaps   = cms.bool(False), # overlaps don't cause the electron to be discared
          ),
      electrons = cms.PSet(
          src       = cms.InputTag("cleanPatElectrons"),
          algorithm = cms.string("byDeltaR"),
          preselection        = cms.string(""),
          deltaR              = cms.double(0.3),
          checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
          pairCut             = cms.string(""),
          requireNoOverlaps   = cms.bool(False), # overlaps don't cause the electron to be discared
          ),
      ),
        # finalCut (any string-based cut on pat::Tau)
        finalCut = cms.string(' '),
)

NominalTESCorrection=-1#in percent

process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   discriminator = cms.string(TAUDISCRIMINATOR),
   NominalUpOrDown = cms.string("Nominal"),
   NominalTESCorrection = cms.double(NominalTESCorrection),                            
   ApplyTESCorrection = cms.bool(APPLYTESCORRECTION),
   flags = cms.PSet(
        isGood = cms.string("")
        )
   )

process.softTausTauUp = process.softTaus.clone(
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   NominalUpOrDown = cms.string("Up"),
   NominalTESCorrection = cms.double(NominalTESCorrection),    
   ApplyTESCorrection = cms.bool(APPLYTESCORRECTION),                        
   discriminator = cms.string(TAUDISCRIMINATOR),
   flags = cms.PSet(
        isGood = cms.string("")
        )
   )

process.softTausTauDown = process.softTaus.clone(
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   NominalUpOrDown = cms.string("Down"),
   NominalTESCorrection = cms.double(NominalTESCorrection), 
   ApplyTESCorrection = cms.bool(APPLYTESCORRECTION),                       
   discriminator = cms.string(TAUDISCRIMINATOR),
   flags = cms.PSet(
        isGood = cms.string("")
        )
   )

process.tauMatch = cms.EDProducer("MCMatcher",
    src = cms.InputTag("softTaus"),
    maxDPtRel = cms.double(999.9),
    mcPdgId = cms.vint32(15),
    mcStatus = cms.vint32(2),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(999.9),
    checkCharge = cms.bool(True),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("prunedGenParticles")
    )


process.taus=cms.Sequence(process.bareTaus + process.softTaus + process.softTausTauUp + process.softTausTauDown)

process.tausMerged = cms.EDProducer("MergeTauCollections",
    src        = cms.InputTag("softTaus"),
    srcTauUp   = cms.InputTag("softTausTauUp"),
    srcTauDown = cms.InputTag("softTausTauDown")
)

#process.tausDummy=cms.Sequence(process.

#process.tausMerged=cms.Sequence(process.softTaus + process.softTausTauUp + process.softTausTauDown + process.tausMerging)

# ### ----------------------------------------------------------------------
# ### b quarks, only from MC
# ### ----------------------------------------------------------------------
# process.bQuarks = cms.EDProducer("bFiller",
#          src = cms.InputTag("prunedGenParticles"),
#          cut = cms.string(BCUT),
#          flags = cms.PSet(
#             isGood = cms.string("")
#         )
#  )                
# if IsMC : process.bquarks = cms.Sequence(process.bQuarks)
# else : process.bquarks = cms.Sequence()

### ----------------------------------------------------------------------
### gen info, only from MC
### ----------------------------------------------------------------------
process.genInfo = cms.EDProducer("GenFiller",
         src = cms.InputTag("prunedGenParticles"),
 )                
if IsMC : process.geninfo = cms.Sequence(process.genInfo)
else : process.geninfo = cms.Sequence()


### ----------------------------------------------------------------------
### Search for FSR candidates
### ----------------------------------------------------------------------
process.load("UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff")
process.appendPhotons = cms.EDProducer("LeptonPhotonMatcher",
    muonSrc = cms.InputTag("softMuons"),
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),#cms.InputTag("cmgPhotonSel"),
    matchFSR = cms.bool(True)
    )

process.fsrSequence = cms.Sequence(process.fsrPhotonSequence + process.appendPhotons)
muString = "appendPhotons:muons"
eleString = "appendPhotons:electrons"
if not APPLYFSR : 
    process.fsrSequence = cms.Sequence()
    muString = "softMuons"
    eleString = "softElectrons"
    tauString = "tausMerged"
#    tauString = "softTaus"
#Leptons
process.softLeptons = cms.EDProducer("CandViewMerger",
    #src = cms.VInputTag(cms.InputTag("slimmedMuons"), cms.InputTag("slimmedElectrons"),cms.InputTag("slimmedTaus"))
    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauString))
)
#print: "lepton collection built"

#tauStringTauUp = "softTausTauUp" 
#process.softLeptonsTauUp = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauStringTauUp))
#)

#tauStringTauDown = "softTausTauDown" 
#process.softLeptonsTauDown = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauStringTauDown))
#)

#
#Jets
#

# add latest pileup jet ID
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
  jets = cms.InputTag("slimmedJets"),
  inputIsCorrected = True,
  applyJec = True,
  vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
)
#print process.pileupJetIdUpdated.dumpConfig()

# apply new jet energy corrections
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
jecLevels = None
if IsMC:
    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
else:
    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJets"),
    levels = jecLevels,
    payload = 'AK4PFchs' # Make sure to choose the appropriate levels and payload here!
)
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)
process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']

process.jets = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("patJetsReapplyJEC"),
    cut = cms.string(JETCUT)
)
##
## Build ll candidates (here OS)
##
decayString="softLeptons softLeptons"
checkcharge=False
if BUILDONLYOS:
    decayString="softLeptons@+ softLeptons@-"
    checkcharge=True
process.barellCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string(decayString),
                                    cut = cms.string(LLCUT),
                                    checkCharge = cms.bool(checkcharge)
)

#decayStringTauUp="softLeptonsTauUp softLeptonsTauUp"
#if BUILDONLYOS:
#    decayStringTauUp="softLeptonsTauUp@+ softLeptonsTauUp@-"
#    checkcharge=True
#process.barellCandTauUp = cms.EDProducer("CandViewShallowCloneCombiner",
#                                    decay = cms.string(decayStringTauUp),
#                                    cut = cms.string(LLCUT),
#                                    checkCharge = cms.bool(checkcharge)
#)

#decayStringTauDown="softLeptonsTauDown softLeptonsTauDown"
#if BUILDONLYOS:
#    decayStringTauDown="softLeptonsTauDown@+ softLeptonsTauDown@-"
#    checkcharge=True
#process.barellCandTauDown = cms.EDProducer("CandViewShallowCloneCombiner",
#                                    decay = cms.string(decayStringTauDown),
#                                    cut = cms.string(LLCUT),
#                                    checkCharge = cms.bool(checkcharge)
#)

## ----------------------------------------------------------------------
## MVA MET
## ----------------------------------------------------------------------

if USEPAIRMET:
    print "Using pair MET (MVA MET)"
    from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
    runMVAMET(process, jetCollectionPF = "patJetsReapplyJEC")
    process.MVAMET.srcLeptons = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
    process.MVAMET.requireOS = cms.bool(False)

    process.MVAMETInputs = cms.Sequence(
        process.slimmedElectronsTight + process.slimmedMuonsTight + process.slimmedTausLoose + process.slimmedTausLooseCleaned + process.patJetsReapplyJECCleaned +
        process.pfCHS + process.pfChargedPV + process.pfChargedPU + process.pfNeutrals + process.neutralInJets +
        process.pfMETCands + process.pfTrackMETCands + process.pfNoPUMETCands + process.pfPUCorrectedMETCands + process.pfPUMETCands +
        process.pfChargedPUMETCands + process.pfNeutralPUMETCands + process.pfNeutralPVMETCands + process.pfNeutralUnclusteredMETCands +
        process.pfChs +
        process.ak4PFCHSL1FastjetCorrector + process.ak4PFCHSL2RelativeCorrector + process.ak4PFCHSL3AbsoluteCorrector + process.ak4PFCHSResidualCorrector +
        process.ak4PFCHSL1FastL2L3Corrector + process.ak4PFCHSL1FastL2L3ResidualCorrector +
        process.tauDecayProducts + process.tauPFMET + process.tauMET + process.tausSignificance
    )
    for met in ["pfMET", "pfTrackMET", "pfNoPUMET", "pfPUCorrectedMET", "pfPUMET", "pfChargedPUMET", "pfNeutralPUMET", "pfNeutralPVMET", "pfNeutralUnclusteredMET"]:
        process.MVAMETInputs += getattr(process, met)
        process.MVAMETInputs += getattr(process, "ak4JetsFor"+met)
        process.MVAMETInputs += getattr(process, "corr"+met)
        process.MVAMETInputs += getattr(process, met+"T1")
        process.MVAMETInputs += getattr(process, "pat"+met)
        process.MVAMETInputs += getattr(process, "pat"+met+"T1")        

    process.METSequence = cms.Sequence(process.MVAMETInputs)

    # python trick: loop on all pairs for pair MET computation
    UnpackerTemplate = cms.EDProducer("PairUnpacker",
        src = cms.InputTag("barellCand")
    )

    MVAPairMET = []
    for index in range(100):
        UnpackerName = "PairUnpacker%i" % index
        UnpackerModule = UnpackerTemplate.clone( pairIndex = cms.int32(index) )
        setattr(process, UnpackerName, UnpackerModule)   #equiv to process.<UnpackerName> = <UnpackerModule>
        process.METSequence += UnpackerModule

        MVAMETName = "patMETMVA%i" % index
        MVAModule = process.MVAMET.clone( srcLeptons = cms.VInputTag (cms.InputTag(UnpackerName) ) )
        setattr(process, MVAMETName, MVAModule)
        process.METSequence += MVAModule
   
        MVAPairMET.append(cms.InputTag(MVAMETName, "MVAMET"))

else:
    print "Using event pfMET (same MET for all pairs)"
    process.load("RecoMET.METProducers.METSignificance_cfi")
    process.load("RecoMET.METProducers.METSignificanceParams_cfi")
    process.METSequence = cms.Sequence(process.METSignificance)

## ----------------------------------------------------------------------
## Z-recoil correction
## ----------------------------------------------------------------------

corrMVAPairMET = []
if IsMC and APPLYMETCORR:
    if USEPAIRMET:
        process.selJetsForZrecoilCorrection = cms.EDFilter("PATJetSelector",
            src = cms.InputTag("jets"),                                      
            cut = cms.string("pt > 30. & abs(eta) < 4.7"), 
            filter = cms.bool(False)
        )
        process.METSequence += process.selJetsForZrecoilCorrection        
        for index in range(100):
            corrMVAMETName = "corrMETMVA%i" % index
            corrMVAModule = cms.EDProducer("ZrecoilCorrectionProducer",                                                   
                srcLeptons = cms.InputTag("PairUnpacker%i" % index),
                srcMEt = MVAPairMET[index],
                srcGenParticles = cms.InputTag("prunedGenParticles"),
                srcJets = cms.InputTag("selJetsForZrecoilCorrection"),
                correction = cms.string("HTT-utilities/RecoilCorrections/data/recoilMvaMEt_76X_newTraining_MG5.root")
            )
            setattr(process, corrMVAMETName, corrMVAModule)
            process.METSequence += corrMVAModule
            corrMVAPairMET.append(cms.InputTag(corrMVAMETName))

    else:
        raise ValueError("Z-recoil corrections for PFMET not implemented yet !!")

## ----------------------------------------------------------------------
## SV fit
## ----------------------------------------------------------------------
process.SVllCand = cms.EDProducer("SVfitInterface",
                                  srcPairs   = cms.InputTag("barellCand"),
                                  usePairMET = cms.bool(USEPAIRMET),
                                  srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                  srcCov     = cms.InputTag("METSignificance", "METCovariance")

)


srcMETTag = None
if USEPAIRMET:
  srcMETTag = cms.VInputTag(corrMVAPairMET) if (IsMC and APPLYMETCORR) else cms.VInputTag(MVAPairMET)
else:
  srcMETTag = cms.VInputTag(PFMetName)

process.SVllCand.srcMET = srcMETTag

# if USEPAIRMET:
#     if IsMC and APPLYMETCORR:
#         process.SVllCand.srcMET    = cms.VInputTag(corrMVAPairMET)
#     else:
#         process.SVllCand.srcMET    = cms.VInputTag(MVAPairMET)
# #   process.SVllCandTauUp.srcMET    = cms.VInputTag(MVAPairMET)
# #   process.SVllCandTauDown.srcMET    = cms.VInputTag(MVAPairMET)
# else:
#    process.SVllCand.srcMET    = cms.VInputTag(PFMetName)

## ----------------------------------------------------------------------
## SV fit BYPASS (skip SVfit, don't compute SVfit pair mass and don't get MET userfloats
## ----------------------------------------------------------------------
process.SVbypass = cms.EDProducer ("SVfitBypass",
                                    srcPairs   = cms.InputTag("barellCand"),
                                    usePairMET = cms.bool(USEPAIRMET),
                                    srcMET     = srcMETTag,
                                    srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                    srcCov     = cms.InputTag("METSignificance", "METCovariance")
)


## ----------------------------------------------------------------------
## Ntuplizer
## ----------------------------------------------------------------------
process.HTauTauTree = cms.EDAnalyzer("HTauTauNtuplizer",
                      JetCollection = cms.InputTag("jets"),    
                      fileName = cms.untracked.string ("CosaACaso"),
                      skipEmptyEvents = cms.bool(True),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      puCollection = cms.InputTag("slimmedAddPileupInfo"),
                      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoMiniRelIsoCollection = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                      PFCandCollection = cms.InputTag("packedPFCandidates"),
                      jetCollection = cms.InputTag("jets"),
                      lepCollection = cms.InputTag("softLeptons"),
                      lheCollection = cms.InputTag("LHEEventProduct"),
                      genCollection = cms.InputTag("generator"),
                      genericCollection = cms.InputTag("genInfo"),
                      genjetCollection = cms.InputTag("slimmedGenJets"),
                      totCollection = cms.InputTag("nEventsTotal"),
                      passCollection = cms.InputTag("nEventsPassTrigger"),
                      lhepCollection = cms.InputTag("externalLHEProducer"),
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
                      triggerSet = cms.InputTag("selectedPatTrigger"),
                      triggerList = HLTLIST,
                      metFilters = cms.InputTag ("TriggerResults","",METfiltersProcess),
                      l1extraIsoTau = cms.InputTag("l1extraParticles", "IsoTau"),
                      HT = cms.InputTag("externalLHEProducer"),
                      useNOHFMet = cms.bool(USE_NOHFMET) # true to run on silver json, false to use ful MET including HF
                      )
if USE_NOHFMET:
    process.HTauTauTree.metCollection = cms.InputTag("slimmedMETsNoHF")
else: 
    process.HTauTauTree.metCollection = cms.InputTag("slimmedMETs")

if SVFITBYPASS:
    process.HTauTauTree.CandCollection = cms.untracked.string("SVbypass")
    process.HTauTauTree.candCollection = cms.InputTag("SVbypass")
    process.SVFit = cms.Sequence (process.SVbypass)


else:
    process.HTauTauTree.CandCollection = cms.untracked.string("SVllCand")
    process.HTauTauTree.candCollection = cms.InputTag("SVllCand")
    process.SVFit = cms.Sequence (process.SVllCand)

#print particles gen level - DEBUG purposes
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(10),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("prunedGenParticles")
)

##
## Paths
##
process.PVfilter = cms.Path(process.goodPrimaryVertices)

#process.HTtruth = cms.Path(process.HT)

# Prepare lepton collections
process.Candidates = cms.Sequence(
    #process.printTree         + # just for debug, print MC particles
    process.nEventsTotal      +
    #process.hltFilter         + 
    process.nEventsPassTrigger+
    process.muons             +
    process.electrons         + process.cleanSoftElectrons +
    process.taus              + process.tausMerged +
    process.fsrSequence       +
    process.softLeptons       + process.barellCand +
#    process.softLeptonsTauUp  + process.barellCandTauUp +
#    process.softLeptonsTauDown  + process.barellCandTauDown +
    process.pileupJetIdUpdated + process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.jets +
    process.METSequence       +
    process.geninfo           +
    process.SVFit             
    )
# always run ntuplizer
process.trees = cms.EndPath(process.HTauTauTree)# + process.HTauTauTreeTauUp + process.HTauTauTreeTauDown)

