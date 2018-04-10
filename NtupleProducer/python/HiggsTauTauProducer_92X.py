import FWCore.ParameterSet.Config as cms
execfile(PyFilePath+"python/triggers_92X.py") # contains the list of triggers and filters
#execfile(PyFilePath+"python/triggers_92X_test.py") # contains the list of triggers and filters

process = cms.Process("TEST")

#set this cut in the cfg file
try: APPLYELECORR
except NameError:
    APPLYELECORR="None"
ELECORRTYPE=APPLYELECORR
try: IsMC
except NameError:
    IsMC=True
try: doCPVariables
except NameError:
    doCPVariables=True       
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
    APPLYMETCORR=True

try: HLTProcessName
except NameError:
    HLTProcessName='HLT'

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")    
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
if IsMC:
    #process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6' # FIXME !!!!!!!!
    #process.GlobalTag.globaltag = '92X_dataRun2_2017Repro_v4'               # 12Sept2017 MC
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'                 # 2017 MC
else :
    # process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0' # ICHEP            # FIXME !!!!!!!!
    #process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' # Run B-G               # FIXME !!!!!!!!
    # process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v16' # Run H                      # FIXME !!!!!!!!
    process.GlobalTag.globaltag = '94X_dataRun2_v6'                    # 2017 Data
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
#process.load("RecoMET.METFilters.python.badGlobalMuonTaggersMiniAOD_cff")
#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

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
    TriggerResultsTag = cms.InputTag("TriggerResults","",HLTProcessName),
    HLTPaths = TRIGGERLIST,
    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False) #if True: throws exception if a trigger path is invalid
)

# Trigger Unpacker Module
#process.patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
#                                            patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
#                                            triggerResults = cms.InputTag("TriggerResults","",HLTProcessName),
#                                            unpackFilterLabels = cms.bool(True)
#)

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
  filter = cms.bool(False), # if True, rejects events . if False, produce emtpy vtx collection
)

#Re-Reco2016 fix from G. Petrucciani
#process.load("RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff")
process.badGlobalMuonTagger = cms.EDFilter("BadGlobalMuonTagger",
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonPtCut = cms.double(20),
    selectClones = cms.bool(False),
    taggingMode = cms.bool(True),
    verbose     = cms.untracked.bool(False),
)
process.cloneGlobalMuonTagger = process.badGlobalMuonTagger.clone(
    selectClones = cms.bool(True)
)

process.removeBadAndCloneGlobalMuons = cms.EDProducer("MuonRefPruner",
    input = cms.InputTag("slimmedMuons"),
    toremove = cms.InputTag("badGlobalMuonTagger", "bad"),
    toremove2 = cms.InputTag("cloneGlobalMuonTagger", "bad")
)

# process.removeCloneGlobalMuons = cms.EDProducer("MuonRefPruner",
#     input = cms.InputTag("removeBadGlobalMuons"),
#     toremove = cms.InputTag("cloneGlobalMuonTagger")
# )

# process.noBadGlobalMuons = cms.Sequence(~process.cloneGlobalMuonTagger + ~process.badGlobalMuonTagger)
process.noBadGlobalMuons = cms.Sequence(process.cloneGlobalMuonTagger + process.badGlobalMuonTagger + process.removeBadAndCloneGlobalMuons) # in tagging mode, these modules return always "true"

#process.softLeptons = cms.EDProducer("CandViewMerger",
#    #src = cms.VInputTag(cms.InputTag("slimmedMuons"), cms.InputTag("slimmedElectrons"),cms.InputTag("slimmedTaus"))
#    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauString))
#)

### Mu Ghost cleaning
# process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
#                                    src = cms.InputTag("cloneGlobalMuonTagger"),
#                                    preselection = cms.string("track.isNonnull"),
#                                    passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
#                                    fractionOfSharedSegments = cms.double(0.499))


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("removeBadAndCloneGlobalMuons"),
    cut = cms.string(MUCUT),
#    Lowering pT cuts
#    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) &&" +
#                     "pt>3 && p>3.5 && abs(eta)<2.4")
)



# # MC matching. As the genParticles are no more available in cmg, we re-match with prunedGenParticles.
# process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
#                                    src     = cms.InputTag("softMuons"), # RECO objects to match  
#                                    matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
#                                    mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
#                                    checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
#                                    mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
#                                    maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
#                                    maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
#                                    resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
#                                    resolveByMatchQuality = cms.bool(False), # False = just match input in order; True = pick lowest deltaR pair first
#                                    )

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

# process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons+ process.softMuons)
process.muons =  cms.Sequence(process.noBadGlobalMuons + process.bareSoftMuons+ process.softMuons)    

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
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',   #Spring16
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',   #Spring16 HZZ
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',     #Fall17 iso
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',   #Fall17 noIso
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
   #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
   #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
   #mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
   #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
   #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
   #eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
   #mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
   #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
   #HZZmvaValuesMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
   eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
   eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90"),
   eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80"),
   mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
   mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
   HZZmvaValuesMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
#  cut = cms.string("userFloat('SIP')<100"),
#  cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),
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

# process.electronMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
#                                        src         = cms.InputTag("bareSoftElectrons"), # RECO objects to match
#                                        matched     = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
#                                        mcPdgId     = cms.vint32(11),               # one or more PDG ID (11 = electron); absolute values (see below)
#                                        checkCharge = cms.bool(True),               # True = require RECO and MC objects to have the same charge
#                                        mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
#                                        maxDeltaR   = cms.double(0.5),              # Minimum deltaR for the match
#                                        maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
#                                        resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
#                                        resolveByMatchQuality = cms.bool(False),    # False = just match input in order; True = pick lowest deltaR pair first
#                                        )


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

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *

tauIdDiscrMVA_trainings_run2_2017 = {
  'tauIdMVAIsoDBoldDMwLT2017' : "tauIdMVAIsoDBoldDMwLT2017",
  }
tauIdDiscrMVA_WPs_run2_2017 = {
  'tauIdMVAIsoDBoldDMwLT2017' : {
    'Eff95' : "DBoldDMwLTEff95",
    'Eff90' : "DBoldDMwLTEff90",
    'Eff80' : "DBoldDMwLTEff80",
    'Eff70' : "DBoldDMwLTEff70",
    'Eff60' : "DBoldDMwLTEff60",
    'Eff50' : "DBoldDMwLTEff50",
    'Eff40' : "DBoldDMwLTEff40"
    }
  }
tauIdDiscrMVA_2017_version = "v1"

def loadMVA_WPs_run2_2017(process):
                
      for training, gbrForestName in tauIdDiscrMVA_trainings_run2_2017.items():

         process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
            cms.PSet(
               record = cms.string('GBRWrapperRcd'),
               tag = cms.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version)),
               label = cms.untracked.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version))
            )
         )

         for WP in tauIdDiscrMVA_WPs_run2_2017[training].keys():
            process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
               cms.PSet(
                  record = cms.string('PhysicsTGraphPayloadRcd'),
                  tag = cms.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP)),
                  label = cms.untracked.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP))
               )
            )

         process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
            cms.PSet(
               record = cms.string('PhysicsTFormulaPayloadRcd'),
               tag = cms.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version)),
               label = cms.untracked.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version))
            )
)

loadMVA_WPs_run2_2017(process)

process.rerunDiscriminationByIsolationMVArun2v1raw = patDiscriminationByIsolationMVArun2v1raw.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   loadMVAfromDB = cms.bool(True),
   mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1"),  # name of the training you want to use: training with 2017 MC_v1 for oldDM
   mvaOpt = cms.string("DBoldDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)
   requireDecayMode = cms.bool(True),
   verbosity = cms.int32(0)
)

process.rerunDiscriminationByIsolationMVArun2v1VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),    
   Prediscriminants = noPrediscriminants,
   toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
   key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw:category'),
   loadMVAfromDB = cms.bool(True),
   mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_mvaOutput_normalization"), # normalization fo the training you want to use
   mapping = cms.VPSet(
      cms.PSet(
         category = cms.uint32(0),
         cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff90"), # this is the name of the working point you want to use
         variable = cms.string("pt"),
      )
   )
)

# here we produce all the other working points for the training
process.rerunDiscriminationByIsolationMVArun2v1VVLoose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff95")
process.rerunDiscriminationByIsolationMVArun2v1Loose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff80")
process.rerunDiscriminationByIsolationMVArun2v1Medium = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff70")
process.rerunDiscriminationByIsolationMVArun2v1Tight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff60")
process.rerunDiscriminationByIsolationMVArun2v1VTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff50")
process.rerunDiscriminationByIsolationMVArun2v1VVTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff40")

# this sequence has to be included in your cms.Path() before your analyzer which accesses the new variables is called.
process.rerunMvaIsolation2SeqRun2 = cms.Sequence(
   process.rerunDiscriminationByIsolationMVArun2v1raw
   *process.rerunDiscriminationByIsolationMVArun2v1VLoose
   *process.rerunDiscriminationByIsolationMVArun2v1VVLoose
   *process.rerunDiscriminationByIsolationMVArun2v1Loose
   *process.rerunDiscriminationByIsolationMVArun2v1Medium
   *process.rerunDiscriminationByIsolationMVArun2v1Tight
   *process.rerunDiscriminationByIsolationMVArun2v1VTight
   *process.rerunDiscriminationByIsolationMVArun2v1VVTight
)

# embed new id's into new tau collection
embedID = cms.EDProducer("PATTauIDEmbedder",
   src = cms.InputTag('slimmedTaus'),
   tauIDSources = cms.PSet(
      byIsolationMVArun2v1DBoldDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
      byVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VLoose'),
      byVVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVLoose'),
      byLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Loose'),
      byMediumIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Medium'),
      byTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Tight'),
      byVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VTight'),
      byVVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVTight'),
      ),
   )

setattr(process, "NewTauIDsEmbedded", embedID)

# old sequence starts here
process.bareTaus = cms.EDFilter("PATTauRefSelector",
   #src = cms.InputTag("slimmedTaus"),
   src = cms.InputTag("NewTauIDsEmbedded"),
   cut = cms.string(TAUCUT),
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

# NominalTESCorrection=-1#in percent\
APPLYTESCORRECTION = APPLYTESCORRECTION if IsMC else False # always false if data
process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   discriminator = cms.string(TAUDISCRIMINATOR),
   NominalTESCorrection = cms.double(-1), #in percent , shift of central value of TES
   ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
   # ApplyTESUpDown = cms.bool(True if IsMC else False), # no shift computation when data
   flags = cms.PSet(
        isGood = cms.string("")
        )
   )

# process.softTausTauUp = process.softTaus.clone(
#    src = cms.InputTag("bareTaus"),
#    genCollection = cms.InputTag("prunedGenParticles"),
#    vtxCollection = cms.InputTag("goodPrimaryVertices"),
#    cut = cms.string(TAUCUT),
#    NominalUpOrDown = cms.string("Up"),
#    NominalTESCorrection = cms.double(NominalTESCorrection),    
#    ApplyTESCorrection = cms.bool(APPLYTESCORRECTION),                        
#    discriminator = cms.string(TAUDISCRIMINATOR),
#    flags = cms.PSet(
#         isGood = cms.string("")
#         )
#    )

# process.softTausTauDown = process.softTaus.clone(
#    src = cms.InputTag("bareTaus"),
#    genCollection = cms.InputTag("prunedGenParticles"),
#    vtxCollection = cms.InputTag("goodPrimaryVertices"),
#    cut = cms.string(TAUCUT),
#    NominalUpOrDown = cms.string("Down"),
#    NominalTESCorrection = cms.double(NominalTESCorrection), 
#    ApplyTESCorrection = cms.bool(APPLYTESCORRECTION),                       
#    discriminator = cms.string(TAUDISCRIMINATOR),
#    flags = cms.PSet(
#         isGood = cms.string("")
#         )
#    )

# process.tauMatch = cms.EDProducer("MCMatcher",
#     src = cms.InputTag("softTaus"),
#     maxDPtRel = cms.double(999.9),
#     mcPdgId = cms.vint32(15),
#     mcStatus = cms.vint32(2),
#     resolveByMatchQuality = cms.bool(False),
#     maxDeltaR = cms.double(999.9),
#     checkCharge = cms.bool(True),
#     resolveAmbiguities = cms.bool(True),
#     matched = cms.InputTag("prunedGenParticles")
#     )


process.taus=cms.Sequence(process.rerunMvaIsolation2SeqRun2 + process.NewTauIDsEmbedded + process.bareTaus + process.softTaus)

# process.tausMerged = cms.EDProducer("MergeTauCollections",
#     src        = cms.InputTag("softTaus"),
#     srcTauUp   = cms.InputTag("softTausTauUp"),
#     srcTauDown = cms.InputTag("softTausTauDown")
# )

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
         storeLightFlavAndGlu = cms.bool(True) # if True, store also udcs and gluons (first copy)
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
    tauString = "softTaus"
#Leptons
process.softLeptons = cms.EDProducer("CandViewMerger",
    #src = cms.VInputTag(cms.InputTag("slimmedMuons"), cms.InputTag("slimmedElectrons"),cms.InputTag("slimmedTaus"))
    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauString))
)


#
#Jets
#

# # add latest pileup jet ID
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
   jets = cms.InputTag("slimmedJets"),
   inputIsCorrected = True,
   applyJec = True,
   vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
)
#print process.pileupJetIdUpdated.dumpConfig()

# apply new jet energy corrections
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

jecLevels = None
if IsMC:
    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
else:
    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None'),
)

process.updatedPatJetsUpdatedJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
process.updatedPatJetsUpdatedJEC.userData.userInts.src    += ['pileupJetIdUpdated:fullId']
process.jecSequence = cms.Sequence(process.pileupJetIdUpdated + process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)
#process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

process.jets = cms.EDFilter("PATJetRefSelector",
                            #src = cms.InputTag("slimmedJets"),
                            src = cms.InputTag("updatedPatJetsUpdatedJEC"),
                            cut = cms.string(JETCUT),
)


##
## QG tagging for jets
##
if COMPUTEQGVAR:
    qgDatabaseVersion = 'v1' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

    from CondCore.DBCommon.CondDBSetup_cfi import *
    QGPoolDBESSource = cms.ESSource("PoolDBESSource",
                                    CondDBSetup,
                                    toGet = cms.VPSet(),
                                    connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
    )
    
    for type in ['AK4PFchs']:
        QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
                    record = cms.string('QGLikelihoodRcd'),
                    tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
                    label  = cms.untracked.string('QGL_'+type)
                    )))


    process.load('RecoJets.JetProducers.QGTagger_cfi')
    process.QGTagger.srcJets          = cms.InputTag("jets")    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
    process.jetSequence = cms.Sequence(process.jets * process.QGTagger)

else:
    process.jetSequence = cms.Sequence(process.jets)




# il primo legge la collezione dei leptoni e stampa quali sono
#process.beforeLLcombiner = cms.EDFilter("beforeCombiner",
#    src = cms.InputTag("softLeptons")
#)

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

#il seconod legge le pairs e stampa quali sono e da chi sono composti
#process.afterLLcombiner = cms.EDFilter("afterCombiner",
#    srcPairs = cms.InputTag("barellCand")
#)

## ----------------------------------------------------------------------
## MVA MET
## ----------------------------------------------------------------------

process.METSequence = cms.Sequence()
if USEPAIRMET:
    print "Using pair MET (MVA MET)"
    from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
    runMVAMET(process, jetCollectionPF = "patJetsReapplyJEC")
    process.MVAMET.srcLeptons = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
    process.MVAMET.requireOS = cms.bool(False)
    process.MVAMET.permuteLeptonsWithinPlugin = cms.bool(False)
    process.MVAMET.leptonPermutations = cms.InputTag("barellCand")

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

    process.METSequence += cms.Sequence(process.MVAMETInputs + process.MVAMET)

    # # python trick: loop on all pairs for pair MET computation
    # UnpackerTemplate = cms.EDProducer("PairUnpacker",
    #     src = cms.InputTag("barellCand")
    # )

    # MVAPairMET = []
    # for index in range(210):
    #     UnpackerName = "PairUnpacker%i" % index
    #     UnpackerModule = UnpackerTemplate.clone( pairIndex = cms.int32(index) )
    #     setattr(process, UnpackerName, UnpackerModule)   #equiv to process.<UnpackerName> = <UnpackerModule>
    #     process.METSequence += UnpackerModule

    #     MVAMETName = "patMETMVA%i" % index
    #     MVAModule = process.MVAMET.clone( srcLeptons = cms.VInputTag (cms.InputTag(UnpackerName) ) )
    #     setattr(process, MVAMETName, MVAModule)
    #     process.METSequence += MVAModule
   
    #     MVAPairMET.append(cms.InputTag(MVAMETName, "MVAMET"))

else:
    print "Using event pfMET (same MET for all pairs)"

    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(
      process,
      isData= (not IsMC),
    )
    # patch to get a standalone MET significance collection
    process.METSignificance = cms.EDProducer ("ExtractMETSignificance",
                                                  srcMET=cms.InputTag("slimmedMETs","","TEST")
                                                  )
    process.METSequence += process.fullPatMetSequence
    process.METSequence += process.METSignificance

# ## always compute met significance
# process.load("RecoMET.METProducers.METSignificance_cfi")
# process.load("RecoMET.METProducers.METSignificanceParams_cfi")
# process.METSequence += cms.Sequence(process.METSignificance)

## ----------------------------------------------------------------------
## Z-recoil correction
## ----------------------------------------------------------------------

# corrMVAPairMET = []
if IsMC and APPLYMETCORR:
    if USEPAIRMET:
        process.selJetsForZrecoilCorrection = cms.EDFilter("PATJetSelector",
            src = cms.InputTag("jets"),                                      
            cut = cms.string("pt > 30. & abs(eta) < 4.7"), 
            filter = cms.bool(False)
        )
        process.corrMVAMET = cms.EDProducer("ZrecoilCorrectionProducer",                                                   
            srcPairs = cms.InputTag("barellCand"),
            srcMEt = cms.InputTag("MVAMET", "MVAMET"),
            srcGenParticles = cms.InputTag("prunedGenParticles"),
            srcJets = cms.InputTag("selJetsForZrecoilCorrection"),
            correction = cms.string("HTT-utilities/RecoilCorrections/data/MvaMET_MG_2016BCD.root")
        )
        process.METSequence += process.selJetsForZrecoilCorrection        
        process.METSequence += process.corrMVAMET

    else:
        raise ValueError("Z-recoil corrections for PFMET not implemented yet !!")


srcMETTag = None
if USEPAIRMET:
  srcMETTag = cms.InputTag("corrMVAMET") if (IsMC and APPLYMETCORR) else cms.InputTag("MVAMET", "MVAMET")
else:
  srcMETTag = cms.InputTag(PFMetName, "", "TEST")

## ----------------------------------------------------------------------
## SV fit
## ----------------------------------------------------------------------
if USECLASSICSVFIT:
    print "Using CLASSIC_SV_FIT"
    process.SVllCand = cms.EDProducer("ClassicSVfitInterface",
                                      srcPairs   = cms.InputTag("barellCand"),
                                      srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                      srcCov     = cms.InputTag("METSignificance", "METCovariance"),
                                      usePairMET = cms.bool(USEPAIRMET),
                                      srcMET     = srcMETTag,
                                      computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False)
    )
else:
    print "Using STANDALONE_SV_FIT"
    process.SVllCand = cms.EDProducer("SVfitInterface",
                                      srcPairs   = cms.InputTag("barellCand"),
                                      srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                      srcCov     = cms.InputTag("METSignificance", "METCovariance"),
                                      usePairMET = cms.bool(USEPAIRMET),
                                      srcMET     = srcMETTag,
                                      computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False)
    )

## ----------------------------------------------------------------------
## SV fit BYPASS (skip SVfit, don't compute SVfit pair mass)
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
                      fileName = cms.untracked.string ("CosaACaso"),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      doCPVariables = cms.bool(doCPVariables),               
                      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      secVtxCollection = cms.InputTag("slimmedSecondaryVertices"), # FRA
                      puCollection = cms.InputTag("slimmedAddPileupInfo"),
                      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoMiniRelIsoCollection = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                      rhoForJER = cms.InputTag("fixedGridRhoAll"), # FRA
                      PFCandCollection = cms.InputTag("packedPFCandidates"),
                      jetCollection = cms.InputTag("jets"),
                      #JECset = cms.untracked.string("patJetCorrFactors"),
                      JECset = cms.untracked.string("patJetCorrFactorsUpdatedJEC"),
                      computeQGVar = cms.bool(COMPUTEQGVAR),
                      QGTagger = cms.InputTag("QGTagger", "qgLikelihood"),
                      ak8jetCollection = cms.InputTag("slimmedJetsAK8"),
                      lepCollection = cms.InputTag("softLeptons"),
                      lheCollection = cms.InputTag("LHEEventProduct"),
                      genCollection = cms.InputTag("generator"),
                      genericCollection = cms.InputTag("genInfo"),
                      genjetCollection = cms.InputTag("slimmedGenJets"),
                      totCollection = cms.InputTag("nEventsTotal"),
                      passCollection = cms.InputTag("nEventsPassTrigger"),
                      lhepCollection = cms.InputTag("externalLHEProducer"),
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", HLTProcessName), #Different names for MiniAODv2 at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.                      
                      #triggerSet = cms.InputTag("selectedPatTrigger"), # FRA
                      triggerSet = cms.InputTag("slimmedPatTrigger"),    # FRA
                      triggerList = HLTLIST,
                      metFilters = cms.InputTag ("TriggerResults","",METfiltersProcess),
                      PUPPImetCollection = cms.InputTag("slimmedMETsPuppi"),
                      srcPFMETCov = cms.InputTag("METSignificance", "METCovariance"),
                      srcPFMETSignificance = cms.InputTag("METSignificance", "METSignificance"),
                      l1extraIsoTau = cms.InputTag("l1extraParticles", "IsoTau"),
                      HT = cms.InputTag("externalLHEProducer"),
                      beamSpot = cms.InputTag("offlineBeamSpot"),
                      nBadMu = cms.InputTag("removeBadAndCloneGlobalMuons"),
                      genLumiHeaderTag = cms.InputTag("generator")
                      )
if USE_NOHFMET:
    process.HTauTauTree.metCollection = cms.InputTag("slimmedMETsNoHF")
else: 
    process.HTauTauTree.metCollection = cms.InputTag("slimmedMETs", "", "TEST") # use TEST so that I get the corrected one

if SVFITBYPASS:
    process.HTauTauTree.candCollection = cms.InputTag("SVbypass")
    process.SVFit = cms.Sequence (process.SVbypass)


else:
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

# Prepare lepton collections
process.Candidates = cms.Sequence(
    # process.printTree         + # just for debug, print MC particles
    process.nEventsTotal       +
    #process.hltFilter         + 
    process.nEventsPassTrigger +
    process.muons              +
    process.electrons          + process.cleanSoftElectrons +
    process.taus               +
    process.fsrSequence        +
    process.softLeptons        + process.barellCand +
    #process.jets              +
    process.jecSequence + process.jetSequence + #process.jets + 
    process.METSequence        +
    process.geninfo            +
    process.SVFit
    )
# always run ntuplizer
process.trees = cms.EndPath(process.HTauTauTree)# + process.HTauTauTreeTauUp + process.HTauTauTreeTauDown)

