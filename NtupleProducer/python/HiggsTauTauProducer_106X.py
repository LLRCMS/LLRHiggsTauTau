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
try: YEAR
except NameError:
    YEAR = 2018
try: PERIOD
except:
    PERIOD =""
print ('Year+Period:', str(YEAR)+PERIOD)
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
### Trigger list
### ----------------------------------------------------------------------
if YEAR == 2016:
  print ('Using HLT trigger 2016')
  exec(open(PyFilePath+"python/triggers_80X.py").read())  # 2016 triggers and filters
if YEAR == 2017:
  print ('Using HLT trigger 2017')
  exec(open(PyFilePath+"python/triggers_92X.py").read())  # 2017 triggers and filters
if YEAR == 2018:
  print ('Using HLT trigger 2018')
  exec(open(PyFilePath+"python/triggers_102X.py").read()) # 2018 triggers and filters


### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
# From PPD table https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis updated 12-05-2021
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")    
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

if IsMC:
  if YEAR == 2016:
    if PERIOD=="postVFP":
        process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v17'       # 2016 postVFP
    else:
        process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_preVFP_v11' # 2016 preVFP
  if YEAR == 2017:
    process.GlobalTag.globaltag = '106X_mc2017_realistic_v10'             # 2017 MC
  if YEAR == 2018:
    process.GlobalTag.globaltag = '106X_upgrade2018_realistic_v16_L1v1'  # 2018 MC
else :
    process.GlobalTag.globaltag = '106X_dataRun2_v37'                    # Data

print ("GT: ",process.GlobalTag.globaltag)

nanosec="25"
if not Is25ns: nanosec="50"

METfiltersProcess = "PAT" # LP: both for data and MC with UL
#METfiltersProcess = "PAT" if IsMC else "RECO" # NB! this is not guaranteed to be true! the following is valid on 2015 Run C + Run D data. Check:
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

### ----------------------------------------------------------------------
### L1ECALPrefiringWeightRecipe (for 2016 and 2017 MC only)
### https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
### ----------------------------------------------------------------------
prefireEraECAL = "UL2016preVFP"
if PERIOD=="postVFP" : prefireEraECAL = "UL2016postVFP"
if YEAR==2017: prefireEraECAL = "UL2017BtoF"
if YEAR==2018: prefireEraECAL = "None"

prefireEraMuon = "2016preVFP"
if PERIOD=="postVFP" : prefireEraMuon = "2016postVFP"
if YEAR==2017: prefireEraMuon = "20172018"
if YEAR==2018: prefireEraMuon = "20172018"


from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"), 
DataEraECAL = cms.string(prefireEraECAL), 
DataEraMuon = cms.string(prefireEraMuon), 
UseJetEMPt = cms.bool(False),
PrefiringRateSystematicUnctyECAL = cms.double(0.2),
PrefiringRateSystematicUnctyMuon = cms.double(0.2)
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


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string(MUCUT),
)

process.softMuons = cms.EDProducer("MuFiller",
    src = cms.InputTag("bareSoftMuons"),
    genCollection = cms.InputTag("prunedGenParticles"),
    rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
    vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    sampleType = cms.int32(LEPTON_SETUP),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    cut = cms.string(""),
    flags = cms.PSet(
        ID = cms.string("userFloat('isPFMuon')" ), # PF ID
        isGood = cms.string(MUCUT)
    )
)

process.muons =  cms.Sequence(process.bareSoftMuons+ process.softMuons)

###
### Electrons
###


# START ELECTRON CUT BASED ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#**********************
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
#**********************

EgammaPostRecoSeq_ERA = '2016preVFP-UL'  # 2016 data preVFP
if PERIOD=='postVFP':
  EgammaPostRecoSeq_ERA = '2016postVFP-UL' # 2016 data postVFP
if YEAR==2017:
  EgammaPostRecoSeq_ERA = '2017-UL'       # 2017 data
if YEAR == 2018:
  EgammaPostRecoSeq_ERA = '2018-UL'       # 2018 data

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era=EgammaPostRecoSeq_ERA)    
		       
process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("slimmedElectrons"),
   rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
   vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   genCollection = cms.InputTag("prunedGenParticles"),
   sampleType = cms.int32(LEPTON_SETUP),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string(ELECUT)
   )

process.electrons = cms.Sequence(process.softElectrons+process.egammaPostRecoSeq)


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

updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig

tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, 
                    updatedTauName = updatedTauName,
                    toKeep = [
                            "deepTau2017v2p1", "deepTau2018v2p5",
                            ])

tauIdEmbedder.runTauID()

process.bareTaus = cms.EDFilter("PATTauRefSelector",
   src = cms.InputTag("slimmedTausNewID"), 
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

# TES: https://github.com/cms-tau-pog/TauIDSFs
# NominalTESCorrection=-1#in percent\
APPLYTESCORRECTION = APPLYTESCORRECTION if IsMC else False # always false if data

TESyear = "UL2016_preVFP"

if PERIOD=='postVFP':
    TESyear = 'UL2016_postVFP'

if YEAR == 2017:
    TESyear = "UL2017"

if YEAR == 2018:
    TESyear = "UL2018"

process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   discriminator = cms.string(TAUDISCRIMINATOR),

   ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
   flags = cms.PSet(
        isGood = cms.string("")
        ),

   year = cms.string(TESyear)
   )

process.taus=cms.Sequence(process.rerunMvaIsolationSequence + process.slimmedTausNewID + process.bareTaus + process.softTaus)

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
    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauString))
)


#
#Jets
#

# apply new jet energy corrections 
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

# JEC corrections
jecLevels = None
if IsMC:
    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
else:
    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]

## B-regression
# following https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/BJetRegression#How_to_treat_the_training
# Y value from json files (could be read directly from there, but this is for another time)

MODEL_DIR = "LLRHiggsTauTau/NtupleProducer/data/DNNmodels/"

if YEAR == 2018: # just one option for weight file
    Y_MEAN = 1.0545977354049683
    Y_STD  = 0.27912887930870056
    WEIGHT_FILE = "breg_training_2018_JECv8.pb"
elif YEAR == 2017: # Weight file to be checked
    Y_MEAN = 1.055067777633667
    Y_STD  = 0.28225210309028625
    WEIGHT_FILE = "breg_training_2017_jecV32.pb"
elif YEAR == 2016: # weight file marked as recommended
    Y_MEAN = 1.047176718711853
    Y_STD  = 0.31976690888404846
    WEIGHT_FILE = "breg_training_2016_JECv11_Oct_2019.pb"

process.bregJets = cms.EDProducer("bRegressionProducer",
                                  JetTag = cms.InputTag("slimmedJets"),
                                  rhoFixedGridCollection = cms.InputTag("fixedGridRhoAll"),
                                  bRegressionWeightfile = cms.FileInPath(MODEL_DIR+WEIGHT_FILE),
                                  pvsrc = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                  svsrc = cms.InputTag('slimmedSecondaryVertices'),
                                  y_mean = cms.double(Y_MEAN),
                                  y_std = cms.double(Y_STD),
)
					  
# Update jet collection
updateJetCollection(
   process,
   jetSource = cms.InputTag('bregJets'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None'),
   labelName = 'UpdatedJEC'
)

# Update the jet sequences
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC *
                                   process.updatedPatJetsUpdatedJEC *
                                   process.selectedUpdatedPatJetsUpdatedJEC)

# Jet Selector after JEC and bTagging
process.jets = cms.EDFilter("PATJetRefSelector",
                            src = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
                            cut = cms.string(JETCUT))

##
## QG tagging for jets
##

if COMPUTEQGVAR:

    QGlikelihood_tag = 'QGLikelihoodObject_v1_AK4PFchs'
    if YEAR == 2017 or YEAR == 2018:
      QGlikelihood_tag = 'QGLikelihoodObject_v1_AK4PFchs_2017'

    from CondCore.CondDB.CondDB_cfi import CondDB
     
    process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDB.clone(
        connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
      ),
      toGet = cms.VPSet(
        cms.PSet(
          record = cms.string('QGLikelihoodRcd'),
          tag    = cms.string(QGlikelihood_tag),
          label  = cms.untracked.string('QGL_AK4PFchs'),
        ),
      ),
    )
    process.es_prefer_qg = cms.ESPrefer("PoolDBESSource", "QGPoolDBESSource")

    process.load('RecoJets.JetProducers.QGTagger_cfi')
    process.QGTagger.srcJets          = cms.InputTag("jets")    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
    process.jetSequence = cms.Sequence(process.jets * process.QGTagger)

else:
    process.jetSequence = cms.Sequence(process.jets)

# Add latest pileup jet ID
process.load("RecoJets.JetProducers.PileupJetID_cfi")
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16, _chsalgos_106X_UL16APV, _chsalgos_106X_UL17, _chsalgos_106X_UL18

if YEAR == 2016:
    
    PUalgo = _chsalgos_106X_UL16APV
    
    if PERIOD=='postVFP':
        PUalgo = _chsalgos_106X_UL16

if YEAR == 2017:
    PUalgo = _chsalgos_106X_UL17

if YEAR == 2018:
    PUalgo = _chsalgos_106X_UL18

process.pileupJetIdUpdated = process.pileupJetId.clone(
   jets = cms.InputTag("jets"),
   inputIsCorrected = True,
   applyJec = False,
   vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
   algos = PUalgo
)
process.jetSequence += cms.Sequence(process.pileupJetIdUpdated)

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

## ----------------------------------------------------------------------
## MVA MET
## ----------------------------------------------------------------------

process.METSequence = cms.Sequence()
if USEPAIRMET:
    print ("Using pair MET (MVA MET)")
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


else:
    print ("Using event pfMET (same MET for all pairs)")

    PFMetName = "slimmedMETs"
    uncorrPFMetTag = cms.InputTag(PFMetName)


    # patch to get a standalone MET significance collection
    process.METSignificance = cms.EDProducer ("ExtractMETSignificance",
                                                  #srcMET=cms.InputTag(PFMetName,"","TEST")
                                                  srcMET=uncorrPFMetTag
                                                  )

    # add variables with MET shifted for TES corrections
    process.ShiftMETforTES = cms.EDProducer ("ShiftMETforTES",
                                             #srcMET  = cms.InputTag(PFMetName,"","TEST"),
                                             srcMET  = uncorrPFMetTag,
                                             tauCollection = cms.InputTag("softTaus")
                                             )

    # add variables with MET shifted for EES corrections (E->tau ES)
    process.ShiftMETforEES = cms.EDProducer ("ShiftMETforEES",
                                             #srcMET  = cms.InputTag(PFMetName,"","TEST"),
                                             srcMET  = uncorrPFMetTag,
                                             tauCollection = cms.InputTag("softTaus")
                                             )

    # Shift met due to central corrections of TES and EES
    process.ShiftMETcentral = cms.EDProducer ("ShiftMETcentral",
                                              srcMET = uncorrPFMetTag,
                                              tauUncorrected = cms.InputTag("bareTaus"),
                                              tauCorrected = cms.InputTag("softTaus")
                                              )


    # Get a standalone Puppi MET significance collection
    process.PuppiMETSignificance = cms.EDProducer ("ExtractMETSignificance",
                                                   srcMET=cms.InputTag("slimmedMETsPuppi")
                                                  )

    # Shift PUPPI met due to central corrections of TES and EES
    process.ShiftPuppiMETcentral = cms.EDProducer ("ShiftMETcentral",
                                                   srcMET = cms.InputTag("slimmedMETsPuppi"),
                                                   tauUncorrected = cms.InputTag("bareTaus"),
                                                   tauCorrected = cms.InputTag("softTaus")
                                                  )

    process.METSequence += process.METSignificance
    process.METSequence += process.ShiftMETforTES
    process.METSequence += process.ShiftMETforEES
    process.METSequence += process.ShiftMETcentral
    process.METSequence += process.PuppiMETSignificance
    process.METSequence += process.ShiftPuppiMETcentral

    # Add both DeepMET tunes from the METCorrectionLevels:
    # https://github.com/cms-sw/cmssw/blob/ef3151a362db68a97ef10c4ff993af637cd163ef/DataFormats/PatCandidates/interface/MET.h#L173-L190
    for tuneName, tuneIdx in zip(["RawDeepResponseTune", "RawDeepResolutionTune"],[13, 14]):

        # Add standalone DeepMET collections
        DeepMET = cms.EDProducer ("CorrectedMETCollectionProducer",
                                  srcMET = cms.InputTag("slimmedMETs"),
                                  correctionLevel = cms.int32(tuneIdx),
                                  )

        setattr(process, "DeepMET" + tuneName, DeepMET)

        # Shift DeepMET due to central corrections of TES and EES
        ShiftDeepMETcentral = cms.EDProducer ("ShiftMETcentral",
                                              srcMET = cms.InputTag("DeepMET" + tuneName),
                                              tauUncorrected = cms.InputTag("bareTaus"),
                                              tauCorrected = cms.InputTag("softTaus"),
                                              )

        setattr(process, "ShiftDeepMETcentral" + tuneName, ShiftDeepMETcentral)

        process.METSequence += getattr(process, "DeepMET" + tuneName)
        process.METSequence += getattr(process, "ShiftDeepMETcentral" + tuneName)

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
  # MET corrected for central TES and EES shifts of the taus
  srcMETTag = cms.InputTag("ShiftMETcentral")

## ----------------------------------------------------------------------
## SV fit
## ----------------------------------------------------------------------
#if USECLASSICSVFIT:
#    print "Using CLASSIC_SV_FIT"
process.SVllCand = cms.EDProducer("ClassicSVfitInterface",
                                  srcPairs   = cms.InputTag("barellCand"),
                                  srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                  srcCov     = cms.InputTag("METSignificance", "METCovariance"),
                                  usePairMET = cms.bool(USEPAIRMET),
                                  srcMET     = srcMETTag,
                                  computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False),
                                  computeForUpDownMET = cms.bool(COMPUTEMETUPDOWNSVFIT if IsMC else False),
                                  METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                  METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                  METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                  METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN"),
                                  METdxUP_EES   = cms.InputTag("ShiftMETforEES", "METdxUPEES"),
                                  METdyUP_EES   = cms.InputTag("ShiftMETforEES", "METdyUPEES"),
                                  METdxDOWN_EES = cms.InputTag("ShiftMETforEES", "METdxDOWNEES"),
                                  METdyDOWN_EES = cms.InputTag("ShiftMETforEES", "METdyDOWNEES")
)
#else:
#    print "Using STANDALONE_SV_FIT"
#    process.SVllCand = cms.EDProducer("SVfitInterface",
#                                      srcPairs   = cms.InputTag("barellCand"),
#                                      srcSig     = cms.InputTag("METSignificance", "METSignificance"),
#                                      srcCov     = cms.InputTag("METSignificance", "METCovariance"),
#                                      usePairMET = cms.bool(USEPAIRMET),
#                                      srcMET     = srcMETTag,
#                                      computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False)
#    )

## ----------------------------------------------------------------------
## SV fit BYPASS (skip SVfit, don't compute SVfit pair mass)
## ----------------------------------------------------------------------
process.SVbypass = cms.EDProducer ("SVfitBypass",
                                    srcPairs   = cms.InputTag("barellCand"),
                                    usePairMET = cms.bool(USEPAIRMET),
                                    srcMET     = srcMETTag,
                                    srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                    srcCov     = cms.InputTag("METSignificance", "METCovariance"),
                                    METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                    METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                    METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                    METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN"),
                                    METdxUP_EES   = cms.InputTag("ShiftMETforEES", "METdxUPEES"),
                                    METdyUP_EES   = cms.InputTag("ShiftMETforEES", "METdyUPEES"),
                                    METdxDOWN_EES = cms.InputTag("ShiftMETforEES", "METdxDOWNEES"),
                                    METdyDOWN_EES = cms.InputTag("ShiftMETforEES", "METdyDOWNEES")
)


## ----------------------------------------------------------------------
## Ntuplizer
## ----------------------------------------------------------------------
process.HTauTauTree = cms.EDAnalyzer("HTauTauNtuplizer",
                      fileName = cms.untracked.string ("CosaACaso"),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      year = cms.int32(YEAR),
                      period = cms.string(PERIOD),
                      doCPVariables = cms.bool(doCPVariables),               
                      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      secVtxCollection = cms.InputTag("slimmedSecondaryVertices"), # FRA
                      puCollection = cms.InputTag("slimmedAddPileupInfo"),
                      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoMiniRelIsoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoForJER = cms.InputTag("fixedGridRhoAll"), # FRA
                      PFCandCollection = cms.InputTag("packedPFCandidates"),
                      jetCollection = cms.InputTag("jets"),
                      JECset = cms.untracked.string(""),   # specified later
                      pileupJetIdUpdatedDiscr = cms.InputTag("pileupJetIdUpdated", "fullDiscriminant"),
                      pileupJetIdUpdatedWP    = cms.InputTag("pileupJetIdUpdated", "fullId"),
                      computeQGVar = cms.bool(COMPUTEQGVAR),
                      QGTagger = cms.InputTag("QGTagger", "qgLikelihood"),
                      stage2TauCollection = cms.InputTag("caloStage2Digis","Tau"),
                      stage2JetCollection = cms.InputTag("caloStage2Digis","Jet"),
                      ak8jetCollection = cms.InputTag("slimmedJetsAK8"),
                      lepCollection = cms.InputTag("softLeptons"),
                      lheCollection = cms.InputTag("LHEEventProduct"),
                      genCollection = cms.InputTag("generator"),
                      genericCollection = cms.InputTag("genInfo"),
                      genjetCollection = cms.InputTag("slimmedGenJets"),
                      totCollection = cms.InputTag("nEventsTotal"),
                      passCollection = cms.InputTag("nEventsPassTrigger"),
                      lhepCollection = cms.InputTag("externalLHEProducer"),
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", HLTProcessName), #Different names for MiniAODv2 at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.                                        #triggerSet = cms.InputTag("selectedPatTrigger"), # FRA
                      triggerSet = cms.InputTag("slimmedPatTrigger"),    # FRA
                      triggerList = HLTLIST,
                      metFilters = cms.InputTag ("TriggerResults","",METfiltersProcess),
                      PUPPImetCollection = cms.InputTag("slimmedMETsPuppi"),
                      metPuppiShiftedCollection = cms.InputTag("ShiftPuppiMETcentral"),
                      srcPuppiMETCov = cms.InputTag("PuppiMETSignificance", "METCovariance"),
                      srcPuppiMETSignificance = cms.InputTag("PuppiMETSignificance", "METSignificance"),                      
                      srcPFMETCov = cms.InputTag("METSignificance", "METCovariance"),
                      srcPFMETSignificance = cms.InputTag("METSignificance", "METSignificance"),
                      metERCollection = uncorrPFMetTag, # save the uncorrected MET (for TES and EES central shifts) just for reference
                      DeepMETResponseTuneCollection = cms.InputTag("DeepMETRawDeepResponseTune"),
                      DeepMETResolutionTuneCollection = cms.InputTag("DeepMETRawDeepResolutionTune"),
                      ShiftedDeepMETResponseTuneCollection = cms.InputTag("ShiftDeepMETcentralRawDeepResponseTune"),
                      ShiftedDeepMETResolutionTuneCollection = cms.InputTag("ShiftDeepMETcentralRawDeepResolutionTune"),
                      HT = cms.InputTag("externalLHEProducer"),
                      beamSpot = cms.InputTag("offlineBeamSpot"),
                      genLumiHeaderTag = cms.InputTag("generator"),
                      L1prefireProb     = cms.InputTag("prefiringweight:nonPrefiringProb"),
                      L1prefireProbUp   = cms.InputTag("prefiringweight:nonPrefiringProbUp"),
                      L1prefireProbDown = cms.InputTag("prefiringweight:nonPrefiringProbDown")
)
if USE_NOHFMET:
    process.HTauTauTree.metCollection = cms.InputTag("slimmedMETsNoHF")
else:
    # MET corrected for central TES and EES shifts of the taus
    process.HTauTauTree.metCollection = srcMETTag

process.HTauTauTree.JECset = cms.untracked.string("patJetCorrFactorsUpdatedJEC")

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
process.PVfilter    = cms.Path(process.goodPrimaryVertices)
process.l1ECALPref  = cms.Path(process.prefiringweight)

# Prepare lepton collections
process.Candidates = cms.Sequence(
    # process.printTree         + # just for debug, print MC particles
    process.nEventsTotal       +
    #process.hltFilter         + 
    process.nEventsPassTrigger +
    process.egammaPostRecoSeq  +
    process.muons              +
    process.electrons          + process.cleanSoftElectrons +
    process.taus               +
    process.fsrSequence        +
    process.softLeptons        + process.barellCand         +
    process.bregJets           + process.jecSequence        + process.jetSequence +
    process.METSequence        +
    process.geninfo            +
    process.SVFit
    )
# always run ntuplizer
process.trees = cms.EndPath(process.HTauTauTree)

