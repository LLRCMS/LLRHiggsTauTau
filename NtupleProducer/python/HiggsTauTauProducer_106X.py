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
print 'Year+Period:', str(YEAR)+PERIOD
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
  print 'Using HLT trigger 2016'
  execfile(PyFilePath+"python/triggers_80X.py")  # 2016 triggers and filters
if YEAR == 2017:
  print 'Using HLT trigger 2017'
  execfile(PyFilePath+"python/triggers_92X.py")  # 2017 triggers and filters
if YEAR == 2018:
  print 'Using HLT trigger 2018'
  execfile(PyFilePath+"python/triggers_102X.py") # 2018 triggers and filters


### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
# From PPD table https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis updated 12-05-2021
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")    
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

if IsMC:
  if YEAR == 2016:
    if PERIOD=="preVFP":
        process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_preVFP_v9' # 2016 preVFP
    else:
        process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v15'       # 2016 postVFP
  if YEAR == 2017:
    process.GlobalTag.globaltag = '106X_mc2017_realistic_v8'             # 2017 MC
  if YEAR == 2018:
    process.GlobalTag.globaltag = '106X_upgrade2018_realistic_v15_L1v1'  # 2018 MC
else :
    process.GlobalTag.globaltag = '106X_dataRun2_v32'                    # Data

print "GT: ",process.GlobalTag.globaltag

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
prefireEra = "2016BtoH"
if YEAR==2017: prefireEra = "2017BtoF"

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
  DataEra = cms.string(prefireEra),
  UseJetEMPt = cms.bool(False),
  PrefiringRateSystematicUncty = cms.double(0.2),
  SkipWarnings = False
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

EgammaPostRecoSeq_ERA = '2016preVFP-UL'  # 2016 data postVFP
if PERIOD=='postVFP':
  EgammaPostRecoSeq_ERA = '2016postVFP-UL' # 2016 data preVFP
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

process.bareTaus = cms.EDFilter("PATTauRefSelector",
   src = cms.InputTag("slimmedTaus"), 
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

# TES corrections: https://indico.cern.ch/event/887196/contributions/3743090/attachments/1984772/3306737/TauPOG_TES_20200210.pdf

# EES corrections: https://indico.cern.ch/event/868279/contributions/3665970/attachments/1959265/3267731/FES_9Dec_explained.pdf

# NominalTESCorrection=-1#in percent\
APPLYTESCORRECTION = APPLYTESCORRECTION if IsMC else False # always false if data

# 2016 data - MVAoldDM2017v2
#NomTESUncDM0   = cms.double(1.0)  # in percent, up/down uncertainty of TES
#NomTESUncDM1   = cms.double(0.9)  # in percent, up/down uncertainty of TES
#NomTESUncDM10  = cms.double(1.1)  # in percent, up/down uncertainty of TES
#NomTESUncDM11  = --> Missing <--  # in percent, up/down uncertainty of TES
#NomTESCorDM0   = cms.double(-0.6) # DecayMode==0
#NomTESCorDM1   = cms.double(-0.5) # DecayMode==1
#NomTESCorDM10  = cms.double(0.0)  # DecayMode==10
#NomTESCorDM11  = --> Missing <--  # DecayMode==11

# 2017 data - MVAoldDM2017v2
#if YEAR == 2017:
#    NomTESUncDM0   = cms.double(0.8)  # in percent, up/down uncertainty of TES
#    NomTESUncDM1   = cms.double(0.8)  # in percent, up/down uncertainty of TES
#    NomTESUncDM10  = cms.double(0.9)  # in percent, up/down uncertainty of TES
#    NomTESUncDM11  = cms.double(1.0)  # in percent, up/down uncertainty of TES
#    NomTESCorDM0   = cms.double(0.7)  # DecayMode==0
#    NomTESCorDM1   = cms.double(-0.2) # DecayMode==1
#    NomTESCorDM10  = cms.double(0.1)  # DecayMode==10
#    NomTESCorDM11  = cms.double(-0.1) # DecayMode==11

# 2018 data - MVAoldDM2017v2
#if YEAR == 2018:
#    NomTESUncDM0   = cms.double(1.1)  # in percent, up/down uncertainty of TES
#    NomTESUncDM1   = cms.double(0.9)  # in percent, up/down uncertainty of TES
#    NomTESUncDM10  = cms.double(0.8)  # in percent, up/down uncertainty of TES
#    NomTESUncDM11  = --> Missing <--  # in percent, up/down uncertainty of TES
#    NomTESCorDM0   = cms.double(-1.3) # DecayMode==0
#    NomTESCorDM1   = cms.double(-0.5) # DecayMode==1
#    NomTESCorDM10  = cms.double(-1.2) # DecayMode==10
#    NomTESCorDM11  = --> Missing <--  # DecayMode==11

# 2016 data - DeepTau2017v2p1
NomTESUncDM0      = cms.double(0.8)  # in percent, up/down uncertainty of TES
NomTESUncDM1      = cms.double(0.6)  # in percent, up/down uncertainty of TES
NomTESUncDM10     = cms.double(0.8)  # in percent, up/down uncertainty of TES
NomTESUncDM11     = cms.double(1.1)  # in percent, up/down uncertainty of TES
NomTESCorDM0      = cms.double(-0.9) # DecayMode==0
NomTESCorDM1      = cms.double(-0.1) # DecayMode==1
NomTESCorDM10     = cms.double(0.3)  # DecayMode==10
NomTESCorDM11     = cms.double(-0.2) # DecayMode==11

#EES BARREL
NomEFakeESCorDM0B     = cms.double(0.679) #DecayMode==0
NomEFakeESUncDM0BUp    = cms.double(0.806) #DecayMode==0
NomEFakeESUncDM0BDown  = cms.double(0.982) #DecayMode==0
NomEFakeESCorDM1B      = cms.double(3.389) #DecayMode==1
NomEFakeESUncDM1BUp    = cms.double(1.168) #DecayMode==1
NomEFakeESUncDM1BDown  = cms.double(2.475) #DecayMode==1
#EES ENDCAP
NomEFakeESCorDM0E      = cms.double(-3.5)   #DecayMode==0
NomEFakeESUncDM0EUp    = cms.double(1.808)  #DecayMode==0
NomEFakeESUncDM0EDown  = cms.double(1.102)  #DecayMode==0
NomEFakeESCorDM1E      = cms.double(5.)      #DecayMode==1
NomEFakeESUncDM1EUp    = cms.double(6.57)   #DecayMode==1
NomEFakeESUncDM1EDown  = cms.double(5.694)  #DecayMode==1

TESyear = "2016Legacy"

# 2017 data - DeepTau2017v2p1
if YEAR == 2017:
    NomTESUncDM0      = cms.double(1.0)  # in percent, up/down uncertainty of TES
    NomTESUncDM1      = cms.double(0.6)  # in percent, up/down uncertainty of TES
    NomTESUncDM10     = cms.double(0.7)  # in percent, up/down uncertainty of TES
    NomTESUncDM11     = cms.double(1.4)  # in percent, up/down uncertainty of TES
    NomTESCorDM0      = cms.double(0.4)  # DecayMode==0
    NomTESCorDM1      = cms.double(0.2)  # DecayMode==1
    NomTESCorDM10     = cms.double(0.1)  # DecayMode==10
    NomTESCorDM11     = cms.double(-1.3) # DecayMode==1

    #EES BARREL
    NomEFakeESCorDM0B      = cms.double(0.911) #DecayMode==0
    NomEFakeESUncDM0BUp    = cms.double(1.343) #DecayMode==0
    NomEFakeESUncDM0BDown  = cms.double(0.882) #DecayMode==0
    NomEFakeESCorDM1B      = cms.double(1.154) #DecayMode==1
    NomEFakeESUncDM1BUp    = cms.double(2.162) #DecayMode==1
    NomEFakeESUncDM1BDown  = cms.double(0.973) #DecayMode==1
    #EES ENDCAP
    NomEFakeESCorDM0E      = cms.double(-2.604)   #DecayMode==0
    NomEFakeESUncDM0EUp    = cms.double(2.249)    #DecayMode==0
    NomEFakeESUncDM0EDown  = cms.double(1.43)     #DecayMode==0
    NomEFakeESCorDM1E      = cms.double(1.5)    #DecayMode==1
    NomEFakeESUncDM1EUp    = cms.double(6.461)      #DecayMode==1
    NomEFakeESUncDM1EDown  = cms.double(4.969)    #DecayMode==1

    TESyear = "2017ReReco"

# 2018 data - DeepTau2017v2p1
if YEAR == 2018:
    NomTESUncDM0          = cms.double(0.9)  # in percent, up/down uncertainty of TES
    NomTESUncDM1          = cms.double(0.5)  # in percent, up/down uncertainty of TES
    NomTESUncDM10         = cms.double(0.7)  # in percent, up/down uncertainty of TES
    NomTESUncDM11         = cms.double(1.2)  # in percent, up/down uncertainty of TES
    NomTESCorDM0          = cms.double(-1.6) # DecayMode==0
    NomTESCorDM1          = cms.double(-0.5) # DecayMode==1
    NomTESCorDM10         = cms.double(-1.2) # DecayMode==10
    NomTESCorDM11         = cms.double(-0.4) # DecayMode==11

    #EES BARREL
    NomEFakeESCorDM0B      = cms.double(1.362)    #DecayMode==0
    NomEFakeESUncDM0BUp    = cms.double(0.904)    #DecayMode==0
    NomEFakeESUncDM0BDown  = cms.double(0.474)    #DecayMode==0
    NomEFakeESCorDM1B      = cms.double(1.954)    #DecayMode==1
    NomEFakeESUncDM1BUp    = cms.double(1.226)    #DecayMode==1
    NomEFakeESUncDM1BDown  = cms.double(1.598)    #DecayMode==1
    #EES ENDCAP
    NomEFakeESCorDM0E      = cms.double(-3.097)   #DecayMode==0
    NomEFakeESUncDM0EUp    = cms.double(3.404)    #DecayMode==0
    NomEFakeESUncDM0EDown  = cms.double(1.25)     #DecayMode==0
    NomEFakeESCorDM1E      = cms.double(-1.5)     #DecayMode==1
    NomEFakeESUncDM1EUp    = cms.double(5.499)    #DecayMode==1
    NomEFakeESUncDM1EDown  = cms.double(4.309)    #DecayMode==1

    TESyear = "2018ReReco"

process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   discriminator = cms.string(TAUDISCRIMINATOR),

   NominalTESUncertaintyDM0         = NomTESUncDM0,
   NominalTESUncertaintyDM1         = NomTESUncDM1,
   NominalTESUncertaintyDM10        = NomTESUncDM10,
   NominalTESUncertaintyDM11        = NomTESUncDM11,
   NominalTESCorrectionDM0          = NomTESCorDM0,
   NominalTESCorrectionDM1          = NomTESCorDM1,
   NominalTESCorrectionDM10         = NomTESCorDM10,
   NominalTESCorrectionDM11         = NomTESCorDM11,

   NominalEFakeESCorrectionDM0B      = NomEFakeESCorDM0B,
   NominalEFakeESUncertaintyDM0BUp   = NomEFakeESUncDM0BUp, 
   NominalEFakeESUncertaintyDM0BDown = NomEFakeESUncDM0BDown, 
   NominalEFakeESCorrectionDM1B      = NomEFakeESCorDM1B,
   NominalEFakeESUncertaintyDM1BUp   = NomEFakeESUncDM1BUp, 
   NominalEFakeESUncertaintyDM1BDown = NomEFakeESUncDM1BDown, 
   NominalEFakeESCorrectionDM0E      = NomEFakeESCorDM0E,
   NominalEFakeESUncertaintyDM0EUp   = NomEFakeESUncDM0EUp, 
   NominalEFakeESUncertaintyDM0EDown = NomEFakeESUncDM0EDown, 
   NominalEFakeESCorrectionDM1E      = NomEFakeESCorDM1E,
   NominalEFakeESUncertaintyDM1EUp   = NomEFakeESUncDM1EUp, 
   NominalEFakeESUncertaintyDM1EDown = NomEFakeESUncDM1EDown, 

   ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
   # ApplyTESUpDown = cms.bool(True if IsMC else False), # no shift computation when data
   flags = cms.PSet(
        isGood = cms.string("")
        ),

   year = cms.string(TESyear)
   )

process.taus=cms.Sequence(process.bareTaus + process.softTaus)

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


# Update jet collection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
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


else:
    print "Using event pfMET (same MET for all pairs)"

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
process.PVfilter = cms.Path(process.goodPrimaryVertices)
process.l1ECALPref = cms.Path(process.prefiringweight)

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
    process.softLeptons        + process.barellCand +
    process.jecSequence        + process.jetSequence +
    process.METSequence        +
    process.geninfo            +
    process.SVFit
    )
# always run ntuplizer
process.trees = cms.EndPath(process.HTauTauTree)

