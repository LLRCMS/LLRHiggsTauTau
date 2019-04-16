import FWCore.ParameterSet.Config as cms
execfile(PyFilePath+"python/triggers_92X.py") # contains the list of triggers and filters

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
    YEAR=2016
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
    if YEAR==2016: process.GlobalTag.globaltag = '102X_mcRun2_asymptotic_v6'
    #if YEAR==2016: process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'
    elif YEAR==2017: process.GlobalTag.globaltag = '102X_mc2017_realistic_v6'
    #elif YEAR==2017: process.GlobalTag.globaltag = '94X_mc2017_realistic_v17'
    elif YEAR==2018: process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v12'
else :
    if YEAR==2016: process.GlobalTag.globaltag = '94X_dataRun2_v10'
    elif YEAR==2017: process.GlobalTag.globaltag = '94X_dataRun2_v11'
    elif YEAR==2018: process.GlobalTag.globaltag = '102X_dataRun2_Sep2018Rereco_v1'
print "GT: ",process.GlobalTag.globaltag

nanosec="25"
if not Is25ns: nanosec="50"

LEPTON_SETUP_LEGACY = YEAR 
print "Lepton setup: ", LEPTON_SETUP_LEGACY

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


### ----------------------------------------------------------------------
### MC stuff
### ----------------------------------------------------------------------

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

### ----------------------------------------------------------------------
### ECAL bad calibration filter
### ----------------------------------------------------------------------

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
)


### ----------------------------------------------------------------------
### Muons
### ----------------------------------------------------------------------

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

process.noBadGlobalMuons = cms.Sequence(process.cloneGlobalMuonTagger + process.badGlobalMuonTagger + process.removeBadAndCloneGlobalMuons) # in tagging mode, these modules return always "true"

# mu isolation from nanoAOD
# different EA for the different years
# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/electrons_cff.py#L114-L120
if YEAR==2016:
        eafileminiisomu = "PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt"
elif YEAR == 2017 or YEAR == 2018:
        eafileminiisomu = "PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt"

process.isoForMu = cms.EDProducer("MuonIsoValueMapProducer",
    src = cms.InputTag("bareSoftMuons"),
    relative = cms.bool(True),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    EAFile_MiniIso = cms.FileInPath(eafileminiisomu),
)

#mu ptRatioRel from nanoAOD
process.ptRatioRelForMu = cms.EDProducer("MuonJetVarProducer",
    srcJet = cms.InputTag("updatedJets"),
    srcLep = cms.InputTag("bareSoftMuons"),
    srcVtx = cms.InputTag("offlineSlimmedPrimaryVertices"),
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
    miniRelIsoChgCollection = cms.InputTag("isoForMu:miniIsoChg"), #nanoAOD
    miniRelIsoAllCollection = cms.InputTag("isoForMu:miniIsoAll"), #nanoAOD
    ptRatioCollection = cms.InputTag("ptRatioRelForMu:ptRatio"), #nanoAOD
    ptRelCollection = cms.InputTag("ptRatioRelForMu:ptRel"), #nanoAOD
    jetNDauChargedMVASelCollection = cms.InputTag("ptRatioRelForMu:jetNDauChargedMVASel"), #nanoAOD
    sampleType = cms.int32(LEPTON_SETUP),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    lep_setup = cms.int32(LEPTON_SETUP_LEGACY),
    cut = cms.string(""),
    flags = cms.PSet(
        ID = cms.string("userFloat('isPFMuon')" ), # PF ID
        isGood = cms.string(MUCUT)
    )
)

# process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons+ process.softMuons)
process.muons =  cms.Sequence(process.noBadGlobalMuons + process.bareSoftMuons + process.isoForMu + process.ptRatioRelForMu + process.softMuons)    

### ----------------------------------------------------------------------
### Electrons
### ----------------------------------------------------------------------

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#**********************
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
#**********************

# Define which electron IDs we want to produce
my_id_modules =[
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',   #Spring16
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',   #Spring16 HZZ
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',     #Fall17 V1 iso
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',   #Fall17 V1 noIso
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',     #Fall17 V2 iso
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',   #Fall17 V2 noIso
    ] 

# Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#electron scale/smear

if (YEAR == 2016):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                          runEnergyCorrections=True,
                          runVID=True,
                          era='2016-Legacy')

if (YEAR == 2017):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                          runEnergyCorrections=True,
                          runVID=True,
                          era='2017-Nov17ReReco')

if (YEAR == 2018):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                          runEnergyCorrections=True,
                          runVID=True,
			  era='2018-Prompt')

# ele iso from nanoAOD
# different EA for the different years
# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/electrons_cff.py#L114-L120

if (YEAR==2016):
        eafileminiisoele = "RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"
        eafilepfisoele = "RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"

if (YEAR == 2017 or YEAR == 2018):
        eafileminiisoele = "RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"
        eafilepfisoele = "RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"

process.isoForEle = cms.EDProducer("EleIsoValueMapProducer",
    src = cms.InputTag("slimmedElectrons"),
    relative = cms.bool(True),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    rho_PFIso = cms.InputTag("fixedGridRhoFastjetAll"),
    EAFile_MiniIso = cms.FileInPath(eafileminiisoele),
    EAFile_PFIso = cms.FileInPath(eafilepfisoele),
)

#ele ptRatioRel from nanoAOD
process.ptRatioRelForEle = cms.EDProducer("ElectronJetVarProducer",
    srcJet = cms.InputTag("updatedJets"),
    srcLep = cms.InputTag("slimmedElectrons"),
    srcVtx = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("slimmedElectrons"),
   rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
   vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   genCollection = cms.InputTag("prunedGenParticles"),
   miniRelIsoChgCollection = cms.InputTag("isoForEle:miniIsoChg"), #nanoAOD
   miniRelIsoAllCollection = cms.InputTag("isoForEle:miniIsoAll"), #nanoAOD
   PFRelIsoChgCollection = cms.InputTag("isoForEle:PFIsoChg"), #nanoAOD
   PFRelIsoAllCollection = cms.InputTag("isoForEle:PFIsoAll"), #nanoAOD
   PFRelIsoAll04Collection = cms.InputTag("isoForEle:PFIsoAll04"), #nanoAOD
   ptRatioCollection = cms.InputTag("ptRatioRelForEle:ptRatio"), #nanoAOD
   ptRelCollection = cms.InputTag("ptRatioRelForEle:ptRel"), #nanoAOD
   jetNDauChargedMVASelCollection = cms.InputTag("ptRatioRelForEle:jetNDauChargedMVASel"), #nanoAOD
   sampleType = cms.int32(LEPTON_SETUP),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   lep_setup = cms.int32(LEPTON_SETUP_LEGACY),

   #MVA ELE ID
   eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpLoose"),
   eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp90"),
   eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp80"),
   eleLooseNoIsoIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wpLoose"),
   eleMediumNoIsoIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90"),
   eleTightNoIsoIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp80"),
   mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
   mvaNoIsoValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values"),
   mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Categories"),
   mvaNoIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Categories"),
   HZZmvaValuesMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
   cut = cms.string(ELECUT),
   flags = cms.PSet(
        ID = cms.string("userInt('isBDT')"), # BDT MVA ID
        isGood = cms.string("")
        )
   )


process.electrons = cms.Sequence(process.egammaPostRecoSeq + process.isoForEle + process.ptRatioRelForEle + process.egmGsfElectronIDSequence * process.softElectrons)

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

### ----------------------------------------------------------------------
### Taus
### ----------------------------------------------------------------------

from LLRHiggsTauTau.NtupleProducer.runTauIdMVA import *
na = TauIDEmbedder(process, cms, # pass tour process object
    debug=True,
    toKeep = ["2017v1", "2017v2", "dR0p32017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)
na.runTauID()

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


#process.taus=cms.Sequence(process.rerunMvaIsolation2SeqRun2 + process.NewTauIDsEmbedded + process.bareTaus + process.softTaus)
process.taus=cms.Sequence(process.rerunMvaIsolationSequence + process.NewTauIDsEmbedded + process.bareTaus + process.softTaus)


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


### ----------------------------------------------------------------------
### Jets
### ----------------------------------------------------------------------

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

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   jetCorrections = ('AK4PFchs', jecLevels, 'None'),
   btagDiscriminators = [
      	'pfDeepFlavourJetTags:probb',
      	'pfDeepFlavourJetTags:probbb',
      	'pfDeepFlavourJetTags:problepb',
      	'pfDeepFlavourJetTags:probc',
      	'pfDeepFlavourJetTags:probuds',
      	'pfDeepFlavourJetTags:probg'
      ],
   postfix='NewDFTraining'
)

#process.updatedPatJetsUpdatedJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
#process.updatedPatJetsUpdatedJEC.userData.userInts.src    += ['pileupJetIdUpdated:fullId']
#process.jecSequence = cms.Sequence(process.pileupJetIdUpdated + process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)
#process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

jetsNameAK4="selectedUpdatedPatJetsNewDFTraining"

#needed for ele/mu ptRatioRel with nanoAOD
process.tightJetId = cms.EDProducer("PatJetIDValueMapProducer",
                          filterParams=cms.PSet(
                            version = cms.string('WINTER17'),
                            quality = cms.string('TIGHT'),
                          ),
                          src = cms.InputTag("slimmedJets")
)

#needed for ele/mu ptRatioRel with nanoAOD
process.tightJetIdLepVeto = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
          version = cms.string('WINTER17'),
          quality = cms.string('TIGHTLEPVETO'),
        ),
                          src = cms.InputTag("slimmedJets")
)

#needed for ele/mu ptRatioRel with nanoAOD
process.slimmedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
     src = cms.InputTag("slimmedJets"),
     userFloats = cms.PSet(),
     userInts = cms.PSet(
        tightId = cms.InputTag("tightJetId"),
        tightIdLepVeto = cms.InputTag("tightJetIdLepVeto"),
     ),
)

#needed for ele/mu ptRatioRel with nanoAOD
from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
## Note: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
##      (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CMSSW_7_6_4_and_above )
process.jetCorrFactors = patJetCorrFactors.clone(src='slimmedJetsWithUserData',
    levels = cms.vstring('L1FastJet',
        'L2Relative',
        'L3Absolute',
        'L2L3Residual'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

#needed for ele/mu ptRatioRel with nanoAOD
from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *
process.updatedJets = updatedPatJets.clone(
  addBTagInfo=False,
  jetSource='slimmedJetsWithUserData',
  jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactors") ),
)

process.jets = cms.EDFilter("PATJetRefSelector",
                            #src = cms.InputTag("slimmedJets"),
                            src = cms.InputTag("updatedPatJetsNewDFTraining"),
                            cut = cms.string(JETCUT),
)


### ----------------------------------------------------------------------
### QG discriminator
### ----------------------------------------------------------------------

if COMPUTEQGVAR:

    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDBSetup,
      toGet = cms.VPSet(
        cms.PSet(
          record = cms.string('QGLikelihoodRcd'),
          tag    = cms.string('QGLikelihoodObject_v1_AK4PFchs_2017'),
          label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    )

    process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')

    process.load('RecoJets.JetProducers.QGTagger_cfi')
    process.QGTagger.srcJets          = cms.InputTag(jetsNameAK4)    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
    process.jetSequence = cms.Sequence(process.jets * process.QGTagger + process.tightJetId + process.tightJetIdLepVeto + process.slimmedJetsWithUserData + process.jetCorrFactors + process.updatedJets)

else:
    process.jetSequence = cms.Sequence(process.jets + process.tightJetId + process.tightJetIdLepVeto + process.slimmedJetsWithUserData + process.jetCorrFactors + process.updatedJets)


## ----------------------------------------------------------------------
## Build ll candidates (here OS)
## ----------------------------------------------------------------------

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
    
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(
      process,
      isData= (not IsMC),
    )
    # patch to get a standalone MET significance collection
    process.METSignificance = cms.EDProducer ("ExtractMETSignificance",
                                                  srcMET=cms.InputTag("slimmedMETs","","TEST")
                                                  )

    # add variables with MET shifted for TES corrections
    process.ShiftMETforTES = cms.EDProducer ("ShiftMETforTES",
                                             srcMET  = cms.InputTag("slimmedMETs","","TEST"),
                                             tauCollection = cms.InputTag("softTaus")
                                             )

    process.METSequence += process.fullPatMetSequence
    process.METSequence += process.METSignificance
    process.METSequence += process.ShiftMETforTES

    # EE noise mitigation 
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(
      process,
      isData= (not IsMC),
      fixEE2017 = True,
      fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
      postfix = "ModifiedMET"
    )
    process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)


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
            src = cms.InputTag(jetsNameAK4),                                      
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
                                      computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False),
                                      computeForUpDownMET = cms.bool(COMPUTEMETUPDOWNSVFIT if IsMC else False),
                                      METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                      METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                      METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                      METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN")
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
                                    srcCov     = cms.InputTag("METSignificance", "METCovariance"),
                                    METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                    METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                    METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                    METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN")
)



## ----------------------------------------------------------------------
## Ntuplizer
## ----------------------------------------------------------------------

if (YEAR==2016):
	rhocol = "fixedGridRhoFastjetCentralNeutral"

if (YEAR==2017 or YEAR==2018):
	rhocol = "fixedGridRhoFastjetAll"

process.HTauTauTree = cms.EDAnalyzer("HTauTauNtuplizer",
                      fileName = cms.untracked.string ("CosaACaso"),
                      applyFSR = cms.bool(APPLYFSR),
                      year = cms.int32(YEAR),
                      IsMC = cms.bool(IsMC),
                      doCPVariables = cms.bool(doCPVariables),               
                      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      secVtxCollection = cms.InputTag("slimmedSecondaryVertices"), # FRA
                      puCollection = cms.InputTag("slimmedAddPileupInfo"),
                      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      #rhoMiniRelIsoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoMiniRelIsoCollection = cms.InputTag(rhocol), #for non-nanoAOD
                      rhoForJER = cms.InputTag("fixedGridRhoAll"), # FRA
                      PFCandCollection = cms.InputTag("packedPFCandidates"),
                      jetCollection = cms.InputTag(jetsNameAK4),
                      #JECset = cms.untracked.string("patJetCorrFactors"),
                      JECset = cms.untracked.string("patJetCorrFactorsTransientCorrectedNewDFTraining"),
                      #JECset = cms.untracked.string("patJetCorrFactorsNewDFTraining"),
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
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", HLTProcessName), #Different names for MiniAODv2 at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.                      
                      #triggerSet = cms.InputTag("selectedPatTrigger"), # FRA
                      triggerSet = cms.InputTag("slimmedPatTrigger"),    # FRA
                      triggerList = HLTLIST,
                      metFilters = cms.InputTag ("TriggerResults","",METfiltersProcess),
                      PUPPImetCollection = cms.InputTag("slimmedMETsPuppi"),
                      srcPFMETCov = cms.InputTag("METSignificance", "METCovariance"),
                      srcPFMETSignificance = cms.InputTag("METSignificance", "METSignificance"),
                      #l1extraIsoTau = cms.InputTag("l1extraParticles", "IsoTau"),
                      HT = cms.InputTag("externalLHEProducer"),
                      beamSpot = cms.InputTag("offlineBeamSpot"),
                      nBadMu = cms.InputTag("removeBadAndCloneGlobalMuons"),
                      genLumiHeaderTag = cms.InputTag("generator"),
                      metERCollection = cms.InputTag("slimmedMETsModifiedMET"),
                      ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")
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


## ----------------------------------------------------------------------
## Paths
## ----------------------------------------------------------------------
process.PVfilter = cms.Path(process.goodPrimaryVertices)
process.ecalBadCalib = cms.Path(process.ecalBadCalibReducedMINIAODFilter)

# Prepare lepton collections
process.Candidates = cms.Sequence(
    process.egammaPostRecoSeq  +
    # process.printTree         + # just for debug, print MC particles
    process.nEventsTotal       +
    #process.hltFilter         + 
    process.nEventsPassTrigger +
    #process.egammaPostRecoSeq  +
    #process.jecSequence + 
    process.patJetCorrFactorsNewDFTraining+
    process.updatedPatJetsNewDFTraining+
    process.pfImpactParameterTagInfosNewDFTraining+
    process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining+
    process.pfDeepCSVTagInfosNewDFTraining+
    process.pfDeepFlavourTagInfosNewDFTraining+
    process.pfDeepFlavourJetTagsNewDFTraining+
    process.patJetCorrFactorsTransientCorrectedNewDFTraining+
    process.updatedPatJetsTransientCorrectedNewDFTraining+
    process.selectedUpdatedPatJetsNewDFTraining+
    #process.jecSequence+
    process.jetSequence + #process.jets + 
    process.muons              +
    #process.egammaPostRecoSeq  +
    process.electrons          + process.cleanSoftElectrons +
    #process.egammaPostRecoSeq  +
    process.taus               +
    process.fsrSequence        +
    process.softLeptons        + process.barellCand +
    #process.jets              +
    #process.jecSequence + process.jetSequence + #process.jets + 
    process.METSequence        +
    process.geninfo            +
    process.SVFit
    )

# always run ntuplizer
process.trees = cms.EndPath(process.HTauTauTree)# + process.HTauTauTreeTauUp + process.HTauTauTreeTauDown)

