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
    #process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'                 # 2017 MC
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'                 # 12Apr2018 2017 MC
else :
    # process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0' # ICHEP            # FIXME !!!!!!!!
    #process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' # Run B-G               # FIXME !!!!!!!!
    # process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v16' # Run H                      # FIXME !!!!!!!!
    process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6'                    # 2017 Data
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

#2017 ECAL bad calibration filter to be rerun, fix from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
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



process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    #src = cms.InputTag("removeBadAndCloneGlobalMuons"), # input for 2016 data only
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string(MUCUT),
#    Lowering pT cuts
#    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) &&" +
#                     "pt>3 && p>3.5 && abs(eta)<2.4")
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

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("slimmedElectrons"),
   rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
   vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   genCollection = cms.InputTag("prunedGenParticles"),
   sampleType = cms.int32(LEPTON_SETUP),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string(ELECUT)
   )

process.electrons = cms.Sequence(process.softElectrons)


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
   # --> Correct values for 2017 data - Uncertainty 3.0%
   # NominalTESUncertainty      = cms.double(3.0) , # in percent, up/down uncertainty of TES
   # NominalTESCorrection1Pr    = cms.double(-3.0), #DecayMode==0
   # NominalTESCorrection1PrPi0 = cms.double(-2.0), #DecayMode==1
   # NominalTESCorrection3Pr    = cms.double(-1.0), #DecayMode==10

   # --> Correct values for 2017 data - Chiara update Jan2019
   NominalTESUncertainty1Pr         = cms.double(0.8) , # in percent, up/down uncertainty of TES      
   NominalTESUncertainty1PrPi0      = cms.double(0.8) , # in percent, up/down uncertainty of TES      
   NominalTESUncertainty3Pr         = cms.double(0.9) , # in percent, up/down uncertainty of TES      
   NominalTESUncertainty3PrPi0      = cms.double(1.0) , # in percent, up/down uncertainty of TES      
   NominalTESCorrection1Pr          = cms.double(0.7), #DecayMode==0                                 
   NominalTESCorrection1PrPi0       = cms.double(-0.2), #DecayMode==1                                 
   NominalTESCorrection3Pr          = cms.double(0.1), #DecayMode==10                                
   NominalTESCorrection3PrPi0       = cms.double(-0.1), #DecayMode==11                                

   ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
   # ApplyTESUpDown = cms.bool(True if IsMC else False), # no shift computation when data
   flags = cms.PSet(
        isGood = cms.string("")
        )
   )


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
    qgDatabaseVersion = 'v2b' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

    from CondCore.CondDB.CondDB_cfi import *
    QGPoolDBESSource = cms.ESSource("PoolDBESSource",
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

#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#process.testCands = cms.EDFilter("CandPtrSelector",
#                                 src=cms.InputTag("packedPFCandidates"),
#                                 cut=cms.string("abs(eta)<2.5")
#)
#
#runMetCorAndUncFromMiniAOD(process,
#                           isData=(not IsMC),
#                           pfCandColl=cms.InputTag("testCands"),
#                           reclusterJets=True,
#                           recoMetFromPFCs=True,
#                           postfix="Test",
#)
##metERTag = cms.InputTag("slimmedMETsTest","","TEST")
#
#### somehow MET recorrection gets this lost again...
#process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
#process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']
#
#
#
#process.MET = cms.Path(process.testCands + process.fullPatMetSequenceTest)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD (
        process,
        isData = (not IsMC),
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = "ModifiedMET"
)

# if running in schedule mode add this to your path
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
process.HTauTauTree = cms.EDAnalyzer("HTauTauNtuplizer",
                      fileName = cms.untracked.string ("CosaACaso"),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      doCPVariables = cms.bool(doCPVariables),               
                      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      secVtxCollection = cms.InputTag("slimmedSecondaryVertices"), # FRA
                      puCollection = cms.InputTag("slimmedAddPileupInfo"),
                      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoMiniRelIsoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoForJER = cms.InputTag("fixedGridRhoAll"), # FRA
                      PFCandCollection = cms.InputTag("packedPFCandidates"),
                      jetCollection = cms.InputTag("jets"),
                      JECset = cms.untracked.string("patJetCorrFactorsUpdatedJEC"),
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
                      srcPFMETCov = cms.InputTag("METSignificance", "METCovariance"),
                      srcPFMETSignificance = cms.InputTag("METSignificance", "METSignificance"),
                      HT = cms.InputTag("externalLHEProducer"),
                      beamSpot = cms.InputTag("offlineBeamSpot"),
                      genLumiHeaderTag = cms.InputTag("generator"),
                      #metERCollection = cms.InputTag("slimmedMETsTest","","TEST"),
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

##
## Paths
##
process.PVfilter = cms.Path(process.goodPrimaryVertices)
process.ecalBadCalib = cms.Path(process.ecalBadCalibReducedMINIAODFilter)

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

