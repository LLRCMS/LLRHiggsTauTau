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
### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if IsMC:
    #process.GlobalTag.globaltag = 'PLS170_V6AN1::All'#'GR_70_V2_AN1::All'   #MC in 70X, cf https://twiki.cern.ch/twiki/bin/view/CMS/MiniAOD
    process.GlobalTag.globaltag = 'PHYS14_25_V1::All' #MC in PHYS14
    #process.GlobalTag.globaltag = 'PHYS14_ST_V1::All'
else :
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


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

### ----------------------------------------------------------------------
### Trigger bit Requests 
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi 
# !!!!!!!!!!!!

#process.L1GtStableParametersRcdSource = cms.ESSource("EmptyESSource",
#    iovIsRunNotTime = cms.bool(True),
#    recordName = cms.string('L1GtStableParametersRcd'),
#    firstValid = cms.vuint32(1)
#)

process.hltFilterDiMu = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiMu.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiMu.HLTPaths = TRIGGERLIST
#process.hltCsc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
#process.hltCsc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")
process.hltFilterDiMu.throw = cms.bool(False)
# !!!!!!!!!!!!

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
  #src = cms.InputTag("offlinePrimaryVertices"),
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string(PVERTEXCUT),
  filter = cms.bool(True),
)

##
## Taus
##
process.bareTaus = cms.EDFilter("PATTauRefSelector",
   src = cms.InputTag("slimmedTaus"),
   cut = cms.string(TAUCUT)
   )

process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   cut = cms.string(TAUCUT),
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


process.taus=cms.Sequence(process.bareTaus + process.softTaus)

### ----------------------------------------------------------------------
### b quarks, only from MC
### ----------------------------------------------------------------------
process.bQuarks = cms.EDProducer("bFiller",
         src = cms.InputTag("prunedGenParticles"),
         cut = cms.string(BCUT),
         flags = cms.PSet(
            isGood = cms.string("")
        )
 )                
if IsMC : process.bquarks = cms.Sequence(process.bQuarks)
else : process.bquarks = cms.Sequence()

#
#Jets
#
process.jets = cms.EDFilter("PATJetRefSelector",
                           src = cms.InputTag("slimmedJets"),
                           cut = cms.string(JETCUT)
)
##
## Build ll candidates (here OS)
##
process.barellCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("softTaus@+ softTaus@-"),
                                    cut = cms.string(LLCUT),
                                    checkCharge = cms.bool(True)
)




## ----------------------------------------------------------------------
## MVA MET
## ----------------------------------------------------------------------

# plugin initialization

process.load("RecoJets.JetProducers.ak4PFJets_cfi")
process.ak4PFJets.src = cms.InputTag("packedPFCandidates")

from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

# output collection is defined here, its name is pfMVAMEt
# access this through pfMVAMEtSequence

process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
#process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons") # srcLeptons contains all the hard scatter products, is set later
process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfMVAMEt.minNumLeptons = cms.int32(2) # this is important to skip void collections in the loop on pairs

process.puJetIdForPFMVAMEt.jec = cms.string('AK4PF')
#process.puJetIdForPFMVAMEt.jets = cms.InputTag("ak4PFJets")
process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")



# python trick: loop on all pairs for pair MET computation

if USEPAIRMET:
   print "Using pair MET"

   # template of unpacker
   UnpackerTemplate = cms.EDProducer ("PairUnpacker",
                                   src = cms.InputTag("barellCand"))

   process.METSequence = cms.Sequence(process.ak4PFJets + process.calibratedAK4PFJetsForPFMVAMEt + process.puJetIdForPFMVAMEt)

   MVAPairMET = ()
   for index in range(100):
      UnpackerName = "PairUnpacker%i" % index
      UnpackerModule = UnpackerTemplate.clone( pairIndex = cms.int32(index) )
      setattr(process, UnpackerName, UnpackerModule)   #equiv to process.<UnpackerName> = <UnpackerModule>
      process.METSequence += UnpackerModule

      MVAMETName = "pfMETMVA%i" % index
      MVAModule = process.pfMVAMEt.clone( srcLeptons = cms.VInputTag (cms.InputTag(UnpackerName) ) )
      setattr(process, MVAMETName, MVAModule)
      process.METSequence += MVAModule
    
      MVAPairMET += (cms.InputTag(MVAMETName),)
    
else:
   print "Using event MET (same MET for all pairs)"
   process.pfMVAMEt.minNumLeptons = cms.int32(0) # ONLY FOR DEBUG PURPOSE, OTHERWISE DOES NOT COMPUTE MET AND SVFIT CRASHES DUE TO A SINGULAR MATRIX
   process.METSequence = cms.Sequence(
       process.ak4PFJets         +
       process.pfMVAMEtSequence
   )



## ----------------------------------------------------------------------
## SV fit
## ----------------------------------------------------------------------
process.SVllCand = cms.EDProducer("SVfitInterface",
                                  srcPairs   = cms.InputTag("barellCand"),
                                  #srcMET     = cms.InputTag("pfMVAMEt"),
                                  #srcMET     = cms.InputTag("slimmedMETs"),
                                  usePairMET = cms.untracked.bool(USEPAIRMET),
								  useMVAMET  = cms.untracked.bool(True),
								  triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT")								
)

if USEPAIRMET:
   process.SVllCand.srcMET    = cms.VInputTag(MVAPairMET)
else:
   process.SVllCand.srcMET    = cms.VInputTag("pfMVAMEt")



#Global configuration
TreeSetup = cms.EDAnalyzer("HHTauTauNtuplizer",
                      CandCollection = cms.untracked.string("SVllCand"),
                      fileName = cms.untracked.string ("CosaACaso"),
                      skipEmptyEvents = cms.bool(True),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
                      )

process.HHTauTauTree = TreeSetup.clone()

##
## Paths
##
process.PVfilter = cms.Path(process.goodPrimaryVertices)
#process.triggerDiMu = cms.Path(process.hltFilterDiMu)
#SkimPaths = cms.vstring('triggerDiMu') #Do not apply skim 

# Prepare lepton collections
process.Candidates = cms.Sequence(
    #process.hltFilterDiMu     + 
    #process.muons             +
    #process.electrons         + process.cleanSoftElectrons +
    process.taus              +
    #process.fsrSequence       +
    #process.softLeptons       + 
    process.barellCand        +
    process.jets              +
    process.METSequence       +
    process.bquarks           +
    process.SVllCand          + 
    process.HHTauTauTree
    )
