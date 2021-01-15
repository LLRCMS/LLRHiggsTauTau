### Description

This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analyses.

The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

### Instructions for older releases:
<details>
### Instructions for 7_2_X

```
cmsrel CMSSW_7_2_3_patch1
cd CMSSW_7_2_3_patch1/src
cmsenv
git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg FWCore/Version
git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720 
git cms-merge-topic ikrav:egm_id_phys14
git-cms-merge-topic -u cms-met:72X-MetSig-150311
git-cms-merge-topic -u cms-met:72X-mvaMETForMiniAOD
cd RecoMET/METPUSubtraction/ ; git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 72X-13TeV-Phys14_25_V4-26Mar15 
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout for72X ; cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h ; cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h ; cd -
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons ; cd -
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone

THEN EDIT THE RecoMET/METPUSubtraction/python/mvaPFMET_cff.py at LINE 75 (could change)
```

### Instructions for 7_4_7 (miniAOD_v1)

```
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src
cmsenv
git cms-merge-topic ikrav:egm_id_747_v2
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout master
cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h
cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons
cd -
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone
scram b -j 4
```

### Instructions for 7_4_12 (miniAODv2)

```
cmsrel CMSSW_7_4_12
cd CMSSW_7_4_12/src
cmsenv
git cms-merge-topic ikrav:egm_id_7.4.12_v1
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout master
cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h
cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons
cd -
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone
scram b -j 4
```

### Instructions for 7_6_3 (miniAODv2)

```
cmsrel CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src
cmsenv
# MVA MET
git cms-merge-topic --unsafe l-cadamuro:MVAMETExtCombinatorics
# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout master
cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h
cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -
# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons
cd -
# updated pileup jet ID (see JetMET twiki https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Information_for_13_TeV_data_anal)
git cms-merge-topic --unsafe jbrands:pileupJetId76X
cd RecoJets/JetProducers/data/
wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta0to2p5_BDT.weights.xml.gz
wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta2p5to2p75_BDT.weights.xml.gz
wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta2p75to3_BDT.weights.xml.gz
wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta3to5_BDT.weights.xml.gz
cd -
# SVfit
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd -
scram b -j 4
```

### Instructions for 8_0_6 (miniAOD 2016)

```
cmsrel CMSSW_8_0_6
cd CMSSW_8_0_6/src
cmsenv
# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout master
cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h
cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -
# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons
cd -
# SVfit
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd -
scram b -j 4
```

### Legacy Instructions (2016 data)
```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src/
cmsenv

# MET Recipe
git cms-merge-topic cms-met:METRecipe_8020 -u
git cms-merge-topic cms-met:METRecipe_80X_part2 -u

#ReReco muons fix
git cms-merge-topic gpetruc:badMuonFilters_80X_v2

# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections

git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau
git checkout master
cd -

git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools
git checkout master -- interface/MuonEffectiveArea.h
cd -

git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools
git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -

# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2
git checkout master FSRPhotons
cd -

# SVfit
git clone git@github.com:veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout d115239192d3eb7531e213767e02ef4777b3fbfe

cd $CMSSW_BASE/src
scram b -j 8

```


### Instructions for 8_0_25

```
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src
cmsenv
# MET Recipe for ICHEP dataset
git cms-merge-topic cms-met:METRecipe_8020
# Spring-16 Electron MVA ID
git cms-merge-topic ikrav:egm_id_80X_v2
#ReReco muons fix
git cms-merge-topic gpetruc:badMuonFilters_80X_v2
# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau; git checkout master
cd -
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h
cd -
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -
# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons
cd -
# SVfit
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd $CMSSW_BASE/src
scram b -j 4
cd $CMSSW_BASE/external/$SCRAM_ARCH
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout egm_id_80X_v1
cd $CMSSW_BASE/src
scram b -j 4
```
### Instructions for 92X

```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src/
cmsenv

# MVA EleID Fall 2017
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
scram b -j 8
cd $CMSSW_BASE/external
# below, you may have a different architecture, this is just one example from lxplus (same on polui)
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/src

# Remove some of the unused weights (otherwise the crab tarball is too big for submission)
cd $CMSSW_BASE/external/slc6_amd64_gcc630/data/RecoEgamma/ElectronIdentification/data
rm -r PHYS14 Spring15 
cd $CMSSW_BASE/external/slc6_amd64_gcc630/data/RecoEgamma/PhotonIdentification/data
rm -r PHYS14 Spring15 Spring16
cd $CMSSW_BASE/src

# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections

# LLRHiggsTauTau framework
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau
git checkout 92X
cd -

git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools
git checkout master -- interface/MuonEffectiveArea.h
cd -

git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools
git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -

# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2
git checkout master FSRPhotons
# need to fix: - FSRPhotons/plugins/FSRPhotonProducer.cc
#              - FSRPhotons/plugins/PhotonPFIsoCalculator.cc
# replace 'std::auto_ptr' with 'std::unique_ptr' 
# search for 'iEvent.put( XXXX );' and replace with 'iEvent.put( std::move(XXXX) );'
cd -

# SVfit
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF
git clone git@github.com:veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout HIG-16-006
# need to fix: TauAnalysis/SVfitStandalone/src/SVfitStandaloneQuantities.cc 
# add '#include <numeric>'

cd $CMSSW_BASE/src
scram b -j 8
```


### Instructions for 94X

```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv

# MVA EleID Fall 2017
git cms-init
git cms-merge-topic cms-egamma:EgammaID_949

# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections

# LLRHiggsTauTau framework
git clone https://github.com/LLRCMS/LLRHiggsTauTau
cd LLRHiggsTauTau
git checkout 94X_HH
cd -

git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools
git checkout master -- interface/MuonEffectiveArea.h
cd -

git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools
git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -

# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2
git checkout master FSRPhotons
# need to fix: - FSRPhotons/plugins/FSRPhotonProducer.cc
#              - FSRPhotons/plugins/PhotonPFIsoCalculator.cc
# replace 'std::auto_ptr' with 'std::unique_ptr'
# search for 'iEvent.put( XXXX );' and replace with 'iEvent.put( std::move(XXXX) );'
cd -

# bad MET filter fix
git cms-addpkg RecoMET/METFilters


# MET - EE noise mitigation
git cms-merge-topic cms-met:METFixEE2017_949_v2

# SVfit
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF
git clone git@github.com:veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout HIG-16-006
# need to fix: TauAnalysis/SVfitStandalone/src/SVfitStandaloneQuantities.cc 
# add '#include <numeric>'

cd $CMSSW_BASE/src
scram b -j 8
```

</details>


### Instructions for 102X

```
## IMPORTANT ##
# SLC6 is not supported anymore on lxplus machines, SLC7 must be used
# Make sure your architecture is slc7_amd64_gcc700
# (you can set the architecture with: export SCRAM_ARCH=slc7_amd64_gcc700 )

cmsrel CMSSW_10_2_23
cd CMSSW_10_2_23/src
cmsenv

git cms-init

# MVA EleID Fall 2018
git cms-merge-topic cms-egamma:EgammaPostRecoTools  #if you want the V2 IDs, otherwise skip

# PU jet ID
git cms-addpkg RecoJets/JetProducers
git clone -b 94X_weights_DYJets_inc_v2 git@github.com:cms-jet/PUjetID.git PUJetIDweights/
cp PUJetIDweights/weights/pileupJetId_{94,102}X_Eta* $CMSSW_BASE/src/RecoJets/JetProducers/data/
rm -rf PUJetIDweights/
git cms-merge-topic alefisico:PUID_102X

# Z-recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections

# LLRHiggsTauTau framework
git clone git@github.com:LLRCMS/LLRHiggsTauTau.git
cd LLRHiggsTauTau
git checkout 102X_HH
cd -

git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools
git checkout master -- interface/MuonEffectiveArea.h
cd -

git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools
git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -

# FSR corrections
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2
git checkout origin/102X FSRPhotons
cd -

# bad MET filter fix
git cms-addpkg RecoMET/METFilters

# SVfit
git clone https://github.com/LLRCMS/ClassicSVfit.git TauAnalysis/ClassicSVfit -b bbtautau_LegacyRun2
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF

#Add TauPOG corrections (TES and EES)
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs

#Add DeepTau code from Tau POG repository (note "-u" option preventing checkout of unnecessary stuff)
git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_DeepTau2017v2p1_nanoAOD

cd $CMSSW_BASE/src
scram b -j 8
```

### Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py

