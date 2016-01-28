### Description

This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analyses.

The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

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
cmsrel CMSSW_7_6_3
cd CMSSW_7_6_3/src
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



### Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py
Please note that since 7_4_7 we switched to a new eleID recipe and we are not using anymore git cms-merge-topic sregnard:Phys14ElectronMvaIdFor745

