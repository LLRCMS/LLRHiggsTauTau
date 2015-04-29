This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analysis.

The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

Instruction to get and run the package:

git cms-addpkg PhysicsTools/PatAlgos

git cms-addpkg FWCore/Version

git-cms-merge-topic -u cms-met:72X-MetSig-150311

git-cms-merge-topic -u cms-met:72X-mvaMETForMiniAOD
( cd RecoMET/METPUSubtraction/ ; git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 72X-13TeV-Phys14_25_V4-26Mar15 )

git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720 
git cms-merge-topic ikrav:egm_id_phys14

git clone https://github.com/LLRCMS/LLRHiggsTauTau
(cd LLRHiggsTauTau; git checkout master)

git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
(cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h )

git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
(cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h) 

git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
(cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons)

git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone

THEN EDIT THE RecoMET/METPUSubtraction/python/mvaPFMET_cff.py at LINE 75 (could change)

Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py

