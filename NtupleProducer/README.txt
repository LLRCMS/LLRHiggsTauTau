This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analysis.

The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

Instruction to get and run the package (under CMSSW_7_2_3_patch1):

git cms-addpkg PhysicsTools/PatAlgos

git cms-addpkg FWCore/Version

git-cms-merge-topic -u cms-met:72X-MetSig-150311
git-cms-merge-topic -u cms-met:72X-mvaMETForMiniAOD
#revert to older MVAMET commit, to be fixed later
(cd RecoMET/METPUSubtraction/ ; git reset --hard 177aa533286cb567678e0c7fde41757055377e72)

git clone https://github.com/LLRCMS/LLRHiggsTauTau.git
(cd LLRHiggsTauTau; git checkout master)

#effective areas (to be updated)
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
(cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h )

#Note we did not eventually adopt the latest EA update, we kept V00-00-30-00 of UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h .  c0db796 corresponds to it.
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
(cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h) 

#MuScleFit: probably tbf
#git clone https://github.com/scasasso/usercode MuScleFit

git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
(cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons)

# SVfit standalone algorithm
git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone

Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py

