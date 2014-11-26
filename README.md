This package provides tools to READ miniAOD events and store the results in an TNtuple optimized for the HTauTau analysis.

The package workflow is the following:
Main flags/cut/configuration tools can be set in analyzer.py
Modules and plugins are called by python/HiggsTauTauProducer.py 
	Cuts on leptons/pairs are also set here
	This module creates AF/OS and AF/SS pairs, runs the SVfit on the pairs and store the useful variables in a TNtuple
The stored variables are set in plugins/HTauTauNtupleMaker.cc

Instruction to get and run the package:

git cms-addpkg PhysicsTools/PatAlgos

git-cms-merge-topic -u cms-met:72X-MVANoPUMet-20140718

git clone https://github.com/LLRCMS/LLRHiggsTauTau
(cd LLRHiggsTauTau; git checkout master)

git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
(cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h )

git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
(cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h) 

git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
(cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons)

git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone

Quick usage:
Define the files you want to run in analyzer.py and run cmsRun analyzer.py

