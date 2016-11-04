#!/usr/bin/python

#json = eval( open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt').read() ) #for Run2015B DCS only JSON
#json = eval( open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt').read() ) #upfate Run2015B
#json = eval( open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt').read() ) #1.56/fb
# json = eval( open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt').read() ) #2.6/fb in 2016 
json = eval( open('Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt').read())

print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
 for ls in lumisections:
   print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'
