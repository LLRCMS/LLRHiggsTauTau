#!/usr/bin/env python

#TO DO: source cmssw grid env

import os,sys
import optparse
import commands
import time


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                                        , default='cms')
parser.add_option('-s', '--sample'     ,    dest='samplefile'         , help='python file names'                                  , default='TTJets_files.py')
parser.add_option('-n', '--njobs'      ,    dest='njobs'              , help='number of jobs'                                     , default=1,  type=int)
parser.add_option('-m', '--maxev'      ,    dest='maxevents'          , help='max total number of events (-1 for all)'            , default=-1, type=int)
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'                                   , default='/data_CMS/cms/cadamuro/test_submit_to_tier3/')
parser.add_option('-c', '--cfg'        ,    dest='cfg'                , help='cfg file'                                           , default='analyzer_forTier3.py')
parser.add_option('-r', '--rep'        ,    dest='customReplacements' , help='sed replacements for cfg  key1:val1,key2:val2,...'  , default=None)
parser.add_option('-t', '--tag'        ,    dest='nameTag'            , help='tag defining the sample and production'             , default='')
parser.add_option('-z', '--offset'     ,    dest='offset'             , help='overall offset (events to be skipped)'              , default=0, type=int)
(opt, args) = parser.parse_args()


# compute nEvts per job and check if something missing
# WHAT TO DO WITH -1?? --> I cannot know the tot number of data
if opt.maxevents == -1:
    sys.exit("maxevents = -1 feature not implemented yet")
EvPerJob = opt.maxevents / opt.njobs
print "Events per job: %d" % EvPerJob
wasted = opt.maxevents - (EvPerJob * opt.njobs)
if wasted != 0:
    print "NOTE: %d events will NOT be analyzed (job number is not an integer divisor)" % wasted

if opt.offset != 0:
    print "Skipping the first %d events" % opt.offset

#prepare output
timeTag = time.time()
cmsswBase=os.environ['CMSSW_BASE']
jobsDir=cmsswBase+'/src/LLRHiggsTauTau/NtupleProducer/test/JobLauncher_%s_%dEvents_%dSkipped_%s'%(opt.nameTag, opt.maxevents, opt.offset, timeTag)
os.system('mkdir -p %s'%jobsDir)
print '** INFO ** Jobs folder is: %s'%(jobsDir)

# create output folder
outFullPath = opt.output
if not outFullPath.endswith('/'):
    outFullPath += '/'
outFullPath += 'HiggsTauTauOutput_%s_%dEvents_%dSkipped_%s/'%(opt.nameTag, opt.maxevents, opt.offset, timeTag)
os.system('mkdir -p %s' % outFullPath)
print "** INFO ** Saving output files in: %s" % outFullPath

def replfunc(match):
    return repldict[match.group(0)]

#initialize t3 for submission
os.system('source /opt/exp_soft/cms//t3/t3setup')

#loop over the required number of jobs
for n in xrange(0,opt.njobs):

    #sed the cfg template 
    inCfg = open(opt.cfg).read()
    outCfg = open('%s/cmssw_%d_cfg.py'%(jobsDir,n), 'w')
    outFullLFN = outFullPath
    outFullLFN += 'output_%d.root' % n
    
    replacements = {
                    'XXX_MAXEVENTS_XXX':str(EvPerJob),
                    'XXX_SKIPEVENTS_XXX':str(n*EvPerJob + opt.offset),
                    'XXX_SAMPLEFILENAME_XXX':opt.samplefile,
                    'XXX_OUTPUTFILE_XXX':outFullLFN
                   }
    if opt.customReplacements is not None:
        for rep in opt.customReplacements.split(','):
            repKeys=rep.split(':')
            replacements[repKeys[0]]=repKeys[1]
    
    for i in replacements.keys():
        inCfg = inCfg.replace(i, replacements[i])

    outCfg.write(inCfg)
    outCfg.close()
    
    
    #create a wrapper for standalone cmssw job
    scriptFile = open('%s/runJob_%d.sh'%(jobsDir,n), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('export X509_USER_PROXY=~/.t3/proxy.cert\n')
    scriptFile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    scriptFile.write('export SCRAM_ARCH=slc6_amd64_gcc472\n')
    scriptFile.write('cd {}\n'.format(cmsswBase + '/src/'))
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd %s\n'%jobsDir)
    scriptFile.write('cmsRun cmssw_%d_cfg.py &> log_%d_job.txt \n' % (n,n) )
    #scriptFile.write('cp SingleElectronPt35_PU0_GEN_SIM_%d.root %s\n'%(n,opt.output))
    #scriptFile.write('rm SingleElectronPt35_PU0_GEN_SIM_%d.root\n'%n)    
    scriptFile.write('echo "All done for job %d" \n'%n)
    scriptFile.close()
    os.system('chmod u+rwx %s/runJob_%d.sh'%(jobsDir,n))

    #submit it to the batch or run it locally
    if opt.queue=='':
        os.system('%s/runJob_%d.sh'%(jobsDir.jobSeed))
    else:
        os.system("/opt/exp_soft/cms/t3/t3submit -k eo -q cms \'%s/runJob_%d.sh\'"%(jobsDir,n))
    
