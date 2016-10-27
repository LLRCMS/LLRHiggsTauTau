#!/usr/bin/env python

import os,sys
import optparse
import fileinput
import commands
import time
import subprocess


def trimInputFileList (originalFileList) :
    lines = [line[:-1] for line in fileinput.input(originalFileList) if not '#' in line and '.root' in line]
    return lines        


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def makeInputFileList (filesList, jobsDir, index) :
    fullfilename = jobsDir + '/inputFiles_' + str (index) + '.py'
    localFilesList = open (fullfilename, 'w')
    localFilesList.write ('FILELIST = cms.untracked.vstring()\n')
    localFilesList.write ('FILELIST.extend ([\n')
    for i in range (len (filesList)) :
        localFilesList.write (filesList[i] + '\n')
    localFilesList.write ('])\n')
    return fullfilename


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == "__main__":

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-q', '--queue'         ,    dest='queue'              , help='batch queue'                                        , default='short')
    parser.add_option('-s', '--sample'        ,    dest='samplefile'         , help='python file names'                                  , default='TTJets_files.py')
    parser.add_option('-n', '--njobs'         ,    dest='njobs'              , help='number of jobs'                                     , default=1,  type=int)
    parser.add_option('-w', '--wholefile'     ,    dest='wholefile'          , help='use whole files [false]'                            , default=False)
    parser.add_option('-m', '--maxev'         ,    dest='maxevents'          , help='max total number of events (-1 for all)'            , default=-1, type=int)
    parser.add_option('-o', '--out'           ,    dest='output'             , help='output directory'                                   , default='/data_CMS/cms/cadamuro/test_submit_to_tier3/')
    parser.add_option('-c', '--cfg'           ,    dest='cfg'                , help='cfg file'                                           , default='analyzer_forTier3.py')
    parser.add_option('-r', '--rep'           ,    dest='customReplacements' , help='sed replacements for cfg  key1:val1,key2:val2,...'  , default=None)
    parser.add_option('-t', '--tag'           ,    dest='nameTag'            , help='tag defining the sample and production'             , default='')
    parser.add_option('-z', '--offset'        ,    dest='offset'             , help='overall offset (events to be skipped)'              , default=0, type=int)
    parser.add_option('-i', '--ismc'          ,    dest='isMC'               , help='is an MC sample'                                    , default=True)
    (opt, args) = parser.parse_args()

    if not opt.wholefile: 
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
    else:
        # prepare a list of lists of files to run on subsets of the total files list
        fileslist = trimInputFileList (opt.samplefile)
        if opt.njobs > len (fileslist) : opt.njobs = len (fileslist)
        nfiles = (len (fileslist) + len (fileslist) % opt.njobs) / opt.njobs
        chunks = [fileslist[x:x+nfiles] for x in xrange (0, len (fileslist), nfiles)]

    #prepare output
    timeTag = time.time()
    cmsswBase=os.environ['CMSSW_BASE']
    relativeJobDir = 'JobLauncher_%s_%dEvents_%dSkipped_%s'%(opt.nameTag, opt.maxevents, opt.offset, timeTag)
    jobsDir=cmsswBase+'/src/LLRHiggsTauTau/NtupleProducer/test/' + relativeJobDir
    # NB samplefile should stay in test and the script has to be run from test
    os.system('mkdir -p %s'%jobsDir)
    # print ('cp ' + opt.samplefile + ' ' + jobsDir)
    os.system ('cp ' + opt.samplefile + ' ' + jobsDir)
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
    # usually I call it externally
    # proc = subprocess.Popen ('voms-proxy-info', stdout=subprocess.PIPE)
    # tmp = [word for word in proc.stdout.read ().split ('\n') if 'timeleft' in word]
    # if len (tmp) == 0 or int (tmp[0].split (':')[1]) < 24 : # hours
    #     os.system ('source /opt/exp_soft/cms/t3/t3setup')

    #loop over the required number of jobs
    nloop = len (chunks) if opt.wholefile else opt.njobs
    for n in xrange (0, nloop) :

        #sed the cfg template 
        inCfg = open(opt.cfg).read()
        outCfg = open('%s/cmssw_%d_cfg.py'%(jobsDir,n), 'w')
        outFullLFN = outFullPath
        outFullLFNNtuples  = outFullLFN + 'HTauTauAnalysis_%d.root' % n
        outFullLFNEnriched = outFullLFN + 'Enriched_miniAOD_%d.root' % n
    
        if not opt.wholefile: 
            replacements = {
                            'XXX_MAXEVENTS_XXX':str(EvPerJob),
                            'XXX_SKIPEVENTS_XXX':str(n*EvPerJob + opt.offset),
                            'XXX_SAMPLEFILENAME_XXX':opt.samplefile,
                            'XXX_OUTPUTFILENTUPLE_XXX':outFullLFNNtuples,
                            'XXX_OUTPUTFILEENRICHED_XXX':outFullLFNEnriched,
                            'XXX_ISMC_XXX':str(opt.isMC)
                           }
        else:
            # scrivi un file name con dentro la lista giusta
            chunkfilename = makeInputFileList (chunks[n], relativeJobDir, n)
            replacements = {
                            'XXX_MAXEVENTS_XXX':'-1',
                            'XXX_SKIPEVENTS_XXX':str (opt.offset),
                            'XXX_SAMPLEFILENAME_XXX':chunkfilename,
                            'XXX_OUTPUTFILENTUPLE_XXX':outFullLFNNtuples,
                            'XXX_OUTPUTFILEENRICHED_XXX':outFullLFNEnriched,
                            'XXX_ISMC_XXX':str(opt.isMC)
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
            os.system("/opt/exp_soft/cms/t3/t3submit_new -%s \'%s/runJob_%d.sh\'"%(opt.queue,jobsDir,n))
    
