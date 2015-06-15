import os
import glob
import ROOT
files = glob.glob('./output_*.root')
if not os.path.isdir("./zombies"): os.system ('mkdir zombies')
if not os.path.isdir("./recovered"): os.system ('mkdir recovered')

nTot = 0
nZombies = 0
nRec = 0

for fname in files:
    nTot += 1
    ff = ROOT.TFile (fname)
    isZombie = ff.IsZombie()
    isRecovered = ff.TestBit(ROOT.TFile.kRecovered)
    if isZombie:
        os.system ('mv %s zombies' % fname)
        nZombies += 1
    if isRecovered:
        os.system ('mv %s recovered' % fname)
        nRec += 1
print "Integrity checked, %i zombies, %i recovered" % (nZombies, nRec)
