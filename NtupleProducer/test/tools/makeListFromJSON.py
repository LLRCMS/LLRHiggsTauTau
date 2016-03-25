import os

## call several time the das_cllient script to retrieve files satisfying requirements in a certain lumi+run range
## useful to select subset of files to run locally what failed on CRAB with a lumi mask

dataset = "/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD"

# copy here JSON from missingLumiSection
JSON = {"254231": [[1, 24]],
 "254232": [[1, 81]],
 "254790": [[90, 90], [93, 630], [633, 697], [701, 715], [719, 784]],
 "254852": [[47, 94]],
 "254879": [[52, 52], [54, 140]],
 "254906": [[1, 75]],
 "254907": [[1, 52]],
 "254914": [[32, 32], [34, 78]],
 "256676": [[42, 42]],
 "257490": [[106, 106]],
 "257613": [[341, 341]],
 "257816": [[120, 120]],
 "257822": [[494, 497], [1093, 1093]],
 "258158": [[147, 147], [1195, 1195], [1199, 1199]],
 "258214": [[170, 170]],
 "258448": [[24, 24]],
 "258702": [[225, 225], [402, 402]],
 "259626": [[109, 109]],
 "259683": [[25, 25]],
 "259686": [[45, 45], [49, 49], [81, 81], [100, 100], [163, 163]],
 "259822": [[132, 132]],
 "260424": [[492, 492]],
 "260627": [[511, 511]]}


for run in JSON:
    for lumirange in JSON[run]:
        command = "python das_client.py --query=\"file dataset=%s run=%s" % (dataset, run)
        low = lumirange[0]
        up = lumirange[1]
        if (low == up):
            command += " lumi=%i" % low
        else:
            command += " lumi between [%i,%i]" % (low, up)
        command += "\" --limit=0"
        print command
        print ""
        os.system (command)
