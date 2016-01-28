#!/bin/bash
# Use > 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).

# create file list for CMSSW compatible with ntuplizer for tier3

# input: use DATASET!
#usage: source MakeFileListDAS.sh -t "Data_MuMu" -o Data_MuMu.py -p /DoubleMuon/Run2015B-PromptReco-v1/MINIAOD

PRINTHEADER=true
PHYS03=false #instance of phys03 DBS, look here for published datasets

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -p|--path)
    DATASET="$2"
    shift # past argument
    ;;
    -o|--out)
    OUTFILE="$2" #output file name + path eventually
    shift # past argument
    ;;
    -t|--tag)
    MESSAGE="$2"
    shift # past argument
    ;;
    -n|--noheader)
    PRINTHEADER=false
    ;;
    -d|--dbsphys3)
    PHYS03=true
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

## create file
echo "#File list is: " $MESSAGE > $OUTFILE
#echo "\n" >> $OUTFILE

if [ $PRINTHEADER = true ] ; then
    echo "FILELIST = cms.untracked.vstring()" >> $OUTFILE
    echo "FILELIST.extend ([" >> $OUTFILE
fi

if [[ $DATASET = *"/" ]] ; then DATASET=${DATASET%?} ; fi

# file list formatted using awk
#ls $FILEPATH/*.root | awk -v pathto=$FILEPATH '{print "'\''file:" pathto "/"$1 "'\''," }' >> $OUTFILE
if [ $PHYS03 = true ] ; then
    python das_client.py --query="file dataset=$DATASET instance=prod/phys03" --limit=0 --verbose=1 | awk '{print "'\''"$1"'\'',"}' >> $OUTFILE
else
    python das_client.py --query="file dataset=$DATASET" --limit=0 --verbose=1 | awk '{print "'\''"$1"'\'',"}' >> $OUTFILE    
fi
#ls $FILEPATH/*.root | awk '{print "'\''file:"$1 "'\''," }' >> $OUTFILE

if [ $PRINTHEADER = true ] ; then
    echo "])" >> $OUTFILE
fi
