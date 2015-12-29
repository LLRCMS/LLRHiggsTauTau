#!/bin/bash
# Use > 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).

# create file list for CMSSW compatible with ntuplizer for tier3

#usage: source MakeFileList.sh -t "HH 74X lambda-4 gen info fix " -o HH_Lambdam4.txt -p /data_CMS/cms/salerno/Lambdam4_74x_step2/miniAOD_lambdam4_2_300000Events_0Skipped_1435089918.88

PRINTHEADER=true
PRINTCMSSWFORMAT=true # lines begin with 'file: ...
PRINTMESSAGE=true
ONLYNTUPLES=false #select only ntuples and not EMiniAOD

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -p|--path)
    FILEPATH="$2"
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
    -c|--cmsswformatoff)
    PRINTCMSSWFORMAT=false
    ;;
    -m|--messageoff)
    PRINTMESSAGE=false
    ;;
    -g|--getonlyntuples)
    ONLYNTUPLES=true
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

## create file
if [ -f $OUTFILE ] ; then
    rm $OUTFILE
fi
touch $OUTFILE

if $PRINTMESSAGE = true ] ; then
    echo "#File list is: " $MESSAGE > $OUTFILE
fi
#echo "\n" >> $OUTFILE

if [ $PRINTHEADER = true ] ; then
    echo "FILELIST = cms.untracked.vstring()" >> $OUTFILE
    echo "FILELIST.extend ([" >> $OUTFILE
fi

# file list formatted using awk
#ls $FILEPATH/*.root | awk -v pathto=$FILEPATH '{print "'\''file:" pathto "/"$1 "'\''," }' >> $OUTFILE
if [ $PRINTCMSSWFORMAT = true ] ; then
    ls $FILEPATH/*.root | awk '{print "'\''file:"$1 "'\''," }' >> $OUTFILE
else
    if [ $ONLYNTUPLES = true ] ; then
       ls $FILEPATH/*.root | grep HTauTauAnalysis | awk '{print $1}' >> $OUTFILE
    else
       ls $FILEPATH/*.root | awk '{print $1}' >> $OUTFILE   
    fi
fi

if [ $PRINTHEADER = true ] ; then
    echo "])" >> $OUTFILE
fi
