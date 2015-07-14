#!/bin/bash

RECOM_SHARLIB=true

# read from command line
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -n|--notrecompile)
    RECOM_SHARLIB=false
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

if [ $RECOM_SHARLIB = true ] ;  then
   echo "removing old files"
   rm -f libHHKinFit.so run

   echo "creating shared library"
   g++ -fPIC -shared HHKinFit/src/*.cpp `root-config --cflags --glibs` -D HHKINFIT_STANDALONE -I ./include -o libHHKinFit.so
fi

echo "creating executable"
g++ ComputeHHKinFit.C `root-config --cflags --glibs` -D HHKINFIT_STANDALONE -I ./include -L . -lHHKinFit  -o ComputeHHKinFit
