#!/bin/bash
echo "removing old files"
rm -f libHHKinFit.so run

echo "creating shared library"
g++ -fPIC -shared HHKinFit/src/*.cpp `root-config --cflags --glibs` -D HHKINFIT_STANDALONE -I ./include -o libHHKinFit.so

echo "creating executable"
g++ ComputeHHKinFit.C `root-config --cflags --glibs` -D HHKINFIT_STANDALONE -I ./include -L . -lHHKinFit  -o ComputeHHKinFit
