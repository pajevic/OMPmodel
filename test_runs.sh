#!/bin/bash
echo '*********************** Example 1:'
python3 ompsimulator.py 0 dpname=omp1 lamm=0.1 lama=0.01 lamh=0.00001 nax=10 nol=2 taug=30 fixedel=nrn_5 tref=0 taus=200 jitter=1 spksig=sync saver=5 nwarmup=1 nreps=2 nepochs=20 Tesec=5 name=testrun1omp randseed=12345
echo '----------------------------------------------------------------------------------------------'
diff -q results/testrun1omp-results.npy results/testrun1orig-results.npy && echo SUCCESS! || echo "ERROR: Files are different!"
echo '----------------------------------------------------------------------------------------------'
echo '*********************** Example 2:'
python3 ompsimulator.py 0 dpname=omp1 lamm=0.05 lama=0.01 lamh=0.00001 nax=10 nol=3 taug=20 fixedel=nrn_10 tref=0 taus=500 jitter=3 spksig=sync saver=5 nwarmup=1 nreps=2 nepochs=20 Tesec=5 name=testrun2omp randseed=12345
echo '----------------------------------------------------------------------------------------------'
diff -q results/testrun2omp-results.npy results/testrun2orig-results.npy && echo SUCCESS! || echo "ERROR: Files are different!"
echo '----------------------------------------------------------------------------------------------'
echo '*********************** Example 3:'
python3 ompsimulator.py 0 dpname=iomp1 lamm=0.02 lama=0.01 lamh=0.00001 nax=10 nol=3 taug=30 fixedel=nrn_5 tref=0 taus=500 jitter=1 spksig=sync saver=5 nwarmup=1 nreps=2 nepochs=20 Tesec=5 name=testrun3omp randseed=12345
echo '----------------------------------------------------------------------------------------------'
diff -q results/testrun3omp-results.npy results/testrun3orig-results.npy && echo SUCCESS! || echo "ERROR: Files are different!"
echo '----------------------------------------------------------------------------------------------'
