# Example 1
python3 ompsimulator.py 0 dpname=omp1 lamm=0.1 lama=0.01 lamh=0.00001 nax=10 nol=2 taug=30 taudA=30 taumax=100 taumin=3 taunom=50 fixedel=nrn_5 tref=0 tr=200 jitter=1 spksig=sync saver=5 nwarmup=1 nreps=2 nepochs=20 Tesec=5 name=testrun1omp randseed=12345
# Example 2
python3 ompsimulator.py 0 dpname=omp1 lamm=0.05 lama=0.01 lamh=0.00001 nax=10 nol=3 taug=30 taudA=30 taumax=100 taumin=3 taunom=50 fixedel=nrn_10 tref=0 tr=500 jitter=1 spksig=sync saver=5 nwarmup=1 nreps=2 nepochs=20 Tesec=5 name=testrun2omp randseed=12345
# Example 3
python3 ompsimulator.py 0 dpname=iomp1 lamm=0.02 lama=0.01 lamh=0.00001 nax=10 nol=3 taug=30 taudA=30 taumax=100 taumin=3 taunom=50 fixedel=nrn_7 tref=0 tr=500 jitter=1 spksig=sync saver=5 nwarmup=1 nreps=2 nepochs=20 Tesec=5 name=testrun3omp randseed=12345

