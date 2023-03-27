#!/usr/bin/env python

from sputils import keyboard
import os, sys, time, glob, re
import numpy as np
import scipy as sp
from itertools import product
from analresultsutils import parse_params_file, dict2prmlist, valstr2prmeqs, remap_resultfiles_submissions

sys.path.append("./")

pycodefile='ompsimulator'
paramdictnames='dpmodel_params, signal_params, run_params, osimul_params, predelay_defprms'
exec('from %s import %s' % (pycodefile, paramdictnames))
exec('paramdicts=[%s]' % paramdictnames)
#from ompsimulator import dpmodel_params, signal_params, run_params, osimul_params, predelay_defprms
#paramdicts=[dpmodel_params, signal_params, run_params, osimul_params, predelay_defprms]

# in each working directory edit get_paramdicts.py that will define the default parameter dicts and their names

pythoncommand='python3'

nargs=len(sys.argv)

if nargs<2:
  print('Usage: '+sys.argv[0]+' projname [genfile]')
  sys.exit(-1)

# directory where information about the submission will be stored
subminfodir='scripts/'
# directory where generator/text files are stored
genfilesdir='scripts/'
# directory where scripts to be executed are stored (mainly via swarm command on a parallel computer)
swfinfodir='scripts/'

if nargs==2:
  if genfilesdir in sys.argv[1]:
     inprmfile=sys.argv[1]
     globrunname=inprmfile.replace(genfilesdir,'').replace('.txt','')
  else:
     globrunname=sys.argv[1]
     inprmfile=genfilesdir+globrunname+'.txt'
else:
  globrunname=sys.argv[1]
  inprmfile=sys.argv[2]

print("Generating file=", inprmfile)
print("Project Name =", globrunname)

subminfofile = subminfodir + globrunname + '-subminfo.npy'

if os.path.exists(subminfofile):
  swffiles=glob.glob(swfinfodir+'%s*.swarm' % globrunname)
  rmfiles=swffiles + [subminfofile, subminfodir + globrunname + '-params.npy']
  print(rmfiles)
  ans=input('The project %s already exists! The files above will be removed! Do you want to proceed anyway? ("Yes" to proceed) default (n)' % globrunname)
  if ans == 'Yes':
    print("Removing files:")
    rmcommand='rm '+' '.join(rmfiles)
    os.system(rmcommand)
  else:
     print('Skipping the creation! Change the name or answer "Yes" at the prompt!')
     sys.exit(-1)

submlines,[pycommand,commandargs,refcodefile,prmsldict,submissiontext]=parse_params_file(inprmfile)
keyboard('check here')

if pycommand is None:
  print('Using the defauld pythoncommand!')
else:
  pythoncommand=pycommand
  
print('Using the pythoncommand=', pythoncommand)
  
theexeccode=open(pycodefile+'.py').read()

if True:
  from diff_match_patch import diff_match_patch
  refexeccode=open(refcodefile).read()
  dmp = diff_match_patch()
  patches = dmp.patch_make(refexeccode, theexeccode)
  execcodeinfo=[refcodefile,patches]
#  diff = dmp.patch_toText(patches)
#  keyboard('check diff')
else:
  execcodeinfo=['', theexeccode]

codeargscommand='%s.py %s' % (pycodefile,' '.join(commandargs))
  
print("submlines=", submlines)
print("len(submlines)=", len(submlines))
for slin in submlines:
  print(' '.join(slin))

subminfo={}

irun=0

allswarmcommands=[]
for irun, argline in enumerate(submlines):
   runname='%s-%d' % (globrunname, irun+1)
   argstring=' '.join(argline)
   command='%s %s %s name=%s' % (pythoncommand, codeargscommand, argstring, runname)
   subminfo[runname]=' '+argstring+' '
   print(command)
   allswarmcommands.append(command)

# save the submission run information
print('Saving submission info!')
np.save(subminfofile, subminfo)

commandinfo=[pythoncommand, codeargscommand]
print('Saving params/code info!')
np.save(subminfodir+ globrunname+'-params.npy', np.array([prmsldict, submissiontext, execcodeinfo, commandinfo, paramdicts, paramdictnames],dtype=object))

maxperfile=1000
iswf=0
fswm=0

for iswm, swmcommand in enumerate(allswarmcommands):
  if iswm%maxperfile==0:
     iswf+=1
     if fswm:
       fswm.close()
     swarmfile= swfinfodir + '%s-%d.swarm' % (globrunname, iswf)
     fswm=open(swarmfile, 'w')
  print(swmcommand, file=fswm)

fswm.close()

  
