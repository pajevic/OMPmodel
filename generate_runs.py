#!/usr/bin/env python

import os, sys, time, glob, re
import numpy as np
import scipy as sp
from itertools import product

sys.path.append("./")

pycodefile='ompsimulator'
paramdictnames='dpmodel_params, signal_params, run_params'
exec('from %s import %s' % (pycodefile, paramdictnames))
exec('paramdicts=[%s]' % paramdictnames)
#from ompsimulator import dpmodel_params, signal_params, run_params, osimul_params, predelay_defprms
# predelay_defprms changed to fixedel_defprms
#paramdicts=[dpmodel_params, signal_params, run_params, osimul_params, predelay_defprms]
# in each working directory edit get_paramdicts.py that will define the default parameter dicts and their names

pythoncommand='python3'


def valstr2prmeqs(prmname, prmvalstring):
     if ':' in prmvalstring:
       ncol=prmvalstring.count(':')
       if ncol==1:
         st,end=map(tryany, prmvalstring.split(':'))
         step=1
       elif ncol==2:
         st,step,end=map(tryany, prmvalstring.split(':'))
       return ["%s=%s" % (prmname,val) for val in np.arange(st,end,step)]
     else:
       vals=prmvalstring.split(',')
       return ["%s=%s" % (prmname,val.strip()) for val in vals]


def parse_params_file(prmfile):
   prmsdict={}#not used currently
   prmeqslist=[]
   prmslist=[]
   allines=[]
   commandargs=[]
   reference_code_file=''
   pycommand=None
   
   for rprmline in open(prmfile).readlines():
     prmline=rprmline.strip()
     allines.append(prmline)
     if (not prmline) or (prmline[0]=='#'): continue
     if '#' in prmline:
        prmline=prmline[:prmline.index('#')].strip()
     prmname,prmvalstring=prmline.split('=')
     if prmname=='command_args':
       commandargs=[cpstr.replace('.eq.','=') for cpstr in prmvalstring.split(',')]
     elif prmname=='pycommand':
       pycommand=prmvalstring
     elif prmname=='reference_code':
        reference_code_file=prmvalstring
     elif prmname=='defparams':
          os.system(prmvalstring)
     elif prmname=='savefile':
        pass
     else: # found the parameter variable
       if prmname in prmslist:
         iprm=prmslist.index(prmname)
         print('Already found the same parameter %s in the list with the values %s!' % (prmname, prmeqslist[iprm]))
         print('OVERWRITTING with:', prmvalstring)
         prmeqslist[iprm]=valstr2prmeqs(prmname, prmvalstring)
       else:     
         if prmvalstring in prmslist:
            iprm=prmslist.index(prmvalstring)
            print('NOT IMPLEMENTED ASSIGN')
            newstrings=[eqp+ ' ' + eqp.replace(prmvalstring,prmname) for eqp in prmeqslist[iprm]]
            prmeqslist[iprm]=newstrings
         else:
           prmslist.append(prmname)
           prmeqslist.append(valstr2prmeqs(prmname, prmvalstring))
   return list(product(*prmeqslist)), [pycommand,commandargs, reference_code_file, dict(zip(prmslist, prmeqslist)), '\n'.join(allines)]

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

if pycommand is None:
  print('Using the defauld pythoncommand!')
else:
  pythoncommand=pycommand
  
print('Using the pythoncommand=', pythoncommand)
  
theexeccode=open(pycodefile+'.py').read()

if refcodefile.lower() != 'none':
  from diff_match_patch import diff_match_patch
  refexeccode=open(refcodefile).read()
  dmp = diff_match_patch()
  patches = dmp.patch_make(refexeccode, theexeccode)
  execcodeinfo=[refcodefile,patches]
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

  
