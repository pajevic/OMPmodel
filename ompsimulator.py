#!/usr/bin/env python

from omputils import *

try:
  import pyspike as pyspk  #  to be used in future versions for feeding spikes into OMP model
  have_pyspike=True
except:
  have_pyspike=False
  
try:
  import elephant as elph  # to be used in future versions for feeding spikes into OMP model
  have_elephant=True
except:
  have_elephant=False

verbose = 0 # non-zero values incrementaly increase verbosity of print-out reports
vreport=False

def_dpmodel_params = {'dpname': 'oliga2', 'tau0': 'taunom', 'taunom': 20, 'taumin': 3, 'taumax':40, 'drreps': 1.e-5, 'Qa': 1, 'taurA': 10, 'taudA': 20, 'Qb': 10, 'taurB':5, 'taudB':10, 'OEamp': 40, 'OEtau': 500., 'zthresh': 5, 'noligos': 10,  'omat': None, 'lrnrate':1, 'adelta': 1/100}

def_signal_params = {'spdata': 'pois', 'tspan': [0, 10000], 'tr': 100, 'ntrains': 10, 'tunit': 'ms', 'reft': 50, 'naxons': 10, 'jitter': 0, 'pnsync': 0, 'kg' : 0, 'predelay':'nrn_10'}

def_run_params={'name': 'test0p1runs', 'solver':'ivp', 'tres':0.01, 'method': 'RK45'}

def_osimul_params={'vtspan': -500, 'saveframes': True}

if verbose>3: ncalls=0

defsaveresults=5  # saveresults=0: Nothing saved!  saveresults=1 (or True): save most of things 2: omit individual oligo timings 3: save only spread information and the used parameters 4: save only spread info (use only for large parameter explorations)

savegmodelcallshistory=False
saveindividualspreadcurves=False

nrecordaver=20
savematrixhistory=True

yyhistrecord=[]
yysehistrecord=[]
yymean=[]
yyse=[]
yyfirst=[]
yyfirstse=[]
tthistrecord=[]

results_dir, logs_dir = get_results_dir(return_log=True)

if not os.path.exists(logs_dir): os.mkdir(logs_dir)

codefile=sys.argv[0]
refcodefile='vrefs/vref_9May2022_oligosimulator.py'

if not os.path.exists(refcodefile):
  if not os.path.exists('vrefs/'): os.mkdir('vrefs')
  refcodefile='vrefs/vref_'+currentdate()+codefile
  if not os.path.exists(refcodefile):
    print('UPDATE REFFILE in the python code : refcodefile : Make it automatic in the future: ')
    os.system('cp %s %s' % (codefile, refcodefile))
    
execcodeinfo=get_codeinfo(codefile, refcodefile)

if not (defsaveresults or savegmodelcallshistory or saveindividualspreadcurves):
  print('WARNING: You are not saving any results!!!')
  ans=input("Do you want to exit? ")
  if ans.lower()[0]=='y':
     sys.exit(-1)
     
nargs=len(sys.argv)

isigen2=np.random.exponential(2, size=10000000)

if savegmodelcallshistory:
  nmaxstorg=10000000
  gmodelcallshistory=np.zeros((nmaxstorg,5))

Tsecs=3
normalize_delays=0
sprrfact=20  # factor to scaled the spread, to be plotted with other similarity measures  
gather_all_spikes=False
gprogvars={}
gprogvars['nbinsisi']=30
gprogvars['mintisi']=0
gprogvars['maxtisi']='tau0'
gprogvars['maxtisi']='tr_9.5'
gprogvars['displayfig']=0
#gprogvars['gax0']
#gprogvars['gah0']

ddrprevtemp=-1
grunhelpers=[0,0,0,0,0,0]

dmratemin=1e-12
dmratemin_assign=1e-11

dpmodel_params = {'dpname': 'oliga2', 'tau0': 'taunom', 'taunom': 20, 'taumin': 3, 'taumax':40, 'dmrateinit': None, 'drreps': 1.e-6, 'Qa': 1, 'taurA': 10, 'taudA': 20, 'Qb': 10, 'taurB':5, 'taudB':10, 'OEamp': 40, 'OEtau': 500., 'zthresh': 5, 'naxons': 10, 'noligos': 5,  'omat': None, 'lrnrate':1, 'myelrate': 0.1, 'adelta': 1/100}

Tms=1000*Tsecs
# if ntrains=0 assign it to be the same as naxons
signal_params = {'spdata': 'pois', 'tspan': [0,Tms], 'tr': 100, 'ntrains': 0, 'tunit': 'ms', 'reft': 50, 'naxons': 10, 'jitter': 0, 'trpjitt' : 0, 'pnsync': 0, 'kg' : 0, 'predelay':'nrn_10'}

run_params={'name': 'test0p1runs', 'solver':'ivp', 'Tsecs': Tsecs, 'tres':0.01, 'method': 'RK45', 'nreps': 1, 'nsreps': 10, 'nwarmup':1, 'randseed':-1, 'saver':  defsaveresults}

predelay_defprms={'rnd': [20], 'rng': [3], 'none':[] }

osimul_params={'vtspan': -3000, 'saveframes': False}

if __name__ == '__main__':

  if nargs<2:
     print('Usage: '+sys.argv[0]+' verbose[=%d]' % (verbose))
     sys.exit(0)

  print(' '.join(sys.argv))
  
  eqargs, oargs= filter_eq_args(sys.argv)

  noargs=len(oargs)
  
  iarg=1
  if noargs>iarg: verbose=int(oargs[iarg])
  iarg+=1
  
  for ieq in eqargs:
    var,val=ieq.split('=')
    if  var in dpmodel_params:
        dpmodel_params[var]=tryeval(val)
    elif  var in signal_params:
        signal_params[var]=tryeval(val)
    elif  var in osimul_params:
        osimul_params[var]=tryeval(val)
    elif  var in run_params:
        run_params[var]=tryeval(val)
    elif  var in gprogvars:
        gprogvars[var]=tryeval(val)
    else:
       print("ieq=", ieq)
       print("***** WARNING var=%s is not recognized" % var)

  if run_params['randseed']>2:
    np.random.seed(run_params['randseed'])
    rnd.seed(run_params['randseed']*2)
#    keyboard('setting random seed')
  else:            
    rndseed=9057
    rndseed=int(time.time()%10000)
    if verbose>1: print("rndseed=", rndseed)
    np.random.seed(rndseed)
    rnd.seed(rndseed+17)

  Tsecs=run_params['Tsecs']
  Tms=1000*Tsecs
  signal_params['tspan']=[0,Tms]

  if isinstance(gprogvars['maxtisi'], str):
    if 'tr' in gprogvars['maxtisi']:
       if '_' in gprogvars['maxtisi']:
         itms= gprogvars['maxtisi'].split('_')
         gprogvars['maxtisi']=signal_params['tr']/float(itms[1])
         if verbose>5: print("maxtisi=", gprogvars['maxtisi'])
       else:
         gprogvars['maxtisi']=signal_params['tr']
       
  def update_signal_params(sigp):
      if 'tr' in sigp:
        visi=sigp['tr']
        if type(visi)==str:
            if visi[:3]=='rnd':
               v1,v2,v3=visi.split('_')
               mnv=float(v2)
               mxv=float(v3)
               sigp['tr']=np.random.rand(sigp['ntrains'])*(mxv-mnv)+mnv
  #    return sigp
  
  def update_dpmodel_params(dpmp, sgpm):
      if 'omat' in dpmp:
        vomat=dpmp['omat']
        if vomat in ['all', 'ones', 'full']:
          dpmp['omat']=np.ones((dpmp['noligos'], dpmp['naxons']))
        elif type(vomat)==str:
            if vomat[:3]=='rnd':
               items = vomat.split('_')
               nax = int(items[1])
               noligo=1
               prob=0.5
               if len(items)>2: noligo = int(items[2])
               if len(items)>3: prob=float(items[3])
               dpmp['omat']=(np.random.rand((noligo,nax))<prob).astype('int')
               dpmp['naxons']=nax
               dpmp['noligos']=noligo
        elif type(vomat)==tuple:
            if len(vomat)==2:
               dpmp['omat']=np.ones(vomat)
               dpmp['noligos']=vomat[0]
               dpmp['naxons']=vomat[1]
  #    return dpmp
  
  if True:
     paramabbrevs={}
     paramabbrevs['ntrains']=('nax','%d')
     paramabbrevs['truedelay']=('td','%d'), 
     paramabbrevs['nruns']=('nruns' ,'%d')
     paramabbrevs['Tsecs']=('T' ,'%d')
     paramabbrevs['tr1']=('tra','%d')
     paramabbrevs['tr2']=('trb','%d')
     paramabbrevs['reft']=('rf', '%d')
     paramabbrevs['actrise']=('ar','%d')
     paramabbrevs['actdecay']=('ad','%d')
  
  nreps=run_params['nreps']
  nsreps=run_params['nsreps']
  nwarmup=run_params['nwarmup']

  # adjust saveresults to see what things need to be saved
  saveresults,nrepsave,nsrepsave,modelhistories=get_save_params(run_params['saver'])
  if nsrepsave: ioligrecord = 0
  else: ioligrecord = nrepsave - 1

  flog=open(logs_dir + '%s-runlogs.log' % run_params['name'], 'a')
  logmaxs=[1000,50,5000,100]
  logcounts=[0,0,0,0]
  
  #remove later meantstdarr=np.zeros((nreps,nsreps))
  #remove later stdtstdarr=np.zeros((nreps,nsreps))
     
#  print("signal_params=", signal_params)

  update_signal_params(signal_params)
  update_dpmodel_params(dpmodel_params, signal_params)
  if verbose>3: 
    print("dpmodel_params=", signal_params)
    print("signal_params=", signal_params)
  
  omat=dpmodel_params['omat']
  taumin=dpmodel_params['taumin']
  taumax=dpmodel_params['taumax']
  taunom=dpmodel_params['taunom']
  tau0=dpmodel_params['tau0']
  drreps=dpmodel_params['drreps']
  if drreps<(-floateps[4]):
    drreps=10.**drreps
  
  if tau0=='taunom':
  #  tau0=dpmodel_params['taunom']
    tau0=taunom
    
  lrnrate=dpmodel_params['lrnrate']
  myelrate=dpmodel_params['myelrate']
  
  adelta=dpmodel_params['adelta']

  noligos=dpmodel_params['noligos']
  naxons=dpmodel_params['naxons']
  
  if omat is not None:
    assert omat.shape == (noligos1,naxons1)
  
  # set of global variables to store results
  kg=signal_params['kg']
  pnsync = signal_params['pnsync']
  spdata = signal_params['spdata']
  trpjitt = signal_params['trpjitt']

  if pnsync>0 or 'psync' in spdata: kg=2
  if 'ksync_' in spdata: kg=int(spdata.split('_')[1])
  if 'msync_' in spdata:
    kg=len(spdata.split('_')[1])
    signal_params['kg']=kg

  if kg>1:
    tstdarr=np.zeros((nreps,nsreps,noligos,kg+5))
    inittstdarr=np.zeros((nreps, kg+5))
  else:
    tstdarr=np.zeros((nreps,nsreps,noligos,1))
    inittstdarr=np.zeros((nreps,1))
  gaxondelays=np.zeros((nreps,nsreps,naxons))
  gpredelays=np.zeros((nreps,naxons))
  goligoinfo={'oevents': [[[] for i in range(nsreps)] for j in range(nreps)], 'orates': np.zeros((nreps,nsreps, noligos)), 'odelays': None}
  
  print("nsreps=", nsreps)
  print("nreps=", nreps)

  genexp2=lambda x: np.random.exponential(2)
  
  otaumin=taumin/noligos
  otaumax=taumax/noligos
  otaunom=taunom/noligos
  
  zthresh=dpmodel_params['zthresh']
  taudA=dpmodel_params['taudA']
  taurA=dpmodel_params['taurA']
  Qa=dpmodel_params['Qa']
  
  taudB=dpmodel_params['taudB']
  taurB=dpmodel_params['taurB']
  Qb=dpmodel_params['Qb']
  
  aicA, bicA, cicA, maxfaktA, tmaxA = coeffs_from_release_params(taurA, taudA, Qa)
  aicB, bicB, cicB, maxfaktB, tmaxB = coeffs_from_release_params(taurB, taudB, Qb)
  
  krecov=1./dpmodel_params['OEtau']
  ampReset=dpmodel_params['OEamp']*maxfaktA
  
  maxvfunc=factor_impresp(tmaxA, Qa, taudA, taurA)
    
  tr1=signal_params['tr']
  gjitter=signal_params['jitter']
  
  texpolig=tr1/naxons  # spike rate/period seen by oligo dendrocyte
  totArate = 2*texpolig/(taurA+taudA)
# dmrateinit = \lambda_M N_A/\tau_s^2$   
  procrate=1/tr1
  Aaver1=maxfaktA*totArate
  dmrateinitpass1= lrnrate * Aaver1 * procrate
  xtr=np.arange(1,50)*texpolig
  Aaver2=sum(factor_impresp(xtr, Qa, taudA, taurA))
  dmrateinitpass= lrnrate * Aaver2 * procrate
  dmrateinitpass= lrnrate * Aaver2 * procrate
  dmrateinitpass=lrnrate*naxons/tr1**2
  dmrateinitgiven=dpmodel_params['dmrateinit']
  if dmrateinitgiven is not None:
    dmrateinitpass = dmrateinitgiven
  print("dmrateinitpass=", dmrateinitpass)
  
  if verbose>1 and dpmodel_params['dpname'] in ['iomp1','omp1']: 
    print("Aaver1=", Aaver1)
    print("Aaver2=", Aaver2)
    print("dmrateinitpass1=", dmrateinitpass1)
    print("dmrateinitpass=", dmrateinitpass)
  
  oeventrate=0.001
  dmrateinit= lrnrate * maxfaktB * oeventrate
  crit_lrnrate=(otaumax-otaumin)/maxfaktB
  
  fsfunc=lambda x: fslin(x, otaumin, otaumax, otaunom)
  fsfuncprim=lambda x: fslinprim(x, otaumin, otaumax, otaunom)
  
  def fsfuncmyel(x):
    if x > otaumin: return (x-otaumin)/(otaumax-otaumin)
    else: return 0
    
  def fsfuncdemyel(x):
    if x < otaumax: return (otaumax - x)/(otaumax-otaumin)
    else: return 0

#  vec_fsfuncdemyel = np.vectorize(fsfuncdemyel)
#  vec_fsfuncmyel = np.vectorize(fsfuncmyel)

  gnaxons=naxons
  if  dpmodel_params['dpname'] in ['iomp1', 'omp1']: itau=3
  if dpmodel_params['dpname'] in ['iomp1', 'oliga1']: gdydts=np.zeros(itau+gnaxons)
  elif dpmodel_params['dpname'] in ['omp1', 'oliga2']: gdydts=np.zeros(itau+2*gnaxons)
  
  def dydt_oligo_passive_diracmyelination(t,y):
  # passive Dirac change
  # passive
      global logcounts, ddrprevtemp, grunhelpers, ncalls
      ncalls+=1
      du = y[1]
      dv = - aicA * bicA * y[0] - (aicA + bicA)*y[1] # + c \delta(t)
      dmrate = y[2]
      if omat is None:
        ddr = dmrate*drreps*(otaunom-y[itau:itau+gnaxons].mean())
        if dmrate<dmratemin:
           if logcounts[0]< logmaxs[0]:
              if verbose>2: print('#%d dmrate problem!!!!!  dmrate=%g ddrprevtemp=%g ddr=%g' % (logcounts[0], dmrate, ddrprevtemp, ddr))
              logcounts[0] +=1
              print('#%d: dmrate problem dmrate=%g ddrprevtemp=%g ddr=%g' % (logcounts[0], dmrate, ddrprevtemp,ddr), file=flog)
           dmrate=dmratemin_assign
        for iax in range(gnaxons):
          ctau=y[iax+itau]
          fsvald=fsfuncdemyel(ctau)
          dtau = dmrate * fsvald
          gdydts[itau + iax]=dtau
      else:
        print('Running OL-axon connectivity studies is left a faster implemention! Not implemented for diracmyelination')
        pass

          
      if verbose>5: print("dmrate=", dmrate, " ddr=", ddr, '\b'*300, end='')
      
      ddrprevtemp=ddr
      if savegmodelcallshistory:
        if ncalls<=nmaxstorg:
           gmodelcallshistory[ncalls-1,:]=[t, dmrate, ddr, fsvald, ctau]
      gdydts[:itau]=[du, dv, ddr]
      return gdydts
  
  def dydt_oligo_passive(t,y):
  # passive OMP model
      if myelrate>1.e4: return dydt_oligo_passive_diracmyelination(t,y)
      global logcounts, logmaxs, ddrprevtemp, grunhelpers, ncalls
      ncalls+=1
      du = y[1]
      dv = - aicA * bicA * y[0] - (aicA + bicA)*y[1] # + c \delta(t)
      dmrate = y[2]
      
      if omat is None:
        meantau=y[itau:itau+gnaxons].mean()
        ddr = dmrate*drreps*(otaunom-meantau)
        if dmrate<dmratemin:
           if logcounts[0]< logmaxs[0]:
              print('#%d dmrate problem!!!!!  dmrate=%g ddrprevtemp=%g ddr=%g' % (logcounts[0], dmrate, ddrprevtemp, ddr))
              logcounts[0] +=1
              print('#%d: dmrate problem dmrate=%g ddrprevtemp=%g ddr=%g' % (logcounts[0], dmrate, ddrprevtemp,ddr), file=flog)
           dmrate=dmratemin_assign
        for iax in range(gnaxons):
          ctau=y[iax+itau]
          mprom=y[itau + gnaxons+ iax]
          fsvald=fsfuncdemyel(ctau)
          fsvalm = fsfuncmyel(ctau)
          dtau = dmrate * fsvald - myelrate * mprom * fsvalm# - xxx \delta(t)
          dmprom = - myelrate * mprom
          gdydts[itau + gnaxons +  iax]= dmprom
          gdydts[itau + iax] = dtau
          
        if savegmodelcallshistory:
          if ncalls<=nmaxstorg:
             gmodelcallshistory[ncalls-1,:]=[t, dmrate, ddr, fsvald, ctau]
        gdydts[:itau]=[du, dv, ddr]
      else:
        pass
      
      return gdydts
  
  evthresh=zthresh*maxfaktA
  
  def dydt_model3(t,y):
      itau=0
      dtau = dmrate*fsfunc(y[itau])
      du = y[2]
      dv = - aic * bic * y[1] - (aic + bic)*y[2] # + c \delta(t)
      return [ dtau, du, dv]
  
  def event_func3(t, y):
      return y[1] - evthresh
  event_func3.terminal = True
  event_func3.direction = 1
  
  def jac_model3(t,y):
      return [[ dmrate*fsfuncprim(y[0]), 0, 0], [ 0, 0, 1], [ 0, -aic*bic, -aic-bic]]
  
  def event_none(t, y):
      print('THIS SHOULD NOT EVER HAPPEN!\n'*20)
      return -1
  event_none.terminal = True
  event_none.direction = 1
  
  oligoinfo={'odelays': None, 'oevents': None, 'orate': 0}
  
  def run_dp_learning(dp_params, sig_params, runparams={}, plotit=False, params={}):
    
    global axondelays, oligodelays, oligomfactors, gaxondelays, gpredelays, dmrateinit, gprogvars, logcounts, ncalls, axmathistrecord

    if verbose>3: print("runparams=", runparams)
  
    csolver=runparams['solver']
    tres=runparams['tres']
    odemethod=runparams['method']
    dpmname=dp_params['dpname']
    
    try:
      tstart = sig_params['tspan'][0]
      tend = sig_params['tspan'][1]
    except:
      tstart=0
      tend = sig_params['tspan']
  
  # INITIALIZATION OF MODELS
    if dpmname=='iomp1':
      y0A_nontau=[0, cicA, dmrateinitpass]
      if csolver in ['odeint']:
        keyboard('is this happening')
        DPmodel=lambda y,t: dydt_oligo_passive_diracmyelination(t,y)
      else:
        DPmodel=dydt_oligo_passive_diracmyelination
        event_func=event_none
#        mjac=jac_model3
    elif dpmname in ['omp1']: # passive model with myelrate myelfactor
      print("Assigned initial dmrateinitpass=", dmrateinitpass)
      y0A_nontau=[0, cicA, dmrateinitpass]
      if csolver in ['odeint']:
        keyboard('is this happening')
        DPmodel=lambda y,t: dydt_oligo_passive(t,y)
      else:
        DPmodel=dydt_oligo_passive
        event_func=event_none
#        event_func=event_func5
#        mjac=jac_model3
    elif dpmname=='simp1':
      y0A=[tau0]
      if csolver in ['odeint']:
        DPmodel=lambda y,t: dydt_model1(t,y)
      else:
        DPmodel=dydt_model1

    retparams={}
    tr1=sig_params['tr']
    ntr=signal_params['ntrains']
    reft=signal_params['reft']

    lpnsync=pnsync
    if 'psync' in spdata:
        psyncitms=spdata.split('_')
        try:
          lpnsync=int(psyncitms[1])
        except:
          if psyncitms[1]=='half':
            lpnsync=naxons//2
          elif psyncitms[1]=='third':
            lpnsync=naxons//3
        print("lpnsync=", lpnsync)
           
    predelayspec=signal_params['predelay']
    
  #  signal_params = {'spdata': 'pois', 'tspan': [0,Tms], 'tr': 50, 'ntrains': 10, 'tunit': 'ms', 'reft': 0}
  #  tr1=100
  #  spta=generate_spike_train(tr1, Tms, reft=0)
  #  spta=SpikeTrains('pois', tr=tr1, ntrains=naxons, tspan=Tms, reft=reft)
  #  if type(spdata)==str and len(spdata)>3 and spdata[:4]=='case':
    labels=None
    prandtau0=0.05 # randomize initial taus for all axons around mean tau0
    if nrepsave:
      axmathistrecord=np.zeros((2,nreps,nsreps,noligos,naxons))
    for irep in range(nreps):  # LOOP 1 (Offs 2)
       ncalls=0
       logcounts[0] = 0
       tthistrecord.append([])
       yyhistrecord.append([])
       yysehistrecord.append([])
       yymean.append([])
       yyse.append([])
       yyfirst.append([])
       yyfirstse.append([])
       print('Running main loop %d of %d' % (irep+1, nreps))
       sys.stdout.flush()
       oligorates=np.zeros(noligos)
       oligodelays=tau0*(1. + prandtau0*randn(noligos, naxons))/noligos
       
       oligomfactors=np.zeros((noligos, naxons))
       axondelays=oligodelays.sum(0)
       y0B_nonax=list(y0A_nontau)
       if verbose: print("y0A_nontau=", y0A_nontau)
       if isinstance(predelayspec, str):
          if '_' in predelayspec:
             itms=predelayspec.split('_')
             pdname=itms[0]
             prms=[tryany(pstr) for pstr in itms[1:]]
          else:
             pdname=predelayspec
             prms=predelay_defprms[pdname]
          if pdname == 'rnn':
             predelays=prms[0]*(10 + randn(naxons))
          elif pdname == 'nrn':
             print("kg=", kg)
             predelays = assign_normalized_predelays(naxons, kg, sigma=prms[0], normall=False)
          elif pdname == 'nrna':
             predelays = assign_normalized_predelays(naxons, kg, sigma=prms[0], normall=True)
          elif pdname == 'rnd':
             predelays=prms[0]*rand(naxons)
          elif pdname == 'rng':
             predelays=prms[0]*np.arange(naxons)
          elif pdname == 'none':
             predelays=np.zeros(naxons)
       elif np.isscalar(predelayspec):
            predelays=predelayspec*np.ones(naxons)
         
       if verbose>-1: print("PREDELAYS: irep=%d predelays=" % irep, predelays)
       gpredelays[irep,:]=predelays
       if lpnsync: sspreads=get_atime_spread2(predelays+axondelays, lpnsync)
       else: sspreads=get_atime_spread_kgroups(predelays+axondelays, kg=kg)
       
       pdelspread=get_atime_spread_kgroups(predelays, kg=kg)
       print("predelay spread=", pdelspread)
       print("Init spread =", sspreads)

       if verbose>5:
         print("irep=", irep)
         print("pdelspread=", pdelspread)
         print("sspreads=", sspreads)
         print("pdelspread[:kg-4]=", pdelspread[:kg+2])
         print("sspreads-4=", sspreads[:kg+2])
#       keyboard('check the spread')

       inittstdarr[irep,:]=sspreads

       for isrepw in range(nsreps+nwarmup):  # LOOP 2 (creates a new set of signals, but keeps already learned delays and dmrate) Offs 5
            isrep=isrepw-nwarmup
            if vreport:
              print("\n\nStarting new isrep with ncalls=",  ncalls)
            if  isrep>=0 and isrep%10==9:
               print( "tr1=", tr1, "taudA=", taudA, "naxons=", naxons, "noligos=", noligos )
               print("\nisrep=%d of %d" % (isrep+1, nsreps))

            oligoevents=[[] for _ in range(noligos)]
# GET SPIKE TRAINS TO FEED into the chain
# get the spike train for this repetition
#            print("tr1=", tr1)
            spta = quick_generate_spiketrain(spdata, naxons, tspan=Tms, tr=tr1, reft=reft, trpjitt=trpjitt)
            if False:
              spta.plot_spikes()
#              keyboard('Check generated spikes spta.plot_spikes()')
#            keyboard('check train')

            if gjitter:
              spta.jitterspikes(gjitter)
          #    plt.figure()

            ntotspikes = len(spta.spikes)
            retparams['ntotspikes']=ntotspikes
            retparams['ospiketrain']=spta.sptimes
            
            prespikes1=spta.spikes
            prespikes2=shiftspikes(prespikes1, predelays)
            prespikes=post2prespikes(prespikes2)
            
            if gather_all_spikes:
               allspikes=[prespikes]

#            first_calculate_distances_section_in_the_loop
            for iolig in range(noligos):  # LOOP 3  loop along the OLs in the chain Offs10
                postspikes=[0 for _ in prespikes]
                tprev=tstart
                tprevolig=tstart
                # initialize the state for the next oligdendrocite in the chain; start where it stopped in the last nsrep iteration
                oligorate=0
                if dpmname in ['iomp1','oliga1']:
                    y0=y0B_nonax + list(oligodelays[iolig,:])
                elif dpmname in ['omp1', 'oliga2']:
                    y0=y0B_nonax + list(oligodelays[iolig,:]) + list(oligomfactors[iolig,:])
                if verbose>2:
                   print("y0=", y0)
          
          #      axondelays=oligodelays.sum(0)
                if normalize_delays: # the old and naive way of normalizing delays (OBSOLETE: Will be removed)
                   alltaus=axondelays
                   axondelays=tau0*alltaus/alltaus.mean()
                   if (min(axondelays)<taumin) or (max(axondelays)>taumax):
                      raise ValueError('exceeed tau range')
                      keyboard('exceeed tau range')
                      
                # initialize the lists for collecting solutions for a given oligodendrocyte in the chain
                solst=[]
                solsy=[]
                for ispk, spk1  in enumerate(prespikes): # @START OF SPIKE LOOP
                  tsp=spk1[0]
                  inode=spk1[1]
                  if verbose>6: print("Starting with y0=", y0[itau+inode])
              #    nptint=int(np.round((tsp-tprev)/tres))
              #    t1 = np.linspace(tprev, tsp, nptint)
              #    y1 = odeint(DPmodel, y0, t1)
                  t1spike=[]
                  if verbose>5: print('*** NEW SPIKE ***'*5)
                  if (tprev + floateps[1]) >= tsp:
                    if verbose>3:
                       print('******* REPEATED SPIKE !!! '*50, tsp, tprev, tsp-tprev)
                    odestatus=2
                  else:
                    odestatus=1
          
                  while odestatus==1: # loop for integrating over all oligo-events between two spikes, or over all simultaneous spikes 
              #       sol1 = solve_ivp(DPmodel, (tprev, tsp), y0, method=odemethod, jac=mjac, dense_output=True, events=event_func)
                     if verbose>4: print("INTEGRATING NOW between prev=%g and tsp=%g" % (tprev, tsp))
#                     sol1 = solve_ivp(DPmodel, (tprev, tsp), y0, method=odemethod, dense_output=True, events=event_func)
                     sol1 = solve_ivp(DPmodel, (tprev, tsp), y0, method=odemethod, dense_output=True)
#                     try:
 #                    except:
 #                      print('sol1 = solve_ivp(DPmodel, (tprev, tsp), y0, method=odemethod, dense_output=True, events=event_func)')
 #                      print('FAILED solve_ivp')
 #                      keyboard('FAILED solve_ivp: check DPmodel y0 ')
                     odestatus=sol1.status
                     
                     if odestatus == 1: # OLIGO-SPIKE : Oligoevent was hit
                       passivemodels=['iomp1', 'omp1', 'sync1']
                       if dpmname in passivemodels:
                          print('There should never be an oligoevent in passive models', passivemodels)
                          keyboard('Inspect! GOT oligoevent in passive ')

                       tev=sol1.t_events[0][0]
                       if verbose>3: print('Oligoevent HAPPENED tev=', tev, tsp, iolig)
                       oligoevents[iolig].append(tev)
                       oligorate=len(oligoevents[iolig])/tev
                       if verbose>5:
                          print('OLIGO OLIGO OLIGO --->:\n'*20)
                          print("oligorate=", oligorate)
                          print("oligoevents[iolig']=", oligoevents[iolig])
          
                       if verbose==5:
                          print('- - '*20)
                          print("tprev=", tprev)
                          print("tsp=", tsp)
                          print("tev=", tev)
                       if tprev>=tev:
                         print("tprev=", tprev)
                         print("tev=", tev)
                         keyboard('Deal breaker: tprev >= tev')
          
                       nptsint=int(np.round((tev-tprev)/tres))+2
                       if verbose>5:
                          print("ASSIGNING AT EVENT 2 tprev=%g tev=%g" % (tprev, tev))
                       t1 = np.linspace(tprev, tev, nptsint)
                       t1spike.extend(t1)
                       y1=sol1.sol(t1).T
                       tprev=tev
                       if verbose>3:
                         print("BEFORE OLIG isrep=", isrep)
                         print("Appending results iolig=", iolig)
                         print("len(t1)=", len(t1))
                         print("len(y1)=", len(y1))
                       solst.append(t1)
                       solsy.append(y1)
                       # Prcoesses Oligo-event changes  # integrate delta due to Oligo-spikes
                       y0 = y1[-1] # record the last values from integration
                       t0last=t1[-1]
                     else: # odestatus != 1 
                       if tprev>=tsp:
                           print('Safety Check: THIS SHOULD NEVER HAPPEN!!!')
                           print("Second tprev=", tprev)
                           print("Second tsp=", tsp)
                           tims = [sp1[0] for sp1 in prespikes]
                           tdiffs=np.diff(tims)
                           iltz=np.where(tdiffs<0)[0]
                           print("iltz=", iltz)
                           print('Second tprev >= tsp')
                           sys.exit(0)
                       nptsint=int(np.round((tsp-tprev)/tres))+2
                       t1 = np.linspace(tprev, tsp, nptsint)
                       t1spike.extend(t1)
                       y1=sol1.sol(t1).T
                       y0=y1[-1]  # record the last values for fully integrating
                       t0last=t1[-1]
                       solst.append(t1)
                       solsy.append(y1)

          ####### integrate over NEURONAL SPIKES
          ####### integrate over the delta function (DIRAC handle neuronal spikes NEURON SPIKES
                  tprev=tsp
                  # integrating delta function(s) for the current spike
                  if verbose>3:
                    print('Integrating delta function(s) for inode=%d spike' % inode)
                  if dpmname=='oliga2':  # integrate delta function(s) ; over "nspikes"
          #            thlast,u0,lastv,uB0,lastvB,dmratelast=y0[:itau] # y1[-1] or y0
                      y0[2] += cicA #   v0 = lastvA + cicA
                      uB0=y0[3]
                      mpromlast=y0[itau+gnaxons+inode]
                      taulast=y0[itau+inode]
#  Adjust myelination promoter after the neuronal spikes (or taus)          
#                      jmpval1=lrnrate*uB0*fsfunc(taulast)
#                      jmpval2=lrnrate*uB0*fsfunc(taulast)*(taulast-otaumin)
                      jmpval=lrnrate*uB0*fsfuncmyel(taulast) # QQQ2 should this be here also, or just for directmyelination. Use lrnrate only for the last step
                      newmprom = mpromlast + jmpval  # adjustment at nxonal spike
                      y0[ itau + naxons + inode ] = newmprom
                      
                      if verbose>5:
                         thlast,u0,lastv,uB0, lastvB, dmratelast=y0[:itau] # y1[-1] or y0
                         print('---- '*10)
                         print("lmrate=", lrnrate)
                         print("lastvB=", lastvB)
                         print("uB0=", uB0)
                         print("tau0c=", tau0c)
                         print("fsfunc(tau0c)=", fsfunc(tau0c))
                         print("dmratelast=", dmratelast)
              #         a0=goligo[0][0][-1]
              #         goligo[0][0].extend(
              #         goligo[0][0] -= 0.2*adelta*fsfunc(lastauv)   # integrate delta function
          #             tau0c = taulast - 0.2*adelta*fsfunc(taulast)  
          #            oligodelays[iolig,inode] = newtau
          #            axondelays[inode]=newtau
          #            y0[:itau]=[thlast,u0,v0,uB0, lastvB,dmratelast]  # don't name varaibles -- just changed them directly
          #########################################################  sudden-change MODELS ############################################################
                  elif dpmname=='iomp1':  # integrate delta function(s) ; over "nspikes" for PASSIVE-sudden change model
                      y0[1] += cicA #   v0 = lastvA + cicA
#                      mpromlast=y0[itau+gnaxons+inode]
                      Gfact=y0[0]
                      dmratepick=y0[2]
                      taulast=y0[itau+inode]
                      taujmpval=lrnrate*Gfact*fsfuncmyel(taulast)
                      y0[itau+inode] = taulast - taujmpval # integrate delta function(s) ; over "nspikes"
                  elif dpmname=='omp1':  # integrate delta function(s) ; over "nspikes" for PASSIVE- slow change model; no saturation
                      Gfact=y0[0]
                      y0[1] += cicA #   v0 = lastvA + cicA
                      dmratepick=y0[2]
                      mpromlast=y0[itau+gnaxons+inode]
#                      taulast=y0[itau+inode]
#                      jmpval=lrnrate*Gfact*fsfuncmyel(taulast) # QQQ2 should this be here also, or just for directmyelination. Use lrnrate only for the last step
                      jmpval=lrnrate*Gfact  # QQQ2 should this be here also, or just for directmyelination. Use lrnrate only for the last step
                      y0[ itau + naxons + inode ] = mpromlast + jmpval  # adjustment at nxonal spike
                  elif dpmname=='simp1':
                    lastauv=y0
                    tau0c = lastauv - 0.2*adelta*fsfunc(lastauv)  # integrate delta function
                    y0=tau0c
                  else:
                    print("Unknown model=",dpmname)
              #    postspikes.append(tsp+y0[5])
          #        postspikes.append([tsp+taudelays[iolig,inode], inode])
          #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ sudden-change MODELS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
                  
                  if verbose>1:
                     print("Shifting the delay on postspikes delay=%g  inode=%d" % (y0[itau+inode], inode))
          
          # add delay for the next oligo:  # technically it should happen on subseqent spikes, but it doesn't really matter
                  if verbose>6: print("Ending with y0=", y0[itau+inode])
          #        postspikes.append([tsp+y0[itau+inode]-otaunom, inode])
                  postspikes[ispk]=[tsp+y0[itau+inode]-otaunom, inode]
          #      ^^^^^^^^ @END: LOOP OVER SPIKES ^^^^^^
          
  # FINISHED PROCESSING OLIGODENDROCYTES in THE CHAIN
  # USE NOW y0B in pre-loop ; save what you learned in this loop for: dmrateinit=y0[itau-1]
                if verbose: print('Finished processing olig#%d (r#=%d sr#=%d)' % (iolig+1, irep+1, isrep+1))
                y0B_nonax[itau-1]=y0[itau-1] # preserved learned dmrate for the next step (even for warmup run)
                if verbose>2: print("dmrateinit=%.7g dmrateinitpass=%.7g dmratepick=%.7g dmratelast=%.7g" % (dmrateinit, dmrateinitpass, dmratepick, y0[itau-1]))
                if isrep>=0: # save oligodelays, oligomfactors state for the next sub-repetiton; update axondelays and ttimes/not used anymore, sync-spread measure, as well as prepare postspikes as the next prespikes
                  oligodelays[iolig,:]=y0[itau:(naxons+itau)]
                  if dpmname in ['oliga2', 'omp1']:
                      oligomfactors[iolig,:]=y0[(naxons + itau):]
                  axondelays=oligodelays.sum(0)
#                  ttimes=predelays+axondelays
                  if verbose>1:
                     print("Inside loop predelays=", predelays)

                  if pnsync: sspreads=get_atime_spread2(predelays+axondelays, pnsync)
                  else: sspreads=get_atime_spread_kgroups(predelays+axondelays, kg=kg)
                  tstdarr[irep,isrep,iolig,:]=sspreads
                  
                  prespikes=post2prespikes(postspikes)

                oligorates[iolig]=oligorate
                if gather_all_spikes:
                   allspikes.append(prespikes)

                retparams['rspiketrain']=prespikes
                if verbose>1:
                    print("Concatenating the model solutions at iolig=%d irep=%d isrep=%d" % (iolig, irep, isrep))
                    if verbose>3:
                      print("The length should be equal to the number of spikes as it integrates till the last spike!!!")
                      print("len(solst)=", len(solst))
                      print("len(solsy)=", len(solsy))
                    
                if nrepsave>0:
                   tt=np.concatenate(solst)
                   yy=np.concatenate(solsy)
                if nrepsave:
                  if ioligrecord>0 and iolig==0:
                     yyfirst[irep].append(yy.mean(0))
                     yyfirstse[irep].append(yy.std(0)/np.sqrt(yy.shape[0]))
                  if iolig == ioligrecord:
                   yyrecord=yy[-nrecordaver:,:].mean(0)
                   yyrecse=yy[-nrecordaver:].std(0)/np.sqrt(nrecordaver)
                   tthistrecord[irep].append(tt[-nrecordaver:].mean())
                   yyhistrecord[irep].append(yyrecord)
                   yysehistrecord[irep].append(yyrecse)
                   yymean[irep].append(yy.mean(0))
                   yyse[irep].append(yy.std(0)/np.sqrt(yy.shape[0]))
                   if modelhistories and (irep<nrepsave) and (isrep>=0) and (isrep<nsrepsave):
                        print('Adding solution to model histories (irep=%d isrep=%d)' % (irep, isrep))
                        modelhistories[irep].append([solst, solsy])
                        if verbose>3:
                          print("len(solst)=", len(solst))
                          print("len(solsy)=", len(solsy))
                if verbose>4:
                  print("Time for solving csolver %s =" % csolver, time.time()-t0)
                  print("yy.max(0)=", yy.max(0))
          
            # END OF LOOP 3 ( loop along the OLs in the chain )
            if verbose>4: print('END OF LOOP3! Just finished the oligo chain loop for irep=%d isrep=%d ' % (irep, isrep))

            if False:
                retparams['prespikes']=prespikes
                retparams['orates']=oligorates
                retparams['t']=tt
                retparams['y']=yy

            if isrep>=0: # save oligoevents and orates
               goligoinfo['oevents'][irep][isrep]=oligoevents
               goligoinfo['orates'][irep, isrep, :]=oligorates
  
            spread2=np.std(axondelays+predelays)
            if False:
#            keyboard('check tt yy spread2 spreadc1')
               print('spreadc1=[sprrfact*spikedistances[0][ic][5] for ic in range(noligos+1)]')
               spreadc1=[sprrfact*spikedistances[0][ic][-1] for ic in range(noligos+1)]
               if not np.allclose(spreadc1[-1], spread2):
                 print("spreadc1=", spreadc1)
                 print("spread2=", spread2)
                 keyboard('what?')


#   allcurves not used anymore #            allcurves[irep,:,isrep]=spreadc1
            if nrepsave and savematrixhistory:
               axmathistrecord[0, irep, isrep, :, :] = oligodelays
               axmathistrecord[1, irep, isrep, :, :] = oligomfactors
            if isrep>=0: # save oligoevents and orates
               gaxondelays[irep,isrep,:]=axondelays
            if saveindividualspreadcurves:
               np.savetxt(results_dir + 'spreadcrvs/%s-distsLoop%d-%d.txt' % (run_params['name'], irep,isrep), spreadc1)
  
# END OF LOOP 2 ( loop along sub-repeats! Only the signal changes statistically)
       if verbose>3: print('END OF LOOP2! Just finished isreps for irep=', irep)
       if savegmodelcallshistory:
         np.save('modelcallshistory-%d.npy' % irep, gmodelcallshistory)
         gmodelcallshistory[:]=0
         
       if verbose>4: print('^^^^^^^^^^^^ Finished the loop \n'*verbose)
# finished all the loops
# END OF LOOP 1 ( loop along repeats! The FINAL LOOP)

    if verbose>3: print('END OF LOOP1, the Final Loop : Done!')
    if modelhistories:
        print("Saving the history1")
        np.save(results_dir + 'modelhistories-%s.npy' % run_params['name'], np.array(modelhistories, dtype=object))
    return retparams
  
  rparams=run_dp_learning(dpmodel_params, signal_params, run_params, plotit=False, params={})
  
  dummy1=[]

  print("nrepsave=", nrepsave)
  if nrepsave>0: # and nsrepsave==0
     if ioligrecord==0:
        hrecordarr=np.array([tthistrecord, yymean, yyse, yyhistrecord, yysehistrecord], dtype=object)
     else:
        hrecordarr=np.array([tthistrecord, yymean, yyse, yyfirst, yyfirstse], dtype=object)
        
     hrsavefile='hrecs/'+'%s-dynsnapshot.npy' % run_params['name']
     np.save(hrsavefile, hrecordarr)
     if savematrixhistory:
        np.save(hrsavefile.replace('dynsnapshot','mathist'), axmathistrecord)

  if saveresults:
     savefile=results_dir+'%s-results.npy' % run_params['name']
     goligoinfo['odelays']=oligodelays
     goligoinfo['omfs']=oligomfactors
     if saveresults==1:
       infotext='tstdarr, inittstdarr, gpredelays, gaxondelays, goligoinfo, dpmodel_params, signal_params, run_params, infotext'
       np.save(savefile, [tstdarr, inittstdarr, gpredelays, gaxondelays, goligoinfo, dpmodel_params, signal_params, run_params, infotext])
     elif  saveresults==2:
       infotext='tstdarr, inittstdarr, gpredelays, gaxondelays, NULL-goligoinfo, dpmodel_params, signal_params, run_params, infotext'
       np.save(savefile, np.array([tstdarr, inittstdarr, gpredelays, gaxondelays, [], dpmodel_params, signal_params, run_params, infotext],dtype=object))
     elif  saveresults==3:
       infotext='tstdarr, inittstdarr, NULL-gpredelays, NULL-gaxondelays, NULL-goligoinfo, dpmodel_params, signal_params, run_params, infotext'
       np.save(savefile, np.array([tstdarr, inittstdarr, [], [], [], dpmodel_params, signal_params, run_params, infotext],dtype=object))
     elif saveresults>3: # 4
       savedatadict={}
       savedatadict['saveresults'] = saveresults
       savedatadict['tstdarr'] = tstdarr
       savedatadict['tstdinit'] = inittstdarr
       savedatadict['predelays'] = gpredelays
       if saveresults>4:
         savedatadict['dpmodel_params'] = dpmodel_params
         savedatadict['signal_params'] = signal_params
         savedatadict['run_params'] = run_params
       if saveresults>5:
         savedatadict['codeinfo'] = execcodeinfo
         savedatadict['gaxondelays'] = gaxondelays
         savedatadict['goligoinfo'] = goligoinfo
       np.save(savefile, savedatadict)

  else:
    print("""This was used to save the following when saveresults=False (changed the format, so now nothing is saved):
    np.save(results_dir + '\%s-allcurves.npy' \% run_params['name'], allcurves)""")
  
  print('Done!')
  sys.exit(0)
  
  spto=arr(rparams['ospiketrain'])
  sptr=arr([rspk[0] for rspk in rparams['rspiketrain']])
  #oligoevents=rparams['oevents'][0]
  
  defalpha=0.9
  plt.figure()
  if rparams['ntotspikes'] < 5000:
    plt.plot(tt,yy[:, 1], alpha=defalpha)
    nevs=len(oligoevents)
    plt.plot(arr(oligoevents),np.ones(nevs)*zthresh*maxfaktA*0.98, 'om')
    plt.plot([tt[0],tt[-1]], 1.001*arr([maxfaktA, maxfaktA]),'g', alpha=defalpha)
    plt.plot([tt[0],tt[-1]], 1.001*arr([zthresh*maxfaktA, zthresh*maxfaktA]),'r', alpha=defalpha)
    plt.ylim(bottom=-0.01)
    plt.eventplot([spto], colors=['C2'], lineoffsets=[0.],linelengths=[0.7*maxfaktA])
    plt.eventplot([sptr-tau0], colors=['C4'], lineoffsets=[0.],linelengths=[0.7*maxfaktA])
    plt.show()
  
  plt.plot(sptr-spto)
  plt.show()
  plt.eventplot([oligoevents], colors=['C2'], lineoffsets=[0.],linelengths=[0.7])
  plt.show()
  oevisi=np.diff(oligoevents)
  plt.hist(oevisi,30)
  plt.show()
