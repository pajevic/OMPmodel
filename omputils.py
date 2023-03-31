import os, sys, time, re
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random as rnd
import numpy.random as npr
from numpy.random import rand,randn,randint,choice
from scipy.integrate import solve_ivp, odeint

verbose=0
mplotlibcolors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
floateps=[10.**(fexp)*np.finfo(float).eps for fexp in range(10)]

if 1:
  import matplotlib
  matplotlib.rcParams['font.size'] = 16
  matplotlib.rcParams['legend.fontsize'] = 12
  matplotlib.rc('xtick', labelsize=14)
  matplotlib.rc('ytick', labelsize=14) 

def tryany(x):
   try:
      rv=int(x)
   except:
     try:
        rv=float(x)
     except:
        rv=x
   return rv

def tryeval(x):
   try:
      rv=eval(x)
   except:
      rv=x
   return rv

def currentdate():
  from datetime import date, datetime
  return date.today().strftime("%d%b%Y")

def filter_eq_args(args):
    eqargs=[]
    oargs=[]
    for carg in args:
       if '=' in carg:
           eqargs.append(carg)
       else:
           oargs.append(carg)
    return eqargs, oargs
  
def lkeyboard(messg):
  print('*** Local Keyboard not used (not needed anymore). Just broadcasting the message:')
  print(messg)
  
def get_results_dir(username='pajevic',rzenv='RZHOME', return_log=False, incsep=True, verbose=0):
  rzenv=os.getenv(rzenv)
  if rzenv is None:
    if verbose: print('results directory is not defined with RZHOME! Using current!')
    rzenv='current'
  if verbose: print("RZHOME=", rzenv)
  if rzenv == 'current':
    if incsep:  rzdir = './results/'
    else: rzdir = './results'
  else:
    curdir=os.getcwd()
    lastuind=[m.start() for m in re.finditer(username,curdir)][-1]
    reldir=curdir[lastuind+7:]
    if incsep: rzdir = rzenv + reldir + '/results/'
    else: rzdir = rzenv + reldir + '/results'
  
  if not os.path.exists(rzdir): os.mkdir(rzdir)
  
  if return_log:
    logdir = rzdir.replace('/results','/logs')
    if not os.path.exists(logdir): os.mkdir(logdir)
    return rzdir, logdir
  else:
    return rzdir

def get_codeinfo(codefile, refcodefile):
  theexeccode=open(codefile).read()
  if not refcodefile: return ['', theexeccode]
  try:
    from diff_match_patch import diff_match_patch
    refexeccode=open(refcodefile).read()
    dmp = diff_match_patch()
    patches = dmp.patch_make(refexeccode, theexeccode)
    return [refcodefile, patches]
  except:
    return ['', theexeccode]

usage_demos=[
"""
print('Demo: for how to combine synced trains with time-delay vector:')
ntrains=10
nsync=5
ndelays=3
print('When having %d trains and %d delays:' % (ntrains, ndelays))
nspksig=['train%d' % (itr+1) for itr in range(nsync)]
nremtaus=ntrains-nsync
td=['td%d' % (itr+1) for itr in range(ndelays)]
for irem in range(nremtaus):
    print('Creating additional train using synctrain %s with delay %s' % (nspksig[irem%nsync], td[irem%ndelays]))
""",
"""
print('Demo: Basic Usage for ompsimulator:
Example 1:
python3 oligosimulator_pruned.py 0 dpname=omp1 lamm=0.1 ... ????? ADD THE EXAMPLE
Example 2:
"""
]

msyncdict={'s':'sync', 'i': 'indep', 'p':'pois', 'r': 'regsync', 'j': 'regind'}

pyvar2texdict={'lamh' : r'$\lambda_H$', 'lamm': r'$\lambda_M$', 'lama': r'$\lambda_A$', 'taug': r'$\tau_G$', 'taud': r'$\tau_d$', 'taud': r'$\tau_d$', 'nol':r'$N_O$', 'jitter': r'$\sigma_j$', 'nax': r'$N_A$', 'taus': r'$\tau_s$', 'fixedel' : r'$\sigma_D$', 'taumax' : r'$\tau_\mathrm{max}$', 'taumin' : r'$\tau_\mathrm{min}$', 'taunom' : r'$\tau_\mathrm{nom}$', 'tref' : '$t_R$', 'Tesec':'$T_e$', 'nreps':'$n_r$', 'nepochs':'$n_e$', 'pfjitt' : r'$\sigma_s$'}

pvar2texdict={'lamr' : r'$\lambda_R$', 'lamr' : r'$\lambda_R$', 'gav' : r'$G_{av}$', 'mps': r'$M_a$', 'sigtau': r'$\sigma_\tau$', 'taus': r'$\tau_a$', 'fdtaus': r'$D_a+\tau_a$' }

pvar_name_change_dict={'lame':'lamr'}

unitstr= {'taud':' ms', 'taug': 'ms', 'lamm':'', 'jitter': ' ms', 'lamh':'', 'nol':'', 'nax':'', 'taus':' ms', 'lama':'', 'fixedel':' ms', 'taumin':' ms', 'taumax':' ms', 'taunom':' ms', 'Tesec': ' s'}

def pvar2tex(vstr):
   if vstr in pvar2texdict: return pvar2texdict[vstr]
   elif vstr in pvar_name_change_dict: return pvar2tex(pvar_name_change_dict[vstr])
   elif vstr in pyvar2texdict: return pyvar2texdict[vstr]
   return vstr
   
def pyvar2tex(vstr, unitvar=''):
   if vstr in pyvar2texdict:
     return pyvar2texdict[vstr]
   else:
     if vstr in  pvar_name_change_dict: return pyvar2tex(pvar_name_change_dict[vstr], unitvar=unitvar)
     if '=' in vstr:
       itms=vstr.split('=')
       if unitvar:
         return pyvar2tex(itms[0]) + '=' + ''.join(itms[1:])+ unitstr[unitvar]
       else:
         if itms[0] in unitstr:
           return pyvar2tex(itms[0]) + '=' + ''.join(itms[1:]) + unitstr[itms[0]]
         else:
           return pyvar2tex(itms[0]) + '=' + ''.join(itms[1:])
     else:
       return vstr

def varstr2val(vstr, vname='fixedel'):
  if isinstance(vstr, list):
    return [varstr2val(vs,vname=vname) for vs in vstr]
  if vname=='fixedel':
    if 'nrn_' in vstr: return float(vstr[4:])
    else: return tryany(vstr)
  else:
    print('Variable %s not implemented in varstr2val' % vname)
    return None
  
# Saturation Functions

def fsquad(tau, taumin, taumax):
    vfunc=2*(tau-taumin)/(taumax-tau)
    if vfunc<0: vfunc=0
    return vfunc

# first derivative for a jacobian
def fsquadprim(tau, taumin, taumax):
    vfunc=2*(tau-taumin)/(taumax-tau)
    if vfunc<0: vfunc=0
    return vfunc

def fslinmid(tau, taumin, taumax, taunom=None):
    taunom=(taumin+taumax)/2.
    if tau>taumin and tau<taunom:
       return (tau-taumin)/(taunom-taumin)
    elif tau<taumax:
       return (taumax-tau)/(taumax-taunom)
    else:
      return 0
  
def fslin(tau, taumin, taumax, taunom):
#    if taunom is None:
#       taunom=(taumin+taumax)/2.
    if tau>taumin and tau<taunom:
       return (tau-taumin)/(taunom-taumin)
    elif tau<taumax and tau>=taunom:
       return (taumax-tau)/(taumax-taunom)
    else:
      return 0
   
# first derivative for a jacobian
def fslinprim(tau, taumin, taumax, taunom):
    if tau>taumin and tau<taunom:
       return 1/(taunom-taumin)
    elif tau<taumax and tau>=taunom:
       return -1/(taumax-taunom)
    else:
      return 0

def factor_impresp(tt, Qv, taud, taur):
     taua=(taud+taur)
     aa=taua/(taud*taur)
     bb=1./taud
     cc = Qv*taua/(taur*taud**2)
     return cc*(np.exp(-tt*bb)-np.exp(-tt*aa))/(aa-bb)
  
def coeffs_from_release_params(taur, taud, Qv):
  tausum=taud+taur
  vpeak=Qv*(tausum/taur)**(-taur/taud)/taud #(qq ((td + taus)/taus)^(-(taus/td)))/td
  tpeak=taur*np.log(tausum/taur)
  aic=tausum/(taud*taur)
  bic=1./taud
  cic = Qv*tausum/(taur*taud**2)
  return (aic,bic,cic,vpeak,tpeak)

def get_atime_spread_kgroups(atimes, kg=2, get_labels=False):
    """ characterize multi-group mixed signals (within each group, between groups, overall within group: run get_atime_spread_kgroups(kg) for help"""
    if isinstance(atimes, (str,int,np.integer)):
       if isinstance(atimes, (int,np.integer)): kg=atimes
       return [r'$\sigma_{tot}$']+["$\\sigma_%d$" % (ck+1) for ck in range(kg)]+[r"$\mu_{RS1}$",r"$\sigma_{RS1}$" ,r"$\sigma_W$",r"$\sigma_B$"]
    nax=len(atimes)
    sspread= [ np.std(atimes) ]
    if kg>1:
      if ( nax<(kg*2) ) or ( nax%kg != 0 ):
         return sspread + [None]*(kg+4)
      k1=nax//kg
      subspreads = [np.std(atimes[ik*k1:(ik+1)*k1]) for ik in range(kg)]
      subspreads += get_rndsubsample_estimate_atime(atimes, k1)
      subspreads += get_wbgroup_spread(atimes, kg)
      return sspread + subspreads
    else:
      return sspread

# OMDP model class

class OMDPmodel:
  def __init__(self, dparams=None, global_params=None, nax=None):
      """ taus inverse spike rate / ISI time interval in "tunit"s;
      """
      if dparams is None:
        dparams=def_dpmodel_params
        
      if 'ode_order' in global_params:
          odeord=global_params['ode_order']
      else: odeord=1
      if nax is None:
        if 'nax' in dparams:
           nax=dparams['nax']
        else:
           nax=10
           
      self.dparams=dparams
           
      self.init_model_params()
        
  def init_model_params(self, dparams=None):
      if dparams is None:
         dparams=self.dparams
      self.taumins=dparams['taumin']*np.ones(nax)
      self.taumaxs=dparams['taumax']*np.ones(nax)
      self.taunoms=dparams['taunom']*np.ones(nax)
      self.dtaus=np.ones(nax)
      self.drate=drate
      self.spiketrains=[]
      self.sptimes=[]
      self.itrain=[]
      self.nvars=5+nax
      self.dydts=np.zeros(self.nvars)
      
      if odeord:      
        self.dydt_func = lambda t, y: self.dydt(t,y)
      else:
        self.dydt_func = lambda y, t: self.dydt(t,y)

      self.fsfunc=lambda x: fslin(x, self.taumins[0], self.taumaxs[0], self.taunoms[0])
      self.fsfuncprim=lambda x: fslinprim(x, self.taumins[0], self.taumaxs[0], self.taunoms[0])
      self.itau=5

  def init_model_params(self):
      pass
      
  def dydt_simp1(self, t, y):
    dbasea = - self.krecov * y[0]
    du = y[2]
    dv = - self.aA * self.bA * y[1] - (self.aA + self.bA)*y[2] # + c \delta(t)
    dub = y[4]
    dvb = - self.aB * self.bB * y[3] - (self.aB + self.bB)*y[4] # + c \delta(t)
    itau=5
    self.dydts[:itau]=[dbasea, du, dv, dub, dvb, dtau]
    for iax in range(nax):
       dtau = self.drate * self.fsfunc(y[iax+itau]) # - xxx \delta(t)
       self.dydts[itau+iax]=dtau
    return dydts

  def dydt_oligo1(self, t,y):
      dbasea = - krecov * y[0]
      du = y[2]
      dv = - self.aA * self.bA * y[1] - (self.aA + self.bA)*y[2] # + c \delta(t)
      dub = y[4]
      dvb = - self.aB * self.bB * y[3] - (self.aB + self.bB)*y[4] # + c \delta(t)
      dv = - aicA * bicA * y[1] - (aicA + bicA)*y[2] # + c \delta(t)
      dub = y[4]
      dvb = - aicB * bicB * y[3] - (aicB + bicB)*y[4] # + c \delta(t)
      dratec = y[5]
      if dratec<1.e-35:
         lkeyboard('dratec problem')
      ddr=0
      for iax in range(gnax):
         ctau=y[iax+itau]
         fsval=fsfuncdemyel(ctau)
  #       ddr += dratec* 1.e-3*(otaunom-ctau)/gnax # /(1-fsval)
         ddr += dratec* lamh *(otaunom-ctau)/gnax # /(1-fsval)
         dtau = dratec * fsval # - xxx \delta(t)
         gdydts[itau+iax]=dtau
         if verbose>6:
  #          print("iax=", iax,ctau, otaumin, otaumax, ctau-otaumin, otaumax-ctau, fsval, dtau)
            print("iax=%d dratec=%g ctau=%g ntau=%g fsval=%g" % (iax, dratec, ctau, normtau(ctau), fsval))
      if verbose>5:
        print("dratec=", dratec)
        print("ddr=", ddr)
      gdydts[:itau]=[dbasea, du, dv, dub, dvb,ddr]
      return gdydts
  
  
# Results file management

def params2filename(paramabbrevs=None, prefix='results/', **params):
   if paramabbrevs is None:
      paramabbrevs={}
      paramabbrevs['ntotsig']=('nsig','%d')
      paramabbrevs['truedelay']=('td','%d'), 
      paramabbrevs['nruns']=('nruns' ,'%d')
      paramabbrevs['Tesec']=('T' ,'%d')
      paramabbrevs['taus1']=('tausa','%d')
      paramabbrevs['taus2']=('tausb','%d')
      paramabbrevs['tref']=('rf', '%d')
      paramabbrevs['actrise']=('ar','%d')
      paramabbrevs['actdecay']=('ad','%d')

   print("paramabbrevs=", paramabbrevs)
   print("params=", params)
   
   savefilename=prefix+params['rname']
   if type(params)==dict:
     for prm in paramabbrevs:
        if prm in params:
           val=params[prm]
           prmabrs=paramabbrevs[0]
           if np.isscalar(prmabrs):
             cstr='_%s%d' % (prmabrs,val)
           else:
             frmt='_%%s%s' % prmabrs[1]
             cstr=frmt % (prmabrs[0],val)
        savefilename+=cstr
   elif type(params)==tuple:
      print('Not using this anymore')
      sys.exit(-1)
   return savefilename

def Tfunc(taus1,tref=0):
    Ttaus=75*(taus1+tref) + 750
    return Ttaus
      
def generate_runs(runparams=None, defparams=None, prefix='results/', **params):
   """python oligolearn.py 20 25 actrise=1 actdecay=12 Tesec=2100 taus1=20 truedelay=5 tref=5"""
   if defparams is None:
      defparams={}
      defparams['rname']='allrunsC'
      defparams['ntotsig']=('nsig','%d', 20)
      defparams['nruns']=('nruns' ,'%d', 25)
      defparams['truedelay']=('td','%d', 2), 
      defparams['Tesec']=('T' ,'%d', Tfunc)
      defparams['taus1']=('tausa','%d', 10)
      defparams['taus2']=('tausb','%d', 10)
      defparams['tref']=('rf', '%d', 0)
      defparams['actrise']=('ar','%d', 10)
      defparams['actdecay']=('ad','%d', 10)

   if runparams is None:
      runparams={}
      runparams['truedelay']=[1,2,5,10]
      runparams['Tesec']='Tfunc'
      runparams['taus1']=[5,10,20, 50]
      runparams['tref']=[0]
      runparams['actrise']=[1,2,5,10,20,30,50]
      runparams['actdecay']=[1,2,5,10,20,30,50]

   print("runparams=", runparams)
   rnstr='python oligolearn.py 20 25'
#  for rprm in runparams:
#         for 
#         rnstr+=' %s=%s
#   """python oligolearn.py 20 25 actrise=1 actdecay=12 Tesec=2100 taus1=20 truedelay=5 tref=5"""
   fsw=open('subinewF.swarm','w')
   for ar in runparams['actrise']:
     for ad in runparams['actdecay']:
       for taus1 in runparams['taus1']:
         for td in runparams['truedelay']:
           for rf in runparams['tref']:
              Tc=Tfunc(taus1,rf)
              print("python oligolearn.py 20 25 actrise=%d actdecay=%d Tesec=%d taus1=%d truedelay=%d tref=%d" % (ar,ad,Tc,taus1,td,rf), file=fsw)
   fsw.close()     
  
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

def plot_lrn_curve(actdecay, actrise, T=50, npoints=300, showplot=True, colors='crgbkmy', **args):
    if type(actdecay)==list:
      if actrise is None:
        actrise=actdecay
        for iplt, ad in enumerate(actdecay):
           plot_lrn_curve(ad, ad, T=T, npoints=npoints, showplot=False, color=colors[iplt%7],  label='$\\tau_d$=$\\tau_r$=%d ms' % ad, **args)
      else:
          assert len(actdecay) == len(actrise)
          for iplt, (ad,ar) in enumerate(zip(actdecay,actrise)):
              plot_lrn_curve(ad, ar, T=T, npoints=npoints, showplot=False, color=colors[iplt%7],  label='$\\tau_d$=%d, $\\tau_r$=%d ms' % (ad, ar), **args)
      plt.xlabel('Time [ms]')
      plt.ylabel('$c_A$ impulse response curve')
      plt.plot([-5,0],[0,0],'k-')
      plt.xlim([-5, T])
      plt.yticks([])
      plt.legend()
      plt.show()
    else:
       xx=np.linspace(0,T,npoints)
       yy=(actdecay+actrise)*np.exp(-xx/actdecay)*(1.-np.exp(-xx/actrise))/actdecay**2
       plt.plot(xx,yy, **args)
       if showplot:
         plt.show()
         
def plot_roft(taug, T=50, npoints=300, savefile='', showplot=True, colors='crgbkmy', **args):
    if type(taug)==list:
      for iplt, tg in enumerate(taug):
        plot_roft(tg, T=T, npoints=npoints, showplot=False, color=colors[iplt%7],  label='$\\tau_G$=%d ms' % tg, **args)
      plt.xlabel('Time [ms]')
      plt.ylabel('$R(t)$ impulse response curve')
      plt.plot([-5,0],[0,0],'k-')
      plt.xlim([-5, T])
      plt.yticks([])
      plt.legend()
      if savefile:
        plt.savefig(savefile)
      plt.show()
    else:
       xx=np.linspace(0,T,npoints)
       yy=2*np.exp(-xx/taug)*(1.-np.exp(-xx/taug))/taug
       plt.plot(xx,yy, **args)
       if showplot:
         plt.show()

def plot_lrn_curve_lactives(tcur, lastactivslist, actdecay, actrise, T=50, npoints=300, message='', plotall=False):
       if plotall:
          tmin=min(lastactivslist)
       else:
          tmin=tcur-1
#         plt.eventplot([spt1,spt2], colors=['r','g'], lineoffsets=[-0.1,-0.2],linelengths=[0.1,0.1])
       plt.eventplot([lastactivslist, [tcur]], colors=['b','r'], lineoffsets=[0.1,0.09],linelengths=[0.3, 0.32])
       tt=np.linspace(tmin, tcur+T, npoints)
       ampfact=0.
       ampfacts=0.
       for lactiv in lastactivslist:
          tsnc=tcur-lactiv
          tsncs=tt-lactiv
          ampa1=np.exp(-tsnc/actdecay)*(1.-np.exp(-tsnc/actrise))
          print("ampa1=", ampa1)
          ampfact+=ampa1
          ampas=np.exp(-tsncs/actdecay)*(1.-np.exp(-tsncs/actrise))
          ampas[tsncs<0.]=0.
#          lkeyboard('tsncs ampas')
          if plotall:
             plt.plot(tt, ampas,'k--', linewidth=1)
          ampfacts += ampas
       plt.plot(tt,ampfacts,'k', linewidth=2)
       plt.plot(tcur,ampfact, 'ro', markersize=10)
       print(lastactivslist)
       plt.ylim(bottom=-0.02)
       plt.xlabel('Time t [ms]')
       plt.ylabel('$c_A(t)$')
       plt.yticks([])
       plt.title('Spike train response (last spike at  t=%g ms)' % tcur, fontsize=14)
       plt.show()

def force_tspan(tspan,spikes=None, toffs=1):
  if spikes is None:
    tmin=0
    tmax=tspan
  else:
     tmin=spikes[0][0]
     tmax=spikes[-1][0]
    
  if np.isscalar(tspan):
     if tspan<0:
        tmin=tmax+tspan
     else:
        tmax=tmin+tspan
  return [tmin-toffs,tmax+toffs] 
       

def spikes2sptrains(spikes, tspan=None, nax=None):
      if nax is None:
           iax=np.array([spk[1] for spk in spikes])
           nax=max(iax)+1
#      if tspan is None:
#        times=np.array([spk[0] for spk in spikes})
#        tspan=(min(times)-1,max(times)+1)

      if tspan is None:
         sptrains=[[spk[0] for spk in spikes if spk[1]==iax] for iax in range(nax)]
      else:
         if np.isscalar(tspan):
             if tspan<0:
               tmax=spikes[-1][0]+1
               tmin=tmax+tspan
               sptrains=[[spk[0] for spk in spikes if (spk[1]==iax and spk[0]>tmin)] for iax in range(nax)]
             else:
               tmin=spikes[0][0]-1
               tspan=(tmin, tmin+tspan)
               sptrains=[[spk[0] for spk in spikes if (spk[1]==iax and spk[0]<tmax)] for iax in range(nax)]
         else:
            sptrains=[[spk[0] for spk in spikes if (spk[1]==iax and spk[0]>tspan[0] and spk[0]<tspan[1])] for iax in range(nax)]
      return sptrains
      

def spikes2pyspktrains(spikes, tspan=None, nax=None):
      if nax is None:
           iax=np.array([spk[1] for spk in spikes])
           nax=max(iax)+1
#      if tspan is None:
#        times=np.array([spk[0] for spk in spikes})
#        tspan=(min(times)-1,max(times)+1)
      if tspan is None:      
         tmin=spikes[0][0]
         tmax=spikes[-1][0]
         tspan=(tmin-1, tmax+1)
         if verbose>2:
             print("Assigned here tspan=", tspan)
         
      if verbose>2:
          print("Final tspan=", tspan, spikes[0][0], spikes[-1][0])
      sptrains=[[spk[0] for spk in spikes if spk[1]==iax] for iax in range(nax)]
      pyspktrains=[pyspk.SpikeTrain(sptr,edges=tspan) for sptr in sptrains]
      return pyspktrains
      
def shiftspikes(spikes, tshift):
      if np.isscalar(tshift):
        return [[spk[0]+tshift, spk[1]] for spk in spikes]
      else:
        return [[spk[0]+tshift[spk[1]], spk[1]] for spk in spikes]
      
def post2prespikes(spikes, toffs=0.5):
      sptimes=np.array([spk[0] for spk in spikes])
      itrain=np.array([spk[1] for spk in spikes])
      mintime=sptimes.min()
      sptimes=sptimes-mintime+toffs
      isrt = np.argsort(sptimes)
      sptimes=sptimes[isrt]
      itrain=itrain[isrt]
      newspikes=[[spt1,it1] for spt1,it1 in zip(sptimes,itrain)]
      return newspikes
      
def plot_sptrains_v1(sptrains, tspan=None, show=True, toffset=0):
      nax=len(sptrains)
      spklen=1/nax
      for ispt,sptr in enumerate(sptrains):
          spkoffs=ispt*(spklen*1.05)
          plt.eventplot([toffset+arr(sptr)], colors=['C%d' % (ispt%5 + 1)], lineoffsets=[spkoffs], linelengths=[spklen])
      if tspan is not None:
         plt.xlim(tspan)
      plt.xlabel('time [ms]')
      if show:
        plt.show()

def plot_sptrains(sptrains, tspan=None, show=True, toffset=None, labels=None, axplt=None, hplot=None, color=None, title=''):
    if np.isscalar(tspan): tspan=force_tspan(tspan,spikes,toffs=5)
    nax=len(sptrains)
    if labels is None:
       labels=['A%d' % (iax+1) for iax in range(nax)]
    if color is None:
      colors=['C%d' % (iax%9 + 1) for iax in range(nax)]
    else:
      if type(color)==list: colors=color
      else: colors=[color for iax in range(nax)]
    if toffset:
      if toffset=='zero':
         sptrains=[arr(sptr)-tspan[0] for sptr in sptrains]
         tspan=[0, tspan[1]-tspan[0]]
      elif toffset=='neg':
         sptrains=[arr(sptr)-tspan[1] for sptr in sptrains]
         tspan=[tspan[0]-tspan[1], 0]
      else:
        sptrains=[toffset+arr(sptr) for sptr in sptrains]
      
    if hplot is None:
       if axplt is None:
         hplot=plt.eventplot(sptrains, linelengths=0.75, color=colors)
         plt.xlabel('Time [ms]', fontsize=16)
         plt.xlim(tspan)
         plt.yticks(range(nax), labels=labels, fontsize=16)
         if title: plt.title(title)
       else:
         hplot=axplt.eventplot(sptrains, linelengths=0.75, color=colors)
         axplt.set_xlabel('Time [ms]', fontsize=16)
         axplt.set_xlim(tspan)
         axplt.set_yticks(range(nax))
         axplt.set_yticklabels(labels=labels, fontsize=16)
         if title: axplt.set_title(title)
    else:
        for hplt,spos in zip(hplot,sptrains):
            if len(spos):
               hplt.set_positions(spos)
            else:
              hplt._is_horizontal=False
              hplt.set_positions(spos)
    return hplot

  
def plotspikes(spikes, tspan=None, nax=None, show=True, toffset=0, labels=None, axplt=None, hplot=None, color=None):
      if np.isscalar(tspan): tspan=force_tspan(tspan,spikes,toffs=5)
      sptrains=spikes2sptrains(spikes, tspan=tspan, nax=nax)
      hplot=plot_sptrains(sptrains, tspan=tspan, show=show, toffset=toffset, labels=labels, axplt=axplt, hplot=hplot, color=color)
      return hplot

def plot_prepost_spread(fixedels, axdelays, pltparams={}, ax=None, title='Axon delays', legends=[], **pltargs):
      if 'axloffs' in pltparams:
         axloffs=pltparams['axloffs']
      else:
         axloffs=0.5
      nax=len(fixedels)
      labels=['A%d' % (iax+1) for iax in range(nax)]
      
      totdelays=fixedels+axdelays
      
      if ax is None: fig, pax = plt.subplots()
      else: pax=ax

      mntot=totdelays.mean()
      mnstd=totdelays.std()
       
      if 'taumin' in pltparams:
         taumin=pltparams['taumin']
#         pax.plot([taumin, taumin],[-axloffs, nax+axloffs],'k--', linewidth=3, label=None)
         pax.axvline(taumin, color=[0.2,0,0.], linestyle="--", linewidth=3)
         xmint=taumin*0.25
      else: xmint=0.
      if 'taumax' in pltparams:
         taumax=pltparams['taumax']
#         pax.plot([taumax, taumax],[-axloffs, nax+axloffs],'k--', linewidth=3, label=None)
         pax.axvline(taumax, color=[0,0,0.2], linestyle="--", linewidth=3)
         xmaxt=taumax*1.05
      else: xmaxt=max(axdelays)

#      pax.barh(np.arange(nax)*1, [fixedels, axdelays], align='center', height=0.5, stacked=True, **pltargs)

      xmin=xmint
      mxtot=max(totdelays)
      if mxtot<xmaxt: pdadd=(xmaxt-mxtot)/3.
      else: pdadd=0.

      mntotplt=mntot+pdadd
      pax.axvline(mntotplt, color='k', linestyle="--", linewidth=1)
      pax.axvline(mntotplt-mnstd, color='k', linestyle=":", linewidth=0.75)
      pax.axvline(mntotplt+mnstd, color='k', linestyle=":", linewidth=0.75)
     
      xmax=max([xmaxt, mxtot])*1.03
#      df = pd.DataFrame({'pds':fixedels,'tds': axdelays})
      if isinstance(legends,list) and len(legends)==2: colnames=legends
      else: colnames=['OMMP','Pre']
      df=pd.DataFrame(np.array([axdelays,fixedels+pdadd]).T, columns=colnames)
#      dfpax = df.plot.barh(stacked = True, ax=ax, color=[[0.01,0.2,0.9],[0.8,0.2, 0.]], **pltargs);
      dfpax = df.plot.barh(stacked = True, ax=ax, color=[mplotlibcolors[0],mplotlibcolors[1]], **pltargs);
# alterantively specify color using column names:  ax = df.plot.barh(color={"Pre": "red", "OMMP": "blue"})
      if False:
        if False:
          for rowNum,row in df.iterrows():
            xpos = 0
            for val in row:
              xpos += val
              ax.text(xpos + 1, rowNum-0.05, str(val), color='black')
            xpos = 0
            
      pax.set_xlim([xmin, xmax])
      pax.set_yticks(range(nax))
      pax.set_yticklabels(labels=labels, fontsize=14)
#      dfpax.legend(['Pre-delay','Myel. delay'])
      dfpax.set(xlabel='τ [ms]', ylabel='Axons')
      if legends:
         dfpax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.05))
      else:
         dfpax.get_legend().remove()
      if title:
         pax.set_title(title, fontsize=12)

def jitter_ftaus(taus, pfjitt, list2list=True, mintr=10, minf=0.00001):
    if np.isscalar(taus):
      return 1./np.clip((1./taus)*(1.+pfjitt*randn()/100.), minf, 1./mintaus)
    elif isinstance(taus, list) and list2list:
      return [jitter_taus(taus1, pfjitt, list2list=list2list, mintaus=mintaus, minf=minf) for taus1 in taus]
    else:
      return 1./np.clip(np.array(1./taus)*(1.+pfjitt*randn(len(taus))/100.), minf, 1./mintaus)

def jitter_taus(taus, pfjitt, list2list=True, mintaus=10, minf=0):
    if minf:
      return jitter_ftaus(taus, pfjitt, list2list=list2list, mintaus=mintaus, minf=minf)
    else:
      if np.isscalar(taus):
        return np.clip(taus*(1.+pfjitt*randn()/100.), mintaus, None)
      elif isinstance(taus, list) and list2list:
        return [jitter_taus(taus1, pfjitt, list2list=list2list, mintaus=mintaus) for taus1 in taus]
      else:
        return np.clip(np.array(taus)*(1.+pfjitt*randn(len(taus))/100.), mintaus, None)
  
class SpikeTrains:
  def __init__(self, spksig=[], tspan=1000, taus=50, ntrains=1, tunit='ms', jitter=0, pfjitt=0,  jnodes=None, td=0, tref=0, kg=1, do_sort=True):
      """ taus inverse spike rate / ISI time interval in "tunit"s;
      """
      self.spikes=[] # Main format
      self.sptimes=[] # if not present, create it from .spikes
      self.itrain=[]  # must be present if sptimes is present
      self.sptrains=[] # optional (see also for get_spike_train(s) to directly obtaining
      self.ntrains=None

      if type(spksig) == str:
        nspksig=[]
        if spksig[:2] in ['po', 'in']:
          items=spksig.split('_')
          if len(items)>1: ntrs=int(items[1])
          else: ntrs=ntrains
          if len(items)>2: taus=float(items[2])
          else: taus=taus
          if items[0] in ['pois', 'indep']:
            if np.isscalar(taus):
               nspksig=[generate_rpoisson_spike_train(jitter_taus(taus, pfjitt), tspan, tref=tref, tunit=tunit) for _ in range(ntrs)]
            else:
               nspksig=[generate_rpoisson_spike_train(jitter_taus(taus1, pfjitt), tspan, tref=tref, tunit=tunit) for taus1 in taus]
        elif spksig[:2]=='re':
          if np.isscalar(tspan): tspan=(0,tspan)
          if spksig=='regsync':
            if np.isscalar(taus): nspksig= [np.arange(tspan[0],tspan[1], jitter_taus(taus,pfjitt)) for _ in range(ntrains)]
            else: nspksig= [np.arange(tspan[0],tspan[1], jitter_taus(taus1, pfjitt)) for taus1 in taus]
          else:
            if np.isscalar(taus): nspksig= [rand()*taus + np.arange(tspan[0],tspan[1], jitter_taus(taus, pfjitt)) for _ in range(ntrains)]
            else: nspksig= [rand()*taus1 + np.arange(tspan[0],tspan[1], jitter_taus(taus1, pfjitt)) for taus1 in taus]
          for inspd,nspd in enumerate(nspksig):
              nspd+=inspd*td
        elif spksig[:2]=='sy':
          items=spksig.split('_')
          if len(items)>1: nsync=int(items[1])
          else: nsync=ntrains//2
          if len(items)>2: ntrs=int(items[2])
          else: ntrs=ntrains
          if len(items)>3: td=float(items[3])
          else: td=td
          if len(items)>4: taus=float(items[4])
          else: taus=taus

          self.ntrains=ntrs

          if items[0]=='sync':
               if ntrs<nsync:
                  print('Number of trains is less than the number of synced trains (ntrains=%d nsync=%d)' % (ntrains, nsync))
                  print('Setting nsync to the number of trains!')
                  nsync=ntrs
               nspksig=[generate_rpoisson_spike_train(taus, tspan, tref=tref, tunit=tunit)]*nsync
#               nspksig=[generate_rpoisson_spike_train(jitter_taus(taus, pfjitt), tspan, tref=tref, tunit=tunit) for _ in range(nsync)]
               nremtaus=ntrs-nsync # remaining unsynced trains
               if nremtaus:
                 if np.isscalar(td):
                    td=[td]
                 ndelays=len(td)
                 nspksig.extend([nspksig[irtaus%nsync]+td[irtaus%ndelays] for irtaus in range(nremtaus)])
#                      nspksig.extend([generate_rpoisson_spike_train(taus, tspan, tref=tref, tunit=tunit) for _ in range(ntrd)])
        elif 'std_' in spksig:
          items=spksig.split('_')
          if len(items)>1: nsync=int(items[1])
          else: nsync=ntrains//2
          if len(items)>3: ntrs=int(items[3])
          else: ntrs=ntrains
          if len(items)>2:
             tds=list((np.arange(ntrs-nsync)+1)*float(items[2]))
          else:
             if isinstance(td, list):
               tds=td
             elif np.isscalar(td):
               tds=list((np.arange(ntrains-nsync)+1)*td)
          sptr=SpikeTrains('sync_%d' % nsync, tspan=tspan, taus=taus, ntrains=ntrs, tunit=tunit, jitter=jitter, pfjitt=pfjitt, jnodes=jnodes, td=tds, tref=tref, kg=kg, do_sort=do_sort)
          nspksig=sptr.get_spike_trains()
        elif spksig[:2]=='mi':
          items=spksig.split('__')
          sptr=SpikeTrains(items[1], tspan=tspan, taus=taus, ntrains=ntrains, tunit=tunit, jitter=jitter, pfjitt=pfjitt, jnodes=jnodes, td=td, tref=tref, kg=kg, do_sort=do_sort)
          for spdat1 in items[2:]:
            sptr1=SpikeTrains(spdat1, tspan=tspan, taus=taus, ntrains=ntrains, tunit=tunit, jitter=jitter, pfjitt=pfjitt, jnodes=jnodes, td=td, tref=tref, kg=kg, do_sort=do_sort)
            sptr.add(sptr1)
          nspksig=sptr.get_spike_trains()

        if nspksig:
            self.__init__(nspksig, tspan=tspan, taus=taus, ntrains=ntrains, tunit=tunit, jitter=jitter, pfjitt=pfjitt, jnodes=jnodes, td=td, tref=tref, do_sort=do_sort)

            return
      elif type(spksig) == list:
          if len(spksig)==0:
             self.ntrains=0
             tspan=0
          elif type(spksig[0])==np.ndarray or type(spksig[0])==list:
            self.ntrains=len(spksig)
            for itr,ctrain in enumerate(spksig):
               self.sptrains.append(ctrain)
               self.spikes.extend( [ [sptm,itr] for sptm in ctrain])
#               self.sptimes.extend([sptm for sptm in ctrain])
      elif isinstance(spksig, SpikeTrains):
           self.spikes=list(spksig.spikes)
           self.ntrains=spksig.ntrains
      elif type(spksig) == np.ndarray:
         if spksig.ndim==1:
            nspikes=len(spksig)
            self.ntrains=1
            self.sptrains[itr]=ctrain
            self.spikes=[[spt,0] for isp in spksig]
#            self.sptimes=list(spksig)
         elif spksig.ndim==2:
            ntrains,nspikes=spksig.shape
            self.ntrains=ntrains
            for itr in range(ntrains):
               self.sptrains[itr]=list(spksig[itr])
               self.spikes.extend([[spd,itr] for spd in spksig[itr]])
#               self.sptimes.extend([spd for spd in spksig[itr]])
         else:  # spksig.ndim>=3:
            print('Not implemented yet for nd-arrays.ndim > 2!')
      elif have_pyspike and isinstance(spksig, pyspk.SpikeTrain):
           pass
      else:
        print('UNRECOGNIZED Data type for spksig=', spksig)
        return None

      self.update_spikes()
      
      if jitter:
        print("jitter=", jitter)
        self.jitterspikes(jitter)

      if tspan is None:
          self.tspan=[min(self.sptimes), max(self.sptimes)]
      elif np.isscalar(tspan):
           try: self.tspan=(0,float(tspan))
           except:
              print('Invalid tspan value')
      else:
        self.tspan=tspan
        
  def update_spikes(self):
      self.spikes.sort(key=lambda x: x[0])
      self.ntotspikes=len(self.spikes)
      self.spikes2sptit()
      self.sptrains=[]
      if self.ntrains is None: self.ntrains=max(self.itrain)+1
      
  def shift(self, tshift):
      if np.isscalar(tshift):
        self.spikes = [[spk[0]+tshift, spk[1]] for spk in spikes]
      else:
        self.spikes = [[spk[0]+tshift[spk[1]], spk[1]] for spk in spikes]
      self.update_spikes()
      
  def create_sptrains(self):
      self.sptrains=spikes2sptrains(self.spikes, nax=self.ntrains)
      
  def subtrains(self, itrains, keepinds=False):
    if keepinds:
       spt=SpikeTrains([])
    else:
       if len(self.sptrains)==0:
          self.create_sptrains()
       return SpikeTrains([self.sptrains[itr] for itr in itrains])

  def __repr__(self):
      return "SpikeTrains: nspikes=%d ntrains=%d: #1(%.2f %d) ... #%d(%.2f %d)" % (self.ntotspikes, self.ntrains, self.spikes[0][0], self.spikes[0][1], self.ntotspikes, self.spikes[-1][0], self.spikes[-1][1])

  def __str__(self):
       return "spikes="+str(self.spikes)
    
  def add(self, sptr2, mode=0, inplace=False, spanall=True, toffs=1):
      """ add spike trains
          if mode==0:
              When adding lists of trains  assumes a list of 2-elem lists are [spiketime,itrain] spikes
              When  adding SpikeTrains assumes that the train indices are unrelated (index 0 is not the same as index 0 in the second train)             Basically, all trains in sptr2 are treated as new trains
      """

      if inplace:
         nspt=self
      else:
         nspt=SpikeTrains()
         
      if type(sptr2)==list:
        telem=type(sptr2[0])
        if telem==tuple or  telem==list:
          if telem==tuple or ( len(sptr2[0])==2 and mode==0):
            self.spikes+=sptr
            self.ntrains=None
          else:
            for itr,tr in enumerate(sptr2):
              self.spikes.extend([[tpt,itr+self.ntrains] for tpt in tr])
            self.ntrains + len(sptr2)
        else:
          self.spikes.extend([[tpt,self.ntrains] for tpt in sptr2])
          self.ntrains += 1
      else: # assume you are dealing with another instance of SpikeTrains
         if mode==1:
            # assuming train indices are the same in both
            self.spikes += sptr2.spikes
         else:
            # assuming trains are independent; just add new spike trains
            if len(self.sptimes)==0:
                 self.spikes2sptit()
            self.sptimes = np.concatenate([self.sptimes, sptr2.sptimes])
            self.itrain = np.concatenate([self.itrain, sptr2.itrain+self.ntrains])
#            isrt = np.argsort(self.sptimes)
#            self.sptimes=self.sptimes[isrt]
#            self.itrain=self.itrain[isrt]
            self.sptit2spikes()
            self.ntrains +=sptr2.ntrains
            if len(self.sptrains):
               print('UPDATING sptrains !!!')
               self.create_sptrains()
            
      self.update_spikes()
      
      if spanall or (self.tspan[0]+floateps[1] >= self.tspan[1]):
         self.tspan=[self.spikes[0][0] - toffs, self.spikes[-1][0] + toffs]
         
      return self
      
  def jitterspikes(self, jittsigma):
      if verbose>5:
        print("jittsigma=", jittsigma)
        print("self.ntotspikes=", self.ntotspikes)
        print("type(self.sptimes)=", type(self.sptimes))
      self.sptimes += npr.randn(self.ntotspikes)*jittsigma
      isrt = np.argsort(self.sptimes)
      self.sptimes=self.sptimes[isrt]
      self.itrain=self.itrain[isrt]
      self.sptit2spikes()
                            
  def spikes2sptit(self):
      self.sptimes=np.array([spk[0] for spk in self.spikes])
      self.itrain=np.array([spk[1] for spk in self.spikes])
      
  def sptit2spikes(self):
      self.spikes=[[spt1,it1] for spt1,it1 in zip(self.sptimes,self.itrain)]
    
  def set_multispikes(self, mindiff=1e-7):
        sprev=self.tspan[0]-1
        newspikes=[]
        for spk in self.spikes:
           if (spk[0]-sprev)<mindiff:
             newspikes[-1].append(spk[1])
           else:
             newspikes.append(spk)
             sprev=spk[0]
        self.spikes=newspikes
                            
  def get_spike_train(self, itrain):      
       return [spk[0] for spk in self.spikes if spk[1]==itrain]
     
  def get_spike_trains(self):      
       return spikes2sptrains(self.spikes, tspan=self.tspan, nax=self.ntrains)
     
  def get_spike_times(self, itrain=None):
      if itrain is None:
         return [spk[0] for spk in self.spikes]
      else:
         return [spk[0] for spk in self.spikes if spk[1]==itrain]
  
  def plot_spikes(self, title='', show=True):
      sptrains=spikes2sptrains(self.spikes, tspan=self.tspan, nax=self.ntrains)
      plot_sptrains(sptrains, tspan=self.tspan)
      if title: plt.title(title)
      if show:
        plt.show()
     
  def plot_isi(self, itrain=None, nbins=30):      
        isis=np.diff(self.get_spike_times(itrain=itrain))
        print("mean ISI=", isis.mean())
        plt.hist(isis,nbins, density=True)
        plt.xlabel('ISI [ms]')
        plt.ylabel('p(ISI)')
        plt.show()
       
def generate_rpoisson_spike_train(taus, tspan, tref=0, tunit='ms'):
    try:
      tmin=tspan[0]
      Tmax=tspan[1]
    except:
      tmin=0
      Tmax=tspan
    Tms=Tmax-tmin
    nspikes1=int(Tms/(taus+tref))
    stdnsp=np.sqrt(nspikes1)
    nspikes1 += int(0.01*stdnsp)
    vals=npr.exponential(taus, nspikes1)+tref
    spt=np.cumsum(vals)+tmin
    spmx1=spt.max()
    while spmx1<Tmax:
       if verbose>2:
         print('Generating more spikes!')
       nspikes1 += int(2.*stdnsp)
       vals = np.concatenate((vals,npr.exponential(taus, int(2*stdnsp))+tref))
       spt=np.cumsum(vals)
       spmx1=spt.max()
    return spt

       
def generate_spike_train(taus, Tms, tref=0):
    nspikes1=int(Tms/(taus+tref))
    stdnsp=np.sqrt(nspikes1)
    nspikes1 += int(0.01*stdnsp)
    vals=npr.exponential(taus, nspikes1)+tref
    spt=np.cumsum(vals)
    spmx1=spt.max()
    while spmx1<Tms:
       if verbose>2:
         print('Adding more spikes!')
       nspikes1 += int(2.*stdnsp)
       vals = np.concatenate((vals,npr.exponential(taus, int(2*stdnsp))+tref))
       spt=np.cumsum(vals)
       spmx1=spt.max()
    return spt

def quick_generate_spiketrain(spksig, nax, tspan=1000, taus=50, tref=0,jitter=0, pfjitt=0):
# (self, spksig=[], tspan=None, taus=None, ntrains=1, tunit='ms', jitter=0, jnodes=None, td=0, tref=0, do_sort=True):
     taus1=taus
     if spksig=='sync' or spksig=='syncall':
        nsynced=nax
        tds=[]
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, td=tds, jitter=jitter, pfjitt=pfjitt)
     elif spksig=='psync_half':
        nsynced=nax//2
        spta=SpikeTrains('mix__sync__pois', taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
     elif spksig=='ksynctwo':
        nsynced=nax//2
        spta=SpikeTrains('mix__sync__sync', taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
     elif 'ksync_' in spksig:
        nkg=int(spksig.split('_')[1])
        nsynced=int(np.round(nax/nkg))
        if nax%nkg:
           print('WARNING: The number of axons WILL change since %d axons cannot be split evenly in %d groups' % (nax, nkg))
           print('WARNING: The new number of axons is %d!!!' % (nkg*nsynced))
        spta=SpikeTrains('sync', taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        for ikg in range(1,nkg):
          spt1=SpikeTrains('sync', taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
          spta.add(spt1)
     elif 'msync_' in spksig:
        syncstr=spksig.split('_')[1]
        nkg=len(syncstr)
        nsynced=int(np.round(nax/nkg))
        if nax%nkg:
           print('WARNING: The number of axons WILL change since %d axons cannot be split evenly in %d groups' % (nax, nkg))
           print('WARNING: The new number of axons is %d!!!' % (nkg*nsynced))
        spta=SpikeTrains(msyncdict[syncstr[0]], taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        for ikg in range(1,nkg):
          spt1=SpikeTrains(msyncdict[syncstr[ikg]], taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
          spta.add(spt1)
     elif 'krates_' in spksig:
        nkg=int(spksig.split('_')[1])
        nsynced=int(np.round(nax/nkg))
        if nax%nkg:
           print('WARNING: The number of axons WILL change since %d axons cannot be split evenly in %d groups' % (nax, nkg))
           print('WARNING: The new number of axons is %d!!!' % (nkg*nsynced))
        spta=SpikeTrains('sync', taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        for ikg in range(1,nkg):
          spt1=SpikeTrains('sync', taus=(1+ikg)*taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
          spta.add(spt1)
     elif 'kratespois_' in spksig:
        nkg=int(spksig.split('_')[1])
        nsynced=int(np.round(nax/nkg))
        if nax%nkg:
           print('WARNING: The number of axons WILL change since %d axons cannot be split evenly in %d groups' % (nax, nkg))
           print('WARNING: The new number of axons is %d!!!' % (nkg*nsynced))
        spta=SpikeTrains('pois', taus=taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)

                       
        for ikg in range(1,nkg):
          spt1=SpikeTrains('pois', taus=(1+ikg)*taus1, ntrains=nsynced, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
          spta.add(spt1)
     elif spksig=='case_td1':
        nsynced=nax//2
        tds=list((np.arange(nax-nsynced)+1)*5.)
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, td=tds, pfjitt=pfjitt)
     elif spksig=='case_td2':
        nsynced=nax//5
        tds=list((np.arange(nax-nsynced)+1)*3.)
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, td=tds, pfjitt=pfjitt)
     elif spksig=='case_td3':
        nsynced=nax//2
        tds=list((-np.arange(nax-nsynced)+1)*5.)
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, td=tds, pfjitt=pfjitt)
     elif 'case_syncj' in spksig:
        jttr=float(spksig[10:])
        print("jttr=", jttr)
        nsynced=nax
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, jitter=jttr, pfjitt=pfjitt)
     elif 'case_sync1' in spksig:
        nsynced=nax
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, jitter=0, pfjitt=pfjitt)
     elif spksig=='case_pois1' or spksig=='pois':
         spta=SpikeTrains('pois', taus=taus1, ntrains=nax, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
     elif spksig=='case_trtd1':
        nsynced=5
     elif spksig=='case_trpois1':
        nsynced=nax//2
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax//2, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        sptb=SpikeTrains('pois', taus=(taus1/2), ntrains=nax//2, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        spta.__add__(sptb)
     elif spksig=='case_syncpois1':
        nsynced=nax//2
        spta=SpikeTrains('sync_%d' % nsynced, taus=taus1, ntrains=nax//2, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        sptb=SpikeTrains('pois', taus=taus1, ntrains=nax//2, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
        spta.__add__(sptb)
     elif spksig=='case_pois2':
        tds=[] # tds were missing in this implementation!? Check!
        spta=SpikeTrains('pois_sync_%d' % nsynced, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, td=tds, pfjitt=pfjitt)
     else:
        spta=SpikeTrains(spksig, taus=taus1, ntrains=nax, tspan=tspan, tref=tref, jitter=jitter, pfjitt=pfjitt)
     return spta
   
#font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 16}
#matplotlib.rc('font', **font)

filenamefmt_old='results/allruns-nsig%d_td%d_nruns%d_T%d_ar%d_ad%d.npy'
filenamefmt='results/%s-nsig%d_td%d_nruns%d_T%s_tausa%d_rf%d_ar%d_ad%d.npy'

def plotfig1(rpref='allrunsA', nsig=20, td=2, nruns=25, T='1500', tausa=10, rf=0, ar=[5,10,15,20,25,30,40,50], ad=None):
  matplotlib.rcParams['font.size'] = 16
  matplotlib.rcParams['legend.fontsize'] = 12
  matplotlib.rc('xtick', labelsize=14) 
  matplotlib.rc('ytick', labelsize=14) 
#  matplotlib.rcParams['figure.titlesize'] = 'medium
  plt.figure(figsize=(8,6))
  colrs=list('krgbcm')+[[0.9,0.6,0.1],[0.5,0.5,0.5], [0.1,0.5,0.9],[0.4,0.9,0.2],[0.,0.5,0.4],[0.2,0.9,1.],[1,1,0],[0.8,0.8,0.8]]
  xvals=np.linspace(0,1.,nsig+1)
#  legendstr=[]
  for i1,ar1 in enumerate(ar):
      clr=colrs[i1]
#      legendstr.append(lbl)
      if ad is None:
         ad1=ar1
      else:
        if np.isscalar(ad):
          ad1=ad
        elif len(ad)==len(ar):
          ad1=ad[i1]
        elif len(ad)==1:
          ad1=ad[0]
        else:
          print("ar and ad must have same length!")
          sys.exit(-1)
          
      lbl='$a_r$=%d $a_d$=%d ms' % (ar1,ad1)
      filename= filenamefmt % (rpref, nsig, td, nruns, T, tausa, rf, ar1, ad1)
      lrndata, cprms, stprms, hlprms = mynpload(filename)
      diffs=lrndata[:,0,:] - lrndata[:,1,:]                                             
      mndiff=diffs.mean(0)                                                    
      sddiff=diffs.std(0)                                                     
      plt.plot(xvals, mndiff, color=clr, linewidth=2, label=lbl)
      plt.plot(xvals, mndiff+sddiff, '--', color=clr, linewidth=0.5)
      plt.plot(xvals, mndiff-sddiff, '--', color=clr, linewidth=0.5)

  plt.ylabel('Induced delay per node ['+r'$\mu s$]')
  plt.xlabel('Delayed Signal Fraction')
#  plt.legend(legendstr)
  plt.legend()
  plt.savefig('figure1a.pdf')
  plt.show()

def plotfig_rfs(rpref='allrunsA', nsig=20, td=2, nruns=25, T=1500, tausa=10, rf=[0,1,2,5], ar=[5,10,15,20,25,30,40,50], ad=None):
  matplotlib.rcParams['font.size'] = 16
  matplotlib.rcParams['legend.fontsize'] = 12
  matplotlib.rc('xtick', labelsize=14) 
  matplotlib.rc('ytick', labelsize=14) 
  plt.figure(figsize=(8,6))
  colrs=list('krgbcm')+[[0.9,0.6,0.1],[0.5,0.5,0.5], [0.1,0.5,0.9],[0.4,0.9,0.2],[0.,0.5,0.4],[0.2,0.9,1.],[1,1,0],[0.8,0.8,0.8]]
  xvals=np.linspace(0,1.,nsig+1)
  for i1,rf1 in enumerate(rf):
      clr=colrs[i1]
      lbl=''
      if np.isscalar(ar):
          ar1=ar
      elif len(ar)==1:
          ar1=ar[0]
      else:
          ar1=ar[i1]
          lbl+='ar=%d ms ' % ar1
          
      if ad is None:
         ad1=ar1
      else:
        if np.isscalar(ad):
          ad1=ad
        elif len(ad)==len(ar):
          ad1=ad[i1]
        elif len(ad)==1:
          ad1=ad[0]
        else:
          print("ar and ad must have same length!")
          sys.exit(-1)
      lbl+='rf=%d ' % rf1
      filename= filenamefmt % (rpref, nsig, td, nruns, T, tausa, rf1, ar1, ad1)
      lrndata, cprms, stprms, hlprms = mynpload(filename)
      diffs=lrndata[:,0,:] - lrndata[:,1,:]                                             
      mndiff=diffs.mean(0)                                                    
      sddiff=diffs.std(0)                                                     
      plt.plot(xvals, mndiff, color=clr, linewidth=2, label=lbl)
      plt.plot(xvals, mndiff+sddiff, '--', color=clr, linewidth=0.5)
      plt.plot(xvals, mndiff-sddiff, '--', color=clr, linewidth=0.5)

  plt.ylabel('Induced delay per node ['+r'$\mu s$]')
  plt.xlabel('Delayed Signal Fraction')
  plt.legend()
  plt.savefig('oldfigure1a.pdf')
  plt.show()

def plotfig_any(figparams,  xprm='nsig', plotmode='outer'):
  "rpref='allrunsA', nsig=20, td=2, nruns=25, T=1500, tausa=[10,20,50], rf=[0,1,2,5], ar=[5,10,15,20,25,30,40,50], ad=None):"
  varsymbols={'ar': r'$\tau_r$', 'ad' : r'$\tau_d$', 'tausa': '$r_s$', 'rf': '$t_R$'}
  units={'ar': 'ms', 'ad' : 'ms', 'tausa': 'ms', 'T':'ms', 'td':'ms', 'rf':'ms'}
  runparams=list(figparams.keys())
  matplotlib.rcParams['font.size'] = 16
  matplotlib.rcParams['legend.fontsize'] = 12
  matplotlib.rc('xtick', labelsize=14) 
  matplotlib.rc('ytick', labelsize=14) 
#  matplotlib.rcParams['figure.titlesize'] = 'medium
  plt.figure(figsize=(8,6))
  colrs=list('krgbcm')+[[0.9,0.6,0.1],[0.5,0.5,0.5], [0.1,0.5,0.9],[0.4,0.9,0.2],[0.,0.5,0.4],[0.2,0.9,1.],[1,1,0],[0.8,0.8,0.8], rand(3), rand(3), rand(3)]
  
  rpref=figparams['rpref']
  nsig=figparams['nsig']
  
  xprmval=figparams[xprm]
  if xprm=='nsig':
     xvals=np.linspace(0,1.,xprmval+1)
  else:
     lkeyboard('Unknown xprm %s' % xprm)
     
  fixedparams=[]
  fixedparamvals=[]
  varparams=[]
  varparamvals=[]
  varparamslen=[]
  for rprm in runparams:
    if rprm==xprm: continue
    if rprm=='rpref': continue
    vprm=figparams[rprm]
    print("rprm=", rprm)
    print("vprm=", vprm)
    if np.isscalar(vprm) or vprm is None:
       fixedparams.append(rprm)
       fixedparamvals.append(vprm)
    elif type(vprm)==list:
       lenlist=len(vprm)
       if lenlist==1:
          fixedparams.append(rprm)
          fixedparamvals.append(vprm)
       else:
          varparams.append(rprm)
          varparamvals.append(vprm)
          varparamslen.append(lenlist)
    else:
        lkeyboard('For now vprm has to be a scalar or a list or None')
         
  ncurves=0
  nvparams=len(varparams)
  
  if plotmode=='same':
    if len(np.unique(varparamslen))==1:
       ncurves=varparamslen[0]
    else:
      lkeyboard('For plotmode="same" lengths of vparams should be the same')
    pltvparams=np.array(varparamvals).T
  elif plotmode=='outer':
      ncurves=np.prod(varparamslen)
      pltvparams = list(itertools.product(*varparamvals))

  labelist=['' for _ in range(ncurves)]

  titlestring=''

  for fprm,fprmval in zip(fixedparams, fixedparamvals):
       if (fprmval is not None):
          if type(fprmval)==str and ('*' in fprmval):
            print('For now not including "starred"  (*) parameters in the title string')
          else:
            if fprm in varsymbols:
              vname=varsymbols[fprm]
            else:
              vname=fprm

            if fprm in units:
               titlestring+='%s=%s %s, ' % (vname, fprmval, units[fprm])
            else:
               titlestring+='%s=%s, ' % (vname, fprmval)

  print("titlestring=", titlestring)
  if titlestring[-2:]==', ': titlestring=titlestring[:-2]
  print("titlestring after=", titlestring)
  
  rvnames=[]
  for vprm in varparams:
      if vprm in varsymbols:
         vname=varsymbols[vprm]
      else:
          vname=vprm
      rvnames.append(vname)

  for icrv in range(ncurves):
      clr=colrs[icrv]
      if 'ar' in varparams:
        iwprm=varparams.index('ar')
        vname=rvnames[iwprm]
        ar1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d ms, ' % (vname, ar1)
      elif 'ar' in fixedparams:
        iar=fixedparams.index('ar')
        ar1=fixedparamvals[iar]
      else:
        print("Parameter ar not specified!")
      
      if 'ad' in varparams:
         iwprm=varparams.index('ad')
         vname=rvnames[iwprm]
         ad1=pltvparams[icrv][iwprm]
         labelist[icrv]+='%s=%d ms, ' % (vname, ad1)
      elif 'ad' in fixedparams:
         iad=fixedparams.index('ad')
         ad1=fixedparamvals[iad]
         if ad1 is None:
            ad1=ar1
            labelist[icrv]=labelist[icrv].replace('ar=', 'ar=ad=')
      else:
         print("Parameter ad not specified!")
        
      if 'rf' in varparams:
        iwprm=varparams.index('rf')
        vname=rvnames[iwprm]
        rf1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d ms, ' % (vname, rf1)
      elif 'rf' in fixedparams:
        irf=fixedparams.index('rf')
        rf1=fixedparamvals[irf]
      else:
        print("Parameter rf not specified!")
        
      if 'td' in varparams:
        iwprm=varparams.index('td')
        vname=rvnames[iwprm]
        td1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d ms, ' % (vname, td1)
      elif 'td' in fixedparams:
        itd=fixedparams.index('td')
        td1=fixedparamvals[itd]
      else:
        print("Parameter td not specified!")
        
      if 'nruns' in varparams:
        iwprm=varparams.index('nruns')
        vname=rvnames[iwprm]
        nruns1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d, ' % (vname, nruns1)
      elif 'nruns' in fixedparams:
        inruns=fixedparams.index('nruns')
        nruns1=fixedparamvals[inruns]
      else:
        print("Parameter nruns not specified!")
        
      if 'tausa' in varparams:
        iwprm=varparams.index('tausa')
        vname=rvnames[iwprm]
        tausa1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d Hz, ' % (vname, tausa1)
      elif 'tausa' in fixedparams:
        itausa=fixedparams.index('tausa')
        tausa1=fixedparamvals[itausa]
      else:
        print("Parameter tausa not specified!")
        
      if 'T' in varparams:
        iwprm=varparams.index('T')
        vname=rvnames[iwprm]
        T1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d ms, ' % (vname, T1)
      elif 'T' in fixedparams:
        iT=fixedparams.index('T')
        T1=fixedparamvals[iT]
      else:
        print("Parameter T not specified!")
        
      if 'xxxp' in varparams:
        iwprm=varparams.index('xxxp')
        vname=rvnames[iwprm]
        xxxp1=pltvparams[icrv][iwprm]
        labelist[icrv]+='%s=%d, ' % (vname, xxxp1)
      elif 'xxxp' in fixedparams:
        ixxxp=fixedparams.index('xxxp')
        xxxp1=fixedparamvals[ixxxp]
      else:
        pass
#        print("Parameter xxxp not specified!")
          
#      lbl+='rf=%d ' % rf1
      lbl=labelist[icrv]
      if lbl[-2:]==', ': lbl=lbl[:-2]
      filename= filenamefmt % (rpref, nsig, td1, nruns1, T1, tausa1, rf1, ar1, ad1)
      if '*' in filename:
        files=glob.glob(filename)
        if len(files)==1:
           filename=files[0]
        else:
          print('Found none or multiple files that match criterion:')
          print("Pattern=", filename)
          print("Found files:", files)
          lkeyboard('Inspect files')
          
      lrndata, cprms, stprms, hlprms = mynpload(filename)
      diffs=lrndata[:,0,:] - lrndata[:,1,:]                                             
      mndiff=diffs.mean(0)                                                    
      sddiff=diffs.std(0)
      mndiff = 100.*mndiff/5.
      sddiff  = 100.*sddiff/5.
      plt.plot(xvals, mndiff, color=clr, linewidth=2, label=lbl)
      plt.plot(xvals, mndiff+sddiff, '--', color=clr, linewidth=0.5)
      plt.plot(xvals, mndiff-sddiff, '--', color=clr, linewidth=0.5)

#  plt.ylabel('Induced delay per node ['+r'$\mu s$]')
  plt.ylabel('Induced delay change [%]')
  plt.xlabel('Delayed Signal Fraction')
#  plt.legend(legendstr)
  plt.legend()
  plt.title(titlestring)
  plt.savefig('figure1a.pdf')
  plt.show()
  
def viewresults(filename='allruns-nsig10-T100.npy'):
   lrndata, cprms, stprms, hlprms = mynpload(filename)
#   lrndata = mynpload(filename)
#   lkeyboard('aha')
   diffs=lrndata[:,0,:] - lrndata[:,1,:]                                             
   mndiff=diffs.mean(0)                                                    
   sddiff=diffs.std(0)                                                     
   plt.plot(mndiff,'k')                                                    
   plt.plot(mndiff+sddiff,'g')                                             
   plt.plot(mndiff-sddiff,'r')
   plt.show()

def calculate_spikedistances_spread(prespikes, sspread, Tms, nax, sprrfact=20, ndists=2, ndescomps=3, params=[]):
  
            comptrains=range(max([1,nax-ndescomps]),nax)
            ncomptrains=len(comptrains)
            spikedistances=[[] for _ in range(ndists)]
            spikedistmats=[[] for _ in range(ndists)]
            
            spikedistnames=['SPKD(0,%d)' % ctr for ctr in comptrains]
            spikedistnames.append('spread/%d' % round(sprrfact))

            if have_pyspike: # calculate spike distances at the very start (before signals pass a single oligodendrocytes)
              pyspktrains=spikes2pyspktrains(prespikes)
              Tmx=max([max(pyspktrains[itr1].spikes) for itr1 in range(len(pyspktrains))])
              if Tmx<Tms:
                if params:
                   mstr='First: irep=%d isrep=%d Tmax<Tms Tmax=%g vs Tms=%g' % (params[0], params[1], Tmx, Tms)
                else:
                   mstr='First: Tmax<Tms Tmax=%g vs Tms=%g' % (Tmx, Tms)

#                if logcnt2<flogmax2:
                if params[2][1] < params[3][1]:
                   params[2][1]+=1
                   print(mstr, file=params[4])
                   print("params[2]=", params[2])

                  
              Tmxd=min(Tms,Tmx)
              try:
                 cspkdists = [pyspk.spike_distance(pyspktrains[0], pyspktrains[itrcomp],  interval=(0, Tmxd)) for itrcomp in comptrains]
              except:
                 lkeyboard('Failed calculating cspdists ')

              cspkdists.append(sspread/sprrfact)
              spikedistances[0].append(cspkdists)
              if verbose>2:
                 print("sspread=", sspread)
                 print("cspkdists=", cspkdists)
              try:
                 cspkdists = [pyspk.isi_distance(pyspktrains[0], pyspktrains[itrcomp],  interval=(0, Tmxd)) for itrcomp in comptrains]
              except:
                 cspkdists = [pyspk.isi_distance(pyspktrains[0], pyspktrains[itrcomp],  interval=(0, Tmx)) for itrcomp in comptrains]
                 lkeyboard('pyspk problem 2')
  
              cspkdists.append(0.)
              spikedistances[1].append(cspkdists)

              return [spikedistances,spikedistnames,spikedistmats,comptrains]

def update_spikedistances_spread(spikedistances, spikedistmats, comptrains, prespikes, sspread, Tms, sprrfact=20, params=[]):
                                            
                   pyspktrains=spikes2pyspktrains(prespikes)

                   Tmx=max([max(pyspktrains[itr1].spikes) for itr1 in range(len(pyspktrains))])
                   if Tmx<Tms:
                     if params:
                       mstr='Second: irep=%d isrep=%d Tmax<Tms Tmax=%g vs Tms=%g' % (params[0], params[1], Tmx, Tms)
                     else:
                       mstr='Second: Tmax<Tms Tmax=%g vs Tms=%g' % (Tmx, Tms)
                       # if logcnt2<flogmax2:
                     if params[2][1] < params[3][1]:
                         params[2][1]+=1
                         print(mstr, file=params[4])
#                  print(mstr, file=flog)

                   Tmxd=min(Tms,Tmx)
                   if verbose>4:
                      print("comptrains=", comptrains)
                      print("Tms=", Tms)
                   try:
                      cspkdists = [pyspk.spike_distance(pyspktrains[0], pyspktrains[itrcomp],  interval=(0, Tmxd)) for itrcomp in comptrains]
                   except:
                      lkeyboard('pyspk problem 3')
                      cspkdists = [-0.1 for itrcomp in comptrains]
                   cspkdists.append(sspread/sprrfact)
                   spikedistances[0].append(cspkdists)
                   try:
                      cspkdists = [pyspk.isi_distance(pyspktrains[0], pyspktrains[itrcomp],  interval=(0, Tmxd)) for itrcomp in comptrains]
                   except:
                      lkeyboard('pyspk problem 4')
                      cspkdists = [-0.1 for itrcomp in comptrains]
                   cspkdists.append(0.)
                   spikedistances[1].append(cspkdists)
                   
                   try:
                      spikedistmats[0]=pyspk.spike_distance_matrix(pyspktrains, interval=(0, Tmxd))
                      spikedistmats[1]=pyspk.isi_distance_matrix(pyspktrains, interval=(0, Tmxd))
                   except:
                     spikedistmats[0] =np.zeros((len(pyspktrains), len(pyspktrains)))-1.
                     spikedistmats[1] =np.zeros((len(pyspktrains), len(pyspktrains)))-1.
                     lkeyboard('   Problem with spikedistmatrices 2')

            
def robus_fit(xvals, yvals, params=None):
  pass

def assign_normalized_fixedels(nvals, kg, sigma, mnshift=None, normall=False):
  if kg > 1:
    if nvals%kg:
      raise ValueError('nvals=%d must be congruent with k=%d!' % (nvals, kg))

  if np.isscalar(sigma): sigmas=[sigma for _ in range(kg)]
  else: sigmas=sigma
  
  fixedels=[]
  if verbose>4: print("sigmas=", sigmas)

  if sigmas:
    k1 = nvals//kg
    if mnshift is None:
      mnshift=sigmas[0]*(2 + rand()/100) # make them positive / doesn't really matter how; just for esthetics; the mean is irrelevant
    for ikg in range(kg):
      pred1=randn(k1)
      pred1=sigmas[ikg]*pred1/pred1.std()
      fixedels.append(pred1)
    fixedel=np.concatenate(fixedels)
    if np.isscalar(sigma) and normall:
      fixedel = sigma * fixedel/ fixedel.std()
  else:
    fixedel=randn(nvals)
    fixedel = sigma * fixedel/ fixedel.std()
    if mnshift is None:
      mnshift=sigma*(2 + rand()/100)
      
  return mnshift + fixedel - fixedel.min()

def get_wbgroup_spread(atimes, grpinds, ddof=0):
   nax=len(atimes)
   mtot=atimes.mean()
   Swp = 0.
   Sbp = 0.
   if isinstance(grpinds,int):
      kg=grpinds
      if kg==0:
        return atimes.std(), 0.
      k1=nax//kg
      gmns=np.zeros(kg)
      for ikg in range(kg):
        atis=atimes[ikg*k1:(ikg+1)*k1]
        mi=atis.mean()
        vari  = np.var(atis, ddof=ddof)
        Swp += k1 * vari
        Sbp += k1 *(mi-mtot)**2
   else:
     glabls = np.unique(grpinds)
     kg = len(glabls)
#     mis=np.zeros(kg)
#     nis=np.zeros(kg,dtype=int)
     for  igrp in range(kg):
       gids = np.where(grpinds == glabls[igrp])[0]
       atis = atimes[gids]
       mi=atis.mean()
       ni=len(atis)
#       nis[igrp]=ni
       vari  = np.var(atis, ddof=ddof)
#       Swp += (ni-1) * vari
       Swp += ni * vari
       Sbp += ni *(mi-mtot)**2
   return np.sqrt(Swp/(nax-ddof)), np.sqrt(Sbp/(nax-ddof))

def get_rndsubsample_estimate_atime(atimes, ksub, ntrials=4000):
   stds=np.zeros(ntrials)
   nax=len(atimes)
   for itr in range(ntrials):
      satimes=atimes[np.random.choice(nax, ksub, False)]
      stds[itr]=satimes.std()
   return [stds.mean(), stds.std()/np.sqrt(ntrials)]

def get_atime_spread2(atimes, pnsync=0):
    nax=len(atimes)
    sspread=np.std(atimes)
    if pnsync>1 and pnsync<nax:
      sspread1=np.std(atimes[:pnsync])
      sspread2=np.std(atimes[pnsync:])
      return [sspread, sspread1,sspread2] + get_rndsubsample_estimate_atime(atimes, pnsync)
    else:
      return [sspread, None, None, None, None]

def get_save_params(saverespec, adict=None, verbose=3):
    ''' Specify what results need to be saved with an integer!
       Returns :  saveresults, nrepsave, epochsave, modelhistsavelist
       if saverespec is integer then:
         saverespec 0-9; but IF saverespec > 9:
         savemodelhistory=saverespec//10
         saverespec=saverespec%10
         nrepsave=(savemodelhistory%10)
         nepochsave=savemodelhistory//10
         epochsave=[0, nepochsave]
         Example: saverespec=525
       if saverespec is string then:
         split and try as integer but also use special characters:
         'a' for all reps, or epochs (then needs 'nreps' and 'nepochs' in adict)
         for epochs one can also choose range using 'r' as a separator: ##r##:
         Example:
              saver=5-7525 will set saveresults to 5, will save 2 replicates, and will save epochs 5 through 75, inclusive.
         if string has two characters only then the first character is iolig step in the chain from which to record ... I know, it's getting too complicated!
              saver=a25 ... saves OMP variable history all epochs for the first two replicates with saveresults=5
              saver=aa5 ... saves OMP variable history for all epochs and replicates with saveresults=5
    '''
    if isinstance(saverespec,int):
      if saverespec>9:
        savemodelhistory=saverespec//10
        saveresults=saverespec%10
        nrepsave=(savemodelhistory%10)
        nepochsave=savemodelhistory//10
        epochsave=[0, nepochsave]
        if nepochsave>0:
           if nrepsave: modelhistsavelist=[[] for _ in range(nrepsave)]
           else: modelhistsavelist=[[]]
        else:
          modelhistsavelist=[]
        if verbose>3: print([saveresults,nrepsave,epochsave,modelhistsavelist])
        return saveresults,nrepsave,epochsave,modelhistsavelist
      else:
        return saverespec, 0, [0,0], []
    elif isinstance(saverespec,str):
      try:
        svrspec=int(saverespec)
        return get_save_params(svrspec, adict=adict)
      except:
        strchars=list(saverespec)
        saveresults=int(strchars[-1])
        try:
          nrepsave=int(strchars[-2])
        except:
          if strchars[-2] == 'a':
            nrepsave=int(adict['nreps'])
          else:
            raise ValueError("Only integers < 10 or 'a' are accepted for the nreps character (2nd from the last)")
        if 'r' in strchars:
          irngsep=strchars.index('r')
          epstartstr=''.join(strchars[:irngsep])
          if epstartstr: epstart=int(epstartstr)-1
          else: epstart=0
          ependstr=''.join(strchars[irngsep+1:-2])
          if ependstr: epend=int(ependstr)
          else: epend=int(adict['nepochs'])+1
        else:
          epstart=0
          try:
            epend=int(''.join(strchars[:-2]))
          except:
            if strchars[-2] == 'a':
              epend=int(adict['nepochs'])
            else:
              raise ValueError("Improperly specified epochsave with saver=%s" % saverespec)
        epochsave=[epstart, epend]
        nepochsave=epend-epstart
        if nepochsave>0:
           if nrepsave: modelhistsavelist=[[] for _ in range(nrepsave)]
           else: modelhistsavelist=[[]]
        else:
          modelhistsavelist=[]
        if verbose>3: print([saveresults,nrepsave,epochsave,modelhistsavelist])
        return saveresults,nrepsave,epochsave,modelhistsavelist
    elif isinstance(saverespec,dict):
      if verbose>3: print('Using dictionary with saver=', saverespec['saver'])
      return get_save_params(saverespec['saver'], adict=saverespec)








