#Script to create simulated data

import os, glob
import numpy as np
import libstempo as T2
from libstempo import toasim as LT
from libstempo import plot as LP
import enterprise
from enterprise.pulsar import Pulsar
import random


def create_psrs(parfiles, timfiles):
    psrs = []
    for file in parfiles:
        try:
            print(file)
            tim = glob.glob(file.strip(".par")+"*tim")[0]
            psr = T2.tempopulsar(parfile = file, timfile = tim)
            #LT.make_ideal(psr)
            #LT.add_efac(psr,efac=1.0,seed=1234)
            psrs.append(psr)
        except:
            pass
        #psrs.append(LT.tempopulsar(parfile = T.data + 'B1953+29_NANOGrav_dfg+12.par', timfile = T.data + 'B1953+29_NANOGrav_dfg+12.tim'))
    return psrs

def create_ideal_psrs(parfiles, timfiles):
    psrs = []
    for file in parfiles:
        try:
            print(file)
            tim = glob.glob(file.strip("_tdb.par")+"*tim")[0]
            psrT2 = T2.tempopulsar(parfile = file, timfile = tim)
            #psr = LT.fakepulsar(parfile = file, obstimes = psrT2.toas(), freq=list(psrT2.freqs), toaerr=1e-3, observatory="meerkat")
            LT.make_ideal(psrT2)
            LT.add_efac(psrT2,efac=1.0,seed=1)
            psrT2.fit()
            psrs.append(psrT2)
        except:
            pass
        #psrs.append(LT.tempopulsar(parfile = T.data + 'B1953+29_NANOGrav_dfg+12.par', timfile = T.data + 'B1953+29_NANOGrav_dfg+12.tim'))
    return psrs

def add_red(psrs, amp=1e-13, gamma=3., rand=True):
    red_noise = {}
    for psr in psrs:
        try:
            if rand==True: 
                amp = random.uniform(1e-17,1e-13)
                gamma = random.uniform(1,7)
            LT.add_rednoise(psr,amp,gamma)
            psr.fit()
            red_noise[psr.name] = {'red_amp': amp, 'red_gamma': gamma}
        except:
            pass
    return psrs, red_noise


def add_red_singular(psr, amp=1e-13, gamma=3., rand=True):
    psr_temp = psr
    try:
        if rand==True: 
            amp = random.uniform(1e-17,1e-13)
            gamma = random.uniform(1,7)
        LT.add_rednoise(psr_temp,amp,gamma)
        psr_temp.fit()
        print("worked")
    except:
        pass
    return psr_temp, amp, gamma

def add_dm(psrs, amp=1e-13, gamma=3., rand=True):
    dm_noise = {}
    for psr in psrs:
        try: 
            if rand==True: 
                amp = random.uniform(1e-17,1e-13)
                gamma = random.uniform(1,7)
            LT.add_dm(psr,amp,gamma)
            psr.fit()
            dm_noise[psr.name] = {'dm_amp': amp, 'dm_gamma': gamma}
        except:
            pass
    return psrs, dm_noise

def add_scattering(psrs, amp=1e-13, gamma=3., scatt_idx=4, rand=True):
    chrom_noise = {}
    for psr in psrs:
        try: 
            if rand==True: 
                amp = random.uniform(1e-17,1e-13)
                gamma = random.uniform(1,7)
                scatt_idx = np.random.normal(4, 1)
            LT.add_scattering(psr, amp, gamma, scatt_idx=scatt_idx)
            psr.fit()
            chrom_noise[psr.name] = {'chrom_amp': amp, 'chrom_gamma': gamma, 'chrom_idx': scatt_idx}
        except:
            pass
    return psrs, chrom_noise

def add_gwb(psrs, amp=1e-14, gam=13./3.):
    """ Add a GWB.
        Takes in pulsar objects, amplitude of RN, and spectral index gamma.
    """
    LT.createGWB(psrs, amp, gam) #modifies pulsars in place, no need to return anything

    
def add_cgw(psrs, pdict, tref, iters):
    """ Add a continuous wave signal.
        Takes in pulsar objects, a parameter dictionary for the single source, and
        a reference time for observations.
    """
    for psr in psrs:
        #also modifying pulsars in place
        LT.add_cgw(psr, gwtheta=pdict['gwtheta'], gwphi=pdict['gwphi'], mc=pdict['mc'], dist=pdict['dist'], fgw=pdict['fgw'], phase0=pdict['phase0'], psi=pdict['psi'], inc=pdict['inc'], pdist=1., pphase=None, psrTerm=False, evolve=False, phase_approx=False, tref=tref)
        
        #iterate the timing model fit a few times
        try:
            psr.fit(iters=iters)
            print(psr.name)
        except:
            print(psr.name, 'had timing model fit issue. Excluding from PTA.')

def lt2ent(psrs):
    """ Converts libstempo pulsar objects to enterprise pulsar objects.
    """
    Epsrs = [Pulsar(psr) for psr in psrs]
    return Epsrs
    
def save_sims(psrs, outdir):
    """ Save simulated timing files.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for p in psrs:
        p.savepar(outdir + p.name + '_simulated.par')
        p.savetim(outdir + p.name + '_simulated.tim')
