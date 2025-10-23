import numpy as np
import matplotlib.pyplot as plt
import os, glob
import libstempo as T2
from libstempo import toasim as LT
from libstempo import plot as LP
import enterprise
from enterprise.pulsar import Pulsar
from simFuncs import *
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Fake data maker.")
parser.add_argument("-noise", type = str.lower, nargs="+",dest="noise", help="Noise values to induce", \
    choices={"red", "dm", "chrom", "gwb", "cw"}, required=False)
parser.add_argument('--rand', default=False, help="If this is turned on, randomise a spectral index and amplitude", action=argparse.BooleanOptionalAction)
parser.add_argument('--out', type = str.lower, dest="out", help="full path of savefile", required=True)
args = parser.parse_args()

noise = args.noise
out = str(args.out)

# get parfiles containing pulsar params
datadir = '/home/mattm/projects/HSYMT/partim_real/tdb/'
rand_select = np.random.randint(0, 82, size=10)

parfiles = sorted(glob.glob(datadir + '/*.par'))
timfiles = sorted(glob.glob(datadir + '/*.tim'))

parfiles = np.array(parfiles)[[rand_select]].tolist()[0]
timfiles = np.array(timfiles)[[rand_select]].tolist()[0]

# create pulsar objects from parfiles
# default TOA errors are set to 0.4 microseconds
psrs = []
psrs = create_ideal_psrs(parfiles, timfiles=timfiles)

cw_p_dict = {
    'gwtheta': 1.75,  
    'gwphi': 5.,  
    'mc': 5e9,  
    'dist': 60, 
    'fgw': 2e-8, 
    'phase0': 0.0,
    'psi': np.pi/4.0,
    'inc': 0.0,
}

tref = np.min([p.toas().min() for p in psrs])  # reference time for observations

total_noise = []
#Give pulsars noise
for n in noise:
    # if not args.rand:
    #     if n == "red":
    #         add_red(psrs)
    #     if n == "dm":
    #         add_dm(psrs)
    # else:
    if n == "red":
        psrs, red_noise = add_red(psrs,rand=True)
        print("Red noise added")
        print(red_noise)
    if n == "dm":
        psrs, dm_noise = add_dm(psrs,rand=True)
        print("DM noise added")
        print(dm_noise)
    if n == "chrom":
        psrs, chrom_noise = add_scattering(psrs,rand=True)
        print("Chromatic noise added")
        print(chrom_noise)
    if n == "gwb":  
        add_gwb(psrs)
        print("GWB added")
    if n == "cw":
        add_cgw(psrs, cw_p_dict, tref, iters=10)
        print("Continuous wave added")
        
total_noise.append(red_noise, dm_noise, chrom_noise)

#Save the pulsars as MJD, frequencies residuals

#toa_all, freq_all, res_all = [], [], []
data_all = []
for psr in psrs:
    psrname = psr.name
    toas = psr.toas()
    freq = psr.freqs.astype(np.float128)
    res = psr.residuals()
    psrname_list = [psrname]*len(toas)
    data = np.array([psrname_list, toas,freq,res])
    data_all.append(data.T)

df = pd.DataFrame(np.row_stack(data_all),columns=["psr","toas","freqs","res"])
df.to_pickle(out)