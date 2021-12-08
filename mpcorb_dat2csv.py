#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 08:55:37 2021

@author: iggymatheson
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 09:34:15 2021

@author: iggymatheson
"""
#%%
import numpy as np
import pandas as pd
import time
import os
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
#%% first, move data from DAT to a Pandas csv
# mpcorb_in = 'MPCORB_20161004.DAT'
# mpcorb_out = 'MPCORB_20161004.csv'
# mpcorb_in = '/Users/iggymatheson/Documents_off_iCloud/mm22 mpcorb databases/MPCORB_20161004.DAT'
# mpcorb_out = '/Users/iggymatheson/Documents_off_iCloud/mm22 mpcorb databases/MPCORB_20161004.csv'
mpcorb_in = '/Users/iggymatheson/Documents_off_iCloud/mm22 mpcorb databases/MPCORB_20211202.DAT'
mpcorb_out = '/Users/iggymatheson/Documents_off_iCloud/mm22 mpcorb databases/MPCORB_20211202.csv'
designation_list = []
H_list = []
G_list = []
Epoch_list = []
M_list = []
Argperi_list = []
Node_list = []
Inc_list = []
e_list = []
n_list = []
a_list = []
U_list = []
ref_list = []
Nobs_list = []
Nopp_list = []
yrfirst_list = []
yrlast_list = []
arclength_list = []
rmsres_list = []
coarse_list = []
precise_list = []
computer_list = []
readable_list = []
lastobs_list = []
file = open(mpcorb_in,'r')
lines = file.readlines()
line_ct = 0
Nlines = len(lines)
# for iline in range(1151450,1151461):
# for iline in range(0,50):
for iline in range(Nlines):
    if iline > 42:
        line = lines[iline]
        if len(line) > 3: # protect against blank lines
            designation = line[0:7]
            designation = designation.lstrip()
            designation = designation.rstrip()
            designation_list.append(designation)
            # print(designation)
            magnitude = line[8:13]
            magnitude = magnitude.lstrip()
            magnitude = magnitude.rstrip()
            H_list.append(magnitude)
            # print(designation,magnitude)
            slope = line[14:19]
            slope = slope.lstrip()
            slope = slope.rstrip()
            G_list.append(slope)
            # print(designation,magnitude,slope)
            epoch = line[20:25]
            epoch = epoch.lstrip()
            epoch = epoch.rstrip()
            Epoch_list.append(epoch)
            # print(designation,magnitude,slope,epoch)
            M = line[26:35]
            M = M.lstrip()
            M = M.rstrip()
            M_list.append(M)
            # print(designation,epoch,M)
            argperi = line[37:46]
            argperi = argperi.lstrip()
            argperi = argperi.rstrip()
            Argperi_list.append(argperi)
            # print(designation,epoch,M,argperi)
            node = line[48:57]
            node = node.lstrip()
            node = node.rstrip()
            Node_list.append(node)
            incl = line[59:68]
            incl = incl.lstrip()
            incl = incl.rstrip()
            Inc_list.append(incl)
            ecc = line[70:79]
            ecc = ecc.lstrip()
            ecc = ecc.rstrip()
            e_list.append(ecc)
            meanmotion = line[80:91]
            meanmotion = meanmotion.lstrip()
            meanmotion = meanmotion.rstrip()
            n_list.append(meanmotion)
            sma = line[92:103]
            sma = sma.lstrip()
            sma = sma.rstrip()
            a_list.append(sma)
            U = line[105]
            if U == 'E':
                U = 999
            if U == 'D':
                U = 999
            if U == 'F':
                U = 999
            if U == ' ':
                U = 999
            U = int(U)
            U_list.append(U)
            ref = line[107:116]
            ref = ref.lstrip()
            ref = ref.rstrip()
            ref_list.append(ref)
            Nobs = line[117:122]
            Nobs = Nobs.lstrip()
            Nobs = Nobs.rstrip()
            try:
                Nobs = int(Nobs)
                Nobs_list.append(Nobs)
            except:
                Nobs_list.append(-999)
            # print(Nobs)
            Nopp = line[123:126]
            Nopp = Nopp.lstrip()
            Nopp = Nopp.rstrip()
            Nopp = int(Nopp)
            Nopp_list.append(Nopp)
            if Nopp > 1:
                yrfirst = line[127:131]
                yrfirst = yrfirst.lstrip()
                yrfirst = yrfirst.rstrip()
                yrfirst = int(yrfirst)
                yrfirst_list.append(yrfirst)
                yrlast = line[132:136]
                yrlast = yrlast.lstrip()
                yrlast = yrlast.rstrip()
                yrlast = int(yrlast)
                yrlast_list.append(yrlast)
            else:
                yrfirst_list.append(-999)
                yrlast_list.append(-999)
            if Nopp == 1:
                arclength = line[127:131]
                arclength = arclength.lstrip()
                arclength = arclength.rstrip()
                arclength_list.append(arclength)
            else:
                arclength_list.append(-999)
            rmsres = line[137:141]
            rmsres = rmsres.lstrip()
            rmsres = rmsres.rstrip()
            rmsres_list.append(rmsres)
            coarse = line[142:145]
            coarse = coarse.lstrip()
            coarse = coarse.rstrip()
            coarse_list.append(coarse)
            precise = line[146:149]
            precise = precise.lstrip()
            precise = precise.rstrip()
            precise_list.append(precise)
            computer = line[150:160]
            computer = computer.lstrip()
            computer = computer.rstrip()
            computer_list.append(computer)
            readable = line[166:194]
            readable = readable.lstrip()
            readable = readable.rstrip()
            readable_list.append(readable)
            lastobs = line[194:202]
            lastobs = lastobs.lstrip()
            lastobs = lastobs.rstrip()
            lastobs_list.append(lastobs)
dictionary = {'Number or provisional designation in packed form':designation_list,
              'H absolute magnitude':H_list,
              'G slope parameter':G_list,
              'Epoch in packed form .0 TT':Epoch_list,
              'Mean anomaly at epoch in degrees':M_list,
              'Argument of perihelion J2000.0 degrees':Argperi_list,
              'Longitude of ascending node J2000.0 degrees':Node_list,
              'Inclination to ecliptic J2000.0 degrees':Inc_list,
              'Orbital eccentricity':e_list,
              'Mean daily motion degrees per day':n_list,
              'Semimajor axis AU':a_list,
              'Uncertainty parameter U':U_list,
              'Reference':ref_list,
              'Number of observations':Nobs_list,
              'Number of oppositions':Nopp_list,
              'Year of first observation':yrfirst_list,
              'Year of last observation':yrlast_list,
              'Arc length days':arclength_list,
              'RMS residual arcsec':rmsres_list,
              'Coarse indicator of perturbers blank if unperturbed one opposition object':coarse_list,
              'Precise indicator of perturbers blank if unperturbed one opposition object':precise_list,
              'Computer name':computer_list,
              'Readable designation':readable_list,
              'Date of last observation included in orbit solution YYYYMMDD':lastobs_list
              }
df_out = pd.DataFrame(dictionary)
df_out.to_csv(mpcorb_out,index=False)
