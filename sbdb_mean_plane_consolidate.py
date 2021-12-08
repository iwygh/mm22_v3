#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 09:05:08 2021
@author: iggymatheson
check for maxiter overruns and consolidate mean plane results
"""
#%%
Njob = 300
rootstr = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3/ianmatheson-sbdb_classify_clones-2705728-'
for ijob in range(Njob):
    filestr = rootstr + str(ijob+1) + '.out'
    file1 = open(filestr, 'r')
Lines = file1.readlines()
count = 0
# Strips the newline character
for line in Lines:
    try:
        print(line.index('maxiter'),count+1,ijob+1)
    except:
        2-2
    count = count + 1
#%%
import pandas as pd
import numpy as np
import os
amin_strlist = ['30','34','34.79','34.79','34.79','40','40.524',\
                '42','42','42','43','44','45','45','50','50']
amax_strlist = ['150','36','40.524','50','150','41','42',\
                '43','45','48','44','45','48','50','80','150' ]
objct_strlist = []
Nbin = len(amin_strlist)
Njob = 300
Nrep = 134
for ibin in range(Nbin):
    amin = amin_strlist[ibin]
    amax = amax_strlist[ibin]
    str1 = 'sbdb_query_results_delcols__objct'
    str2 = '_amin' + amin + '_amax' + amax + '_Nrep' + str(Nrep) + '_set1.txt'
    maxobj = 10000
    for iobj in range(maxobj):
        if os.path.exists(str1 + str(iobj+1) + str2):
            str3 = str(iobj+1)
    # print(ibin+1,str3,amin,amax)
    objct_strlist.append(str3)
# df = pd.DataFrame()
for ibin in range(Nbin):
# for ibin in range(1):
    df = pd.DataFrame()
    amin = amin_strlist[ibin]
    amax = amax_strlist[ibin]
    objct = objct_strlist[ibin]
    str1 = 'sbdb_query_results_delcols__objct' + objct + '_amin' + amin + '_amax' + amax + '_Nrep'
    str2 = str(Nrep) + '_set'
    for ijob in range(Njob):
        infile = str1 + str2 + str(ijob+1) + '.txt'
        df2 = pd.read_csv(infile,delim_whitespace=True)
        frames = [df,df2]
        df = pd.concat(frames)
    Nrep_total = df.shape[0]
    # outfile = str1 + str(Nrep_total) + '_meanplanes.txt'
    outfile = 'sbdb_query_results_delcols_meanplanes_objct' + objct + \
                '_amin' + amin + '_amax' + amax + '_Nrep' + str(Nrep_total) + '.txt'
    df.to_csv(outfile,index=False)
