#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 13:37:47 2021

@author: iggymatheson
"""
template_file = 'sbdb_classify_clones.py'
with open(template_file,'r') as file:
    template_data = file.read()
Nsets = 275
for i in range(Nsets):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'sbdb_classify_clones_' + ctstr + '.py'
    filedata = template_data
    filedata = filedata.replace('THIS_INSTANCE = 1','THIS_INSTANCE = '+str(int(i+1)))
    with open(outfile,'w') as file:
        file.write(filedata)
