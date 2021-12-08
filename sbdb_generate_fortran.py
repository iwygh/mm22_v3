#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 13:37:47 2021

@author: iggymatheson
"""
template_file = 'fortran_template.f95'
with open(template_file,'r') as file:
    template_data = file.read()
Nsets = 300
for i in range(Nsets):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'fortran_mean_plane' + ctstr + '.f95'
    filedata = template_data
    filedata = filedata.replace('setstr = "1"','setstr = "'+ctstr+'"')
    with open(outfile,'w') as file:
        file.write(filedata)
