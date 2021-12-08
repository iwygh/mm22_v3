#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 08:37:09 2021
@author: iggymatheson
read in the results of the clone classifications and add those as columns to
the sbdb_query_results_clones_MPCDES.csv files and to the sbdb_query_results_addcov.csv file
"""
import numpy as np
import pandas as pd
import time
import os
import json
import urllib
all_objects_input = 'sbdb_query_results_addcov.csv'
all_objects_output = 'sbdb_query_results_addcloneclass.csv'
sbdb_clones_input_stem = 'sbdb_query_results_clones_'
sbdb_clones_output_stem = 'sbdb_query_results_clones_classified_'
dfall = pd.read_csv(all_objects_input,low_memory=False)
alldes = dfall['Packed MPC designation'].tolist()
Nobj = dfall.shape[0]
Njobs = 275
job_input_stem = 'ianmatheson-sbdb_classify_clones-2689905-'
nominal_class_list = []
most_common_class_list = []
detached_count_list = []
scattering_count_list = []
classical_count_list = []
resonant_count_list = []
most_common_class_percent_list = []
old_des = 'blank'
for ijob in range(Njobs):
# for ijob in range(0,1):
# for ijob in range(0,3):
    jobct = ijob+1
    job_input_file = job_input_stem + str(jobct) + '.out'
    file = open(job_input_file,'r')
    lines = file.readlines()
    line_ct = 0
    line = lines[line_ct] # first line: Classifier is  98.33948339483395 % accurate on testing set
    line_ct = line_ct + 1
    line = lines[line_ct]
    # this line is either the start of the next clone classification, or it's the end of the file
    while line[0:9] == 'sim.N = 6':
        line_ct = line_ct + 1
        line = lines[line_ct] # This object has the following probabilities of class membership:
        line_ct = line_ct + 1
        line = lines[line_ct] # Category : 0.007104932245933555 %
        index1 = line.index(':')
        index2 = line.index('%')
        this_class = line[0:index1-1]
        this_probability = float(line[index1+1:index2])
        if this_class == 'Detached':
            detached_probability = this_probability
        if this_class == 'Resonant':
            resonant_probability = this_probability
        if this_class == 'Scattering':
            scattering_probability = this_probability
        if this_class == 'Classical':
            classical_probability = this_probability
        line_ct = line_ct + 1
        line = lines[line_ct] # Category : 0.007104932245933555 %
        index1 = line.index(':')
        index2 = line.index('%')
        this_class = line[0:index1-1]
        this_probability = float(line[index1+1:index2])
        if this_class == 'Detached':
            detached_probability = this_probability
        if this_class == 'Resonant':
            resonant_probability = this_probability
        if this_class == 'Scattering':
            scattering_probability = this_probability
        if this_class == 'Classical':
            classical_probability = this_probability
        line_ct = line_ct + 1
        line = lines[line_ct] # Category : 0.007104932245933555 %
        index1 = line.index(':')
        index2 = line.index('%')
        this_class = line[0:index1-1]
        this_probability = float(line[index1+1:index2])
        if this_class == 'Detached':
            detached_probability = this_probability
        if this_class == 'Resonant':
            resonant_probability = this_probability
        if this_class == 'Scattering':
            scattering_probability = this_probability
        if this_class == 'Classical':
            classical_probability = this_probability
        line_ct = line_ct + 1
        line = lines[line_ct] # Category : 0.007104932245933555 %
        index1 = line.index(':')
        index2 = line.index('%')
        this_class = line[0:index1-1]
        this_probability = float(line[index1+1:index2])
        if this_class == 'Detached':
            detached_probability = this_probability
        if this_class == 'Resonant':
            resonant_probability = this_probability
        if this_class == 'Scattering':
            scattering_probability = this_probability
        if this_class == 'Classical':
            classical_probability = this_probability
        line_ct = line_ct + 1
        line = lines[line_ct]
        space1 = line.index(' ')
        space2 = line[space1+1:-1].index(' ') + space1 + 1
        space3 = line[space2+1:-1].index(' ') + space2 + 1
        space4 = line[space3+1:-1].index(' ') + space3 + 1
        space5 = line[space4+1:-1].index(' ') + space4 + 1
        int1 = int(line[space2:space3])
        int2 = int(line[space3:space4])
        int3 = int(line[space4:space5])
        int4 = int(line[space5:len(line)])
        # iclone,Nclones,iobj,Nobj
        iclone_here = int1
        Nclones_here = int2
        iobj_here = int3
        Nobj_here = int4
        des = alldes[iobj_here]
        if des != old_des: # open file, don't need to open it anew each time
            clones_input = sbdb_clones_input_stem + des + '.csv'
            clones_output = sbdb_clones_output_stem + des + '.csv'
            dfdes = pd.read_csv(clones_input,low_memory=False)
            old_des = des # so we don't have to reopen this every time
            obj_class_list = []
            obj_classical_list = []
            obj_resonant_list = []
            obj_scattering_list = []
            obj_detached_list = []
            obj_classical_count = 0
            obj_resonant_count = 0
            obj_scattering_count = 0
            obj_detached_count = 0
        probabilities_list = np.array([detached_probability,resonant_probability,classical_probability,scattering_probability])
        if np.max(probabilities_list) == detached_probability:
            category_here = 'Detached'
            obj_class_list.append(category_here)
            obj_detached_count = obj_detached_count + 1
        if np.max(probabilities_list) == resonant_probability:
            category_here = 'Resonant'
            obj_class_list.append(category_here)
            obj_resonant_count = obj_resonant_count + 1
        if np.max(probabilities_list) == classical_probability:
            category_here = 'Classical'
            obj_class_list.append(category_here)
            obj_classical_count = obj_classical_count + 1
        if np.max(probabilities_list) == scattering_probability:
            category_here = 'Scattering'
            obj_class_list.append(category_here)
            obj_scattering_count = obj_scattering_count + 1
        obj_classical_list.append(classical_probability)
        obj_resonant_list.append(resonant_probability)
        obj_scattering_list.append(scattering_probability)
        obj_detached_list.append(detached_probability)
        # print(int1,int2,int3,int4,des)
        # print(category_here,detached_probability,resonant_probability,classical_probability,scattering_probability)
        if iclone_here == Nclones_here:
            dfdes['Category'] = obj_class_list
            dfdes['Resonant probability'] = obj_resonant_list
            dfdes['Classical probability'] = obj_classical_list
            dfdes['Scattering probability'] = obj_scattering_list
            dfdes['Detached probability'] = obj_detached_list
            dfdes.to_csv(clones_output,index=False)
        if iclone_here == 0:
            nominal_class_list.append(category_here)
        if iclone_here == Nclones_here:
            detached_count_list.append(obj_detached_count)
            resonant_count_list.append(obj_resonant_count)
            classical_count_list.append(obj_classical_count)
            scattering_count_list.append(obj_scattering_count)
            countarray = np.array([obj_detached_count,obj_resonant_count,obj_classical_count,obj_scattering_count])
            maxcount = np.max(countarray)
            most_common_class_percent_list.append(maxcount/(Nclones_here+1))
            if np.max(countarray) == obj_detached_count:
                most_common_class_list.append('Detached')
            if np.max(countarray) == obj_resonant_count:
                most_common_class_list.append('Resonant')
            if np.max(countarray) == obj_classical_count:
                most_common_class_list.append('Classical')
            if np.max(countarray) == obj_scattering_count:
                most_common_class_list.append('Scattering')
        # move on to next reading cycle
        line_ct = line_ct + 1
        line = lines[line_ct]
dfall['Nominal class'] = nominal_class_list
dfall['Most common class'] = most_common_class_list
dfall['Most common class percent'] = most_common_class_percent_list
dfall['Detached count'] = detached_count_list
dfall['Scattering count'] = scattering_count_list
dfall['Classical count'] = classical_count_list
dfall['Resonant count'] = resonant_count_list
dfall.to_csv(all_objects_output,index=False)
# nominal_class_list = []
# most_common_class_list = []
# detached_count_list = []
# scattering_count_list = []
# classical_count_list = []
# resonant_count_list = []
# most_common_class_percent_list = []
