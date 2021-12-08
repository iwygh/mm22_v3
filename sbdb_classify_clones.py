#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 12:53:52 2021
classification simulation and data parsing is copied from sv20 code with minimal changes to not have to query Horizons through Rebound when setting up a sim, because Rebound's Horizons query is SOOOO SLOOOOW
@author: iggymatheson
"""
#%%
def parsedata(data,classifier):
    '''
    parse(data) computes the necessary features to classify
    data MUST be a 101 row x 6 column array
    columns are t, a, e, i, Omega, omega
    rows are different time outputs: MUST be 1000yr outputs, ie [0, 1E3, 2E3....99E3,100E3]
    Returns features for classification
    '''
    import numpy as np
    # Take stats of simulations
    initials = data[0,1:] # a, e, i, Omega, omega
    finals = data[-1,1:]
    mins = np.amin(data[:,1:],axis = 0)
    maxes = np.amax(data[:,1:],axis = 0)
    dels = maxes-mins
    means = np.mean(data[:,1:],axis = 0)
    stdev = np.std(data[:,1:],axis = 0)
    # Take time derivatives
    diffs = data[1:,:]-data[:-1,:]
    dxdt = diffs[:,1:]/diffs[:,0, np.newaxis] # add on new axis to time to give same dimensionality as the numerator
    mindxdt = np.amin(dxdt,axis = 0)
    meandxdt = np.mean(dxdt,axis = 0)
    maxdxdt = np.amax(dxdt,axis = 0)
    deldxdt = maxdxdt-mindxdt
    # rearrange data into the order I want
    arrs = [initials,finals,mins,means,maxes,stdev,dels,mindxdt,meandxdt,maxdxdt,deldxdt]
    inds = [0,1,2,3,4] # a, e, i, Omega, omega
    features = []
    ## features contains all x values, then all y, etc: xi, xf, xmin, xmean, xmax, xsigma, Deltax, xdotmin, xdotmean, xdotmax
    for i in inds:
        for a in arrs:
            features += [a[i]]
    features_out = np.array(features).reshape(1,-1) # make sure features is a 2d array
    prediction = classifier.predict_proba(features_out) # Predict the probabilities of class membership for object
    if np.max(prediction) == prediction[0][0]:
        category = int_dict[0]
    elif np.max(prediction) == prediction[0][1]:
        category = int_dict[1]
    elif np.max(prediction) == prediction[0][2]:
        category = int_dict[2]
    elif np.max(prediction) == prediction[0][3]:
        category = int_dict[3]
    print('This object has the following probabilities of class membership:')
    p=prediction[0]
    for i,k in enumerate(list(int_dict.keys())):
        print(int_dict[k],':',p[i]*100,'%')
    return category,prediction
#%%
def handle_warning(message, category, filename, lineno, file=None, line=None):
    import time
    print('A warning occurred:')
    print(message)
    time.sleep(1)
#%%
def runsim(sim_template,iobj,alist,elist,ilist,wlist,Wlist,Mlist):
    # try:
    sim = sim_template
    # sim.integrator = 'mercurius'
    # sim.integrator = 'ias15'
    # sim.ri_ias15.epsilon = 1e-10
    sim.integrator = "whfast"
    sim.dt = 0.01
    # sim.ri_whfast.corrector = 17
    primary = sim.calculate_com()
    ahere = alist[iobj]
    ehere = elist[iobj]
    ihere = ilist[iobj]
    where = wlist[iobj]
    Where = Wlist[iobj]
    Mhere = Mlist[iobj]
    print('a =',round(ahere,2),'au')
    print('e =',round(ehere,2))
    print('i =',round(ihere,2),'deg')
    print('w =',round(where,2),'deg')
    print('W =',round(Where,2),'deg')
    print('M =',round(Mhere,2),'deg')
    sim.add(a=ahere,e=ehere,inc=np.radians(ihere),omega=np.radians(where),\
            Omega=np.radians(Where),M=np.radians(Mhere),primary=primary)
    time_outs = np.linspace(0,100E3,101)*2*np.pi
    data = []
    for i,t in enumerate(time_outs):
        if t>0: sim.integrate(t, exact_finish_time=True) # integrate to next output
        orbits = sim.calculate_orbits(primary=sim.calculate_com())
        o = orbits[-1] # take KBO
        step = np.array([t/2/np.pi, o.a, o.e, np.degrees(o.inc), np.degrees(o.Omega)%360, np.degrees(o.omega)%360]) # save t, a, e, i, Omega, omega - time in data needs to be in years, so divide by 2pi
        # add step to data
        if len(data)==0: data = step
        else: data = np.vstack((data,step))
    return data
#%%
import numpy as np
import pandas as pd
import time
# import os
# import matplotlib.pyplot as plt
import rebound
from astroquery.jplhorizons import Horizons
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
import warnings
#%%
t0 = time.time()
training_file = 'KBO_features.csv'
df_training = pd.read_csv(training_file)
training_designations = df_training['MPC ID'].tolist()
training_classifications = df_training['Class'].tolist()
training_status = df_training['Securely Classified'].tolist()
all_KBOs = pd.read_csv(training_file,skipinitialspace=True)
secure_KBOs = all_KBOs[all_KBOs['Securely Classified']==True]
all_types = list(set(secure_KBOs['Class']))
types_dict = { all_types[i] : i for i in range( len(all_types) ) }
int_dict = { i : all_types[i] for i in range( len(all_types) ) }
classes = secure_KBOs['Class'].map(types_dict)
features_train, features_test, classes_train, classes_test = train_test_split(secure_KBOs, classes, test_size=0.3, random_state=30)
ids_train = features_train['MPC ID'].to_numpy()
features_train.drop(['MPC ID', 'Securely Classified', 'Class'], axis=1, inplace=True)
features_train = features_train.to_numpy()
ids_test = features_test['MPC ID'].to_numpy()
features_test.drop(['MPC ID', 'Securely Classified', 'Class'], axis=1, inplace=True)
features_test = features_test.to_numpy()
classifier = GradientBoostingClassifier( learning_rate=0.1, loss='deviance', max_depth=3, max_features='log2', n_estimators=130, random_state=30 )
classifier.fit(features_train, classes_train)
classes_predict = classifier.predict( features_test )
print('Classifier is ', accuracy_score(classes_test, classes_predict) * 100, '% accurate on testing set' )
# time.sleep(5)
#%%
warnings.showwarning = handle_warning
JD = 2459600.50000 # 2022 January 21 00:00:00
# sim = rebound.Simulation()
# center = '500@10' # center of the sun
# GM_list = [1,1/6023600,1/408523.71,1/328900.56,1/3098708,1/1047.3486,\
#             1/3497.898,1/22902.98,1/19412.24]
# m_sun_2 = GM_list[0] + GM_list[1] + GM_list[2] + GM_list[3] + GM_list[4]
# m_jupiter_2 = GM_list[5]/m_sun_2
# m_saturn_2 = GM_list[6]/m_sun_2
# m_uranus_2 = GM_list[7]/m_sun_2
# m_neptune_2 = GM_list[8]/m_sun_2
# m_sun_2 = m_sun_2/m_sun_2
# m_sun = m_sun_2
# m_jupiter = m_jupiter_2
# m_saturn = m_saturn_2
# m_uranus = m_uranus_2
# m_neptune = m_neptune_2
# sim.add(m = m_sun)
# name_in = '5' # jupiter barycenter
# obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
# el = obj.elements()
# a = float(el['a']) # au
# e = float(el['e'])
# inc = float(el['incl']) # degrees
# w = float(el['w']) # degrees
# W = float(el['Omega']) # degrees
# M = float(el['M']) # degrees
# inc = np.radians(inc)
# w = np.radians(w)
# W = np.radians(W)
# M = np.radians(M)
# w = np.mod(w,2*np.pi)
# W = np.mod(W,2*np.pi)
# M = np.mod(M,2*np.pi)
# sim.add(primary = sim.particles[0],m = m_jupiter,\
#         a = a,e = e,inc = inc,omega = w,Omega = W,M = M)
# name_in = '6' # saturn barycenter
# obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
# el = obj.elements()
# a = float(el['a']) # au
# e = float(el['e'])
# inc = float(el['incl']) # degrees
# w = float(el['w']) # degrees
# W = float(el['Omega']) # degrees
# M = float(el['M']) # degrees
# inc = np.radians(inc)
# w = np.radians(w)
# W = np.radians(W)
# M = np.radians(M)
# w = np.mod(w,2*np.pi)
# W = np.mod(W,2*np.pi)
# M = np.mod(M,2*np.pi)
# sim.add(primary = sim.particles[0],m = m_saturn,\
#         a = a,e = e,inc = inc,omega = w,Omega = W,M = M)
# name_in = '7' # uranus barycenter
# obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
# el = obj.elements()
# a = float(el['a']) # au
# e = float(el['e'])
# inc = float(el['incl']) # degrees
# w = float(el['w']) # degrees
# W = float(el['Omega']) # degrees
# M = float(el['M']) # degrees
# inc = np.radians(inc)
# w = np.radians(w)
# W = np.radians(W)
# M = np.radians(M)
# w = np.mod(w,2*np.pi)
# W = np.mod(W,2*np.pi)
# M = np.mod(M,2*np.pi)
# sim.add(primary = sim.particles[0],m = m_uranus,\
#         a = a,e = e,inc = inc,omega = w,Omega = W,M = M)
# name_in = '8' # neptune barycenter
# obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
# el = obj.elements()
# a = float(el['a']) # au
# e = float(el['e'])
# inc = float(el['incl']) # degrees
# w = float(el['w']) # degrees
# W = float(el['Omega']) # degrees
# M = float(el['M']) # degrees
# inc = np.radians(inc)
# w = np.radians(w)
# W = np.radians(W)
# M = np.radians(M)
# w = np.mod(w,2*np.pi)
# W = np.mod(W,2*np.pi)
# M = np.mod(M,2*np.pi)
# sim.add(primary = sim.particles[0],m = m_neptune,\
#         a = a,e = e,inc = inc,omega = w,Omega = W,M = M)
# sim.N_active = 5
# # sim.units = ('yr', 'AU', 'Msun')
# sim.move_to_com()
# # sim.integrator = 'mercurius'
# sim.integrator = 'ias15'
# # sim.dt = 0.01 #timestep in years
# # sim_template = sim
# # del sim
# sim.save("sim_template_20220122.bin")
# sim = None # Remove reference, allow python to release memory
#%%
sbdb_input = 'sbdb_query_results_addcov.csv'
df = pd.read_csv(sbdb_input,low_memory=False)
packed_designation_list = df['Packed MPC designation'].tolist()
Nobj = df.shape[0]
THIS_INSTANCE = 1 # runs 1 through 275
obj_per_instance = int(np.ceil(Nobj/300)) # 300 jobs at a time on the cluster
# obj_per_instance = 1
start_obj = (THIS_INSTANCE-1) * obj_per_instance
stop_obj = THIS_INSTANCE * obj_per_instance
if stop_obj > Nobj:
    stop_obj = Nobj
for iobj in range(start_obj,stop_obj):
    des = packed_designation_list[iobj]
    clones_input = 'sbdb_query_results_clones_' + des + '.csv'
    df2 = pd.read_csv(clones_input,low_memory=False)
    Nclones = df2.shape[0] - 1 # 301 rows, 1 nominal + 300 clones
    e_list = df2['e'].tolist()
    q_list = df2['q'].tolist() # au
    tp_list = df2['tp'].tolist() # days since pericenter
    node_list = df2['node'].tolist() # degrees
    peri_list = df2['peri'].tolist() # degrees
    i_list = df2['i'].tolist() # degrees
    a_list = []
    M_list = []
    class_list = []
    classical_list = []
    scattering_list = []
    detached_list = []
    resonant_list = []
    # t0 = time.time()
    for iclone in range(Nclones+1):
        e_here = e_list[iclone]
        q_here = q_list[iclone]
        tp_here = tp_list[iclone]
        node_here = node_list[iclone]
        peri_here = peri_list[iclone]
        i_here = i_list[iclone]
        a_here = q_here / (1-e_here)
        a_list.append(a_here)
        tp_here_in = tp_here / 365.25 # Julian days to Julian years
        tp_here_in = tp_here_in * 2*np.pi # yrs to yr2pi
        sim = rebound.Simulation("sim_template_20220122.bin")
        sim.t = 0
        sim.integrator = 'ias15'
        primary = sim.calculate_com()
        sim.add(a=a_here,e=e_here,inc=np.radians(i_here),omega=np.radians(peri_here),\
                Omega=np.radians(node_here),T=tp_here_in,primary=primary)
        M_here = sim.particles[5].M # radians
        M_here = np.degrees(M_here)
        M_list.append(M_here)
        print('sim.N =',sim.N)
        time_outs = np.linspace(0,100E3,101)*2*np.pi
        data = []
        for i,t in enumerate(time_outs):
            if t>0: sim.integrate(t, exact_finish_time=True) # integrate to next output
            orbits = sim.calculate_orbits(primary=sim.calculate_com())
            o = orbits[-1] # take KBO
            step = np.array([t/2/np.pi, o.a, o.e, np.degrees(o.inc), np.degrees(o.Omega)%360, np.degrees(o.omega)%360]) # save t, a, e, i, Omega, omega - time in data needs to be in years, so divide by 2pi
            # add step to data
            if len(data)==0: data = step
            else: data = np.vstack((data,step))
        category,prediction = parsedata(data,classifier)
        class_list.append(category)
        if int_dict[0] == 'Detached':
            detached_list.append(prediction[0][0])
        elif int_dict[0] == 'Resonant':
            resonant_list.append(prediction[0][0])
        elif int_dict[0] == 'Scattering':
            scattering_list.append(prediction[0][0])
        elif int_dict[0] == 'Classical':
            classical_list.append(prediction[0][0])
        else:
            raise Exception('int_dict is all weird',int_dict)
        if int_dict[1] == 'Detached':
            detached_list.append(prediction[0][1])
        elif int_dict[1] == 'Resonant':
            resonant_list.append(prediction[0][1])
        elif int_dict[1] == 'Scattering':
            scattering_list.append(prediction[0][1])
        elif int_dict[1] == 'Classical':
            classical_list.append(prediction[0][1])
        else:
            raise Exception('int_dict is all weird',int_dict)
        if int_dict[2] == 'Detached':
            detached_list.append(prediction[0][2])
        elif int_dict[2] == 'Resonant':
            resonant_list.append(prediction[0][2])
        elif int_dict[2] == 'Scattering':
            scattering_list.append(prediction[0][2])
        elif int_dict[2] == 'Classical':
            classical_list.append(prediction[0][2])
        else:
            raise Exception('int_dict is all weird',int_dict)
        if int_dict[3] == 'Detached':
            detached_list.append(prediction[0][3])
        elif int_dict[3] == 'Resonant':
            resonant_list.append(prediction[0][3])
        elif int_dict[3] == 'Scattering':
            scattering_list.append(prediction[0][3])
        elif int_dict[3] == 'Classical':
            classical_list.append(prediction[0][3])
        else:
            raise Exception('int_dict is all weird',int_dict)
        print('done with',iclone,Nclones,iobj,Nobj)
        # t1 = time.time()
        # print(t1-t0)
    df2['a'] = a_list
    df2['M'] = M_list
    df2['class'] = class_list
    df2['classical_probability'] = classical_list
    df2['resonant_probability'] = resonant_list
    df2['scattering_probability'] = scattering_list
    df2['detached_probability'] = detached_list
    clones_output = 'sbdb_query_results_clones_classified_' + des + '.csv'
    df2.to_csv(clones_output,index=False)
