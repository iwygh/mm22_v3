#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 12:33:22 2021

@author: iggymatheson
delete columns from sbdb_query_results_addcloneclass.csv to make it more Fortran-friendly
"""
import pandas as pd
import numpy as np
import time
sbdb_input = 'sbdb_query_results_addcloneclass.csv'
sbdb_output = 'sbdb_query_results_delcols.csv'
df = pd.read_csv(sbdb_input,low_memory=False)
df = df.drop('spkid', 1)
df = df.drop('full_name', 1)
df = df.drop('pdes', 1)
df = df.drop('name', 1)
df = df.drop('prefix', 1)
df = df.drop('neo', 1)
df = df.drop('pha', 1)
df = df.drop('H', 1)
df = df.drop('G', 1)
df = df.drop('M1', 1)
df = df.drop('M2', 1)
df = df.drop('K1', 1)
df = df.drop('K2', 1)
df = df.drop('PC', 1)
df = df.drop('diameter', 1)
df = df.drop('extent', 1)
df = df.drop('albedo', 1)
df = df.drop('rot_per', 1)
df = df.drop('GM', 1)
df = df.drop('BV', 1)
df = df.drop('UB', 1)
df = df.drop('IR', 1)
df = df.drop('spec_B', 1)
df = df.drop('spec_T', 1)
df = df.drop('H_sigma', 1)
df = df.drop('diameter_sigma', 1)
df = df.drop('orbit_id', 1)
df = df.drop('epoch.mjd', 1)
df = df.drop('epoch.cal', 1)
df = df.drop('equinox', 1)
df = df.drop('tp.cal', 1)
df = df.drop('per.y', 1)
df = df.drop('moid', 1)
df = df.drop('moid.ld', 1)
df = df.drop('moid_jup', 1)
df = df.drop('t_jup', 1)
df = df.drop('class', 1)
df = df.drop('producer', 1)
df = df.drop('data_arc', 1)
df = df.drop('first_obs', 1)
df = df.drop('last_obs', 1)
df = df.drop('n_obs_used', 1)
df = df.drop('n_del_obs_used', 1)
df = df.drop('n_dop_obs_used', 1)
df = df.drop('condition_code', 1)
df = df.drop('rms', 1)
df = df.drop('two_body', 1)
df = df.drop('A1', 1)
df = df.drop('A2', 1)
df = df.drop('A3', 1)
df = df.drop('DT', 1)
df = df.drop('pdes2', 1)
df = df.drop('Unpacked MPC designation', 1)
df = df.drop('e_e', 1)
df = df.drop('e_q', 1)
df = df.drop('e_tp', 1)
df = df.drop('e_node', 1)
df = df.drop('e_peri', 1)
df = df.drop('e_i', 1)
df = df.drop('q_e', 1)
df = df.drop('q_q', 1)
df = df.drop('q_tp', 1)
df = df.drop('q_node', 1)
df = df.drop('q_peri', 1)
df = df.drop('q_i', 1)
df = df.drop('tp_e', 1)
df = df.drop('tp_q', 1)
df = df.drop('tp_tp', 1)
df = df.drop('tp_node', 1)
df = df.drop('tp_peri', 1)
df = df.drop('tp_i', 1)
df = df.drop('node_e', 1)
df = df.drop('node_q', 1)
df = df.drop('node_tp', 1)
df = df.drop('node_node', 1)
df = df.drop('node_peri', 1)
df = df.drop('node_i', 1)
df = df.drop('peri_e', 1)
df = df.drop('peri_q', 1)
df = df.drop('peri_tp', 1)
df = df.drop('peri_node', 1)
df = df.drop('peri_peri', 1)
df = df.drop('peri_i', 1)
df = df.drop('i_e', 1)
df = df.drop('i_q', 1)
df = df.drop('i_tp', 1)
df = df.drop('i_node', 1)
df = df.drop('i_peri', 1)
df = df.drop('i_i', 1)
df.to_csv(sbdb_output,index=False)
 # barycentric 2022-01-21 00:00:00
 # rr,eclon,eclat,xhat,yhat,zhat
Nobj = df.shape[0]
rr_list = []
eclon_list = []
eclat_list = []
xhat_list = []
yhat_list = []
zhat_list = []
M_list = df['Mean anomaly degrees barycentric 2022-01-21 00:00:00'].tolist()
w_list = df['Argument of perihelion J2000.0 degrees barycentric 2022-01-21 00:00:00'].tolist()
W_list = df['Longitude of ascending node J2000.0 degrees barycentric 2022-01-21 00:00:00'].tolist()
i_list = df['Inclination to ecliptic J2000.0 degrees barycentric 2022-01-21 00:00:00'].tolist()
e_list = df['Orbital eccentricity barycentric 2022-01-21 00:00:00'].tolist()
a_list = df['Semimajor axis AU barycentric 2022-01-21 00:00:00'].tolist()
for iobj in range(Nobj):
    M = np.radians(M_list[iobj])
    w = np.radians(w_list[iobj])
    W = np.radians(W_list[iobj])
    i = np.radians(i_list[iobj])
    a = a_list[iobj]
    e = e_list[iobj]
    Etol = 1e-14
    maxiter = 1e7
    convergence = 0
    while convergence == 0:
        E = M
        diff = 1
        iter_ct = 0
        while np.abs(diff) > Etol and iter_ct < maxiter:
            top = E - e*np.sin(E) - M
            bottom = 1 - e*np.cos(E)
            diff = top / bottom
            E = E - diff
            iter_ct = iter_ct + 1
        if iter_ct >= maxiter:
            print('convergence failure at this Etol',np.abs(diff),Etol,e)
            Etol = Etol * 10
        else:
            convergence = 1
            # time.sleep(20)
    f = 2 * np.arctan(np.sqrt(1+e)/np.sqrt(1-e)*np.tan(E/2))
    r = a * (1-e*e) / (1+e*np.cos(f))
    rr_list.append(r)
    xperi = r * np.cos(f)
    yperi = r * np.sin(f)
    zperi = 0
    # using mu = 1
    h = np.sqrt(a*(1-e*e))
    vxperi = 1/h * -np.sin(f)
    vyperi = 1/h * (e+np.cos(f))
    vzperi = 0
    Q11 = -np.sin(W)*np.cos(i)*np.sin(w)+np.cos(W)*np.cos(w)
    Q12 = -np.sin(W)*np.cos(i)*np.cos(w)-np.cos(W)*np.sin(w)
    Q13 = np.sin(W)*np.sin(i)
    Q21 = np.cos(W)*np.cos(i)*np.sin(w)+np.sin(W)*np.cos(w)
    Q22 = np.cos(W)*np.cos(i)*np.cos(w)-np.sin(W)*np.sin(w)
    Q23 = -np.cos(W)*np.sin(i)
    Q31 = np.sin(i)*np.sin(w)
    Q32 = np.sin(i)*np.cos(w)
    Q33 = np.cos(i)
    rperivec = np.array([xperi,yperi,zperi])
    vperivec = np.array([vxperi,vyperi,vzperi])
    Qmat = np.array([[Q11,Q12,Q13],[Q21,Q22,Q23],[Q31,Q32,Q33]])
    rvec = np.dot(Qmat,rperivec)
    vvec = np.dot(Qmat,vperivec)
    x = rvec[0]
    y = rvec[1]
    z = rvec[2]
    xhat = x/r
    yhat = y/r
    zhat = z/r
    xhat_list.append(xhat)
    yhat_list.append(yhat)
    zhat_list.append(zhat)
    eclat = np.arcsin(zhat)
    eclon = np.arctan2(yhat/np.cos(eclat),xhat/np.cos(eclat))
    eclat = np.degrees(eclat)
    eclon = np.degrees(eclon)
    eclon = np.mod(eclon,360)
    eclon_list.append(eclon)
    eclat_list.append(eclat)
    print(iobj,Nobj,Etol,e)
df['Range au barycentric 2022-01-21 00:00:00'] = rr_list
df['xhat barycentric 2022-01-21 00:00:00'] = xhat_list
df['yhat barycentric 2022-01-21 00:00:00'] = yhat_list
df['zhat barycentric 2022-01-21 00:00:00'] = zhat_list
df['Ecliptic latitude beta J2000.0 2022-01-21 00:00:00'] = eclat_list
df['Ecliptic longitude lambda J2000.0 2022-01-21 00:00:00'] = eclon_list
df.to_csv(sbdb_output,index=False)
