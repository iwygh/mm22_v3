#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:35:41 2021

@author: iggymatheson
"""
import numpy as np
import pandas as pd
from astroquery.jplhorizons import Horizons
JD = 2459601.50000 # 2022-01-22 00:00:00
outfile = 'sun_and_planets_Horizons_barycentric_20220122.csv'
center = '500@0' # solar system barycenter
GM_list = [1,1/6023600  ,1/408523.71,1/328900.56,1/3098708,\
             1/1047.3486,1/3497.898 ,1/22902.98 ,1/19412.24]
a_list = []
e_list = []
i_list = []
w_list = []
W_list = []
M_list = []
name_in = '10' # center of the sun
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '1' # earth barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '2' # venus barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '3' # earth-moon barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '4' # mars barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '5' # jupiter barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '6' # saturn barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '7' # uranus barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
name_in = '8' # neptune barycenter
obj = Horizons(id=name_in,location=center,epochs=JD,id_type='majorbody')
el = obj.elements()
a = float(el['a']) # au
e = float(el['e'])
inc = float(el['incl']) # degrees
w = float(el['w']) # degrees
W = float(el['Omega']) # degrees
M = float(el['M']) # degrees
w = np.mod(w,360)
W = np.mod(W,360)
M = np.mod(M,360)
a_list.append(a)
e_list.append(e)
i_list.append(inc)
w_list.append(w)
W_list.append(W)
M_list.append(M)
dictionary = {'GM_solar_masses':GM_list,
              'semimajor_axis_au':a_list,
              'eccentricity':e_list,
              'inclination_degrees':i_list,
              'argument_of_pericenter_degrees':w_list,
              'longitude_of_node_degrees':W_list,
              'mean_anomaly_degrees':M_list
              }
df_out = pd.DataFrame(dictionary)
df_out.to_csv(outfile,index=False)
