#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 14:31:55 2021
@author: iggymatheson
show convergence of 1sigma imin,imax,Omegamin,Omegamax for Nreps up to 40,200
"""
#%%
def covariance_ellipse(q_list,p_list,N_sigma):
    import numpy as np
    from shapely.geometry import Polygon, Point
    q_ellipse,p_ellipse = covariance_ellipse_2(q_list,p_list,N_sigma)
    sin_i_ellipse = np.sqrt(q_ellipse**2+p_ellipse**2)
    i_ellipse = np.arcsin(sin_i_ellipse)
    Omega_ellipse = np.arctan2(p_ellipse/sin_i_ellipse,q_ellipse/sin_i_ellipse)
    i_ellipse_degrees = np.degrees(i_ellipse)
    Omega_ellipse_degrees = np.degrees(Omega_ellipse) # -180 to + 180 degrees
    i_min_degrees = np.min(i_ellipse_degrees)
    i_max_degrees = np.max(i_ellipse_degrees)
    Omega_min_degrees = np.min(Omega_ellipse_degrees)
    Omega_max_degrees = np.max(Omega_ellipse_degrees)
    # if ellipse straddles the second and third quadrants
    if (-180<=Omega_min_degrees<-90) and (90<Omega_max_degrees<=180):
        Omega_ellipse_degrees = np.mod(Omega_ellipse_degrees,360)
        Omega_min_degrees = np.min(Omega_ellipse_degrees)
        Omega_max_degrees = np.max(Omega_ellipse_degrees)
    else:
        Omega_min_degrees = np.mod(Omega_min_degrees,360)
        Omega_max_degrees = np.mod(Omega_max_degrees,360)
    # if ellipse contains origin, Omega runs 0 to 360 degrees
    Ne = len(q_ellipse)
    linestring = []
    for i2 in range(Ne):
        pt = (q_ellipse[i2],p_ellipse[i2])
        linestring.append(pt)
    pt = (q_ellipse[0],p_ellipse[0])
    linestring.append(pt)
    poly = Polygon(linestring)
    pt = Point(0,0)
    checkstatus = pt.within(poly)
    if checkstatus == True:
        Omega_min_degrees = 0
        Omega_max_degrees = 360
    return i_min_degrees,i_max_degrees,Omega_min_degrees,Omega_max_degrees
#%%
def covariance_ellipse_2(q_list,p_list,N_sigma):
    import numpy as np
    from numpy import linalg as LA
    Npts = 5000
    Nqp = len(q_list)
    q_list = np.reshape(q_list,(Nqp,1))
    p_list = np.reshape(p_list,(Nqp,1))
    qp_matrix = np.column_stack((q_list,p_list))
    cov_qp = np.cov(np.transpose(qp_matrix))
    eigenvalues,eigenvectors = LA.eig(cov_qp)
    eigenvalues_list = eigenvalues.tolist()
    max_eigenvalue = np.max(eigenvalues)
    max_eigenvalue_index = eigenvalues_list.index(max_eigenvalue)
    min_eigenvalue = np.min(eigenvalues)
    max_eigenvector = eigenvectors[:][max_eigenvalue_index]
    angle = np.arctan2(max_eigenvector[1],max_eigenvector[0])
    if angle < 0:
        angle = angle + 2*np.pi
    theta_grid = np.linspace(0,2*np.pi,Npts,endpoint=False)
    phi = angle
    a = N_sigma*np.sqrt(max_eigenvalue)
    b = N_sigma*np.sqrt(min_eigenvalue)
    ellipse_x_r = a*np.cos(theta_grid)
    ellipse_y_r = b*np.sin(theta_grid)
    R = np.array([[np.cos(phi),np.sin(phi)],[-np.sin(phi),np.cos(phi)]])
    r_ellipse = np.dot(np.column_stack((ellipse_x_r,ellipse_y_r)),np.transpose(R))
    q_ellipse = r_ellipse[:,0] + np.mean(q_list)
    p_ellipse = r_ellipse[:,1] + np.mean(p_list)
    q_ellipse = np.array(q_ellipse)
    p_ellipse = np.array(p_ellipse)
    return q_ellipse,p_ellipse
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import matplotlib
# Nreps_list = [10,20,30,40,50,60,70,80,90,\
#               100,200,300,400,500,600,700,800,900,\
#               1000,2000,3000,4000,5000,6000,7000,8000,9000,\
#               10000,20000,30000,40000,40200]
Nreps_list = np.geomspace(start=10,stop=40200,num=1000,endpoint=True)
for i in range(len(Nreps_list)):
    Nreps_list[i] = int(Nreps_list[i])
datadir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3/'
datafiles = [\
                  'sbdb_query_results_delcols_meanplanes_objct179_amin34.79_amax40.524_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct149_amin40.524_amax42_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct215_amin42_amax43_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct394_amin43_amax44_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct281_amin44_amax45_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct408_amin45_amax48_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct448_amin45_amax50_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct296_amin50_amax80_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct392_amin50_amax150_Nrep40200.txt',\
                  'sbdb_query_results_delcols_meanplanes_objct2058_amin34.79_amax150_Nrep40200.txt',\
                  ]
objct_list = [179,149,215,394,281,408,448,296,392,2058]
amin_list = (34.79, 40.524,42,43,44,45,45,50,50, 34.79)
amax_list = (40.524,42,    43,44,45,48,50,80,150,150)
Nbin = len(amin_list)
Nreps_options = len(Nreps_list)
imin_array = np.zeros((Nreps_options,Nbin))
imax_array = np.zeros((Nreps_options,Nbin))
Omegamin_array =  np.zeros((Nreps_options,Nbin))
Omegamax_array = np.zeros((Nreps_options,Nbin))
for ibin in range(Nbin):
    objct = objct_list[ibin]
    amin = amin_list[ibin]
    amax = amax_list[ibin]
    datafile = datafiles[ibin]
    df = pd.read_csv(datafile)
    all_q = df['q_vm'].tolist()
    all_p = df['p_vm'].tolist()
    for ireps_options in range(Nreps_options):
        Nreps = Nreps_list[ireps_options]
        Nreps =  int(Nreps)
        q_list = all_q[0:Nreps]
        p_list = all_p[0:Nreps]
        i_min_degrees,i_max_degrees,Omega_min_degrees,Omega_max_degrees = \
            covariance_ellipse(q_list,p_list,N_sigma=1)
        imin_array[ireps_options,ibin] = i_min_degrees
        imax_array[ireps_options,ibin] = i_max_degrees
        Omegamin_array[ireps_options,ibin] = Omega_min_degrees
        Omegamax_array[ireps_options,ibin] = Omega_max_degrees
xlist = Nreps_list
matplotlib.rcParams['pdf.fonttype'] = 42 # makes text editable in svg and pdf
matplotlib.rcParams['ps.fonttype'] = 42 # makes text editable in svg and pdf
for ibin in range(Nbin):
    ylist = []
    objct = objct_list[ibin]
    amin = amin_list[ibin]
    amax = amax_list[ibin]
    # fig,axs = plt.subplots(ncols=1,nrows=1,figsize=(2,1))
    imin_list = imin_array[:,ibin]
    imax_list = imax_array[:,ibin]
    Omegamin_list = Omegamin_array[:,ibin]
    Omegamax_list = Omegamax_array[:,ibin]
    ylist = np.abs(imin_list-imin_list[-1])/np.abs(imin_list[-1])
    plt.loglog(Nreps_list,ylist,'b.',markersize=0.5)
    plt.xlabel('Repetition count')
    plt.ylabel('|$i_{min}$-$i_{min,40200}$| / |$i_{min,40200}$|')
    plt.title('$i_{min}$ convergence, '+str(amin)+'-'+str(amax)+' au')
    # plt.rcParams["figure.figsize"] = (2,1)
    savename = 'convergence_imin_amin'+str(amin)+'_amax'+str(amax)
    # print(savename)
    plt.grid()
    plt.savefig(savename+'.pdf',format='pdf',transparent=True)
    plt.savefig(savename+'.svg',format='svg',transparent=True)
    plt.show()
    ylist = np.abs(imax_list-imax_list[-1])/np.abs(imax_list[-1])
    plt.loglog(Nreps_list,ylist,'r.',markersize=0.5)
    plt.xlabel('Repetition count')
    plt.ylabel('|$i_{max}$-$i_{max,40200}$| / |$i_{max,40200}$|')
    plt.title('$i_{max}$ convergence, '+str(amin)+'-'+str(amax)+' au')
    # plt.rcParams["figure.figsize"] = (2,1)
    savename = 'convergence_imax_amin'+str(amin)+'_amax'+str(amax)
    # print(savename)
    plt.grid()
    plt.savefig(savename+'.pdf',format='pdf',transparent=True)
    plt.savefig(savename+'.svg',format='svg',transparent=True)
    plt.show()
    ylist = np.abs(Omegamin_list-Omegamin_list[-1])/np.abs(Omegamin_list[-1])
    plt.loglog(Nreps_list,ylist,'k.',markersize=0.5)
    plt.xlabel('Repetition count')
    plt.ylabel('|$\Omega_{min}$-$\Omega_{min,40200}$| / |$\Omega_{min,40200}$|')
    plt.title('$\Omega_{min}$ convergence, '+str(amin)+'-'+str(amax)+' au')
    # plt.rcParams["figure.figsize"] = (2,1)
    savename = 'convergence_Omegamin_amin'+str(amin)+'_amax'+str(amax)
    # print(savename)
    plt.grid()
    plt.savefig(savename+'.pdf',format='pdf',transparent=True)
    plt.savefig(savename+'.svg',format='svg',transparent=True)
    plt.show()
    ylist = np.abs(Omegamax_list-Omegamax_list[-1])/np.abs(Omegamax_list[-1])
    plt.loglog(Nreps_list,ylist,'m.',markersize=0.5)
    plt.xlabel('Repetition count')
    plt.ylabel('|$\Omega_{max}$-$\Omega_{max,40200}$| / |$\Omega_{max,40200}$|')
    plt.title('$\Omega_{max}$ convergence, '+str(amin)+'-'+str(amax)+' au')
    # plt.rcParams["figure.figsize"] = (2,1)
    savename = 'convergence_Omegamax_amin'+str(amin)+'_amax'+str(amax)
    # print(savename)
    plt.grid()
    plt.savefig(savename+'.pdf',format='pdf',transparent=True)
    plt.savefig(savename+'.svg',format='svg',transparent=True)
    plt.show()
