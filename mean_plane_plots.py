#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 09:59:29 2021
@author: iggymatheson
"""
#%%
def setup_iOmega_axes():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rcParams['pdf.fonttype'] = 42 # makes text editable in svg and pdf
    matplotlib.rcParams['ps.fonttype'] = 42 # makes text editable in svg and pdf
    fig,axs = plt.subplots(ncols=2,nrows=2,figsize=(3.5,2.5),\
                           gridspec_kw={'width_ratios':[2,1]})
    axs[0,0].set_xlim((34,50))
    axs[0,1].set_xlim((50,150))
    axs[1,0].set_xlim((34,50))
    axs[1,1].set_xlim((50,150))
    axs[0,0].set_ylim((0,18))
    axs[0,1].set_ylim((0,18))
    axs[1,0].set_ylim((0,360))
    axs[1,1].set_ylim((0,360))
    axs[0,0].set_xticks((34,36,38,40,42,44,46,48,50))
    axs[0,1].set_xticks((50,100,150))
    axs[1,0].set_xticks((34,36,38,40,42,44,46,48,50))
    axs[1,1].set_xticks((50,100,150))
    axs[0,0].set_yticks((0,2,4,6,8,10,12,14,16,18))
    axs[0,1].set_yticks((0,2,4,6,8,10,12,14,16,18))
    axs[1,0].set_yticks((0,60,120,180,240,300,360))
    axs[1,1].set_yticks((0,60,120,180,240,300,360))
    axs[0,0].set_xticklabels('')
    axs[0,1].set_xticklabels('')
    axs[0,1].set_yticklabels('')
    axs[1,1].set_yticklabels('')
    axs[0,0].set_ylabel('inclination (deg)',fontsize=8)
    axs[1,0].set_ylabel('Ω (deg)',fontsize=8)
    axs[0,0].plot([34,50],[1.578694,1.578694],color='black',linestyle='-',linewidth=0.3)
    axs[0,1].plot([50,150],[1.578694,1.578694],color='black',linestyle='-',linewidth=0.3)
    axs[1,0].plot([34,50],[107.582222,107.582222],color='black',linestyle='-',linewidth=0.3)
    axs[1,1].plot([50,150],[107.582222,107.582222],color='black',linestyle='-',linewidth=0.3)
    # axs[1,0].set_xlabel('a (au)')
    labels = axs[1,0].get_xticks().tolist()
    labels[-1] = ''
    axs[1,0].set_xticklabels(labels)
    axs[0,0].tick_params(direction="in")
    axs[0,1].tick_params(direction="in")
    axs[1,0].tick_params(direction="in")
    axs[1,1].tick_params(direction="in")
    axs[0,0].tick_params(top=True,right=True)
    axs[0,1].tick_params(top=True,right=True)
    axs[1,0].tick_params(top=True,right=True)
    axs[1,1].tick_params(top=True,right=True)
    axs[0,0].tick_params(labelsize='x-small')
    axs[0,1].tick_params(labelsize='x-small')
    axs[1,0].tick_params(labelsize='x-small')
    axs[1,1].tick_params(labelsize='x-small')
    fig.tight_layout(h_pad=0.5,w_pad=0)
    fig.text(0.5,0.02,'a (au)',ha='center',fontsize=8)
    return fig,axs
#%%
def setup_qp_axes():
    import matplotlib.pyplot as plt
    fig,axs = plt.subplots(ncols=3,nrows=3,figsize=(6.5,6.5),\
                           gridspec_kw={'width_ratios':[1,1,1]})
    ## this set of commands is for all plot versions
    axs[2,1].set_xlabel('$q_0\cdot 1000$')
    axs[1,0].set_ylabel('$p_0\cdot 1000$')
    axs[0,0].tick_params(direction="in")
    axs[0,1].tick_params(direction="in")
    axs[0,2].tick_params(direction="in")
    axs[1,0].tick_params(direction="in")
    axs[1,1].tick_params(direction="in")
    axs[1,2].tick_params(direction="in")
    axs[2,0].tick_params(direction="in")
    axs[2,1].tick_params(direction="in")
    axs[2,2].tick_params(direction="in")
    axs[0,0].tick_params(top=True,right=True)
    axs[0,1].tick_params(top=True,right=True)
    axs[0,2].tick_params(top=True,right=True)
    axs[1,0].tick_params(top=True,right=True)
    axs[1,1].tick_params(top=True,right=True)
    axs[1,2].tick_params(top=True,right=True)
    axs[2,0].tick_params(top=True,right=True)
    axs[2,1].tick_params(top=True,right=True)
    axs[2,2].tick_params(top=True,right=True)
    axs[0,0].tick_params(labelsize='x-small')
    axs[0,1].tick_params(labelsize='x-small')
    axs[0,2].tick_params(labelsize='x-small')
    axs[1,0].tick_params(labelsize='x-small')
    axs[1,1].tick_params(labelsize='x-small')
    axs[1,2].tick_params(labelsize='x-small')
    axs[2,0].tick_params(labelsize='x-small')
    axs[2,1].tick_params(labelsize='x-small')
    axs[2,2].tick_params(labelsize='x-small')
    # axs[0,0].set_title('34.79-40.525 au')
    # axs[0,1].set_title('40.525-42 au')
    # axs[0,2].set_title('42-43 au')
    # axs[1,0].set_title('43-44 au')
    # axs[1,1].set_title('44-45 au')
    # axs[1,2].set_title('45-48 au')
    # axs[2,0].set_title('45-50 au')
    # axs[2,1].set_title('50-80 au')
    # axs[2,2].set_title('50-150 au')
    axs[0,0].text(0.85,0.9,'(a)',fontsize=10,transform=axs[0,0].transAxes)
    axs[0,1].text(0.85,0.9,'(b)',fontsize=10,transform=axs[0,1].transAxes)
    axs[0,2].text(0.85,0.9,'(c)',fontsize=10,transform=axs[0,2].transAxes)
    axs[1,0].text(0.85,0.9,'(d)',fontsize=10,transform=axs[1,0].transAxes)
    axs[1,1].text(0.85,0.9,'(e)',fontsize=10,transform=axs[1,1].transAxes)
    axs[1,2].text(0.85,0.9,'(f)',fontsize=10,transform=axs[1,2].transAxes)
    axs[2,0].text(0.85,0.9,'(g)',fontsize=10,transform=axs[2,0].transAxes)
    axs[2,1].text(0.85,0.9,'(h)',fontsize=10,transform=axs[2,1].transAxes)
    axs[2,2].text(0.85,0.9,'(i)',fontsize=10,transform=axs[2,2].transAxes)
    axs[0,0].grid(True)
    axs[0,1].grid(True)
    axs[0,2].grid(True)
    axs[1,0].grid(True)
    axs[1,1].grid(True)
    axs[1,2].grid(True)
    axs[2,0].grid(True)
    axs[2,1].grid(True)
    axs[2,2].grid(True)
    axs[0,0].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[0,1].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[0,2].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[1,0].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[1,1].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[1,2].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[2,0].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[2,1].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[2,2].plot((-500,500),(0,0),color='dimgray',linewidth=1)
    axs[0,0].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[0,1].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[0,2].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[1,0].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[1,1].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[1,2].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[2,0].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[2,1].plot((0,0),(-500,500),color='dimgray',linewidth=1)
    axs[2,2].plot((0,0),(-500,500),color='dimgray',linewidth=1)

    # this set of commands is for equal grid, both bootstrap and vm17method plotted
    x_span = 400
    x0 = -175
    x1 = x_span + x0
    y_span = 400
    y0 = -275
    y1 = y_span + y0
    axs[0,0].set_xlim((x0,x1))
    axs[0,0].set_ylim((y0,y1))
    axs[0,1].set_xlim((x0,x1))
    axs[0,1].set_ylim((y0,y1))
    axs[0,2].set_xlim((x0,x1))
    axs[0,2].set_ylim((y0,y1))
    axs[1,0].set_xlim((x0,x1))
    axs[1,0].set_ylim((y0,y1))
    axs[1,1].set_xlim((x0,x1))
    axs[1,1].set_ylim((y0,y1))
    axs[1,2].set_xlim((x0,x1))
    axs[1,2].set_ylim((y0,y1))
    axs[2,0].set_xlim((x0,x1))
    axs[2,0].set_ylim((y0,y1))
    axs[2,1].set_xlim((x0,x1))
    axs[2,1].set_ylim((y0,y1))
    axs[2,2].set_xlim((x0,x1))
    axs[2,2].set_ylim((y0,y1))
    axs[0,0].set_xticks((-100,0,100,200))
    axs[0,1].set_xticks((-100,0,100,200))
    axs[0,2].set_xticks((-100,0,100,200))
    axs[1,0].set_xticks((-100,0,100,200))
    axs[1,1].set_xticks((-100,0,100,200))
    axs[1,2].set_xticks((-100,0,100,200))
    axs[2,0].set_xticks((-100,0,100,200))
    axs[2,1].set_xticks((-100,0,100,200))
    axs[2,2].set_xticks((-100,0,100,200))
    axs[0,0].set_yticks((-200,-100,0,100))
    axs[0,1].set_yticks((-200,-100,0,100))
    axs[0,2].set_yticks((-200,-100,0,100))
    axs[1,0].set_yticks((-200,-100,0,100))
    axs[1,1].set_yticks((-200,-100,0,100))
    axs[1,2].set_yticks((-200,-100,0,100))
    axs[2,0].set_yticks((-200,-100,0,100))
    axs[2,1].set_yticks((-200,-100,0,100))
    axs[2,2].set_yticks((-200,-100,0,100))
    axs[0,0].set_xticklabels('')
    axs[0,1].set_xticklabels('')
    axs[0,2].set_xticklabels('')
    axs[1,0].set_xticklabels('')
    axs[1,1].set_xticklabels('')
    axs[1,2].set_xticklabels('')
    axs[2,0].set_xticklabels('')
    axs[2,1].set_xticklabels('')
    axs[2,2].set_xticklabels('')
    axs[0,0].set_yticklabels('')
    axs[0,1].set_yticklabels('')
    axs[0,2].set_yticklabels('')
    axs[1,0].set_yticklabels('')
    axs[1,1].set_yticklabels('')
    axs[1,2].set_yticklabels('')
    axs[2,0].set_yticklabels('')
    axs[2,1].set_yticklabels('')
    axs[2,2].set_yticklabels('')
    labels = axs[0,0].get_yticks().tolist()
    # labels[-1] = ''
    axs[0,0].set_yticklabels(labels)
    labels = axs[1,0].get_yticks().tolist()
    # labels[-1] = ''
    axs[1,0].set_yticklabels(labels)
    labels = axs[2,0].get_yticks().tolist()
    # labels[-1] = ''
    axs[2,0].set_yticklabels(labels)
    labels = axs[2,0].get_xticks().tolist()
    # labels[-1] = ''
    axs[2,0].set_xticklabels(labels)
    labels = axs[2,1].get_xticks().tolist()
    # labels[-1] = ''
    axs[2,1].set_xticklabels(labels)
    labels = axs[2,2].get_xticks().tolist()
    # labels[-1] = ''
    axs[2,2].set_xticklabels(labels)
    fig.tight_layout(h_pad=0,w_pad=0.0)

    # ## this set of commands is for separategrids, both bootstrap and vm17method plotted
    # axs[0,0].set_xlim((-175,125))
    # axs[0,0].set_ylim((-110,190))
    # axs[0,1].set_xlim((-325,150))
    # axs[0,1].set_ylim((-200,275))
    # axs[0,2].set_xlim((-90,85))
    # axs[0,2].set_ylim((-50,125))
    # axs[1,0].set_xlim((-50,45))
    # axs[1,0].set_ylim((-10,85))
    # axs[1,1].set_xlim((-45,45))
    # axs[1,1].set_ylim((-10,80))
    # axs[1,2].set_xlim((-50,70))
    # axs[1,2].set_ylim((-40,80))
    # axs[2,0].set_xlim((-50,70))
    # axs[2,0].set_ylim((-55,65))
    # axs[2,1].set_xlim((-170,155))
    # axs[2,1].set_ylim((-125,200))
    # axs[2,2].set_xlim((-200,175))
    # axs[2,2].set_ylim((-200,175))
    # fig.tight_layout(h_pad=0.2,w_pad=0.2)

    # ## this set of commands is for separategrids, just bootstrap plotted
    # axs[0,0].set_xlim((-175,125))
    # axs[0,0].set_ylim((-110,190))
    # axs[0,1].set_xlim((-325,150))
    # axs[0,1].set_ylim((-200,275))
    # axs[0,2].set_xlim((-30,30))
    # axs[0,2].set_ylim((10,70))
    # axs[1,0].set_xlim((-15,10))
    # axs[1,0].set_ylim((23,48))
    # axs[1,1].set_xlim((-23,18))
    # axs[1,1].set_ylim((15,60))
    # axs[1,2].set_xlim((-15,40))
    # axs[1,2].set_ylim((-20,55))
    # axs[2,0].set_xlim((-15,40))
    # axs[2,0].set_ylim((-25,50))
    # axs[2,1].set_xlim((-170,155))
    # axs[2,1].set_ylim((-125,200))
    # axs[2,2].set_xlim((-200,175))
    # axs[2,2].set_ylim((-200,175))
    # fig.tight_layout(h_pad=0.2,w_pad=0.2)

    # ## this set of commands is for separategrids, just vm17method plotted
    # axs[0,0].set_xlim((-175,125))
    # axs[0,0].set_ylim((-110,190))
    # axs[0,1].set_xlim((-300,50))
    # axs[0,1].set_ylim((-150,200))
    # axs[0,2].set_xlim((-80,80))
    # axs[0,2].set_ylim((-40,120))
    # axs[1,0].set_xlim((-50,40))
    # axs[1,0].set_ylim((-7,83))
    # axs[1,1].set_xlim((-45,45))
    # axs[1,1].set_ylim((-10,80))
    # axs[1,2].set_xlim((-50,70))
    # axs[1,2].set_ylim((-40,80))
    # axs[2,0].set_xlim((-40,70))
    # axs[2,0].set_ylim((-50,60))
    # axs[2,1].set_xlim((-125,125))
    # axs[2,1].set_ylim((-100,150))
    # axs[2,2].set_xlim((-90,120))
    # axs[2,2].set_ylim((-100,110))
    # fig.tight_layout(h_pad=0.2,w_pad=0.2)

    return fig,axs
#%%
def add_origin_invariable_laplace_to_qp(fig,axs,amin_list,amax_list):
    import numpy as np
    import pandas as pd
    i_invariable = 1.578694
    Omega_invariable = 107.582222
    i_invariable = np.radians(i_invariable)
    Omega_invariable = np.radians(Omega_invariable)
    q_invariable = np.sin(i_invariable)*np.cos(Omega_invariable)
    p_invariable = np.sin(i_invariable)*np.sin(i_invariable)
    q_invariable = q_invariable * 1000
    p_invariable = p_invariable * 1000
    axs[0,0].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[0,1].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[0,2].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[1,0].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[1,1].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[1,2].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[2,0].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[2,1].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    axs[2,2].plot(q_invariable,p_invariable,marker='x',markersize=5,color='black')
    laplace_file_0 = 'a_p0q0_i0Omega0_34.79_40.524_20220122.txt'
    laplace_file_1 = 'a_p0q0_i0Omega0_40.524_42_20220122.txt'
    laplace_file_2 = 'a_p0q0_i0Omega0_42_43_20220122.txt'
    laplace_file_3 = 'a_p0q0_i0Omega0_43_44_20220122.txt'
    laplace_file_4 = 'a_p0q0_i0Omega0_44_45_20220122.txt'
    laplace_file_5 = 'a_p0q0_i0Omega0_45_48_20220122.txt'
    laplace_file_6 = 'a_p0q0_i0Omega0_45_50_20220122.txt'
    laplace_file_7 = 'a_p0q0_i0Omega0_50_80_20220122.txt'
    laplace_file_8 = 'a_p0q0_i0Omega0_50_150_20220122.txt'
    laplace_dir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3/'
    laplace_file = laplace_dir + laplace_file_0
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[0,0].scatter(q0*1000,p0*1000,marker='^',s=0.5,color='blue')
    laplace_file = laplace_dir + laplace_file_1
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[0,1].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_2
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[0,2].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_3
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[1,0].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_4
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[1,1].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_5
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[1,2].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_6
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[2,0].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_7
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[2,1].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    laplace_file = laplace_dir + laplace_file_8
    df = pd.read_csv(laplace_file,header=None)
    p0 = df.iloc[:,1].tolist()
    q0 = df.iloc[:,2].tolist()
    p0 = np.array(p0)
    q0 = np.array(q0)
    axs[2,2].scatter(q0*1000,p0*1000,marker='^',s=2,color='blue')
    # midplane_dir = '/Users/iggymatheson/Documents_off_iCloud/school/jpl2021/much_of_a_muchness/'
    # midplane_files = ['midplane_34.79_40.525.txt',\
    #                   'midplane_40.525_42.txt',\
    #                   'midplane_42_43.txt',\
    #                   'midplane_43_44.txt',\
    #                   'midplane_44_45.txt',\
    #                   'midplane_45_48.txt',\
    #                   'midplane_45_50.txt',\
    #                   'midplane_50_80.txt',\
    #                   'midplane_50_150.txt']
    midplane_dir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3/'
    midplane_files = [\
                      'sbdb_query_results_delcols__objct179_amin34.79_amax40.524_nominal.txt',\
                      'sbdb_query_results_delcols__objct149_amin40.524_amax42_nominal.txt',\
                      'sbdb_query_results_delcols__objct215_amin42_amax43_nominal.txt',\
                      'sbdb_query_results_delcols__objct394_amin43_amax44_nominal.txt',\
                      'sbdb_query_results_delcols__objct281_amin44_amax45_nominal.txt',\
                      'sbdb_query_results_delcols__objct408_amin45_amax48_nominal.txt',\
                      'sbdb_query_results_delcols__objct448_amin45_amax50_nominal.txt',\
                      'sbdb_query_results_delcols__objct296_amin50_amax80_nominal.txt',\
                      'sbdb_query_results_delcols__objct392_amin50_amax150_nominal.txt',\
                      ]
    file = midplane_dir + midplane_files[0]
    a0 = 0
    a1 = 0
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[1]
    a0 = 0
    a1 = 1
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[2]
    a0 = 0
    a1 = 2
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[3]
    a0 = 1
    a1 = 0
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[4]
    a0 = 1
    a1 = 1
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[5]
    a0 = 1
    a1 = 2
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[6]
    a0 = 2
    a1 = 0
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[7]
    a0 = 2
    a1 = 1
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    file = midplane_dir + midplane_files[8]
    a0 = 2
    a1 = 2
    df = pd.read_csv(file,delim_whitespace=True)
    q = df['q_vm'].tolist()
    p = df['p_vm'].tolist()
    q = np.array(q)
    p = np.array(p)
    q = q[0]*1000
    p = p[0]*1000
    axs[a0,a1].plot(q,p,marker='+',markersize=5,color='saddlebrown')
    return fig,axs
#%%
def add_vm17method_to_qp(fig,axs):
    import pandas as pd
    import numpy as np
    # datadir = '/Users/iggymatheson/Documents_off_iCloud/school/jpl2021/check_overruns_mpc21sv20opp1_maxiters107/'
    # datafiles = ['bootstrap_vm17method_supercompilation_34.79_40.525_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_40.525_42_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_42_43_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_43_44_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_44_45_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_45_48_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_45_50_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_50_80_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_50_150_mpc21sv20opp1.txt']
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
                      ]
    file = datadir + datafiles[0]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[0,0].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[1]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[0,1].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[2]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[0,2].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[3]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[1,0].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[4]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[1,1].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[5]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[1,2].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[6]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[2,0].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[7]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[2,1].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    file = datadir + datafiles[8]
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    axs[2,2].scatter(q_qp,p_qp,marker='.',s=0.5,color='tan')
    return fig,axs
#%%
def add_ellipses_vm17method_to_qp(fig,axs):
    import pandas as pd
    import numpy as np
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
                      ]
    file = datadir + datafiles[0]
    a0 = 0
    a1 = 0
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[1]
    a0 = 0
    a1 = 1
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[2]
    a0 = 0
    a1 = 2
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[3]
    a0 = 1
    a1 = 0
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[4]
    a0 = 1
    a1 = 1
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[5]
    a0 = 1
    a1 = 2
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[6]
    a0 = 2
    a1 = 0
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[7]
    a0 = 2
    a1 = 1
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    file = datadir + datafiles[8]
    a0 = 2
    a1 = 2
    df = pd.read_csv(file)
    q_qp = df['q_vm'].tolist()
    p_qp = df['p_vm'].tolist()
    q_qp = np.array(q_qp)*1000
    p_qp = np.array(p_qp)*1000
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=1)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=2)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,N_sigma=3)
    axs[a0,a1].plot(q_ellipse,p_ellipse,linewidth=1,color='saddlebrown')
    return fig,axs
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
def add_vm17results_to_iOmega(fig,axs):
    import numpy as np
    # amin_list = (35,  40.3,42,43,44,45,45,50,50)
    # amax_list = (40.3,42,  43,44,45,48,50,80,150)
    # imin_degrees_list = (1.25,5.76,0.24,1.10,1.82,0.19,0.53,5.35,3.74)
    i_degrees_list_qp = (3.8,9.6,1.6,1.8,2.7,1.1,1.6,9.1,6.8)
    # imax_degrees_list = (8.50,17.47,3.07,2.78,3.65,2.59,2.93,15.81,12.38)
    # Omegamin_degrees_list = (0,271,347,51,65,316,327,188,191)
    Omega_degrees_list_qp = (134,292,51,79,81,29,20,227,222)
    # Omegamax_degrees_list = (360,315,95,105,99,81,75,251,264)
    a0 = 0
    a1 = 0
    i0 = i_degrees_list_qp[0]
    Omega0 = Omega_degrees_list_qp[0]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 0
    a1 = 1
    i0 = i_degrees_list_qp[1]
    Omega0 = Omega_degrees_list_qp[1]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 0
    a1 = 2
    i0 = i_degrees_list_qp[2]
    Omega0 = Omega_degrees_list_qp[2]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 1
    a1 = 0
    i0 = i_degrees_list_qp[3]
    Omega0 = Omega_degrees_list_qp[3]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 1
    a1 = 1
    i0 = i_degrees_list_qp[4]
    Omega0 = Omega_degrees_list_qp[4]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 1
    a1 = 2
    i0 = i_degrees_list_qp[5]
    Omega0 = Omega_degrees_list_qp[5]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 2
    a1 = 0
    i0 = i_degrees_list_qp[6]
    Omega0 = Omega_degrees_list_qp[6]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 2
    a1 = 1
    i0 = i_degrees_list_qp[7]
    Omega0 = Omega_degrees_list_qp[7]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 2
    a1 = 2
    i0 = i_degrees_list_qp[8]
    Omega0 = Omega_degrees_list_qp[8]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    return fig,axs
#%%
def save_fig(fig,axs,savename):
    import matplotlib.pyplot as plt
    plt.savefig(savename+'.pdf',format='pdf',transparent=True)
    plt.savefig(savename+'.svg',format='svg',transparent=True)
    plt.show()
    return
#%%
def draw_laplace_on_iOmega(fig,axs):
    import pandas as pd
    import os
    opendir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3'
    os.chdir(opendir)
    openfile = 'a_p0q0_i0Omega0_30_150_20220122.txt'
    df = pd.read_csv(openfile)
    a_list = df.iloc[:,0].tolist()
    i_list = df.iloc[:,3].tolist()
    Omega_list = df.iloc[:,4].tolist()
    axs[0,0].plot(a_list,i_list,linewidth=0.5,color='blue',label='Barycentric')
    axs[0,1].plot(a_list,i_list,linewidth=0.5,color='blue')
    axs[1,0].plot(a_list,Omega_list,linewidth=0.5,color='blue',label='Barycentric')
    axs[1,1].plot(a_list,Omega_list,linewidth=0.5,color='blue')
    return fig,axs
#%%
def draw_bins_vm17method_on_iOmega(fig,axs,amin_list,amax_list):
    # amin_list = (34.79, 40.525,42,43,44,45,45,50,50)
    # amax_list = (40.525,42,    43,44,45,48,50,80,150)
    # vm17method_list = ['bootstrap_vm17method_supercompilation_34.79_40.525_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_40.525_42_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_42_43_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_43_44_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_44_45_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_45_48_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_45_50_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_50_80_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_50_150_mpc21sv20opp1.txt',\
    #              'bootstrap_vm17method_supercompilation_34.79_150_mpc21sv20opp1.txt']
    vm17method_list = [\
                      'sbdb_query_results_delcols_meanplanes_objct179_amin34.79_amax40.524_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct149_amin40.524_amax42_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct215_amin42_amax43_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct394_amin43_amax44_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct281_amin44_amax45_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct408_amin45_amax48_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct448_amin45_amax50_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct296_amin50_amax80_Nrep40200.txt',\
                      'sbdb_query_results_delcols_meanplanes_objct392_amin50_amax150_Nrep40200.txt',\
                      ]
    objct_list = [179,149,215,394,281,408,448,296,392]
    # import matplotlib.pyplot as plt
    import os
    import numpy as np
    import pandas as pd
    os.chdir('/Users/iggymatheson/Documents_off_iCloud/mm22_v3')
    Nbins = len(amin_list)
    i_degrees_list_qp = []
    Omega_degrees_list_qp = []
    amid_list = []
    imin_degrees_list = []
    imax_degrees_list = []
    Omegamin_degrees_list = []
    Omegamax_degrees_list = []
    for i in range(Nbins):
        amin = amin_list[i]
        amax = amax_list[i]
        objct = objct_list[i]
        amid = (amin+amax)/2 # shift a bit for non-overlapping plotting
        if amin == 50:
            amid = amid + 2
        else:
            amid = amid + 0.2
        amid_list.append(amid)
        # midfile = 'midplane_' + str(amin) + '_' + str(amax) + '.txt'
        midfile = 'sbdb_query_results_delcols__objct' + str(objct) + '_amin' + str(amin) + \
            '_amax' + str(amax) + '_nominal.txt'
        df = pd.read_csv(midfile,delim_whitespace=True)
        i_qp_list = df['i_vm'].tolist()
        Omega_qp_list = df['Omega_vm'].tolist()
        i_qp = i_qp_list[0]
        Omega_qp = Omega_qp_list[0]
        i_degrees_list_qp.append(i_qp)
        Omega_degrees_list_qp.append(Omega_qp)
        vm17file = vm17method_list[i]
        df = pd.read_csv(vm17file)
        q_list = df['q_vm'].tolist()
        p_list = df['p_vm'].tolist()
        i_min_degrees,i_max_degrees,Omega_min_degrees,Omega_max_degrees = \
            covariance_ellipse(q_list,p_list,N_sigma=1)
        imin_degrees_list.append(i_min_degrees)
        imax_degrees_list.append(i_max_degrees)
        Omegamin_degrees_list.append(Omega_min_degrees)
        Omegamax_degrees_list.append(Omega_max_degrees)
    i_up_errs = []
    i_down_errs = []
    i_left_errs = []
    i_right_errs = []
    Omega_up_errs = []
    Omega_down_errs = []
    for i in range(Nbins):
        i_up_errs.append(np.abs(imax_degrees_list[i]-i_degrees_list_qp[i]))
        i_down_errs.append(np.abs(imin_degrees_list[i]-i_degrees_list_qp[i]))
        Omax = Omegamax_degrees_list[i]
        Omin = Omegamin_degrees_list[i]
        Obin = Omega_degrees_list_qp[i]
        if (Omin<=Obin<=Omax):
            Omega_up_errs.append(Omax-Obin)
            Omega_down_errs.append(Obin-Omin)
        if (Omax<=Omin<=Obin):
            Omega_up_errs.append(Omax+360-Obin)
            Omega_down_errs.append(Obin-Omin)
        if (Obin<=Omax<=Omin):
            Omega_up_errs.append(Omax-Obin)
            Omega_down_errs.append(Obin-(Omin-360))
        # Omega_up_errs.append(np.abs(Omegamax_degrees_list[i]-Omega_degrees_list_qp[i]))
        # Omega_down_errs.append(np.abs(Omegamin_degrees_list[i]-Omega_degrees_list_qp[i]))
    for i in range(Nbins):
        i_left_errs.append(np.abs(amid_list[i]-amin_list[i]))
        i_right_errs.append(np.abs(amid_list[i]-amax_list[i]))
    Omega_left_errs = i_left_errs
    Omega_right_errs = i_right_errs
    i_x_errs = [i_left_errs,i_right_errs]
    i_y_errs = [i_down_errs,i_up_errs]
    Omega_x_errs = [Omega_left_errs,Omega_right_errs]
    Omega_y_errs = [Omega_down_errs,Omega_up_errs]
    axs[0,0].errorbar(amid_list, i_degrees_list_qp, xerr=i_x_errs, yerr=i_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='saddlebrown')
    axs[0,1].errorbar(amid_list, i_degrees_list_qp, xerr=i_x_errs, yerr=i_y_errs,\
        marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
        color='saddlebrown')
    axs[1,0].errorbar(amid_list, Omega_degrees_list_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='saddlebrown')
    axs[1,1].errorbar(amid_list, Omega_degrees_list_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='saddlebrown')
    Omega2_qp = np.array(Omega_degrees_list_qp) + 360
    axs[1,0].errorbar(amid_list, Omega2_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='saddlebrown')
    axs[1,1].errorbar(amid_list, Omega2_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='saddlebrown')
    return fig,axs,i_degrees_list_qp,Omega_degrees_list_qp,imin_degrees_list,imax_degrees_list,Omegamin_degrees_list,Omegamax_degrees_list
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
def draw_bins_vm17results_on_iOmega(fig,axs):
    amin_list = (35,  40.3,42,43,44,45,45,50,50)
    amax_list = (40.3,42,  43,44,45,48,50,80,150)
    imin_degrees_list = (1.25,5.76,0.24,1.10,1.82,0.19,0.53,5.35,3.74)
    i_degrees_list_qp = (3.8,9.6,1.6,1.8,2.7,1.1,1.6,9.1,6.8)
    imax_degrees_list = (8.50,17.47,3.07,2.78,3.65,2.59,2.93,15.81,12.38)
    Omegamin_degrees_list = (0,271,347,51,65,316,327,188,191)
    Omega_degrees_list_qp = (134,292,51,79,81,29,20,227,222)
    Omegamax_degrees_list = (360,315,95,105,99,81,75,251,264)
    import os
    import numpy as np
    os.chdir('/Users/iggymatheson/Documents_off_iCloud/mm22_v3')
    Nbins = len(amin_list)
    amid_list = []
    for i in range(Nbins):
        amin = amin_list[i]
        amax = amax_list[i]
        amid = (amin+amax)/2 # shift a bit for non-overlapping plotting
        if amin == 50:
            amid = amid - 2
        else:
            amid = amid - 0.2
        amid_list.append(amid)
    i_up_errs = []
    i_down_errs = []
    i_left_errs = []
    i_right_errs = []
    Omega_up_errs = []
    Omega_down_errs = []
    for i in range(Nbins):
        i_up_errs.append(np.abs(imax_degrees_list[i]-i_degrees_list_qp[i]))
        i_down_errs.append(np.abs(imin_degrees_list[i]-i_degrees_list_qp[i]))
        Omax = Omegamax_degrees_list[i]
        Omin = Omegamin_degrees_list[i]
        Obin = Omega_degrees_list_qp[i]
        if (Omin<=Obin<=Omax):
            Omega_up_errs.append(Omax-Obin)
            Omega_down_errs.append(Obin-Omin)
        if (Omax<=Omin<=Obin):
            Omega_up_errs.append(Omax+360-Obin)
            Omega_down_errs.append(Obin-Omin)
        if (Obin<=Omax<=Omin):
            Omega_up_errs.append(Omax-Obin)
            Omega_down_errs.append(Obin-(Omin-360))
        # Omega_up_errs.append(np.abs(Omegamax_degrees_list[i]-Omega_degrees_list_qp[i]))
        # Omega_down_errs.append(np.abs(Omegamin_degrees_list[i]-Omega_degrees_list_qp[i]))
    for i in range(Nbins):
        i_left_errs.append(np.abs(amid_list[i]-amin_list[i]))
        i_right_errs.append(np.abs(amid_list[i]-amax_list[i]))
    Omega_left_errs = i_left_errs
    Omega_right_errs = i_right_errs
    i_x_errs = [i_left_errs,i_right_errs]
    i_y_errs = [i_down_errs,i_up_errs]
    Omega_x_errs = [Omega_left_errs,Omega_right_errs]
    Omega_y_errs = [Omega_down_errs,Omega_up_errs]
    axs[0,0].errorbar(amid_list, i_degrees_list_qp, xerr=i_x_errs, yerr=i_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='magenta')
    axs[0,1].errorbar(amid_list, i_degrees_list_qp, xerr=i_x_errs, yerr=i_y_errs,\
        marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
        color='magenta')
    axs[1,0].errorbar(amid_list, Omega_degrees_list_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='magenta')
    axs[1,1].errorbar(amid_list, Omega_degrees_list_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='magenta')
    Omega2_qp = np.array(Omega_degrees_list_qp) + 360
    axs[1,0].errorbar(amid_list, Omega2_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='magenta')
    axs[1,1].errorbar(amid_list, Omega2_qp, xerr=Omega_x_errs, yerr=Omega_y_errs,\
            marker='o',markersize=0.5,linewidth=0.5,linestyle='none',\
            color='magenta')
    return fig,axs
#%%
def add_vm17results_to_qp(fig,axs):
    import numpy as np
    # amin_list = (35,  40.3,42,43,44,45,45,50,50)
    # amax_list = (40.3,42,  43,44,45,48,50,80,150)
    # imin_degrees_list = (1.25,5.76,0.24,1.10,1.82,0.19,0.53,5.35,3.74)
    i_degrees_list_qp = (3.8,9.6,1.6,1.8,2.7,1.1,1.6,9.1,6.8)
    # imax_degrees_list = (8.50,17.47,3.07,2.78,3.65,2.59,2.93,15.81,12.38)
    # Omegamin_degrees_list = (0,271,347,51,65,316,327,188,191)
    Omega_degrees_list_qp = (134,292,51,79,81,29,20,227,222)
    # Omegamax_degrees_list = (360,315,95,105,99,81,75,251,264)
    a0 = 0
    a1 = 0
    i0 = i_degrees_list_qp[0]
    Omega0 = Omega_degrees_list_qp[0]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 0
    a1 = 1
    i0 = i_degrees_list_qp[1]
    Omega0 = Omega_degrees_list_qp[1]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 0
    a1 = 2
    i0 = i_degrees_list_qp[2]
    Omega0 = Omega_degrees_list_qp[2]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 1
    a1 = 0
    i0 = i_degrees_list_qp[3]
    Omega0 = Omega_degrees_list_qp[3]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 1
    a1 = 1
    i0 = i_degrees_list_qp[4]
    Omega0 = Omega_degrees_list_qp[4]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 1
    a1 = 2
    i0 = i_degrees_list_qp[5]
    Omega0 = Omega_degrees_list_qp[5]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 2
    a1 = 0
    i0 = i_degrees_list_qp[6]
    Omega0 = Omega_degrees_list_qp[6]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 2
    a1 = 1
    i0 = i_degrees_list_qp[7]
    Omega0 = Omega_degrees_list_qp[7]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    a0 = 2
    a1 = 2
    i0 = i_degrees_list_qp[8]
    Omega0 = Omega_degrees_list_qp[8]
    i0 = np.radians(i0)
    Omega0 = np.radians(Omega0)
    q0 = np.sin(i0)*np.cos(Omega0)*1000
    p0 = np.sin(i0)*np.sin(Omega0)*1000
    axs[a0,a1].scatter(q0,p0,marker='D',s=5,color='magenta')
    return fig,axs
#%%
def invariable_rejection_function_vm17method(Nsigma):
    import numpy as np
    import pandas as pd
    from shapely.geometry import Polygon, Point
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
                      ]
    invariable_rejection = []
    i_invariable = 1.578694
    Omega_invariable = 107.582222
    i_invariable = np.radians(i_invariable)
    Omega_invariable = np.radians(Omega_invariable)
    q_invariable = np.sin(i_invariable)*np.cos(Omega_invariable)
    p_invariable = np.sin(i_invariable)*np.sin(i_invariable)
    q_invariable = q_invariable * 1000
    p_invariable = p_invariable * 1000
    pt_invariable = Point(q_invariable,p_invariable)
    Nbin = len(datafiles)
    for ibin in range(Nbin):
        file = datadir + datafiles[ibin]
        df = pd.read_csv(file)
        q_qp = df['q_vm'].tolist()
        p_qp = df['p_vm'].tolist()
        q_qp = np.array(q_qp)*1000
        p_qp = np.array(p_qp)*1000
        q_ellipse,p_ellipse = covariance_ellipse_2(q_qp,p_qp,Nsigma)
        Ne = len(q_ellipse)
        linestring = []
        for i2 in range(Ne):
            pt = (q_ellipse[i2],p_ellipse[i2])
            linestring.append(pt)
        pt = (q_ellipse[0],p_ellipse[0])
        linestring.append(pt)
        poly = Polygon(linestring)
        checkstatus = pt_invariable.within(poly)
        invariable_rejection.append(checkstatus)
    return invariable_rejection
#%% make vm17fig3
do_this = 1
if do_this == 1:
    import numpy as np
    amin_list = (34.79, 40.524,42,43,44,45,45,50,50)
    amax_list = (40.524,42,    43,44,45,48,50,80,150)
    savename = 'iOmega_vm17results'
    fig,axs = setup_iOmega_axes()
    fig,axs = draw_laplace_on_iOmega(fig,axs)
    fig,axs,i_degrees_list,Omega_degrees_list,imin_degrees_list,imax_degrees_list,Omegamin_degrees_list_vm,Omegamax_degrees_list_vm = \
        draw_bins_vm17method_on_iOmega(fig,axs,amin_list,amax_list)
    fig,axs = draw_bins_vm17results_on_iOmega(fig,axs)
    save_fig(fig,axs,savename)
#%% make nine-bin plot suggested by RM, compute invariable plane rejection
do_this = 1
if do_this == 1:
    amin_list = (34.79, 40.524,42,43,44,45,45,50,50)
    amax_list = (40.524,42,    43,44,45,48,50,80,150)
    savename = 'qp_vm17method_equalgrid'
    fig,axs = setup_qp_axes()
    fig,axs = add_vm17method_to_qp(fig,axs)
    fig,axs = add_ellipses_vm17method_to_qp(fig,axs)
    fig,axs = add_origin_invariable_laplace_to_qp(fig,axs,amin_list,amax_list)
    fig,axs = add_vm17results_to_qp(fig,axs)
    save_fig(fig,axs,savename)
    print(amin_list)
    print(amax_list)
    invariable_rejection = invariable_rejection_function_vm17method(Nsigma=1)
    print(invariable_rejection)
    invariable_rejection = invariable_rejection_function_vm17method(Nsigma=2)
    print(invariable_rejection)
    invariable_rejection = invariable_rejection_function_vm17method(Nsigma=3)
    print(invariable_rejection)
# (34.79,  40.524, 42,    43,    44,    45,    45,    50,   50)
# (40.524, 42,     43,    44,    45,    48,    50,    80,   150)
# [False,  False,  False, False, False, False, False, True, True]
# [True,   False,  True,  False, False, True,  True,  True, True]
# [True,   True,   True,  False, True,  True,  True,  True, True]
#%% compute confidence intervals
do_this = 1
if do_this == 1:
    import pandas as pd
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
    for ibin in range(Nbin):
        amin = amin_list[ibin]
        amax = amax_list[ibin]
        objct = objct_list[ibin]
        midfile = 'sbdb_query_results_delcols__objct' + str(objct) + '_amin' + str(amin) + \
            '_amax' + str(amax) + '_nominal.txt'
        df = pd.read_csv(midfile,delim_whitespace=True)
        i_qp_list = df['i_vm'].tolist()
        Omega_qp_list = df['Omega_vm'].tolist()
        i_mid_degrees = i_qp_list[0]
        Omega_mid_degrees = Omega_qp_list[0]
        datafile = datafiles[ibin]
        df = pd.read_csv(datafile)
        q_list = df['q_vm'].tolist()
        p_list = df['p_vm'].tolist()
        i_min_degrees,i_max_degrees,Omega_min_degrees,Omega_max_degrees = \
            covariance_ellipse(q_list,p_list,N_sigma=1)
        print(amin,amax)
        print(round(i_min_degrees,2),round(i_mid_degrees,2),round(i_max_degrees,2))
        print(round(Omega_min_degrees,2),round(Omega_mid_degrees,2),round(Omega_max_degrees,2))
