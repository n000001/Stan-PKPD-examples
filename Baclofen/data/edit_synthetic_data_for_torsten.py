#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 16:04:29 2022

@author: tng

Adaptation of mrgsove-simulated data for Stan input
"""

import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns

def edit_sim_data(ev_data_path, data_dir, sim_data_file, save_path) :

    sim_data_path = data_dir + sim_data_file 
    
    sim_data = pd.read_json(sim_data_path)
    sim_data['TIME'] = np.round(sim_data.TIME.values.astype('float'), 2)
    
    ev_data = pd.read_csv(ev_data_path).drop('Unnamed: 0', axis=1)
    ev_data.rename({'DV': 'DV_emp'}, axis=1, inplace=True)
    ev_data.loc[ev_data['EVID'] == 1, ['DV_emp']] = -1
    ev_data['TIME'] = np.round(ev_data.TIME.values.astype('float'), 2)
    
    ev_data.sort_values(['ID', 'TIME'], inplace=True)
    sim_data.sort_values(['ID', 'TIME'], inplace=True)
    
    if not (np.array_equal(ev_data['ID'], 
                           sim_data['ID']) and np.array_equal(ev_data['TIME'], 
                                                              sim_data['TIME'])) :
                                                              print('discrepancy between data, check')
                                                                
                                                    
    data = ev_data.copy()
    data['DV'] = sim_data['DV']
    
    #For torsten
    data['RATE'] = 0 #bolus
    data['ADDL'] = 0
    data['CMT']= 1
    
    
    nSub = len(data.ID.unique())
    iObs = np.where(data['EVID'] == 0)[0].tolist()
    nObs = len(iObs)
    nObsSub = data.iloc[iObs].groupby('ID')['ID'].count().values.tolist()
    tObs = data.TIME.iloc[iObs].tolist()
    cObs = data.DV[iObs].tolist()
    
    cObs_emp = data.DV_emp[iObs].tolist()
    
    nt = len(data.TIME)
    nDataSub = data.groupby('ID')['ID'].count().tolist()
    nData = len(data)
    
    starts_bool = data.ID.diff().ne(0)
    starts = data.index[starts_bool].tolist()
    ends = data.index[starts_bool.shift(-1, fill_value=True)].tolist()
    
    out_dict = data.loc[:, data.columns != 'DATE'].to_dict(orient='list')
    
    out_dict['nSub'] = nSub
    out_dict['nObsSub'] = nObsSub
    out_dict['iObs'] = np.add(iObs, 1).tolist()
    out_dict['nObs'] = nObs
    out_dict['tObs'] = tObs
    out_dict['cObs'] = cObs
    
    out_dict['nData'] = nData
    out_dict['nDataSub'] = nDataSub
    out_dict['starts'] = np.add(starts, 1).tolist()
    out_dict['ends'] = np.add(ends, 1).tolist()
    
    out_dict['nt'] = nt
    ID_list = data.ID.unique()
    out_dict['ID_list'] = ID_list.tolist()
    
    Obs_df = data[['TIME', 'ID']].iloc[iObs].reset_index(drop=True)
    starts_bool_tObs = Obs_df.ID.diff().ne(0)
    starts_tObs = Obs_df.index[starts_bool_tObs].tolist()
    ends_tObs = Obs_df.index[starts_bool_tObs.shift(-1, fill_value=True)].tolist()
    #Not adding 1 because not used in Stan rather in python for plotting and all
    out_dict['starts_tObs'] = starts_tObs
    out_dict['ends_tObs'] = ends_tObs
    
    out_dict['CL_pop_meanPrior'] = 7
    out_dict['V1_pop_meanPrior'] = 70
    out_dict['ka_pop_meanPrior'] = 3
    
    
    with open(save_path, "w") as outfile:
        json.dump(out_dict, outfile)
    
    fig, ax = plt.subplots(nrows=np.ceil(nSub/6).astype('int'), ncols=6,
                           sharey=True, figsize=(30,60))
    for isub in range(nSub) :
        ax[isub//6][isub%6].scatter(x=data['TIME'][starts[isub]:ends[isub]+1],
                                    y=data['DV'][starts[isub]:ends[isub]+1],
                                    color='blue',
                                    label='dosing event simulated conc.')
        ax[isub//6][isub%6].scatter(x=tObs[starts_tObs[isub]:ends_tObs[isub]+1],
                                    y=cObs[starts_tObs[isub]:ends_tObs[isub]+1], 
                                    color='green',
                                    label='obs event simulated conc.')
        ax[isub//6][isub%6].scatter(x=tObs[starts_tObs[isub]:ends_tObs[isub]+1],
                                    y=cObs_emp[starts_tObs[isub]:ends_tObs[isub]+1],
                                    color='yellow',
                                    label='empirical conc.')
    for a in ax.flatten():
        a.legend()
    plt.savefig(sim_data_file + '.png')
        
        
#Event schedule data file, empirical data
ev_data_path = 'pk_data.csv'
#Simulated data from event schedule 
data_dir = 'pop_data/'
sim_data_file = 'sim_pop.json'
save_path = data_dir + "torsten_" + sim_data_file

edit_sim_data(ev_data_path, data_dir, sim_data_file, save_path)


