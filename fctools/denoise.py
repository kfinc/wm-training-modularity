# -*- coding: utf-8 -*-


"""
Created on Wed Jul 19 2018
Last edit: Sat Sep 01 2018
@author: kfinc

"""

import pandas as pd
import numpy as np 
from sklearn import preprocessing
from nistats.design_matrix import make_first_level_design_matrix


def motion_24_friston(dataframe):
 	
    """Simple function that calculates 24 motion parameters from pandas dataframe. 
    
    Parameters
    ----------
    dataframe: pandas dataframe including 6 movement parameters with headers
    
    Returns
    -------
    motion_24_friston:  pandas dataframe including 24 motion parameters
    
    - the first 6 are the motion parameters
    - the next 6 are the temporal difference of motion parameters ('_td' suffix)
    - the next 12 are the square of the motion parameters and the differenced values ('_sqrt' suffix)
      
    """

    motion_24_friston = dataframe 

    for col in dataframe.columns:
        temp_diff = np.roll(dataframe[col], 1, axis = 0)
        temp_diff[0] = 0
        temp_diff = pd.DataFrame(temp_diff)
        motion_24_friston[col + '_td'] = temp_diff
    
    for col in motion_24_friston.columns:
        quad = motion_24_friston[col] ** 2
        motion_24_friston[col + '_quad'] = quad
    
    return motion_24_friston


def scrubbing(fd, thr = 0.5, before = True, after = True):
    
    """Function that calculates motion outliers (frames with motion above threshold_.
    
    Parameters
    ----------
    fd:      pandas dataframe including frame-wise displacement (FD)
    thr:     threshold (default: 0.5)
    before:  marks frames before outlier datapoint (default: True)
    after:   marks frames after outlier datapoint (default: True)
    
    Returns
    -------
    scrubbing:  pandas dataframe including all ourliers datapoints
          
    """
    
    scrubbing = pd.DataFrame()
    fd.loc[0] = 0
    fd = fd.astype(float)
    
    scrub1 = fd > thr
    scrub1 = scrub1.astype(int)
    scrubbing['scrubbing'] = scrub1
    
    if before == True:
        scrub2 = np.roll(scrubbing['scrubbing'], -1, axis = 0)
        scrub2[0] = 0
        scrubbing['scrubbing_bef'] =  scrub2
        
    if after == True:
        scrub3 = np.roll(scrubbing['scrubbing'], 1, axis = 0)
        scrub3[0] = 0
        scrubbing['scrubbing_aft'] =  scrub3
        
    return scrubbing



def standardize(dataframe):
    """
    Normalizes each column and returns values set to unit variance.
    
    Parameters
    ----------
    dataframe: pandas dataframe including columns of interest
    
    Returns
    -------
    dataframe_stand:  pandas dataframe with standarized values
    
    
    """
    
    dataframe_stand = pd.DataFrame()
    val = dataframe.values
    standardize = preprocessing.StandardScaler()
    val_scaled = standardize.fit_transform(val)
    dataframe_stand = pd.DataFrame(val_scaled, columns = dataframe.columns)
    
    return dataframe_stand


def motion_temp_diff(dataframe):
 	
    """Simple function that calculates 12 motion parameters from pandas dataframe. 
    
    Parameters
    ----------
    dataframe: pandas dataframe including 6 movement parameters with headers
    
    Returns
    -------
    motion_temp_diff:  pandas dataframe including 24 motion parameters
    
    - the first 6 are the motion parameters
    - the next 6 are the temporal difference of motion parameters ('_td' suffix)
      
    """

    motion_temp_diff = dataframe 

    for col in dataframe.columns:
        temp_diff = np.roll(dataframe[col], 1, axis = 0)
        temp_diff[0] = 0
        temp_diff = pd.DataFrame(temp_diff)
        motion_temp_diff[col + '_td'] = temp_diff

    return motion_temp_diff


def temp_deriv(dataframe, quadratic = False):
    """Simple function that calculates temporal derivatives for each column of pandas dataframe. 
    
    Parameters
    ----------
    dataframe: pandas dataframe with variable to calculate temporal derivarives
    
    Returns
    -------
    temp_deriv:  pandas dataframe including original columns and their temporal derivatives ('_td') and (optional) 
    their quadratic terms
          
    """
    
    temp_deriv = dataframe.copy()
    
    for col in dataframe.columns:
        #--- backward difference algorithm
        temp = np.diff(dataframe[col], 1, axis = 0)
        temp = np.insert(temp, 0, 0)
        temp = pd.DataFrame(temp )
        temp_deriv[col + '_td'] = temp
        
    if quadratic == True:
        for col in temp_deriv.columns:
            quad = temp_deriv[col] ** 2
            temp_deriv[col + '_quad'] = quad

    return temp_deriv

def scrubbing(fd, thr = 0.5, before = True, after = True):
    
    """Function that calculates motion outliers (frames with motion above threshold_.
    
    Parameters
    ----------
    fd:      pandas dataframe including frame-wise displacement (FD)
    thr:     threshold (default: 0.5)
    before:  marks frames before outlier datapoint (default: True)
    after:   marks frames after outlier datapoint (default: True)
    
    Returns
    -------
    scrubbing:  pandas dataframe including all outliers datapoints
          
    """
    
    scrubbing = pd.DataFrame()
    fd.loc[0] = 0
    fd = fd.astype(float)
    
    scrub1 = fd > thr
    scrub1 = scrub1.astype(int)
    scrubbing['scrubbing'] = scrub1
    
    if before == True:
        scrub2 = np.roll(scrubbing['scrubbing'], -1, axis = 0)
        scrub2[0] = 0
        scrubbing['scrubbing_bef'] =  scrub2
        
    if after == True:
        scrub3 = np.roll(scrubbing['scrubbing'], 1, axis = 0)
        scrub3[0] = 0
        scrubbing['scrubbing_aft'] =  scrub3
        
    return scrubbing

def outliers_fd_dvars(dataframe, fd=0.5, dvars=3):
    """Function that calculates motion outliers (frames with frame-wise displacement (FD)
    and DVARS above predefined threshold).
    
    Parameters
    ----------
    dataframe: pandas dataframe including columns with DVARS and FD
    fd:        threshold for FD (default: 0.5)
    dvars:     threshold for DVARS (+/-SD, default: 3)   
    
    Returns
    -------
    outliers:  pandas dataframe including all outliers datapoints
          
    """
    
    df = dataframe.copy()
    df.fillna(value=0, inplace=True)

    dvars_out = np.absolute(df[df.columns[0]].astype(float)) > dvars  
    fd_out = df[df.columns[1]].astype(float) > fd

    outliers = (dvars_out == True) | (fd_out == True)
    outliers = pd.DataFrame(outliers.astype('int'))
    outliers.columns = ['scrubbing']

    return outliers

def get_condition_column(events, tr = 2, n_scans = 340):
    """converts events file to pd dataframe with column representing each condition"""
    frame_times = np.arange(n_scans) * tr
    box = make_first_level_design_matrix(frame_times, events, hrf_model = None)
    box = box.reset_index()

    x = box.iloc[:,1:4] > 0.8
    y = x.astype('int')
    col = pd.DataFrame(y.idxmax(axis=1), columns = ['condition'])
    return col
    
    