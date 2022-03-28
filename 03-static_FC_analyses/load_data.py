import numpy as np
import pandas as pd

ex_dual_mot = ['sub-13', 'sub-21', 'sub-23', 'sub-50'] # Higly motion subjects in one of four sessions
ex_rest_mot = ['sub-21', 'sub-46', 'sub-47'] # Higly motion subjects in one of four sessions / missing data(20-44)
ex_rest_mis = ['sub-20', 'sub-44'] # Missing data in resting (2nd session)
ex_firs_mot = ['sub-21'] # Subject with highly motion on first session

# path_fc_pow = 'data/matrices/LB_static_all_conditions_power.npy'
path_fc_pow = 'data/matrices/LB_static_all_conditions_power.npy'
path_fc_sch = 'data/matrices/LB_static_all_conditions_schaefer.npy'
path_ma_pow = 'data/support/Power_modules.txt'
path_ma_sch = 'data/support/Schaefer2018_300Parcels_7Networks_order.csv'
path_bh_tra = 'data/behavioral/LB_training_mean_tidy.csv'
path_bh_mri = 'data/behavioral/LB_fmri_behaviour_mean_tidy.csv'
path_gr_ass = 'data/behavioral/LB_group_assignment.csv'

# Functional connectivity
print('=== Imaging data ===')
mat_pow = np.load(path_fc_pow)
mat_sch = np.load(path_fc_sch)
print('mat_pow is {} having shape {}'.format(type(mat_pow), mat_pow.shape))
print('mat_sch is {} having shape {}'.format(type(mat_sch), mat_sch.shape))
print('\n')

# A priori module assignment
with open(path_ma_pow, 'r') as f:
    ma_pow = f.readlines()
ma_pow = [ma.rstrip('\n') for ma in ma_pow]
ma_sch = pd.read_csv(path_ma_sch, usecols=[6])

# Behavioral
print('=== Behavioral data ===')
beh_tra = pd.read_csv(path_bh_tra)
beh_mri = pd.read_csv(path_bh_mri)
grp_ass = pd.read_csv(path_gr_ass)
print(beh_tra.head(3))
print('\nbeh_tra shape is {}\n'.format(beh_tra.shape))
print(beh_mri.head(3))
print('\nbeh_mri shape is {}\n'.format(beh_mri.shape))
print(grp_ass.head(3))
print('\ngrp_ass shape is {}\n'.format(grp_ass.shape))

# Exclusion
print('=== Excluded subjects ===')
print('Excluded due to motion in first session {}'.format(ex_firs_mot))
print('Excluded due to motion in dual {}'.format(ex_dual_mot))
print('Excluded due to motion in rest {}'.format(ex_rest_mot))
print('Excluded due to missing data in rest {}'.format(ex_rest_mis))
ex_rest = ex_rest_mis + ex_rest_mot
