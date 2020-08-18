#%% # Import needed modules
import mne
import numpy as np
import os
# import pandas as pd
from mne.preprocessing import ICA
from autoreject import (AutoReject, get_rejection_threshold)
#%% Read raw file
file_loc = '/Users/Neil/Desktop/Data/ppn2_s1_DUO_OXY.cnt'
csv_out_loc = '/Users/Neil/Desktop/Clean_Data/ppn01_s1_DUO_OXY.csv'


# odd
raw = mne.io.read_raw_cnt(file_loc,
                          preload=True, eog=('VEOG+', 'VEOG-',
                                             'HEOG+', 'HEOG-'))

ch_drop = raw.ch_names[32:]
raw.drop_channels('RM')

# even
'''raw = mne.io.read_raw_cnt(file_loc,
                          preload=True, eog=('2VEOG+', '2VEOG-',
                                             '2HEOG+', '2HEOG-'))

ch_drop = raw.ch_names[:32]


raw.drop_channels('2RM')


raw.drop_channels(ch_drop)
ch_dict = {'2VEOG+': 'VEOG+',
           '2VEOG-': 'VEOG-',
           '2F7': 'F7',
           '2F3': 'F3',
           '2Fz': 'Fz',
           '2F4': 'F4',
           '2F8': 'F8',
           '2FC5': 'FC5',
           '2FC1': 'FC1',
           '2FC2': 'FC2',
           '2FC6': 'FC6',
           '2T7': 'T7',
           '2C3': 'C3',
           '2Cz': 'Cz',
           '2C4': 'C4',
           '2T8': 'T8',
           '2FCz': 'FCz',
           '2CP5': 'CP5',
           '2CP1': 'CP1',
           '2CP2': 'CP2',
           '2CP6': 'CP6',
           '2P7': 'P7',
           '2P3': 'P3',
           '2Pz': 'Pz',
           '2P4': 'P4',
           '2P8': 'P8',
           '2O1': 'O1',
           '2Oz': 'Oz',
           '2O2': 'O2',
           '2HEOG+': 'HEOG+',
           '2HEOG-': 'HEOG-'}

raw.rename_channels(ch_dict)
'''


# set channel locations
raw = raw.set_montage('standard_1020')

print(raw.info)
raw.plot()



#%% Referencing to average, filtering, preparing for autoreject and ICA

# use the average of all channels as reference
raw_avg_ref = raw.copy().set_eeg_reference(ref_channels='average')
raw_avg_ref.plot()

# filtering to remove slow drifts
raw = raw_avg_ref
raw.load_data().filter(l_freq=1, h_freq=40, phase='zero')

raw.plot()


#%% Prepare to Autoreject Artifacts
# n_interpolates are the p values that we would like autoreject
n_interpolates = np.array([1, 4, 32])

# consensus_percs are the k values that autoreject will try
consensus_percs = np.linspace(0, 1.0, 11)

# pick channels to autoreject
picks = mne.pick_types(raw.info, meg=False, eeg=True, stim=False, eog=False,
                       include=[], exclude=[])

raw.info['projs'] = list()  # remove proj, don't proj while interpolating

# reject bad channels
raw.info['bads'] = ['T7', 'T8', 'F8', 'FC5', 'FC6', 'Cz']
#raw.plot()

# Create fixed length epochs
# select fixed epoch length in seconds
tstep = 1.0

# create fixed length events and epochs
events = mne.make_fixed_length_events(raw, duration=tstep)
epochs = mne.Epochs(raw, events, tmin=0.0, tmax=tstep, baseline=(0, 0), reject=None,
                    verbose=False, detrend=0, preload=True)

# plot epochs
epochs.plot()


#%% Fit ICA on epochs under rejection threshhold trials

# find peak to peak rejection threshhold
reject = get_rejection_threshold(epochs)

# run ICA on all epochs lower than the threshhold
ica = ICA(n_components=15, random_state=42)
ica.fit(epochs, reject=reject, tstep=tstep)

# plot the components
raw.load_data()
ica.plot_sources(raw)

ica.plot_components()

# notify when done
os.system('say "... I am ready for you Neil."')

#%% Select components to exclude and zero them out

rm_select = [0, 3, 8]

# blinks
ica.plot_overlay(raw, exclude=rm_select, picks='eeg')

# indices chosen based on various plots above
ica.exclude = rm_select

# apply ICA
# ica.apply() changes the Raw object in-place:
ica.apply(raw)

# plot raw data excluding ocular artifacts for quality check
raw.plot()

#%% Fit autoreject

events = mne.make_fixed_length_events(raw, duration=tstep)
epochs = mne.Epochs(raw, events, tmin=0.0, tmax=tstep, baseline=(0, 0), reject=None,
                    verbose=False, detrend=0, preload=True)

ar = AutoReject(n_interpolates, consensus_percs, picks=picks,
                thresh_method='random_search', random_state=42)

# Note that fitting and transforming can be done on different compatible
# portions of data if needed.
ar.fit(epochs)

# epochs_ar, reject_log = ar.fit_transform(epochs, return_log=True)

epochs_clean = ar.transform(epochs)
reject_log = ar.get_reject_log(epochs)
evoked_clean = epochs_clean.average()
evoked = epochs.average()


# visualize rejected epochs
scalings = dict(eeg=40e-6)
reject_log.plot_epochs(epochs, scalings=scalings)

# epochs after cleaning
epochs_clean.plot(scalings=scalings)

# notify when done
os.system('say "... I am ready for you Neil."')

# Visualize repaired data
ylim = dict(eeg=(-2, 2))
epochs.average().plot(ylim=ylim, spatial_colors=True)
epochs_clean.average().plot(ylim=ylim, spatial_colors=True)

 #%% Export clean data to dataframe

index, scale_time = ['time'], 1e3

df_clean = epochs_clean.to_data_frame(picks=picks,
                          index=index)
# export as csv to previously specified location
df_clean.to_csv(csv_out_loc, header=True)



