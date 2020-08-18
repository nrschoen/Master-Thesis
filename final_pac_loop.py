# Import needed modules
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pactools import Comodulogram
# from pactools.dar_model import DAR, extract_driver


# Loop through all files in directory

directory = '/Users/Neil/Desktop/Clean_Data/SOLO_OXY'
out_dir = '/Users/Neil/Desktop/Clean_Data/SOLO_OXY/SOLO_OXY_DAR'
out_dir_Z = '/Users/Neil/Desktop/Clean_Data/SOLO_OXY/SOLO_OXY_DAR_Z'

dir_list = os.listdir(directory)
print(dir_list)

# Settings for PAC

fs = 200  # Hz
high_fq = 20.0  # Hz
high_fq_range = np.linspace(14, 30, 50)
low_fq = 2.5  # Hz
low_fq_width = 1.5  # Hz
low_fq_range = np.linspace(1, 4, 50)


#%% Loooop
for filename in dir_list:
    
    if not filename.endswith('.csv'):
        continue
    
    csv_loc = os.path.join(directory, filename)
    out_loc = os.path.join(out_dir, filename)
    out_loc_Z = os.path.join(out_dir_Z, filename)
    print(csv_loc)
    raw = pd.read_csv(csv_loc)
    signal = raw.loc[:, raw.columns.str.contains('Fz')]
    # signal = raw.iloc[:, 4]
    signal = signal.to_numpy()
    signal = signal[signal != 0]
    signal = np.reshape(signal, (len(signal),))
    signal = signal[180000:480000,]



    #  Significance testing with surrogate analysis

    # number of comodulograms computed in the surrogate analysis.
    # rule of thumb is 10 / p_value. Example: 10 / 0.05 = 200.
    n_surrogates = 200
    # how many cores to use
    n_jobs = 2
    # input the ordar and ordriv from DAR model
    method = 'duprelatour'
    
    estimator = Comodulogram(fs=fs, low_fq_range=low_fq_range,
                             low_fq_width=low_fq_width, method=method,
                             n_surrogates=n_surrogates, high_fq_range=high_fq_range,
                             progress_bar=True, n_jobs=n_jobs)
    estimator.fit(signal)
    
    # plot the significance level on top of the comodulogram
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    
    # suffers from multiple testing and assumes normality
    # titles=['With a z-score on each couple of frequency']
    z_score = 4.
    estimator.plot(contour_method='z_score', contour_level=z_score,
                   axs=[axs[0]])
    
    # does not assume normality nor does it suffer from multiple testing
    # titles=['With a p-value on the distribution of maxima']
    p_value = 0.05
    estimator.plot(contour_method='comod_max', contour_level=p_value,
                   axs=[axs[1]])
    
    
    
    # save and show resuting comodulogram
    figure = filename.replace('.csv', '.jpeg')
    out_fig = os.path.join(out_dir, figure)
    plt.tight_layout()
    plt.savefig(out_fig)
    plt.show()
    

    # Save output
    np.savetxt(out_loc, estimator.comod_, delimiter=",")
    np.savetxt(out_loc_Z, estimator.comod_z_score_, delimiter=",")
    
    
    print('done')

# Notify me

os.system('say "... I am ready for you Neil. Much love and kisses"')
