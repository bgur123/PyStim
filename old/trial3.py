import os  # handy system and path functions
import sys  # to get file system encoding
import numpy as np # working with arrays and numbers
import time

from psychopy import visual,core, gui, monitors,event
from datetime import datetime
# Self written packages
from helper import *
import drawings as dr


#%% Initialize 
# Parameters of projection for stimulus creation
use_nidaq = 0 

proj_params={
    "unit" : 'degFlat',
    "size" : 5,
    "win_size_pix" : 500.0,
    "monitorName" : "iMac-bgur",
    "monitorSizePix" : (4023,2304),
    "monitorWidthcm" : 47.57,
    "observerDistancecm" : 57.0 ,
    "monitorRefreshRate" : 60.0
}

# An output structure for aligning with acquired imaging frames and storing stimulus
# meta data
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
meta = {
    "nidaq" : use_nidaq,
    "proj_params" : proj_params,
    "date_time" : dt_string
}
output_dict = {
    "meta" : meta
    "sample_num" : [] # Just and index number for rows
    "imaging_frame" : [] # Current frame of imaging series
    "stimulus_epoch" : [] # Current epoch being presented
    "time_passed" : [] # Time passed since start of logging
    "stim_info1" : [] # Anything can be placed here depending on the stimulus
    "stim_info2" : []
    "stim_info3" : []
}

#%% Reading stimulus information & generating epochs
# Ask user where the stim input file is located
stim_fname = gui.fileOpenDlg(os.getcwd())[0]
# Pre-organizing the routine and epochs
routine = PyStimRoutine(stim_fname)
print('Stimulus routine with {eN} epochs is generated...'.format(\
    eN = routine.total_epoch_n))
#%% We need a monitor to present the stimulus
mon = monitors.Monitor(proj_params["monitorName"])
mon.setSizePix = proj_params["monitorSizePix"]
mon.setWidth = proj_params["monitorWidthcm"]
mon.setDistance = proj_params["observerDistancecm"]
refresh_rate = proj_params["monitorRefreshRate"]

win = visual.Window(
    size=(proj_params['win_size_pix'], 
    proj_params['win_size_pix']), fullscr=False, 
    monitor=mon, units='degFlat')


# %% Main loop for presenting the stimulus
epoch_start = True
stim_start = True
# Generate a sequence of 1000 trials to have a sequence for all presentation which is faster.
trial_epoch_sequence = routine.generateTrialSeries(1000)
current_epoch_idx = -1

# We can keep both global and epoch time
routine_clock = core.Clock()
routine_time = 42
epoch_clock = core.Clock()
stop = False

print("Stimulus started...")

# tictocs = [] # To check epoch timings
# Main Loop: dit diplays the stimulus unless:
        # keyboard key is pressed (manual stop)
        # stop condition becomse "True"    
while not (len(event.getKeys(['escape'])) \
    or routine_time < routine_clock.getTime()):

    if stim_start = True:
        if meta{"nidaq"}:
            # Send pulse to the imaging computer
        stim_start = False

    if epoch_start:
        # tic = time.perf_counter()
        # We need to restart the clock 
        epoch_clock.reset()
        # Set the current epoch
        current_epoch_idx += 1
        current_epoch_num = trial_epoch_sequence[current_epoch_idx]
        currentEpoch = routine.epochs[current_epoch_num]
        print(current_epoch_num)
        
        # Duration is needed to stop the epoch
        epoch_duration = currentEpoch.total_dur_sec
        
        # Create the stimulus with given properties
        currentEpoch = dr.initialize_stimulus(currentEpoch,win,
            proj_params)
        epoch_start = False

    elif epoch_clock.getTime() > epoch_duration:
        epoch_start = True
        # toc = time.perf_counter()
        # tictocs.append(toc-tic)
    else:
        # Update stimulus
        dr.update_stimulus(currentEpoch,refresh_rate,
            win,epoch_clock)
        


print(f"Stimulus finished...")
print(tictocs)