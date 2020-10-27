import os  # handy system and path functions
import sys  # to get file system encoding
import numpy as np # working with arrays and numbers
import time
from psychopy import visual,core, gui, monitors,event
from datetime import datetime
import PyDAQmx as daq

# Self written packages
from helper import *
import psyTools as pt
import nidaqTools as daqT


#%% Initialize 
# Parameters of projection for stimulus creation
use_nidaq = 1 


(proj_params, use_nidaq) = pt.setup_params()


# An output structure for aligning with acquired imaging frames and storing stimulus
# meta data
now = datetime.now()
dt_string = now.strftime("%Y%m%d_%H%M%S")
meta = {
    "nidaq" : use_nidaq,
    "proj_params" : proj_params,
    "date_time" : dt_string
}
outputObj = OutputInfo(meta=meta)



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
routine_time = 12
epoch_clock = core.Clock()
stop = False

print("Stimulus started...")

# tictocs = [] # To check epoch timings
# Main Loop: dit diplays the stimulus unless:
        # keyboard key is pressed (manual stop)
        # stop condition becomse "True"    
while not (len(event.getKeys(['escape'])) \
    or routine_time < routine_clock.getTime()):

    if epoch_start:
        # tic = time.perf_counter()
        # We need to restart the clock 
        epoch_clock.reset()
        # Set the current epoch
        current_epoch_idx += 1
        current_epoch_num = trial_epoch_sequence[current_epoch_idx]
        currentEpoch = routine.epochs[current_epoch_num]
        
        # Duration is needed to stop the epoch
        epoch_duration = currentEpoch.total_dur_sec
        
        if stim_start == True:
            sample_num = 0
            if meta["nidaq"]:
                # Send pulse to the imaging computer and start the counter
                daq_pulse_h, daq_counter_h, daq_data = daqT.configure_daq()
            stim_start = False
        # Create the stimulus with given properties
        stim_frame = 0
        currentEpoch = pt.initialize_stimulus_v2(currentEpoch,win,
            proj_params,outputObj, stim_frame, sample_num)
        epoch_start = False

        (daq_data, last_data_frame, last_data_frameT) = \
                daqT.check_status_daq(daq_counter_h,daq_data,routine_clock)

        sample_num +=1
        stim_frame +=1
        outputObj.sample_num.append(sample_num)
        outputObj.stim_frame.append(stim_frame)
        outputObj.imaging_frame.append(last_data_frame)
        outputObj.stimulus_epoch.append(current_epoch_num)
        outputObj.time_passed.append(last_data_frameT)

        

    elif epoch_clock.getTime() > epoch_duration:
        epoch_start = True
        # toc = time.perf_counter()
        # tictocs.append(toc-tic)
    else:
        while (epoch_clock.getTime() < currentEpoch.total_dur_sec):
            # Update stimulus
            pt.update_stimulus(currentEpoch,refresh_rate,win,epoch_clock,
                outputObj)

            # Check microscope frame status
            (daq_data, last_data_frame, last_data_frameT) = \
                daqT.check_status_daq(daq_counter_h,daq_data,routine_clock)

            sample_num +=1
            stim_frame +=1
            outputObj.sample_num.append(sample_num)
            outputObj.stim_frame.append(stim_frame)
            outputObj.imaging_frame.append(last_data_frame)
            outputObj.stimulus_epoch.append(current_epoch_num)
            outputObj.time_passed.append(last_data_frameT)


            

daqT.clearTask(daq_counter_h)
daqT.clearTask(daq_pulse_h)

save_loc = 'C:\PyStim_outputs'
if not(os.path.exists(save_loc)):
    os.mkdir(save_loc)
outputObj.saveOutput(save_loc)
print(f"Stimulus finished...")
# print(tictocs)