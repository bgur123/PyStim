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

#%% We need a monitor and a window to present the stimulus
mon = monitors.Monitor(proj_params["monitorName"])
mon.setSizePix = proj_params["monitorSizePix"]
mon.setWidth = proj_params["monitorWidthcm"]
mon.setDistance = proj_params["observerDistancecm"]
refresh_rate = proj_params["monitorRefreshRate"]

win = visual.Window(
    size=(proj_params['win_width_pix'] ,
    proj_params['win_height_pix']), useFBO = True,allowGUI=False, 
    viewOri = 0.0, fullscr=False, monitor=mon, units='degFlat')

#Detecting dropped frames if any
win.setRecordFrameIntervals(True)
win._refreshTreshold = 1/proj_params["monitorRefreshRate"]+0.004 # warn if frame is late more than 4 ms
logging.console.setLevel(logging.WARNING)
# %% Main loop for presenting the stimulus
epoch_start = True
stim_start = True
stop = False # For nidaq stop condition
last_data_frame = -1
data_frame = 0
# Generate a sequence of 1000 trials to have a sequence for all presentation which is faster.
trial_epoch_sequence = routine.generateTrialSeries(1000)
current_epoch_idx = -1

# We can keep both global and epoch time
print("Stimulus started...")
routine_clock = core.Clock()
routine_max_time = 10000 * 60 # in seconds - to stop the loop if there is no other stop condition
epoch_clock = core.Clock()


while not (len(event.getKeys(['escape'])) \
    or routine_max_time < routine_clock.getTime() or stop):

    if stim_start == True:
            sample_num = 0
            if meta["nidaq"]:
                # Send pulse to the imaging computer and start the counter
                daq_pulse_h, daq_counter_h, daq_data = daqT.configure_daq()
            stim_start = False
    
    if epoch_start:
        # We need to restart the clock 
        epoch_clock.reset()

        # Set the current epoch
        current_epoch_idx += 1
        current_epoch_num = trial_epoch_sequence[current_epoch_idx]
        currentEpoch = routine.epochs[current_epoch_num]
        currentEpoch.startSignal = True
        # Duration is needed to stop the epoch
        epoch_duration = currentEpoch.total_dur_sec
    
        # Create the stimulus with given properties
        stim_frame = 0
        
        epoch_start = False

    elif epoch_clock.getTime() > epoch_duration:
        epoch_start = True
  
    else:
        while (epoch_clock.getTime() < currentEpoch.total_dur_sec):
            # Update stimulus
            pt.run_stimulus(epochObj,proj_params,screen_refresh_rate,win,epoch_clock,
                outputObj)


            # Add a stop condition of nidaq doesn't receive signal for 2 seconds
            if meta["nidaq"]:
                # Check microscope frame status
                (daq_data, data_frame, data_frameT) = \
                    daqT.check_status_daq(daq_counter_h,daq_data,routine_clock)
                if last_data_frame != data_frame:
                    last_data_frameT = data_frameT
                elif last_data_frameT - data_frameT > 2:
                    stop = True
            
            sample_num +=1
            stim_frame +=1
            outputObj.sample_num.append(sample_num)
            outputObj.stim_frame.append(stim_frame)
            outputObj.epoch_time(f'{epoch_clock.getTime():.4f}')
            outputObj.imaging_frame.append(data_frame)
            outputObj.stimulus_epoch.append(current_epoch_num)
            outputObj.time_passed.append(f'{data_frameT:.4f}')


            
if meta["nidaq"]:
    daqT.clearTask(daq_counter_h)
    daqT.clearTask(daq_pulse_h)

save_loc = 'C:\PyStim_outputs'
if not(os.path.exists(save_loc)):
    os.mkdir(save_loc)
outputObj.saveOutput(save_loc)
print(f"Stimulus finished...")
# print(tictocs)