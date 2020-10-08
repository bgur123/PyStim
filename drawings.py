# -*- coding: utf-8 -*-
"""

@author: bgur123

"""
from psychopy import visual,core,event, gui, monitors
from helper import *

def initialize_stimulus(epochObj,win,proj_params):
    """ For creating the PsychoPy stimulus objects"""
    if epochObj.stim_type == 'gratings-v1':
        # Moving gratings
        unit = proj_params['unit']
        size = proj_params['size']
        orientation = np.mod(epochObj.direction_deg-90,360) # direction & orientation orthogonal
        grating = visual.GratingStim(
            win=win, name='grating',tex=epochObj.type, 
            size=(size, size), ori = orientation,
            sf=1/epochObj.spatial_wavelength,
            units=proj_params['unit'], contrast=epochObj.michelson_contrast)
        
        grating.autoDraw = True  # Automatically draw every frame
        grating.autoLog = False 
        epochObj.grating = grating
    
    return epochObj


def update_stimulus(epochObj,screen_refresh_rate,win,epoch_clock):
    """ For updating the PsychoPy stimulus objects"""

    
    if epochObj.stim_type == 'gratings-v1':
        # Moving gratings
        while (epoch_clock.getTime() < epochObj.duration_seconds):
            # For moving the grating we will need to advance the phase
            # according to the desired velocity and screen refresh rate
            #| v (degree/s) / refresh rate (update/s) |
            # Since phase advances are related to 1 cycle we should also
            # use spatial wavelength of the grating.
            epochObj.grating.phase += \
                (epochObj.velocity/screen_refresh_rate)/epochObj.spatial_wavelength
            win.flip()

        
        
