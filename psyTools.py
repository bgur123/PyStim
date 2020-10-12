# -*- coding: utf-8 -*-
"""

@author: bgur123

"""
import time

from psychopy import visual,core,event, gui, monitors
from helper import *
def setup_params():
    options_panel = gui.Dlg(title="PyStim options")
    options_panel.addText('Enter desired run options')
    options_panel.addField('Monitor:', choices=monitors.getAllMonitors())
    options_panel.addField('use NIDAQ:',initial=True)
    user_entry = options_panel.show()  

    mon = monitors.Monitor(user_entry[0])
    
    proj_params={
        "unit" : 'degFlat',
        "size" : 20,
        "win_size_pix" : 500.0,
        "monitorName" : "testMonitor",
        "monitorSizePix" : (1920,1080),
        "monitorWidthcm" : mon.getWidth(),
        "observerDistancecm" : mon.getDistance() ,
        "monitorRefreshRate" : 60.0
    }
    return proj_params, user_entry[1]

def initialize_stimulus(epochObj,win,proj_params):
    """ For creating the PsychoPy stimulus objects"""
    if epochObj.stim_type == 'gratings-v1':
        # Moving gratings
        dimension = 1024 # It needs to be square power-of-two (e.g. 256 x 256) for PsychoPy
        unit = proj_params['unit']
        size = proj_params['size']
        # Calculation of BG and FG depending on the michelson contrast definition
        # fg-bg = c * 2 * l
        # fg+bg = 2 * l                    
        contrast = epochObj.michelson_contrast
        luminance = epochObj.mean_luminance
        fg = contrast * 2 * luminance
        bg = 2*luminance - fg
        f = 1# generate a single cycle
        # Generate 1D wave and modulate the luminance and contrast
        if epochObj.type == "sin":
            
            x = np.arange(dimension)
            # Wave needs to be scaled to 0-1 so we can modulate it easier later
            sine_signal = (np.sin(2 * np.pi * f * x / dimension)/2 +0.5)
            # We need to scale and shift the wave to match the fg bg values
            oneD_wave = sine_signal * 2*(fg - bg) # Scaling the signal
            oneD_wave = oneD_wave -1 + bg 
           
        else:
            raise ValueError("Grating type {s} doesn't exist".format(\
                s=epochObj.type))
        grating_texture = np.tile(oneD_wave, [dimension,1])
        # We need to extend the wave to 2D
        
        
        orientation = np.mod(epochObj.direction_deg-90,360) # direction & orientation orthogonal
        grating = visual.GratingStim(
            win=win, name='grating',tex=grating_texture, 
            size=(size, size), ori = orientation,
            sf=1/epochObj.spatial_wavelength,
            units=proj_params['unit'])
        
        grating.autoDraw = True  # Automatically draw every frame
        grating.autoLog = False 
        epochObj.grating = grating

    
    return epochObj


def update_stimulus(epochObj,screen_refresh_rate,win,epoch_clock,outputObj):
    """ For updating the PsychoPy stimulus objects"""

    
    if epochObj.stim_type == 'gratings-v1':
        # Moving gratings
        
        # For moving the grating we will need to advance the phase
        # according to the desired velocity and screen refresh rate
        #| v (degree/s) / refresh rate (update/s) |
        # Since phase advances are related to 1 cycle we should also
        # use spatial wavelength of the grating.
        epochObj.grating.phase += \
            (epochObj.velocity/screen_refresh_rate)/epochObj.spatial_wavelength
        win.flip()

        outputObj.stim_info1.append(epochObj.grating.phase)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)

            
            


        
        
