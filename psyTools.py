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
    options_panel.addField('Refresh rate:',100)
    options_panel.addText("Refresh rate has to match with the monitor's")
    options_panel.addField('Stim win width pix:',500)
    options_panel.addField('Stim win height pix:',500)
    options_panel.addField('Stim size deg:',20)
    user_entry = options_panel.show()  

    mon = monitors.Monitor(user_entry[0])
    
    proj_params={
        "unit" : 'degFlat',
        "size" : float(user_entry[6]),
        "win_width_pix" : float(user_entry[4]),
        "win_height_pix" : float(user_entry[5]),
        "monitorName" : mon.name,
        "monitorSizePix" : mon.getSizePix(),
        "monitorWidthcm" : mon.getWidth(),
        "observerDistancecm" : mon.getDistance() ,
        "monitorRefreshRate" : float(user_entry[2])
    }
    return proj_params, user_entry[1]

def initialize_stimulus(epochObj,win,proj_params,outputObj):
    """ For creating the PsychoPy stimulus objects"""
    if epochObj.stim_type == 'gratings-v1':
        # Moving gratings
        dimension = 1024 # It needs to be square power-of-two (e.g. 256 x 256) for PsychoPy

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
            size=(proj_params['size'], proj_params['size']), ori = orientation,
            sf=1/epochObj.spatial_wavelength,
            units=proj_params['unit'])
        
        grating.autoDraw = True  # Automatically draw every frame
        grating.autoLog = False 
        epochObj.grating = grating
        
        outputObj.stim_info1.append(epochObj.grating.phase)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)

    
    return epochObj
def initialize_stimulus_v2(epochObj,win,proj_params):

    """ For creating the PsychoPy stimulus objects"""
    if epochObj.stim_type == 'gratings-v1':
        orientation = np.mod(epochObj.direction_deg-90,360) # direction & orientation orthogonal
        grating = visual.GratingStim(
            win=win, name='grating',tex=epochObj.grating_texture, 
            size=(proj_params['size'], proj_params['size']), ori = orientation,
            sf=1/epochObj.spatial_wavelength,
            units=proj_params['unit'])
        
        grating.autoDraw = True  # Automatically draw every frame
        grating.autoLog = False 
        epochObj.grating = grating
    else:
        raise ValueError("Grating type {s} doesn't exist".format(\
                    s=self.type))
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

def run_stimulus(epochObj,proj_params,screen_refresh_rate,win,epoch_clock,
                    outputObj):
    """ For drawing and updating the PsychoPy stimulus objects"""

    
    if epochObj.stim_type == 'gratings-v1':
        # Moving gratings
        if epochObj.startSignal:
            orientation = np.mod(epochObj.direction_deg-90,360) # direction & orientation orthogonal
            grating = visual.GratingStim(
                win=win, name='grating',tex=epochObj.grating_texture, 
                size=(proj_params['size'], proj_params['size']), ori = orientation,
                sf=1/epochObj.spatial_wavelength,
                units=proj_params['unit'])
            grating.autoLog = False 
            epochObj.grating = grating
            epochObj.startSignal = False
        
        
        # For moving the grating we will need to advance the phase
        # according to the desired velocity and screen refresh rate
        #| v (degree/s) / refresh rate (update/s) |
        # Since phase advances are related to 1 cycle we should also
        # use spatial wavelength of the grating.
        epochObj.grating.phase += \
            (epochObj.velocity/screen_refresh_rate)/epochObj.spatial_wavelength
        epochObj.grating.draw()
        win.flip()

        outputObj.stim_info1.append(epochObj.grating.phase)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)


            
            


        
        
