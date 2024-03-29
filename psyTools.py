# -*- coding: utf-8 -*-
"""

@author: bgur123

"""
import time

from psychopy import visual,core,event, gui, monitors
from helper import *
def setup_params(stim_fname):
    """ 
        Configures the parameters of stimulus presentation.
    """

    # Take view settings
    view_settings_file = os.path.join(os.path.dirname(os.path.dirname(stim_fname)),
                                        'view_settings.txt')
    f_handle = open(view_settings_file, 'r')
    view_settings = {}
    for line in f_handle:
        line = re.sub('\n', '', line)
        line = re.sub('\r', '', line)
        line = line.split('\t')
        if line[0]: # If there is input
            view_settings[line[0]] = float(line[1])
            
    
    # Put them in a gui so they can be modified and seen before stimulation
    options_panel = gui.Dlg(title="PyStim options")
    options_panel.addText('Enter desired run options')
    options_panel.addText("Below are default options determined from view_setting.txt")

    options_panel.addField('Monitor settings:', choices=monitors.getAllMonitors()) # idx: 0
    options_panel.addText("Check monitor settings from PsychoPy app and make sure that \n they are correct.")
    # DLP and NIDAQ
    options_panel.addField('use DLP:',initial=True) # idx: 1
    options_panel.addField('use NIDAQ:',initial=True) # idx: 2

    # Screen settings
    options_panel.addText("Settings for the screen") 
    options_panel.addField('Refresh rate:',view_settings['screen_refresh_rate']) # idx: 3
    options_panel.addText("Refresh rate has to match with the screen's") 

    # Stimulus window settings
    options_panel.addText("Settings for the stimulus window") 
    options_panel.addField('Stim win pos X:',view_settings['stim_pos_x']) # idx: 4
    options_panel.addField('Stim win pos Y:',view_settings['stim_pos_y']) # idx: 5
    
    # Presentation settings
    options_panel.addText("Settings for the stimulus presentation") 
    options_panel.addField('Stim units:',choices=["deg","degFlat","degFlatPos"]) # idx: 6
    options_panel.addField('Bit depth:',6) # idx: 7
    options_panel.addField('use warper:',initial=True) # idx: 8
    options_panel.addField('view scale X:',view_settings['view_scale_x']) # idx: 9
    options_panel.addField('view scale Y:',view_settings['view_scale_y']) # idx: 10
    

    user_entry = options_panel.show()  
    mon = monitors.Monitor(user_entry[0])

    stim_size_X = np.arctan((mon.getWidth()/2) / mon.getDistance())
    stim_size_X = abs(np.degrees(stim_size_X)) * 2 # we need the full extend

    stim_size_Y = np.arctan((mon.getWidth()/2) / mon.getDistance())
    stim_size_Y = abs(np.degrees(stim_size_Y)) * 2

    proj_params={
        "unit" : user_entry[6],
        "bit_depth" : float(user_entry[7]),
        "warper" : user_entry[8],
        "sizeX" : stim_size_X,
        "sizeY" : stim_size_Y,
        "posX" : float(user_entry[4]),
        "posY" : float(user_entry[5]),
        "win_width_pix" : mon.getSizePix()[0], # Pixel size from the monitor center
        "win_height_pix" : mon.getSizePix()[1],
        "onDLP" : int(user_entry[1]), # 1 is for DLP 0 is for PC monitors
        "monitorName" : mon.name,
        "monitorSizePix" : mon.getSizePix(),
        "monitorWidthcm" : mon.getWidth(),
        "observerDistancecm" : mon.getDistance(),
        "monitorRefreshRate" : float(user_entry[3]),
        "viewScale" : [float(user_entry[9]),float(user_entry[10])]
    }
    print(proj_params)
    return proj_params, user_entry[2]

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

def run_stimulus(epochObj,cur_time,screen_refresh_rate,win,outputObj):
    """ For drawing and updating the PsychoPy stimulus objects"""
    if epochObj.stim_type in ['gratings-v1', 'centered-gratings-v1','center-grat-with-circle-v1']:
       
        # Moving gratings

        if epochObj.stim_type == 'centered-gratings-v1':
            # Draw the background rectangle with the desired luminance
              epochObj.bg_rect.draw()
        if epochObj.stim_type == 'center-grat-with-circle-v1':
            # Draw the desired circle and the background
            epochObj.bg_rect.draw()
            epochObj.circle.draw()
            
            


        # We will need to advance the phase
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
    elif epochObj.stim_type == 'movingStripe-v1':
        # Moving stripes

        # Draw the background rectangle with the desired luminance
        epochObj.bg_rect.draw()

        # Reset stripe position if epoch is starting new
        if epochObj.startSignal:
            epochObj.stripe.pos = epochObj.initial_pos
            epochObj.startSignal = False
        else:
            movX = (epochObj.velocity/screen_refresh_rate) * np.sin(np.deg2rad(epochObj.direction_deg))
            movY = (epochObj.velocity/screen_refresh_rate) * np.cos(np.deg2rad(epochObj.direction_deg))
            epochObj.stripe.pos += (movX,movY)
        epochObj.stripe.draw()
        win.flip()
        
        outputObj.stim_info1.append(epochObj.stripe.pos)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0) 

    elif epochObj.stim_type == 'colored2stripes-v1':
        # Moving stripes

        # Draw the background rectangle with the desired luminance
        epochObj.bg_rect.draw()

        # Draw stripes
        epochObj.stripe1.draw()
        epochObj.stripe2.draw()
        win.flip()
        
        outputObj.stim_info1.append(epochObj.stripe1.pos)
        outputObj.stim_info2.append(epochObj.stripe2.pos)
        outputObj.stim_info3.append(0) 
        
    elif epochObj.stim_type == 'edges-v1':
        # Moving edges
        
        # Draw the background rectangle with the desired luminance
        epochObj.bg_rect.draw()

        # Reset rectangle if epoch is starting new
        if epochObj.startSignal:
            epochObj.rectangle.width = 0 
            epochObj.startSignal = False

        if cur_time > epochObj.pre_dur_sec:
            # Moved by increasing the width (2x)
            epochObj.rectangle.width += \
                2 * (epochObj.velocity/screen_refresh_rate)
        epochObj.rectangle.draw()
        win.flip()

        outputObj.stim_info1.append(epochObj.rectangle.width)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)
    elif epochObj.stim_type == 'centered-circle-v1':
        # Circles with adjustable size and center that 
        # have an option to change luminance with a certain Weber contrast
        # after a pre-defined duration
        

        if cur_time > epochObj.initial_dur_sec:
            epochObj.circle.fillColor = epochObj.second_lum
            epochObj.circle.lineColor = epochObj.second_lum
        else:
            epochObj.circle.fillColor = epochObj.first_lum
            epochObj.circle.lineColor = epochObj.first_lum
        epochObj.bg_rect.draw()
        epochObj.circle.draw()
        win.flip()

        outputObj.stim_info1.append(0)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)
        
    elif epochObj.stim_type == 'stepFlash-v1':
        # Full field flashes that occur with a certain Weber contrast
        # after a pre-defined duration
        

        if cur_time > epochObj.initial_dur_sec:
            epochObj.rectangle.fillColor = epochObj.second_lum
        else:
            epochObj.rectangle.fillColor = epochObj.first_lum
        epochObj.rectangle.draw()
        win.flip()

        outputObj.stim_info1.append(0)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)

    elif epochObj.stim_type == 'fff-v1':
        # Full field flashes

        # Do nothing just run loop
        epochObj.rectangle.draw()
        win.flip()

        outputObj.stim_info1.append(0)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)
        
    elif epochObj.stim_type == 'test-coordinates':
        # Test the coordinates

        # Do nothing just run loop
        epochObj.test_stim.draw()
        win.flip()

        outputObj.stim_info1.append(0)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)
        
    elif epochObj.stim_type == 'whiteNoiseRectangles-v1':
        # white noise
        
        # Reset texture if epoch is starting new
        if epochObj.startSignal:
            epochObj.currIdx = 0
            epochObj.currFrameRep = 1
            epochObj.startSignal = False

        # Present the current texture (corresponding to the frame of the stimulus)
        epochObj.white_noise_stripes.setImage(epochObj.noise_texture_scaled[epochObj.currIdx,:,:])
        epochObj.white_noise_stripes.draw()
        win.flip()

        outputObj.stim_info1.append(epochObj.currIdx)
        outputObj.stim_info2.append(0)
        outputObj.stim_info3.append(0)

        # Check if it is time to change the stimulus frame
        # This calcuation is based on the update rate
        if ((1.0/epochObj.update_rate))<=(epochObj.currFrameRep * (1/screen_refresh_rate)):
            epochObj.currFrameRep = 1
            epochObj.currIdx += 1
            # epochObj.toc = time.time()
            # print(epochObj.toc- epochObj.tic)
        else:
            # if epochObj.currFrameRep == 1:
            #     epochObj.tic = time.time()
            epochObj.currFrameRep += 1
            
    else:
        raise NameError(f"Stimulus type {epochObj.stim_type} could not be initialized.")


        
def set_window_color(win, epochObj):
    """ 
        Sets the window color before stimulus presentation so that the previous epoch 
        luminance doesn't produce a flash. This phenomena, I think, is a bug in PsychoPy
    """

    if epochObj.stim_type in ['movingStripe-v1', 'edges-v1']:
        win.color= epochObj.win_lum
    
    return None