# -*- coding: utf-8 -*-
"""

@author: bgur123

"""
import numpy as np
import re
import os
import pickle # To save .pickle files for python post processing
import scipy.io as sio # To save .mat files
import h5py # To import .mat files
from psychopy import visual



class PyStimRoutine(object):
    """ 
        For constructing a particular set of stimulus consisting of pre-defined epochs.
    """
    def __init__(self,stim_fname):

        
        self.stim_fname = stim_fname
        # Find the name of the stimulus routine
        self.stim_name = os.path.basename(self.stim_fname).split('.')[0]
        # We need to parse the .txt file to understand the stimulus routine structure.
        self.stim_input_raw = self.readStimInput(self.stim_fname)
        self.processStimInput(self.stim_input_raw)
        
        
    def readStimInput(self, stim_fname):
        """ 
            Parses the input .txt files 
        """
        f_handle = open(stim_fname, 'r')
        stim_input_raw = {}
        for line in f_handle:
            line = re.sub('\n', '', line)
            line = re.sub('\r', '', line)
            line = line.split('\t')
            stim_input_raw[line[0]] = line[1:]
        return stim_input_raw

    def processStimInput(self, stim_input_raw):
        """
        """
    
        self.randomization_condition = \
            int(stim_input_raw['randomization_condition'][0])

        self.total_epoch_n = \
            int(stim_input_raw['total_epoch_num'][0])

        self.epoch_nums = \
            list(map(int, stim_input_raw['epoch_nums']))

        
        # Total epoch number should match with given epochs since 
        # there might be a reading mistake 
        if self.total_epoch_n != len(self.epoch_nums):
            raise ValueError("Total epoch number doesn't match with the number of given epochs") 
        
        # Generate the individual epochs
        self.stim_types = \
            list(map(str, stim_input_raw['stim_type']))
        self.epochs = {}
        for idxEpoch, iEpoch in enumerate(self.epoch_nums):
            parameters = {}
            for key in stim_input_raw.keys():
                if "param" in key:
                    parameters[key] = stim_input_raw[key][idxEpoch]
            self.epochs[iEpoch]=(PyStimEpoch(self.stim_types[idxEpoch], 
                                            parameters,self.stim_fname))

    def generateEpochSeries(self):
        """ Generates epoch series for a single trial """
        if self.randomization_condition == 0:
            # No randomization
            series = self.epoch_nums
        
        elif self.randomization_condition == 1:
            # Randomize with background epoch always present
            series = np.ones((self.total_epoch_n - 1)*2)
            epochs = self.epoch_nums[1:]
            np.random.shuffle(epochs)
            series[1::2] = epochs
       
        elif self.randomization_condition == 2:
            # Randomize all epochs
            series = self.epoch_nums
            np.random.shuffle(series)

        elif self.randomization_condition == 3:
            # No randomization with background epoch always present
            series = np.ones((self.total_epoch_n - 1)*2)
            epochs = self.epoch_nums[1:]
            series[1::2] = epochs
        
        else:
            raise ValueError('Randomization condition not defined.')
        
        self.current_epoch_series = series
        return series

    def generateTrialSeries(self,n_trials):
        """ Generates epoch series for a certain number of trials """
        single_trial_len = len(self.generateEpochSeries())

        trial_series = np.zeros(single_trial_len*n_trials)
        for iTrial in range(n_trials):
            trial_series[iTrial*single_trial_len:iTrial*single_trial_len+single_trial_len]=\
                self.generateEpochSeries()
        return trial_series

    def initializeEpochs(self,win,proj_params):
        """ Initializes Epochs that are within this routine """
        for epochId in self.epochs:
            self.epochs[epochId].initializeStimulus(win,proj_params)
        print("Initialized all epochs...")


        

class PyStimEpoch(PyStimRoutine):
    """ 
        For constructing a an epoch with pre-defined parameters.
    """
    def __init__(self, stim_type, params, stim_fname):

        # Stimulus information is stored in an additional .txt file 
        # This file contains parameters that are going to be used for the epoch.
        stim_info_fname = os.path.join(os.path.dirname(os.path.dirname(stim_fname)),
                                        'stim_info', "stim-info_{s}.txt".format(s=stim_type))
        self.stim_type = stim_type
        self.stim_fname = stim_fname
        self.epoch_info = self.readStimInput(stim_info_fname)
        self.processed_params = self.processEpochInfo(self.epoch_info,params)

    def processEpochInfo(self,epoch_info, params):
        """ Processes the epoch parameters"""
        processed_parameters = {}
        for info_param in epoch_info.keys():
            if not('param' in info_param):
                continue
            for param in params.keys():
                if not(info_param == param):
                    continue
                if epoch_info[info_param][1] == 'int':
                    processed_parameters[epoch_info[info_param][0]] =\
                        int(params[info_param])
                    self.__dict__[epoch_info[info_param][0]] =\
                        int(params[info_param])
                elif epoch_info[info_param][1] == 'float':
                    processed_parameters[epoch_info[info_param][0]] =\
                        float(params[info_param])
                    self.__dict__[epoch_info[info_param][0]] =\
                        float(params[info_param])
                elif epoch_info[info_param][1] == 'bool':
                    processed_parameters[epoch_info[info_param][0]] =\
                        bool(params[info_param])
                    self.__dict__[epoch_info[info_param][0]] =\
                        bool(params[info_param])
                elif epoch_info[info_param][1] == 'str':
                    processed_parameters[epoch_info[info_param][0]] =\
                        str(params[info_param])
                    self.__dict__[epoch_info[info_param][0]] =\
                        str(params[info_param])
                else:
                    errormsg = '{s} {p} type not defined correctly.'.format(s=self.stim_type,p=info_param)
                    raise ValueError(errormsg)
        return processed_parameters
    
    def initializeStimulus(self,win,proj_params):
        """ Initializing the stimulus depending on the type.
        """
        # Values must be scaled for matching bit depths which is by default
        # 8 bits in psychopy functions
        bit_depth_scaler = ((2**proj_params["bit_depth"])-1)/float(2**8-1) 
        if self.stim_type == 'gratings-v1':
            grating_texture = self.generateGratingTexture(bit_depth_scaler,dimension = 1024)
            orientation = np.mod(self.direction_deg-90,360) # direction & orientation orthogonal
            
            # Size of the other dimension is bigger to span all the screen
            grating = visual.GratingStim(
                win=win, name='grating',tex=grating_texture, 
                size=(proj_params['sizeX'], proj_params['sizeY']), ori = orientation,
                sf=1/self.spatial_wavelength,
                units=proj_params['unit'])

            grating.autoLog = False 
            self.grating = grating

        elif self.stim_type == 'centered-gratings-v1':
            # Gratings that are centered in a defined position with a certain size
            
            # Background luminance
            self.bg_lum_scaled = ((self.background_lum*2)* bit_depth_scaler -1) 
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            bg_rect = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=self.bg_lum_scaled, fillColorSpace='rgb')
            bg_rect.autoLog = False
            self.bg_rect = bg_rect


            grating_texture = self.generateGratingTexture(bit_depth_scaler,dimension = 1024)
            orientation = np.mod(self.direction_deg-90,360) # direction & orientation orthogonal

            grating = visual.GratingStim(
                win=win, name='grating',tex=grating_texture,mask='circle',
                size=(self.diameter_deg, self.diameter_deg), ori = orientation,
                sf=1/self.spatial_wavelength,pos=(self.x_center,self.y_center),
                units=proj_params['unit'])

            grating.autoLog = False 
            self.grating = grating

        elif self.stim_type == 'centered-circle-v1':
            # Circles that are centered in a defined position with a certain size
            
            # Background luminance
            self.bg_lum_scaled = ((self.background_lum*2)* bit_depth_scaler -1) 
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            bg_rect = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=self.bg_lum_scaled, fillColorSpace='rgb')
            bg_rect.autoLog = False
            self.bg_rect = bg_rect


            # Second luminance is calculated via Weber contrast formula
            second_lum = (self.weber_c * self.initial_lum) + self.initial_lum
             
            self.first_lum = ((self.initial_lum*2)* bit_depth_scaler-1)
            self.second_lum = ((second_lum*2)* bit_depth_scaler-1)
            circle = visual.Circle(
                        win=win, name='circle', units=proj_params['unit'],radius = self.diameter_deg/2,
                        fillColor=self.first_lum, fillColorSpace='rgb',pos=(self.x_center,self.y_center),
                        lineColor = self.first_lum)
            circle.autoLog = False
            self.circle = circle

        elif self.stim_type == 'movingStripe-v1':
            # Moving stripe that covers the whole screen

            
            # Luminances need to be scaled
            stripe_lum = ((self.stripe_lum*2)* bit_depth_scaler -1) 
            self.win_lum = ((self.background_lum*2)* bit_depth_scaler -1) 
            
            # Background is defined by a rectangle that spans whole screen (and also the cut parts)
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            bg_rect = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=self.win_lum, fillColorSpace='rgb')
            bg_rect.autoLog = False
            self.bg_rect = bg_rect

            # Direction of movement
            orientation = np.mod(self.direction_deg-90,360)
        
            # The height should be sufficient to cover all of the extends
            diag = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            height = diag * 2

            self.initial_pos, distance_to_travel = self.defineStimPos(proj_params)
            distance_to_travel += self.width # so that the stripe disappears from the screen
            retractX = (self.width/2) * np.sin(np.deg2rad(self.direction_deg))
            retractY = (self.width/2) * np.cos(np.deg2rad(self.direction_deg))
            self.initial_pos = (self.initial_pos[0]-retractX,\
                self.initial_pos[1]-retractY) # Retracting
                
            # Total time required is calculated automatically
            self.total_dur_sec = distance_to_travel/self.velocity

            stripe = visual.Rect(
                        win=win, name='stripe', units=proj_params['unit'],
                        width=self.width, height=height, pos = self.initial_pos,
                        ori=orientation,lineWidth=0, 
                        fillColor=stripe_lum, fillColorSpace='rgb')
            stripe.autoLog = False
            self.stripe = stripe

        elif self.stim_type == 'edges-v1':
            # Moving edges with a luminance before the edge comes

            # Pre edge window luminance (a.k.a. background)
            # Scaling so it fits to -1 1 range
            # Scaling so it is presented accurately with desired bit depth
            self.win_lum = ((self.pre_lum*2)* bit_depth_scaler -1) 
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            bg_rect = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=self.win_lum, fillColorSpace='rgb')
            bg_rect.autoLog = False
            self.bg_rect = bg_rect

            # Edge luminance
            edge_luminance = ((self.edge_lum*2)* bit_depth_scaler-1)   
            
            # Direction of movement
            orientation = np.mod(self.direction_deg-90,360)

            # The edge should fill the screen - the maximum possible size
            screen_diag = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            height = screen_diag * 2
            stim_pos, _ = self.defineStimPos(proj_params)
            rectangle = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        width=0, height=height, pos = stim_pos,
                        ori=orientation,lineWidth=0, 
                        fillColor=edge_luminance, fillColorSpace='rgb')
            rectangle.autoLog = False
            self.rectangle = rectangle

        elif self.stim_type == 'fff-v1':
            # Full field flashes
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            lum = ((self.lum*2)* bit_depth_scaler-1) 
            self.win_lum = lum
            rectangle = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=lum, fillColorSpace='rgb')
            rectangle.autoLog = False
            self.rectangle = rectangle

        elif self.stim_type == 'stepFlash-v1':
            # Full field flashes that occur with a certain Weber contrast
            # after a pre-defined duration
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
             
            # Second luminance is calculated via Weber contrast formula
            second_lum = (self.weber_c * self.initial_lum) + self.initial_lum
             
            self.first_lum = ((self.initial_lum*2)* bit_depth_scaler-1)
            self.second_lum = ((second_lum*2)* bit_depth_scaler-1)
            rectangle = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=self.first_lum, fillColorSpace='rgb')
            rectangle.autoLog = False
            self.rectangle = rectangle

        elif self.stim_type == 'whiteNoiseRectangles-v1':

            frame_n = int(self.update_rate * self.total_dur_sec) 
            
            # Update rate needs to be a multiple of screen refresh rate otherwise it is not possible to present it
            if not((proj_params['monitorRefreshRate']/self.update_rate).is_integer()):
                raise ValueError(f"Update rate: {self.update_rate} not possible with screen refresh rate: {proj_params['monitorRefreshRate']}.")
            

            # Width of both dimensions have to be a divisor of screen dimensions, this will be based on the
            # rounded value of the screen angular dimensions since otherwise it is not easily possible 
            # to find a divisor
            if self.x_width != 0:
                isNotDivX = bool(np.mod(np.floor(proj_params['sizeX']),self.x_width))
                if isNotDivX:
                    raise ValueError(f"Stripe X width: {self.x_width} not divisible to screen X width {proj_params['sizeX']}.")
                x_dim = int(proj_params['sizeX']/self.x_width)
            else:
                x_dim = 1
            if self.y_width != 0: 
                isNotDivY = bool(np.mod(np.floor(proj_params['sizeY']),self.y_width))
                if isNotDivY:
                    raise ValueError(f"Stripe Y width: {self.y_width} not divisible to screen Y width {proj_params['sizeY']}.")
                y_dim = int(proj_params['sizeY']/self.y_width)
            else:
                y_dim = 1

            

            np.random.seed(self.seed)

            noise_texture= np.random.choice(np.linspace(0,1,self.uniq_vals), 
                size=(frame_n,y_dim,x_dim))
            # noise_texture = np.repeat(noise_texture,stripe_n,axis=1)
            self.processed_params['noise_texture'] = noise_texture
            curr_idx = 0
            # Scale according to bit depth
            noise_texture_scaled = ((noise_texture*2)* bit_depth_scaler-1) 
            self.noise_texture_scaled = noise_texture_scaled

            white_noise_stripes = visual.ImageStim(
                win=win, name='white_noise_stripes',image=noise_texture_scaled[curr_idx,:,:], 
                size=(proj_params['sizeX'], proj_params['sizeY']),
                units=proj_params['unit'])

            white_noise_stripes.autoLog = False 

            self.white_noise_stripes = white_noise_stripes
            self.frame_idx = curr_idx

        elif self.stim_type == 'test-coordinates':

            test_texture = np.zeros(shape=(20,20))
            test_texture[10,10] = 1
            test_texture[1,1] = 1
            test_texture[1,19] = 0.6
            test_texture[19,19] = 0.2
    
            # Scale according to bit depth
            test_texture_scaled = ((test_texture*2)* bit_depth_scaler-1) 
            self.test_texture_scaled = test_texture_scaled

            test_stim = visual.ImageStim(
                win=win, name='test_stim',image=test_texture_scaled, 
                size=(proj_params['sizeX'], proj_params['sizeY']),
                units=proj_params['unit'])

            test_stim.autoLog = False 

            self.test_stim = test_stim
        else:
            raise NameError(f"Stimulus type {self.stim_type} could not be initialized.")
    
    def generateGratingTexture(self, bit_depth_scaler,dimension = 1024):
        # Moving gratings
        dimension = 1024 # It needs to be square power-of-two (e.g. 256 x 256) for PsychoPy
        
        # Calculation of BG and FG depending on the michelson contrast definition
        # fg-bg = c * 2 * l
        # fg+bg = 2 * l                    
        contrast = self.michelson_contrast
        luminance = self.mean_luminance
        fg = contrast * 2 * luminance
        bg = 2*luminance - fg
        f = 1# generate a single cycle

        # Generate 1D wave and modulate the luminance and contrast
        if self.type == "sin":
            
            x = np.arange(dimension)
            # Wave needs to be scaled to 0-1 so we can modulate it easier later
            sine_signal = (np.sin(2 * np.pi * f * x / dimension)/2 +0.5)
            # We need to scale and shift the wave to match the fg bg values
            # We also need to scale for bit depth since DLP has 6 bits
            oneD_wave = sine_signal * 2*(fg - bg) * bit_depth_scaler  # Scaling the signal
            oneD_wave = oneD_wave -1 + bg 
        else:
            raise ValueError("Grating type {s} doesn't exist".format(\
                s=self.type))
        grating_texture = np.tile(oneD_wave, [dimension,1])
        return grating_texture

    def defineStimPos_v2_underConst(self,proj_params):
        """ To determine where the stimulus will be located."""

        
        if (self.stim_type == 'edges-v1') or \
            (self.stim_type == 'movingStripe-v1'):
            
            posX = -(proj_params['sizeX']/2)/ np.cos(np.deg2rad(self.direction_deg))
            posY = -(proj_params['sizeX']/2)/ np.sin(np.deg2rad(self.direction_deg))
            stim_pos = (posX,posY)
            
            diag = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            diag*np.cos(np.deg2rad(self.direction_deg + 45.0))
            angle_y = np.degrees(np.arcsin(proj_params['sizeY'] / diag))
            # Start the stimulus according to the pre-defined directions.
            # Also calculates the distance that the edge or the stripe needs to travel
            # to cover the whole screen
            if self.direction_deg <= 90:
                stim_pos = (-proj_params['sizeX']/2,-proj_params['sizeY']/2)
                distance_to_travel = diag * np.sin(np.deg2rad(angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            elif self.direction_deg <= 180:
                stim_pos = (-proj_params['sizeX']/2,proj_params['sizeY']/2)
                distance_to_travel = diag * np.cos(np.deg2rad(90-angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            elif self.direction_deg <= 270:
                stim_pos = (proj_params['sizeX']/2,proj_params['sizeY']/2)
                distance_to_travel = diag * np.sin(np.deg2rad(angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            elif self.direction_deg <= 360:
                stim_pos = (proj_params['sizeX']/2,-proj_params['sizeY']/2)
                distance_to_travel = diag * np.cos(np.deg2rad(90-angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            else:
                raise NameError(f"Stimulus type {self.stim_type} direction has to be given in degrees.")
            
        else:   
            raise NameError(f"Stimulus type {self.stim_type} positions are not defined.")
        
        return stim_pos, distance_to_travel
        
    def defineStimPos(self,proj_params):
        """ To determine where the stimulus will be located."""

        
        if (self.stim_type == 'edges-v1') or \
            (self.stim_type == 'movingStripe-v1'):
            diag = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            angle_y = np.degrees(np.arcsin(proj_params['sizeY'] / diag))
            # Start the stimulus according to the pre-defined directions.
            # Also calculates the distance that the edge or the stripe needs to travel
            # to cover the whole screen
            if self.direction_deg <= 90:
                stim_pos = (-proj_params['sizeX']/2,-proj_params['sizeY']/2)
                distance_to_travel = diag * np.sin(np.deg2rad(angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            elif self.direction_deg <= 180:
                stim_pos = (-proj_params['sizeX']/2,proj_params['sizeY']/2)
                distance_to_travel = diag * np.cos(np.deg2rad(90-angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            elif self.direction_deg <= 270:
                stim_pos = (proj_params['sizeX']/2,proj_params['sizeY']/2)
                distance_to_travel = diag * np.sin(np.deg2rad(angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            elif self.direction_deg <= 360:
                stim_pos = (proj_params['sizeX']/2,-proj_params['sizeY']/2)
                distance_to_travel = diag * np.cos(np.deg2rad(90-angle_y+self.direction_deg))
                distance_to_travel = np.abs(distance_to_travel)
            else:
                raise NameError(f"Stimulus type {self.stim_type} direction has to be given in degrees.")
            
        else:   
            raise NameError(f"Stimulus type {self.stim_type} positions are not defined.")
        
        return stim_pos, distance_to_travel

class OutputInfo(object):
    """ 
        For constructing an output organization
    """
    def __init__(self,meta=None):

        self.meta = meta
        self.sample_num = []
        self.sample_time = []
        self.sampling_interval = []
        self.imaging_frame = []
        self.frame_time = []
        self.stimulus_epoch = []
        self.stim_frame = []
        self.epoch_time = []
        self.stim_info1 = []
        self.stim_info2 = []
        self.stim_info3 = []
    
    def saveOutput(self,save_loc):
        """ Saves the stimulus output structure in .txt, .pickle and .mat formats"""
        save_str = f"""stimulus_output_NAME-{self.__dict__['meta']['stim_name']}"""\
                    f"""_NIDAQ-{self.__dict__['meta']['nidaq']}"""\
                    f"""_DATETIME-{self.__dict__['meta']['date_time']}"""
        save_strTxt = f"{save_str}.txt"

        file_h = open(os.path.join(save_loc,save_strTxt),'w')

        row_n = len(self.__dict__['sample_num'])
        # Write the column names
        for key in self.__dict__.keys():
            if not(key == "meta"):
                file_h.write(key)
                file_h.write('\t')
        file_h.write('\n')

        for row_num in range(row_n):
            for key in self.__dict__.keys():
                if not(key == "meta") and not(len(self.__dict__[key]) ==0):
                    file_h.write(str(self.__dict__[key][row_num]))
                    file_h.write('\t')
            file_h.write('\n')

        file_h.close()

        # Saving as a .pickle file for easier Python processing
        saveDict = {}
        for key in self.__dict__.keys():
            saveDict[key] = self.__dict__[key]

        savePath = os.path.join(save_loc, f"{save_str}.pickle")
        saveVar = open(savePath, "wb")
        pickle.dump(saveDict, saveVar, protocol=2) # Protocol 2 (and below) is compatible with Python 2.7 downstream analysis
        saveVar.close()

        # Saving as a .mat file
        sio.savemat(os.path.join(save_loc, f"{save_str}.mat"), saveDict)
        
        print(f'Output .txt, .mat and .pickle files saved:\n{os.path.join(save_loc,save_str)}')



    