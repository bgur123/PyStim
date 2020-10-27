# -*- coding: utf-8 -*-
"""

@author: bgur123

"""
import numpy as np
import re
import os
import pickle
import scipy.io as sio # To save .mat files
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
        stim_info_fname = os.path.join(os.path.dirname(stim_fname),'stimInfos',
                                        "stim-info_{s}.txt".format(s=stim_type))
        self.stim_type = stim_type
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
            orientation = np.mod(self.direction_deg-90,360) # direction & orientation orthogonal
            
            # Size of the other dimension is bigger to span all the screen
            grating = visual.GratingStim(
                win=win, name='grating',tex=grating_texture, 
                size=(proj_params['sizeX'], proj_params['sizeY']), ori = orientation,
                sf=1/self.spatial_wavelength,
                units=proj_params['unit'])

            grating.autoLog = False 
            self.grating = grating
        elif self.stim_type == 'edges-v1':
            # Moving edges with a luminance before the edge comes

            # Pre edge window luminance (a.k.a. background)
            # Scaling so it fits to -1 1 range
            # Scaling so it is presented accurately with desired bit depth
            self.win_lum = ((self.pre_lum*2)* bit_depth_scaler -1) 

            # Edge luminance
            edge_luminance = ((self.edge_lum*2)* bit_depth_scaler-1)   
            
            # Direction of movement
            orientation = np.mod(self.direction_deg-90,360)

            # The edge should fill the screen - the maximum possible size
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)

            rectangle = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        width=0, height=span, pos = self.defineStimPos(proj_params),
                        ori=orientation,lineWidth=0, 
                        fillColor=edge_luminance, fillColorSpace='rgb')
            rectangle.autoLog = False
            self.rectangle = rectangle

        elif self.stim_type == 'fff-v1':
            # Full field flashes
            span = np.sqrt(proj_params['sizeY']**2+proj_params['sizeX']**2)
            lum = ((self.lum*2)* bit_depth_scaler-1) 

            rectangle = visual.Rect(
                        win=win, name='rectangle', units=proj_params['unit'],
                        size=span,lineWidth=0, 
                        fillColor=lum, fillColorSpace='rgb')
            rectangle.autoLog = False
            self.rectangle = rectangle
        else:
            raise NameError(f"Stimulus type {self.stim_type} could not be initialized.")
    
    def defineStimPos(self,proj_params):
        """ To determine where the stimulus will be located."""
        if self.stim_type == 'edges-v1':
            orientation = np.mod(self.direction_deg-90,360)
            if orientation == 0:
                stim_pos = (-proj_params['sizeX']//2,0)
            elif orientation == 180:
                stim_pos = (proj_params['sizeX']//2,0)
            elif orientation == 90:
                stim_pos = (0,proj_params['sizeY']//2)
            elif orientation == 270:
                stim_pos = (0,-proj_params['sizeY']//2)
        else:
            raise NameError(f"Stimulus type {self.stim_type} positions are not defined.")
        
        return stim_pos

class OutputInfo(object):
    """ 
        For constructing an output organization
    """
    def __init__(self,meta=None):

        # We need to parse the .txt file to understand the stimulus routine structure.
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
        """ Saves the stimulus output structure in .txt and .pickle formats"""
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
        savePath = os.path.join(save_loc, f"{save_str}.pickle")
        saveVar = open(savePath, "wb")
        saveDict = {'outputObj' : self}
        pickle.dump(saveDict, saveVar, protocol=-1)
        saveVar.close()

        # Saving as a .mat file
        sio.savemat(os.path.join(save_loc, f"{save_str}.mat"), saveDict)
        
        print(f'Output .txt, .mat and .pickle files saved:\n{os.path.join(save_loc,save_str)}')



    