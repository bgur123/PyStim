# -*- coding: utf-8 -*-
"""

@author: bgur123

"""
import numpy as np
import re
import os
from psychopy import visual


class PyStimRoutine(object):
    """ 
        For constructing a particular set of stimulus consisting of pre-defined epochs.
    """
    def __init__(self,stim_fname):

        # We need to parse the .txt file to understand the stimulus routine structure.
        self.stim_fname = stim_fname
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
        stim_info_fname = os.path.join(os.path.dirname(stim_fname),
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
                oneD_wave = sine_signal * 2*(fg - bg) # Scaling the signal
                oneD_wave = oneD_wave -1 + bg 
            
            else:
                raise ValueError("Grating type {s} doesn't exist".format(\
                    s=self.type))
            grating_texture = np.tile(oneD_wave, [dimension,1])
            orientation = np.mod(self.direction_deg-90,360) # direction & orientation orthogonal
            grating = visual.GratingStim(
                win=win, name='grating',tex=grating_texture, 
                size=(proj_params['size'], proj_params['size']), ori = orientation,
                sf=1/self.spatial_wavelength,
                units=proj_params['unit'])
        
            grating.autoLog = False 
            self.grating = grating
        else:
            raise NameError(f"Stimulus type {self.stim_type} could not be initialized.")


class OutputInfo(object):
    """ 
        For constructing an output organization
    """
    def __init__(self,meta=None):

        # We need to parse the .txt file to understand the stimulus routine structure.
        self.meta = meta
        self.sample_num = []
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
        save_str = f"""stimulus_output_NIDAQ-{self.__dict__['meta']['nidaq']}"""\
                    f"""_DATETIME-{self.__dict__['meta']['date_time']}.txt"""
            
        file_h = open(os.path.join(save_loc,save_str),'w')

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
        print(f'Output file saved:\n{os.path.join(save_loc,save_str)}')



    