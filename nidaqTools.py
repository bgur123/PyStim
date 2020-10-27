import PyDAQmx as daq
from psychopy import core
def configure_daq():
    """ Configures NIDAQ communication """
    counterChannel = "Dev2/ctr1"
    pulseChannel = "Dev2/ctr0"

    counterTaskHandle = daq.TaskHandle(0)
    pulseTaskHandle = daq.TaskHandle(0)

    daq.DAQmxCreateTask("2",
        daq.byref(counterTaskHandle))
    daq.DAQmxCreateCICountEdgesChan(\
        counterTaskHandle,counterChannel,"",
        daq.DAQmx_Val_Rising,0,daq.DAQmx_Val_CountUp)
    daq.DAQmxCreateTask("1",\
        daq.byref(pulseTaskHandle))
    daq.DAQmxCreateCOPulseChanTime(\
        pulseTaskHandle,pulseChannel,"",
        daq.DAQmx_Val_Seconds,daq.DAQmx_Val_Low,0,0.05,0.05)

    
    
    return pulseTaskHandle, counterTaskHandle
    
def start_imaging(pulseTaskHandle, counterTaskHandle):

    daq.DAQmxStartTask(counterTaskHandle)
    daq.DAQmxStartTask(pulseTaskHandle)
    
    # Reads incoming signal from microscope computer and stores it to 'data'. A rising edge is send every new frame 
    # the microscope starts to record, thus the 'data' variable is incremented
    daq_data = daq.uInt32(0)
    # daq.DAQmxReadCounterScalarU32(counterTaskHandle,
    #     1.0,daq.byref(daq_data), None)
    return daq_data

def check_status_daq(daq_counter_h,daq_data,global_clock):
    """ """ 

    # check for DAQ Data
    daq.DAQmxReadCounterScalarU32(daq_counter_h,1.0,daq.byref(daq_data), None) 
    
    return (daq_data, daq_data.value, global_clock.getTime())

def clearTask(taskHandle):
    """
    Clears a task from the card.
    """
    daq.DAQmxStopTask(taskHandle)
    daq.DAQmxClearTask(taskHandle)