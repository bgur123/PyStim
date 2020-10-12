import PyDAQmx as daq

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

# DAQmx Start Code
daq.DAQmxStartTask(counterTaskHandle)
daq.DAQmxStartTask(pulseTaskHandle)

data = daq.uInt32(0)
daq.DAQmxReadCounterScalarU32(counterTaskHandle,
    1.0,daq.byref(data), None)

