from psychopy import visual, core, monitors

mon = monitors.Monitor('iMac-bgur')
mon.setSizePix = (4023 ,2304)
mon.setWidth = 47.57
mon.setDistance = 57

refresh_rate = 60.0

pix_cm = mon.getSizePix()[0]/mon.getWidth()


total_degrees = 5

# Setup stimulus
win = visual.Window(
    size=(500 , 500), fullscr=False, 
    screen=0, winType='pyglet',useRetina=True,
    monitor=mon, color=[0,0,0], colorSpace='rgb')
# grating = visual.GratingStim(win, tex='sin', sf=5, name='grating',size=(1,1))


velocity = 2.5 # dps
spatial_wavelength = 2.5

sf = 1/spatial_wavelength
grating = visual.GratingStim(
            win=win, name='grating',tex='sin', 
            size=(total_degrees, total_degrees), 
            sf=sf,units='degFlat')

grating.autoDraw = True  # Automatically draw every frame
grating.autoLog = False  # Or we'll get many messages about phase change

# Let's draw a stimulus for 2s, drifting for middle 0.5s
epoch_clock = core.Clock()
epoch_duration = 10
while epoch_clock.getTime() < epoch_duration:
    # For moving the grating we will need to advance the phase
    # according to the desired velocity and screen refresh rate
    #| v (degree/s) / refresh rate (update/s) |
    grating.phase += (velocity/refresh_rate)/spatial_wavelength
    win.flip()