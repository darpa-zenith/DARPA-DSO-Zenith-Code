# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 23:40:51 2024

@author: jaris
"""

import QuickBasler as qb
import LMTpy as lp
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
var_dict = get_ipython().__dict__['user_module'].__dict__
import pickle
from os.path import join
import os
import time

# =============================================================================
# load the instruments separately, inherits "inst" from kernel if it exists
# note that if inst.close_all() was called then inst.connect_all() must be called to re-establish comms
# =============================================================================
inst = lp.Instruments() if not 'inst' in var_dict.keys() else var_dict['inst']
#inst.motor.set_velocity_params(max_velocity=1,accel=1) #Motor 
path = r'C:\Users\EXAMPLE\Documents\DARPA_Zenith'  #Save path for all files
if not os.path.exists(path): os.mkdir(path)
exp = lp.Experiment(inst=inst,FILEPATH=path)

if not exp.motor.status.IsHomed:
    print("Warning: motor is not homed, position may not be accurate")
    
# =============================================================================
# move stage position and take timed measurements
# =============================================================================
exp.jog_motor_timed(STEP_SIZE=0.5,DURATION=5,IMG_DELAY=0.5,FOLDER_NAME='Test')
    #Data will be saved as pickle files to the folder C:\Users\EXAMPLE\Documents\DARPA_Zenith\Test

# =============================================================================
# take a timed measurement at a single location and return data as dictionary
# note that all of the arguments passed into timed_measurement can also be passed into jog_motor_timed and motor_scan_timed
# =============================================================================
data = exp.timed_measurement(DURATION=20,IMG_DELAY=0.5,START_DELAY=5,SHOW_IMGS=True)

# =============================================================================
# Take a measurement at a list of angles and return data as dictionary
# =============================================================================
degrees = np.linspace(0,10,101)
data = exp.motor_scan_timed(DEGREES=degrees,IMG_DELAY=0.01,DURATION=0.01,START_DELAY=0.01,FOLDER_NAME='TEST1',SHOW_IMGS=False)
    #Data will be saved as pickle files to the folder C:\Users\EXAMPLE\Documents\DARPA_Zenith\TEST1
plt.close('all')
exp.motor.move(exp.deg_to_pos(0))

# =============================================================================
# Load a pickle file
# =============================================================================
fname = 'some_file.pickle'
with open(fname,'rb') as f:
    data = pickle.load(f)
    
    
# =============================================================================
# Reset stage to level
# =============================================================================
exp.motor.home()
exp.motor.move(exp.deg_to_pos(0))


# =============================================================================
# Take an image and display it, but do not save anything
# =============================================================================
img = qb.get_img(exp.camera)
plt.figure()
plt.imshow(img)


# =============================================================================
# Wait for user input, then take picture and save as png
# Does not save raw data. Will take as many pictures as there are elements in the "volts" array
# The exact voltage will be measured and then printed on the picture, while the filename will correspond to the entry in the "volts" array
# =============================================================================
volts = [0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
# volts = [0]
savepath = join(path,r'TEST2')
    #Data will be saved as pickle files to the folder C:\Users\EXAMPLE\Documents\DARPA_Zenith\TEST2
if not os.path.exists(savepath): os.mkdir(savepath)
for v in volts:
    input('Press Enter to take a picture at %i volts'%v)
    img = qb.get_img(exp.camera)
    vol = exp.daq.voltage()
    fig = plt.figure(figsize=(16,9))
    num = fig.number
    plt.tight_layout()
    plt.imshow(img,origin='lower')
    plt.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)
    plt.text(25,50,'~%0.1f volts'%vol,color='white',fontsize='x-large')
    plt.savefig(join(savepath,'%i volts.png'%v),bbox_inches='tight')
    plt.close(num)
    print('Picture taken')
    
    
# =============================================================================
# Move to 10 degrees, then to -8.5 degrees, then back to 0 degrees, waiting 2 seconds at each position
# =============================================================================
exp.motor.move(exp.deg_to_pos(10))
time.sleep(2)
exp.motor.move(exp.deg_to_pos(-8.5))
time.sleep(2)
exp.motor.move(exp.deg_to_pos(0))