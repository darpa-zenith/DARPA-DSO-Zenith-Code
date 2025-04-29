This software suite was developed for control of and data acquisition for the DARPA Zenith experimental setup at General Atomics.
It is written in Python and relies on a ThorLabs BSC201 One-Channel Benchtop Stepper Motor for angle control, an MCC DAQ for measuring voltage, and a Basler camera for image acquisition.
The Python packages for communicating to these various devices are required for this software suite to function.

Examples of code structures that were used to take datasets have been placed in the "example.py"

Below are high level descriptions of commonly used commands.

LMTpy.Instruments()
	-Initiates the Instruments class and connects to all the instruments
	-Example usage:
		> import LMTpy as lp
		> inst = lp.Instruments()
		
LMTpy.Experiment()
	-Initiaties the Experiment class and optionally sets the save path for data acquisition
	-Allows access to more directly control instruments, or allows use of already created commands to automate data acquisition
	-keys:
		-inst: Instruments class object. If none are specified, the script will look for instruments to connect to by default
		-FILEPATH: path to the folder you want to save your data to. If not specified, a new folder will be made in the current working directory named the current date and time
	-Example usage:
		> import LMTpy as lp
		> exp = lp.Experiment(inst=inst,FILEPATH=path)
		
Commands in Experiment Class:
	jog_motor_timed()
		-Moves stage position and repeatedly takes images over a set amount of time
		-Allows user to input direction of step between measurements
		-Data is saved as a pickle file and includes current voltage reading
		-keys:
			-STEP_SIZE: number of degrees to move per step. Default = 0.5 degrees
			-DURATION: how long to spend taking images at each step. Default = 10 seconds
			-START_DELAY: how long to wait at new position before beginning to take images. Default = 5 seconds
			-IMG_DELAY: time between taking two consecutive images. Default = 0.5 seconds
			-FOLDER_NAME: folder name to save data to. If not provided, default will create a new folder named the current date and time
			-SHOW_IMGS: option to display the images as they are taken. Default = False
		
	timed_measurement()
		-Take a series of images at a single location
		-Data is returned as a dictionary that the user can then save
		-keys:
			-DURATION: how long to spend taking images at each step. Default = 10 seconds
			-START_DELAY: how long to wait at new position before beginning to take images. Default = 5 seconds
			-IMG_DELAY: time between taking two consecutive images. Default = 0.5 seconds
			-SHOW_IMGS: option to display the images as they are taken. Default = False
			
	motor_scan_timed()
		-Automatically scan through a list of angles, taking a series of images at each position
		-Data is saved, but can also be returned as a dictionary if a variable is assigned
			example: > data = exp.motor_scan_timed(DEGREES=degrees)
		-keys:
			-DEGREES: list of specific angles to take measurements at.
			-DURATION: how long to spend taking images at each step. Default = 10 seconds
			-START_DELAY: how long to wait at new position before beginning to take images. Default = 5 seconds
			-IMG_DELAY: time between taking two consecutive images. Default = 0.5 seconds
			-FOLDER_NAME: folder name to save data to. If not provided, default will create a new folder named the current date and time
			-SHOW_IMGS: option to display the images as they are taken. Default = False
			-PAUSE_BETWEEN: option to delay data acquisition at a new position until the user inputs a button to begin 
			
	motor.home()
		-Moves the motor to the home position to recalibrate it for angle measurements
	
	motor.move()
		-Moves the motor to a specified position in mm
		-keys:
			-pos_mm: Position within the motor's range of movement, in millimeters from home (0mm)
		
	deg_to_pos()
		-Calibrated conversion from the desired number of degrees to the correlating motor position
		-Used within the motor.move() command to move to a specific stage angle
			example: > exp.motor.move(exp.deg_to_pos(5)) #moves stage to be tilted 5 degrees
		-keys:
			-deg = Number of degrees to convert to motor position
			
	daq.voltage()
		-Returns the current voltage measurements
	
QuickBasler.get_img()
	-Takes an image using a Basler camera
	-keys:
		-camera: ID of camera to take image from
		-convert_to_rgb: Option to convert the image taken to an rgb color image. Default = True
	-example usage:
		> import QuickBasler as qb
		> img = qb.get_img(exp.camera)