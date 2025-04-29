# Controls Files

The code for the Control files was created in MATLAB (R2021b) Update 7. The code is in format **mlx** for the Live Editor. 
Each folder has its own readme_xxx.pdf, which is an export of the **mlx**.

There are 2 folders in the Controls files:

- GSalgo

- Simulation


## Gerchberg-Saxton algorithm (GSalgo folder)
The Gerchberg_Saxton_algorithm.mlx solves the GS algo for one simulated surface image (Zernike poly order 6th).
The files for this are:

* Gerchberg_Saxton_algorithm.mlx

* mksprecon.m

* sprecon.m

* spunwrap.m

* zernike6th.mat

To run this code Run the **mlx** file.

## Simulation (Simulation folder)
A Simulation of the 2 actuators working for this project is created in Simulink. The Constants_for_Simulation.mlx file defines the constants for the motor and PZA, creates the transfer functions and PID controller. The code generates the closed-loop transfer functions for the actuators and for the dual architecture and shows a brief analysis of the Step-response for each.
The code executes the simulation and loads the results of the signals in 3 different scopes.
The files for this are:

* Constants_for_Simulation.mlx

* DualFdbk_PZAwH_MotorDC.slx

To run this code  Run the **mlx** file.
