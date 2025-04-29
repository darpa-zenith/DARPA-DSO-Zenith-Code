This folder contains the files to simulate the surface form for a Mercury and air interface using COMSOL and additional 
files/code for processing files in Matlab and Excel. 

Please note Code V0 Drop material is also included with no required revisions. 
Additional new items are mentioned below.

Code Ver 1:

FILE: 0.5m Pararbola Centered and Symm 4mm equal_4_21_24.mph
FILE: COMSOL_3Dver11.m
Description: The segemented version of the liquid mirror is simulated to 1.0 seconds with 2mm MRS.
2mm wide MRSs provide the lowest PV and separated by 2mm. The file starts w/ a 2mm air gap separation in order to 
tune the widths of each MRS to 3mm and 3.75mm and observe the PV. Electrowetting is applied to each MRS w/ mercury and air interface.
No electrostatic forces were added as the GUI interface is not suited to deterministically add 200 syntax codes for each liquid/wall interface thru a one line only GUI interface.

1. Go to Study  and click COMPUTE. The run time is  2hrs 4min
2. Go to Results and 1D Plot Group 7: Each lineout of the mercury and air interface is selected for export.
3. Right click on Line Graph 1 under 1D Plot Group 7 and click on "Add plot data to export".
4. In the Results section under Export: Plot XX will appear. Select File type .txt and Filename and EXPORT.
5. Open the file and manually remove headers and keep the two data columns.
6. The exported COMSOL lineout contains all sequence frames within a 1 second time interval. After opening the file in EXCEL, copy and paste the last frame and save as .csv file in order
to observe the PV of one time sequence. The last time frame or first time frame contains about 3800 rows that describe the surface form in 2D.
7. Open MATLAB FILE: COMSOL_3Dver11.m and run, a residual surface error plot will be generated with newly fitted coefficents for optical sag of the parabola.
Note:  The interpolated 3D surface will overshoot the true PV values. In order to observe the residual surface map quickly, it is recommended to use 4D Technologies, 4Sight software.
There are steps involved in the transfer process from COMSOL to Matlab to 4Sight( if 4Sight is available to the general public), but it does aid in an accurate and fast analysis to avoid
developement of the anaylsis yourself. 
8. To tune the width of each MRS is done by going into: Component 1-Geometry 1- Rectangle 1,2,3,4... and changing the width from 2mm to 3mm for all MRS.
9. To change the slewing rate please see accompanying word document.
10. To turn off gravity go to Laminar Flow-> Gravity-> Input for gravity the syntax for x and y as zero.
11. Changing the focal length is also possible.
Note: The parabolic fitted equation has been updated:zz=0.7503 + 0.0003334*xtest.^2 + 0.0003334*ytest.^2;
File outputs: 0.5m Full Parabola 4mm equal slew 1 deg per sec BEGIN.csv are provided from the .mph file.

FILE:High Resolution positioning 03_31_24 revised 04_20_24.xlsx
Description: The position and angles have been updated and how to use the excel sheet is explained in the file.

FILE: PS_x10 visc 800V 10um 01_31_24
Description: The continuous Hg film has an applied voltage of 800V across the top surface in a 2D axisymmteric model.
There is a volumetric force applied across mercury that uses two forces, a electrostatic force of attraction and frictional surface force.
Lastly a voltage controlled contact angle determined by the YL equation.
In order to bound the problem the viscosity is adjusted by a factor of x10 in order to overcome the roadblock of long oscillation times.

1. Go to Study 1 and click COMPUTE. The run time is about 35 minutes.
2. Scroll down to Results-> Select: 2D Plot Group 4 -> Surface 1
3. For surface 1, Click on green/red arrow to select the Physical Quantity of interest. Ex: Electric potential (V).
4. Select Animate above.
Note: The initial geometric boundary conditions are sensitive so the voltages selected and volume used have to be near the final equilibrium geometry
or else the simulation will not converge. This is presummably due to the fast changes of the liquid to air interface which produce a nonl

Code Ver 0:

COMSOL files simulates a 2D lineout for a mercury surface in a multi tub or container configuration.
The file incorporates slewing and other adjustable parameters as described in the word file:COMSOL Multiphysics Simulations for the liquid metal surface 04_01_24.docx
The post anaylsis uses the Matlab file: COMSOL_3Dver7.m to generate figures and a dataset for the residual surface form map based on the simulated surface form in
2D done in COMSOL. The excel file is a basic spread sheet that estimates the angle needed for each tub. The file name is: High Resolution positioning 03_31_24.xlsx
---
Demonstration of the Mercury-Air interface with multi containers in a curved geometry with slewing using COMSOL is provided in the file.
FILE: 2D Flat  x 11 2mm Hg and Air w ACDC Pulse FULLY ON  tilted at 10 deg NO Fes Z_ OSCILLA 2sec w SLEWING _PARAB ALIGNED  Center tub 0 deg OFF 4_01_24.mph
FILE: COMSOL_3Dver7.m

1. A word file is provided to describe additional details of the .mph file and export process of the liquid Mercury surface in 2D.
2. The .m file for matlab is used for post analysis and rotates the 2D lineout and provides an interpolated surface and residual surface error.
---
Demonstration of the Mercury-Air interface with a thin dielectric side walls made of quartz and electrostatic surface forces using COMSOL is provided in the file.
FILE: 2D Flat  20mm Hg e is 80 and Air AND DIELECTRIC w ACDC Pulse FULLY ON  _ NEW BC ON RHS _ ADDED ON Fes RHS BOT TOP 4_01_24 CORRECTED Left Side and FORCE_WORKS 200ms.mph

1. Go to Study 1 and click COMPUTE. The run time is about 1 hr and 8 minutes.
2. Scroll down to Results-> Select: 2D Plot Group 3 -> Surface 1
3. For surface 1, Click on green/red arrow to select the Physical Quantity of interest. Ex: Electric Potential.
4. Select Animate above to simulate the electric potential fields over time.

Note: For a metal conductor the vendor provider suggests setting the dielectric constant to one for a metal. I have set the value arbitrarly to 80.
Metals should have a dieletric constant of infinity so further tuning could be done to observe no field lines in the metal conductor.
---
Demonstration of Adhesion between the charge liquid metal and surface interface using COMSOL was performed systematically in the following file.
FILE: x10 WORKs TOP SURFACE FORCE ONLY _ALL SURFACES Inverted TALLER 03_04_24 FINAL CONFIGURATION Floor at 800V surf charge on bottom 02_18_2024.mph

1. Go to Study 1 and click COMPUTE. The run time is about 47 minutes.
2. Scroll down to Results-> Select: 2D Plot Group 4 -> Surface 1
3. For surface 1, Click on green/red arrow to select the Physical Quantity of interest. Ex: Volume Force z-dir.
4. Select Animate above.

Note: Steps for providing a partition or half ellipse shaped Mercury film are listed in the Geometry Tab.
---
The COMSOL files are related to the Continuous Film of Mercury with an applied voltage on a incline. 

Demonstration of electrostatic forces only on the top surface of mercury using Voltage/distance with a dieletric on the floor was performed.
The .mph file assumes a fricitional forces between the fluid and solid interface with a coefficent of 0.1 or 5.5 degrees.
FILE: FINAL CONFIGURATION Floor at 400V surf charge on bottom 02_18_2024.mph

Note: A domain specified dielectric was not integrated until March of 2024 and the threshold of slippage is between 400V and 800V.
Note: The electrostatic repulsion forces are not integrated between the top and bottom charged surfaces as the surfaces deviate away from a parallel plate configuration.
Note: The electrical potential does not integrate a Infinite Element Domain to specifiy the potential at infinite distances away to be zero in a mm scale simulation.
Consultation with COMSOL seems to be indicate that it is an unknown if CFD and AC/DC can be used simultaneously with the Infinite Element Domain feature.

1. Go to Study 1 and click COMPUTE. The run time is about 47 minutes.
2. Scroll down to Results-> Select: 2D Plot Group 4 -> Surface 1
3. For surface 1, Click on green/red arrow to select the Physical Quantity of interest. Ex: Volume Force r-dir.
4. Select Animate above.