This folder contains the files to simulate the surface form for a Mercury and air interface using COMSOL and additional 
files/code for processing files in Matlab and Excel.

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