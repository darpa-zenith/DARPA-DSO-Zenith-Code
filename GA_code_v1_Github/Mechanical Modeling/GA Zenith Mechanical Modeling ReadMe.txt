V1 Notes

---
CHANGES

Included the phase unwrapping .py script intended for simulating piston modes. This code produced a .mat that would be used to test the capability of a Gerch-Berg Saxton Algorithm.
---

This folder contains the ANSYS 2021 R2 Project for the 12 arm Multi-mode Deformable Mirror (MDM) dish used for finite element analysis (FEA), and the post porcessing codes used for subsequent analysis.

The ANSYS simulation is used to model strength of materials and actuator influence functions of the dish while under inertial loads.

ANSYS projects are sensitive to subfolder renaming and relocation, so be careful.

Do not modify the folder structure unless you are certain.

---

1) EXPORTING DATA FROM ANSYS

Unit system; ALWAYS USE MILLIMETERS for these prep/post analysis steps.

When exporting ANSYS nodal results, ensure you are exporting results from a surface only (not a volume).

Users will have to navigate to specific load steps or configure displacments as needed ot produce desired FEA results for influence function, interial loads, etc. They project was saved in a state that enables users to readily extract these results, and also includes a project archive in case they wish to revert to the original model.

a) Open the "Influence Function" system in ANSYS Workbench
b) Select the displacement results for the mirror's top parabolic surface
b) Right click the result and "Export"
c) Save .txt results in the "GA Zenith Mechanical Modeling Outputs" subfolder.

---

2) RESULT OUTPUTS

The "GA Zenith Mechanical Modeling Outputs" Subfolder contains various .txt, .py, .dat, and .mat files.

.txt files are nodal results exported from ANSYS in mesh coordinate x,y,z,USUM format. The .py scripts will process the .txt file into .mat or .dat files for use in subsequent analysis. This is required becuase meshes nodal results are not a regular grid. In short, the .py scripts interpolate the FE nodal data and save it as regular grid data so it can be used by Matlab and Zemax.

---

3) POST - INFLUENCE FUCNTION RETIVAL

Extract influence function by running the "ANSYS_nodes_to_grid_031424.py" script. You will have to modify the script to read the latest .txt nodal result file, or your own defined results file.

You will also have to re-run for every influence function, so it is recommended you name each .txt result in a logical fashion. Examples are already included in the folder.

The script produces a .mat file of the influence funciton, and will name them based on the .txt files name.

Important: Using influence functions to construct a complete matrix for the actuator system is described in the other "READ ME" found within "Controls". This is ultimately used to run a Least Square Error Analysis (LSEA) using matrix math.

After generating .mat files for the center actuator, 1 inner radius actuator, and 1 outer radius actuator, transfer these. mat files to the respect "Controls" folder. Remaining IF's are genereated using a rotatation matrix transformation.

---

4) POST - ZEMAX RAY TRACE ANALYSIS

Zemax's grid-sag analysis uses a .dat file to interpolate a sag onto an ideal optical surface.

Run "ANSYS_nodes_to_DAT_031424.py" to convert a nodal displacment result (sag) into a .dat file that's ready to import.

The .dat is ready to import into Zemax.

---

5) OTHERS

There are other .py codes found within. These generate non-zernike noise parameters that users may find interesting.