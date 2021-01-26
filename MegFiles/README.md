# eVTOL
CFD and experimental code for the eVTOL rig


CFD

Run_Datum defines the block boundary conditions for the meshes and patches together the meshes of the turbomachinery (from Autogrid) and the inlet and outlet (from Pointwise). TURBOSTREAM is used to run RANS simulation with poisson, mixing length and AR solutions.

Thrust_FoM_Calculate analyses the CFD results by creating cuts to extract flow properties. Control volume analysis is used to calculate the thrust components of each ducted fan section.

Experimental

To set up the eVTOL rig, the motor must be clocked and balanced using eVTOL_balance. The eVTOL rig is capable of a 1D radial traverse and 2D XY traverse. These scripts call on general scripts to set the speed using the pxie, control the stepper motors, and record data from the five hole probe. The five hole probe must be calibrated before use. These generic scripts can be found in the Whittle Lab repo.




