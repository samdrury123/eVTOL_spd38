# README.m
# - README for fan design scripts

Two files required are 'Main.m' and 'OptimiseMass.m'. Everything else is 
called from within these two functions.

*** Main.m *** 

- Function controlling fan design.
- Sets blade geometric parameters.
- Change variables in OptimiseMass
    
    Calls: 
        - OptimiseMass (SEE BELOW)
        - Mean_Line_Design (Get velocity triangles for flow)
        - MetalAngles (Get blade angles across span)
        - PrintGeometry (Produce .ibl of blade shapes)


*** OptimiseMass.m ***

- For a given operating point (OP), (phi, psi, omega), will determine size 
    of the blades (casing radius) at a fixed hub radius to satisfy
    thrust = mass.
- Outputs all of the performance metrics at design.
- Uses a mass-model to determine variation in mass and solves numerically
    for thrust = mass.
- Displays required minimum propulsor geometry.
- Displays contour maps of how performance varies with OP, with limits on 
    diffuser length and a limit on thrust = mass. (OPTIONAL)

    Calls:
        - OptimiseMassMESH (OPTIONAL 'mesh = 1') (Calculate critical casing radius (T = M) for
            each OP. Display as a contour map to illustrate design choice.)
        - chicsatdesign (OPTIONAL 'chics = 1') (Display characteristics of various 
            performance metrics at various speeds at design point.)