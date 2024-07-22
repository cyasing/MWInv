=====================
Architectute of MWInv
=====================

This project has the following components, with their associated files:

Physical System Setup
----------------------
#. `MWISystemSetup.cpp` and `MWISystemSetup.h`: 
    #. This class holds all the information about the physical system setup that would be given by the experimenter/user.

Meshing (Based on `Gmsh` package)
---------------------------------
#. `DefaultMesh.cpp` and `DefaultMesh.h`:
    #. This class defines a default mesh that I consider optimal enough to do meshing.
    #. Takes a referenced object from `MWISystemSetup`.

#. `AdaptiveMesh.cpp` and `AdaptiveMesh.h`:
    #. This class defines an adaptive mesh that can be used to refine the mesh based on the forward solver results.
    #. Inherits from `BaseMesh`.
    #. Will be used by `ForwardSolve` only if the user wants to use adaptive meshing.

FEM Forward Solver (Based on `MFEM` package)
--------------------------------------------
#. `ForwardSolve.cpp` and `ForwardSolve.h`:
    #. This class defines the forward solver that would be used to solve the physical system.
    #. Takes a referenced object from `MWISystemSetup`, `BaseMesh` and `AdaptiveMesh`.

Optimization (Based on `dlib` package)
--------------------------------------

#. `OptimizationParams.h`: 
    #. This structure holds all the parameters required for all three types of optimization (Least Squares, Conjugate Gradient, and Newton).
#. `OptimzeBase.cpp` and `OptimizeBase.h`: 
    #. This class defines the base class for optimization. Implementing the functions will require different types of data, so it is left abstract.
    #. Takes a referenced object from the structure `OptimizationParams`.
#. `OptimizeLS.cpp` and `OptimizeLS.h`: 
    #. This class defines the optimization for least squares method.
    #. Inherits from `OptimizeBase`.
#. `OptimizeCG.cpp` and `OptimizeCG.h`: 
    #. This class defines the optimization for conjugate gradient method.
    #. Inherits from `OptimizeBase`.
#. `OptimizeNW.cpp` and `OptimizeNW.h`: 
    #. This class defines the optimization for newton method.
    #. Inherits from `OptimizeBase`.

Inversion
---------

#. `InverseSolve.cpp` and `InverseSolve.h`: 
    #. This class takes the forward solver in an iteratively optimizing loop, with the mesh also iteratively improving.
    #. Inherits from `ForwardSolve` and the optimization classes.

Main Run
--------
#. `run_mwinv.cpp` and `run_mwinv.h`: 
    #. This has the main function that runs the inversion by calling all class objects as required, with a list of all includes in this project.