===============================================
MWInv - Microwave Inversion for Imaging Objects
===============================================

.. contents::

Introduction
------------

This program implements a microwave imaging system using inverse scattering theory based on electromagnetic wave propagation. 
The goal is to develop a computational model that can reconstruct the internal structure of an object from microwave measurements.
The program will divided into the following steps: 

Simulation Setup
----------------
This will be a two-dimensional simulation. The object to be imaged will be represented as a dielectric distribution in the centre of the region of simulation.
The object will be surrounded by an array of point microwave dipole sources and receivers, which model microwave antennas. 
At a time, one point acts as a source and the rest as receivers. This goes on for all points, and the data is recorded for each source-receiver pair.
The number of sources and receivers will be adjustable, and the distance between them and the object will be adjustable as well.
A  Perfectly Matched Layer (PML) boundary is implemented beyond the sources and receivers.

Inverse Scattering Algorithm
----------------------------
The algorithm to be used is the Distorted Born Iterative Method (DBIM). This method is based on the Born approximation for weak scattering. 
It allows for constructing a linear system of equations of the form :math:`Ax = b`, where :math:`x` is the dielectric distribution of the object, :math:`A` is a matrix that depends on the space and sources and receivers, and :math:`b` is the measured data.
It is an iterative solver, which starts with an initial guess of the dielectric distribution, and then updates the guess based on the difference between the measured and simulated data over several iterations.
The algorithm is implemented as an optimization problem, where the objective function is based on the difference between the measured and simulated data, and the dielectric distribution of the object will be the optimization variable.

Program Implementation
-----------------------

The program has the following structure:

#. Model Environment Setup: This involves classes for the object, sources, receivers, simulation region and the boundary.
#. "Measured Data" Generation: For different structures like circle in a region, I will use MEEP (a finite-difference time-domain (FDTD) simulation software) to generate the "measured data"
 (here, I mean measured data as in reference data that corresponds to the actual object - hence it should be accurate). However, a simpler method is to simply get the data from the forward simulation, and add noise to it.
#. Simulation Loop: The simulation loop will have the following steps, based on the DBIM algorithm:
    #. Generating a guess for dielectric distribution, starting from a constant distribution.
    #. Using sources to transmit electromagnetic waves in the simulation region, which then get scattered by the object and are recorded by the receivers. This I do using a Finite Element Library called :code:`MFEM`.
    #. Solving the optimization problem using the data recorded by receivers for this guess dielectric, and the measured data.
    #. Updating the guess for dielectric distribution based on how good the current guess is, then repeating the process from step b.
#. Visualization: This will involve using Python to convert the data generated to images. For now, :code:`MFEM` has a built-in visualization tool :code:`GLViS`, which I am using currently.

Model Validation
---------------

I will begin by using simple objects, such as a circle or a square, as the object to be imaged. Their dielectric distribution will be constant, and I will compare the reconstructed images with the actual objects.
When these are perfectly done, I will then move on to more complex objects, such as layers of the Earth or a human brain model, which can be decided later. The human brain model could be easier, since it is bounded in region of simulation.

Results and Analysis
--------------------

The major results will be images of reconstructed model, error between reconstructed model and actual object, and the time taken for the reconstruction. 
The error will be calculated using the L2 norm of the difference between the reconstructed and actual dielectric distribution.

Conclusion
----------

MWInv is a project that uses inverse theory based on scattering of electromagnetic waves to reconstruct the internal structure of an object.
There are a few challenging aspects, like implementing the DBIM algorithm, and correctly implementing the math that goes into simulation.
The project will be complete if I can reconstruct simple objects with near perfect accuracy, and if I can reconstruct more complex objects with a good enough accuracy (will need test results).

References
----------

There are many research papers, books and other references that I am using for the project. The complete list of references will be added as soon as possible. However, here is one key reference that gives a good idea of the theory and method used (not exact, but good for a fundamental understanding):

#. Lu, P.; Kosmas, P. Three-Dimensional Microwave Head Imaging with GPU-Based FDTD and the DBIM Method. Sensors 2022, 22, 2691. https://doi.org/10.3390/s22072691
#. Lu, Pan, Juan Córcoles, and Panagiotis Kosmas. "Enhanced FEM-based DBIM approach for two-dimensional microwave imaging." IEEE Transactions on Antennas and Propagation 69.8 (2020): 5187-5192.
#. Pastorino, Matteo. Microwave imaging. John Wiley & Sons, 2010.
