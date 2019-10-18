MULTALL-OPEN
Multall-open is a turbomachinery design system based on the 3D, multistage, Navier-Stokes solver Multall.

It can be used to design axial, mixed and radial flow turbines and compressors.

The system is written in FORTRAN77 and can be freely downloaded using the link given on this site.

General Description
Multall-open is a suite of programs for designing turbomachinery blading which has been developed over many years by Prof. John Denton, formerly of the Whittle Lab, Cambridge University Engineering Department, UK, now retired. 

The suite comprises 3 programs, all written in FORTRAN77:

1.    A one-dimensional mean-line program "meangen" which receives input either from the screen or from a file and writes a data input file for "stagen" .

2.    "stagen" a blade geometry generation and manipulation program that produces blade sections on a selection of stream surfaces, stacks the blades, combines them into multiple stages and writes an input file for the 3D solver "multall" .

3.    "multall" a 3D Navier-Stokes solver for multistage turbomachines. The program can perform quasi-3D blade to blade flow and axisymmetric throughflow calculations as well as fully 3D calculations. It predicts the machine overall performance as well as the detailed flow pattern.

The system should work on any computer with a FORTRAN compiler.

**Owner
Former director of the Whittle Lab, and emeritus Professor John Denton (FRS) has released to the public his suite of turbomachinery design tools (Multall).