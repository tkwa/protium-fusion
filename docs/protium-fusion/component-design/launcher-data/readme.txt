#################################################################################
#                                                                               #
# First-Principles Equation of State Database for Warm Dense Matter Computation #
#                                                                               #
# Burkhard Militzer                                      Berkeley, CA, 12-16-20 #
#                                                                               #
#################################################################################

(1) To compile our FPEOS code, simply type 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make

which will call Gnu's C++ compiler g++ with the argument -O3. The output should be

g++ -O3  -c FPEOS_29.C
g++ -O3  -c EOSTableSimple.C
g++ -O3  -c Parser.C
g++ -O3  -c ReadInTable.C
g++ -O3  -c Standard.C
g++ -O3  -c Form.C
g++ -O3  -c Physics.C
g++ -O3  -c Spline.C
g++ -O3  -c FermiGas.C
g++ -O3  -c RelativisticFermiGas.C
g++ -O3  -c Debye.C
g++ -O3  -c GKIntegration.C
g++ -O3  -c RungeKutta.C
g++ FPEOS_29.o EOSTableSimple.o Parser.o ReadInTable.o Standard.o Form.o Physics.o Spline.o FermiGas.o RelativisticFermiGas.o Debye.o GKIntegration.o RungeKutta.o -lm -o fpeos

In case this fails, try editing the Makefile or switch to a different compiler. 
The code has been compiled successfully with these three compilers:

Linux g++ version 4.8.4 (Ubuntu 4.8.4-2ubuntu1~14.04.4) 
Linux g++ version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04)
Apple clang version 12.0.0 (clang-1200.0.32.2)

To run our code, you need a C++ compiler. Python needs to be installed to 
generate the graphs (but C++ code will work without Python). Linux machines 
typically have the g++ compiler and Python already installed.

On Macs, we recommend installing XCode with the command line tools as well as 
installing Python with Anaconda (but there could be more disk space conserving
options to install a C++ compiler and Python).

On Windows, install Python as well as Cygwin with g++ version 10.2.1. 
Then run fpeos.exe and open the resulting PDF files with cygstart. 

(2) Files included
~~~~~~~~~~~~~~~~~~

The EOS tables are stored in files

Al_EOS_09-18-20.txt      C16H24_EOS_09-18-20.txt  LiF_EOS_09-18-20.txt     Ne_EOS_09-18-20.txt
B4C_EOS_09-18-20.txt     C18H18_EOS_09-18-20.txt  MgO_EOS_09-18-20.txt     O_EOS_09-21-20.txt
BN_EOS_09-18-20.txt      C20H10_EOS_09-18-20.txt  MgSiO3_EOS_09-18-20.txt  Si_EOS_10-19-20.txt
B_EOS_09-18-20.txt       C_EOS_09-18-20.txt       Mg_EOS_09-18-20.txt
C10H40_EOS_09-18-20.txt  H_EOS_09-18-20.txt       N_EOS_09-18-20.txt
C14H28_EOS_09-18-20.txt  He_EOS_09-18-20.txt      Na_EOS_09-18-20.txt

The C++ files are

Array.C                            Form.C                             Random.h
Array.h                            Form.h                             ReadInTable.C
DataElement.h                      Function.h                         ReadInTable.h
Debye.C                            GKIntegration.C                    RelativisticCorrectionsFermiGas.h
Debye.h                            GKIntegration.h                    RelativisticFermiGas.C
EOSTableSimple.C                   HeapSort.h                         RelativisticFermiGas.h
EOSTableSimple.h                   IdealMixture.h                     RungeKutta.C
Elements.h                         OrderedArray.h                     RungeKutta.h
FPEOS_29.C                         ParseCommandLineArguments.h        Spline.C
FermiGas.C                         Parser.C                           Spline.h
FermiGas.h                         Parser.h                           Standard.C
FindMinimum.h                      Physics.C                          Standard.h
FindRoot.h                         Physics.h                          Vector.h

The C++ code calls the following Python scripts to write graphics files in .png and .pdf formats.

FPEOS_P_T_plot03.py             FPEOS_comp_P_plot04.py          FPEOS_mixture_rho_T_plot01.py
FPEOS_T_E_Debye_plot03.py       FPEOS_comp_T_plot03.py          FPEOS_mixture_up_us_plot01.py
FPEOS_T_E_Ideal_plot03.py       FPEOS_mixture_P_T_plot01.py     FPEOS_rho_T_plot03.py
FPEOS_T_P_Debye_plot03.py       FPEOS_mixture_comp_P_plot01.py  FPEOS_up_us_plot01.py
FPEOS_T_P_Ideal_plot03.py       FPEOS_mixture_comp_T_plot01.py

These files should work with Python 2.7 and 3.7 as long as the following import commands work:

from matplotlib.pylab import *
import numpy
import sys
import os 

(3) Execution
~~~~~~~~~~~~~

To obtain the latest instructions type

fpeos

which will print the elements and compounds in the database:

 1: Hydrogen     8: Sodium       15: CH4	     
 2: Helium	 9: Magnesium  	 16: CH2	     
 3: Boron	10: Aluminum   	 17: CH1.5    
 4: Carbon	11: Silicon    	 18: CH	     
 5: Nitrogen	12: LiF	       	 19: CH0.5    
 6: Oxygen	13: B4C        	 20: MgO	     
 7: Neon	14: BN         	 21: MgSiO3   

You can run our fpeos code in the following ways:

For helium EOS:         fpeos EOS=2 
For B4C    EOS:         fpeos EOS=13
For MgSiO3 EOS:         fpeos EOS=21
For H2+O mixture:       fpeos binaryMixture EOS1=1  2.0  EOS2=6  1.0 
For He+Ne mixture:      fpeos binaryMixture EOS1=2  1.0  EOS2=7  1.0 
For H+He mixture:       fpeos binaryMixture EOS1=1  0.92 EOS2=2  0.08
For MgO+MgSiO3 mixture: fpeos binaryMixture EOS1=20 1.0  EOS2=21 1.0
For C+O  mixture:       fpeos binaryMixture EOS1=4  1.0  EOS2=6  1.0  rho0=1.0426 E0=-112.9115
For C+O2 mixture:       fpeos binaryMixture EOS1=4  1.0  EOS2=6  2.0  rho0=1.6381 E0=-188.1576

rho0[g/cc] and E0[Ha/FU] overwrite the default values for initial conditions of Hugoniot curves.

The extra numbers are number fractions that define the formula unit (FU),
which is used to normalize extensive quantities like volume and energy.
To normalize per C2H6 unit, write:  fpeos binaryMixture EOS1=15 1.0   EOS2=16 1.0  
To normalize per CH3  unit, write:  fpeos binaryMixture EOS1=15 0.5   EOS2=16 0.5  
To normalize per atom,      write:  fpeos binaryMixture EOS1=15 0.125 EOS2=16 0.125

The FPEOS code will write a number of .txt files, then call Python to make plots, and
finally try to open the .pdf files. To identify all new output files, type "ls -ltr"

Add the options "doNotCallPython" or "doNotOpenPDFFiles" if you prefer. The arguments
"processAllEOSs","processAllMixturesOfElements", "processAllMixtures",
"hugoniotCompressionMaximum", and "noLogInterpolation" are available also.

(4) Resulting files
~~~~~~~~~~~~~~~~~~~

In standard Linux fashion, type "ls -ltr" to figure out, which files have been written. 
For "fpeos EOS=13", the result should be

FPEOS_convex_hull.txt
FPEOS_isobars.txt
FPEOS_isochor_points.txt
FPEOS_isochores_ideal_Debye.txt
FPEOS_isotherms.txt
FPEOS_isotherm_points.txt
FPEOS_adiabats.txt
FPEOS_Hugoniot.txt
FPEOS_Hugoniot_lower_initial_density.txt
FPEOS_Hugoniot_higher_initial_density.txt
FPEOS_Hugoniot_rel.txt
FPEOS_Hugoniot_rad.txt
FPEOS_T_P_Debye_plot03.pdf
FPEOS_T_P_Debye_plot03.png
FPEOS_T_E_Debye_plot03.pdf
FPEOS_T_E_Debye_plot03.png
FPEOS_T_P_Ideal_plot03.pdf
FPEOS_T_P_Ideal_plot03.png
FPEOS_T_E_Ideal_plot03.pdf
FPEOS_T_E_Ideal_plot03.png
FPEOS_P_T_plot03.pdf
FPEOS_P_T_plot03.png

The .txt files were written by the fpeos executable. The .png and .pdf files were 
produced by the Python scripts that fpeos called.

Best wishes and happy FPEOS computation! 

Burkhard Militzer                                       Berkeley, CA, January 7, 2021
