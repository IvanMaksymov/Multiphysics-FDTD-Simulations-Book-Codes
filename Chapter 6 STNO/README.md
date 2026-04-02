
To compile: gfortran -O2 -o aef STNO_main.f90 config_STNO.f90 cdot.f90, where -O2 is the compilator optimisation option and 'aef' is the name of the executable file. To run the code, the text files H0.IN and KVEC.IN should be located in the same directory. The first file contains the value of the static magnetic field in Oe and the second contains the value of the dc current. 
