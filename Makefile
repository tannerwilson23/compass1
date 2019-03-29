FC = gfortran
SRC = allocate.f90 setup.f90 energy.f90 output.f90 spline_kernal.f90 density.f90 ghosts.f90 eos.f90 viscosity.f90 acceleration.f90  smoother.f90 deriven.f90 time_stepping.f90 main.f90
OBJ = ${SRC:.f90 =.o}
FFLAGS = -Wall -fdefault-real-8

fdefault: main

main.o: init.o allocate.o output.o


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

main: ${OBJ}
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm main *.o *.mod
