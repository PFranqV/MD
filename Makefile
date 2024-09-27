MD1.exe: MD1.o MD_mod.o 
	gfortran -o MD1.exe MD_mod.o MD1.o
	./MD1.exe &
MD_mod.o: MD_mod.f90
	gfortran -c MD_mod.f90
MD1.o: MD1.f90 
	gfortran -c MD1.f90 

.PHONY:plot1
plot1: 
	gnuplot MD1_VVERLET_velDISTR.gnu
	gnuplot MD1_VVERLET_EnPmomTinst.gnu
	mv *.jpeg MD1/
	mv *.dat MD1/
