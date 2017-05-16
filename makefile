# call : make "name of program"
FC=gfortran
FFLAG1= -mcmodel=medium -ffixed-line-length-none 
OBJ1 = get_SWZH.o four.o dofilt.o getper.o sort.o get_amp_u.o gaussfilter.o envelope.o hilbert.o correlate.o zfour.o sacio.o getalpha.o
OBJ2 = statistic_swzh.o sort.o getper.o
all: sacio.mod getalpha.mod ../bin/MFT_SWZH ../bin/statistic_swzh
%.o:%.f90
	$(FC) -c $^
sacio.mod:sacio.o
	$(FC) -c $^
getalpha.mod:getalpha.o getalpha.f90
	$(FC) -c $^
../bin/MFT_SWZH: $(OBJ1) 
	$(FC) $(FFLAG1) -o $@ $^
../bin/statistic_swzh: $(OBJ2)
	$(FC) $(FFLAG1) -o $@ $^
uninstall:
	-rm ../bin/MFT_SWZH ../bin/statistic_swzh
clean:
	-rm *.o *.mod
