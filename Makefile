MPIFC = mpif77
FFLAGS = -extend-source
RM = rm
INCLUDE = .
LIB = libaztec.a

.SUFFIXES: .o .f .f90

all: aztest

aztest:mymodule.o az_test.o solvsub.o mylib.o
	$(MPIFC) -o $@ $^ -I $(INCLUDE) $(LIB)

.f.o:
	$(MPIFC) $(FFLAGS) -c $< -o $@

.f90.o:
	$(MPIFC) $(FFLAGS) -c $< -o $@

clean:
	$(RM) *.o *.mod

realclean:
	$(RM) *.o aztest
