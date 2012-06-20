SOURCES = vector.f90 proceed.f90 settings.f90 force.f90 cluster.f90 property.f90
OBJECTS = ${patsubst %.f90, %.o, $(SOURCES)}
FC = gfortran
all: cluster
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@
cluster: $(OBJECTS)
	$(FC) $(OBJECTS) -o $@
depend: 
	perl ./f90makedepend $(SOURCES) > .depend

include .depend
