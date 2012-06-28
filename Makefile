TARGET: main scl
all: $(TARGET)
#General rule to make FOO from FOO.f90.
# % in rule line is the matching pattern for file name.
# $* in the method line is replaced by the matched pattern.
%: %.f90
	gfortran $*.f90 -o $*
depend: 
	perl ./f90makedepend main.f90 scf.f90 > .depend
clean:
	rm $(TARGET) *~ *.o *.mod

include .depend
