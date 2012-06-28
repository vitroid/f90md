TARGET: main scl
all: $(TARGET)
#General rule to make FOO from FOO.f90.
# % in rule line is the matching pattern for file name.
# $* in the method line is replaced by the matched pattern.
%: %.f90
	gfortran $*.f90 -o $*
clean:
	rm $(TARGET) *~ *.o *.mod
