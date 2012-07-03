PREV=100Modular
TARGET=main gen_scl gen_2compo
all: $(TARGET)
#General rule to make FOO from FOO.f90.
# % in rule line is the matching pattern for file name.
# $* in the method line is replaced by the matched pattern.
%: %.f90
	gfortran -Wall $*.f90 -o $*
clean:
	rm $(TARGET) *~ *.o *.mod
diff:
	git diff -U5 --color $(PREV)
