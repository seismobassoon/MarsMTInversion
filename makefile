FC      = gfortran
FFLAGS  = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0
# if ifortran comment the above and use the folowings:
#FC=ifort
#FFLAGS = -O4 -heap-arrays -check nobounds -xAVX -ftz -assume buffered_io -assume byterecl  -implicitnone -warn truncated_source  -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel
# or for the debugging
#FFLAGS= -check all -heap-arrays -debug -g -O0 -fp-stack-check -traceback -ftrapuv -assume byterecl

OBJS0    = module.f90 sgt.f90 others.f90 main.f90 additionalOthers.f90
program0= MarsMTInversion
#program1=tipsvSingle

#all : $(program0) $(program1)
all : $(program0)

$(program0): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS0)  -o $@ 

#$(program1): $(OBJS2)
#                $(FC) $(FFLAGS) $(OBJS2)  -o $@ 
                
clean:
	rm -f MarsMTInversion *.o *.mod
