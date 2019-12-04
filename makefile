# makefile Butterworth 

FC      = gfortran
FFLAGS  = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0
OBJS0    = module.f90 sgt.f90 others.f90 main.f90
program0=MarsMTInversion
#program1=tipsvSingle

#all : $(program0) $(program1)
all : $(program0)

$(program0): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS0)  -o $@ 

#$(program1): $(OBJS2)
#                $(FC) $(FFLAGS) $(OBJS2)  -o $@ 
                
clean:
	rm -f MarsInversion *.o *.mod
