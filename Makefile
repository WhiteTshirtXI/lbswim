#*****************************************************************
#*       Makefile for LBSwim                                     *
#*****************************************************************

COMP_SER = ifort
FLAG_SER =-FR -O3 -ipo -no-prec-div -xHost  #Local machines
FLAG_SER =-FR -O3

COMP_PAR = mpif90
FLAG_PAR = -FR -O3 -ipo -no-prec-div -xHost -DMPI -lmpi -lopen-rte -lopen-pal   # Local machines
FLAG_PAR = -FR -O3 -DMPI 

OBJS_SER=lb_ser.o \
         swimmers_ser.o \
         tracers_ser.o \
         module_ser.o \
         io_ser.o \
         aux_ser.o
          
OBJS_PAR=lb_par.o \
         swimmers_par.o \
         tracers_par.o \
         module_par.o \
         io_par.o \
         aux_par.o
          
serial  : LBswim_ser.o $(OBJS_SER)
	$(COMP_SER) $(FLAG_SER) -o LBswim.exe LBswim_ser.o $(OBJS_SER) 

mpi     : LBswim_par.o $(OBJS_PAR)
	$(COMP_PAR) $(FLAG_PAR) -o LBswim.exe LBswim_par.o $(OBJS_PAR) 


LBswim_ser.o : LBswim.F90 module_ser.o
	$(COMP_SER) $(FLAG_SER) -o LBswim_ser.o -c LBswim.F90
LBswim_par.o : LBswim.F90 module_par.o
	$(COMP_PAR) $(FLAG_PAR) -o LBswim_par.o -c LBswim.F90

lb_ser.o : lb.F90 module_ser.o aux_ser.o
	$(COMP_SER) $(FLAG_SER) -o lb_ser.o -c lb.F90
lb_par.o : lb.F90
	$(COMP_PAR) $(FLAG_PAR) -o lb_par.o -c lb.F90

swimmers_ser.o : swimmers.F90
	$(COMP_SER) $(FLAG_SER) -o swimmers_ser.o -c swimmers.F90
swimmers_par.o : swimmers.F90
	$(COMP_PAR) $(FLAG_PAR) -o swimmers_par.o -c swimmers.F90

tracers_ser.o : tracers.F90
	$(COMP_SER) $(FLAG_SER) -o tracers_ser.o -c tracers.F90
tracers_par.o : tracers.F90
	$(COMP_PAR) $(FLAG_PAR) -o tracers_par.o -c tracers.F90

module_ser.o : module.F90
	$(COMP_SER) $(FLAG_SER) -o module_ser.o -c module.F90
module_par.o : module.F90
	$(COMP_PAR) $(FLAG_PAR) -o module_par.o -c module.F90

io_ser.o : io.F90
	$(COMP_SER) $(FLAG_SER) -o io_ser.o -c io.F90
io_par.o : io.F90
	$(COMP_PAR) $(FLAG_PAR) -o io_par.o -c io.F90

aux_ser.o    : aux.F90
	$(COMP_SER) $(FLAG_SER) -o aux_ser.o -c aux.F90
aux_par.o    : aux.F90
	$(COMP_PAR) $(FLAG_PAR) -o aux_par.o -c aux.F90

tidy :
	rm -f *.o *.mod

clean :
	rm -f *.o *.mod *.exe

