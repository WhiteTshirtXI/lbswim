#*****************************************************************
#*       Makefile for LBSwim                                     *
#*****************************************************************

COMP_SER = ifort
FLAG_SER = -g -traceback -check #Debugging
FLAG_SER =-FR -O3
FLAG_SER =-FR -O3 -ipo -no-prec-div -xHost 

COMP_PAR = mpiifort      #Aurora
COMP_PAR = mpif90        #Local machines
FLAG_PAR = -FR -O3 -ipo -no-prec-div -xHost -DMPI -lmpi # Aurora
FLAG_PAR = -DMPI -g -traceback #Debugging
FLAG_PAR = -FR -O3 -DMPI #Minimal optimization
FLAG_PAR = -FR -O3 -ipo -no-prec-div -xHost -DMPI -lmpi -lopen-rte -lopen-pal   # Local machines

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
          
serial  : lbswim_ser.o $(OBJS_SER)
	$(COMP_SER) $(FLAG_SER) -o lbswim_ser.exe lbswim_ser.o $(OBJS_SER) 

mpi     : lbswim_par.o $(OBJS_PAR)
	$(COMP_PAR) $(FLAG_PAR) -o lbswim_par.exe lbswim_par.o $(OBJS_PAR) 


lbswim_ser.o : lbswim.F90 module_ser.o
	$(COMP_SER) $(FLAG_SER) -o lbswim_ser.o -c lbswim.F90
lbswim_par.o : lbswim.F90 module_par.o
	$(COMP_PAR) $(FLAG_PAR) -o lbswim_par.o -c lbswim.F90

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

