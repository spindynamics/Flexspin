CC = mpicxx -std=c++11 -I. -c


OBJECTS=main.o read_input.o pulse.o simulation_range.o engine.o prints.o path_resolver.o junction_init.o Ndip.o evol_filename.o CalcHeff.o LLG_Heun1.o LLG_Heun2.o LLG_RK4.o LLG_Symplectic.o Rand_dist.o


.SUFFIXES: .cpp .o

Flexspin : 		$(OBJECTS)
	mpicxx -O3  $(OBJECTS) -o Flexspin  -lm -ljson-c



.cpp.o:
		$(CC) $<

.cc.o:
		$(CC) $<

clean: 
		rm -f *.o *~
cleanall: 
		rm -f *.o *~ Flexspin





