
PROG =burger
OBJS = PROG

F90 = ifort

SRCS = 	zfft1d.f\
	int_to_char.f90\
	ci.f90\
	sa.f90\
	constantes.f90\
	burger.f90\
	


all:$(PROG)

$(PROG): $(SRCS) 
	 $(F90) $(SRCS)  -o $(PROG)

$(OBJS): $(SRCS)
	 $(F90)$(SRCS)  -o $(OBJS)

test: both
	./$(PROG) > f1.out
	diff -s c1.out f1.out

clean:
	rm -rf $(OBJS) $(PROG) *out
	rm -rf *.dat
	rm -rf *.err
	rm -rf *.log
	rm -rf *~
	rm -rf *.mod

backup:
	cd ./fftw3; make clean
	cd ..
	make clean
	tar -cf backup.tar  $(SRCS) makefile fftw3
	

	
