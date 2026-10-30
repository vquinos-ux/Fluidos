# =====================================================================
#  Makefile  -  una variante del solver de combustion de HDPE
#  Identico para las 7 variantes: lo unico que cambia es mod_malla.f90.
#  Uso:   make        -> compila
#         make limpiar-> borra objetos, modulos y ejecutable
#         make datos  -> borra las salidas .dat y .png de la corrida
# =====================================================================

FC      = gfortran
FFLAGS  = -O2 -fopenmp
EXE     = combustion_pe

# El orden importa: cada modulo tiene que estar compilado antes de que
# lo usen los que vienen despues.
OBJS = mod_precision.o mod_termo.o mod_cinetica.o mod_malla.o mod_solver.o

all: $(EXE)

$(EXE): $(OBJS) combustion_pe.f90
	$(FC) $(FFLAGS) -o $(EXE) $(OBJS) combustion_pe.f90

mod_precision.o: mod_precision.f90
	$(FC) $(FFLAGS) -c mod_precision.f90

mod_termo.o: mod_termo.f90 mod_precision.o
	$(FC) $(FFLAGS) -c mod_termo.f90

mod_cinetica.o: mod_cinetica.f90 mod_termo.o
	$(FC) $(FFLAGS) -c mod_cinetica.f90

mod_malla.o: mod_malla.f90 mod_termo.o
	$(FC) $(FFLAGS) -c mod_malla.f90

mod_solver.o: mod_solver.f90 mod_malla.o mod_cinetica.o
	$(FC) $(FFLAGS) -c mod_solver.f90

limpiar:
	rm -f *.o *.mod $(EXE)

datos:
	rm -f *.dat *.png

.PHONY: all limpiar datos
