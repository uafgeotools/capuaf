# make clean cap d=WB f=omp
# 
#FFLAGS = -O -no-pie    # needed for gcc version 7.3.0 (Ubuntu 7.3.0-27ubuntu1~18.04)
FFLAGS = -O	
FC = gfortran
CFLAGS = ${FFLAGS}

ifeq (${d}, WB)
	CFLAGS += -D${d}
endif

ifeq (${f}, omp)
	FFLAGS += -fopenmp
	CFLAGS += -DOMP
endif

CAP  = cap cap_dir

SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o sub_tt2cmt.o sub_first_motion_misfit.o sub_fmp_print_params.o sub_misfit.o sub_inversion.o sub_uv2lune.o

all: $(CAP)

# to compile cap with the option of Writing Binary file, run
# 	make cap d=WB
# to compile cap so that it runs in parallel mode, run
# 	make cap f=omp
# for both options, run
#       make cap d=WB f=omp
cap cap_dir: %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

clean:
	rm -f $(CAP) *.o
