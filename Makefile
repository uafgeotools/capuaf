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

#SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o sub_tt2cmt.o sub_first_motion_misfit.o sub_fmp_print_params.o calerr.o caperror.o sub_initSearch.o uv2lune.o
SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o sub_tt2cmt.o sub_first_motion_misfit.o sub_fmp_print_params.o calerr.o sub_inversion.o uv2lune.o

all: $(CAP)

# to compile cap with the option of Writing Binary file, run
# 	make cap d=WB
cap cap_dir: %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

clean:
	rm -f $(CAP) *.o
