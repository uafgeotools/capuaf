FFLAGS = -O
FC = gfortran
#CFLAGS = ${FFLAGS} -I$(NR)/NRC_includes
#NRLIB = $(MY_BIN)
NRLIB = $(NRHOME_C)
CFLAGS = ${FFLAGS} -I$(NRLIB)/includes

CAP  = cap cap_dir

#SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o sub_tt2cmt.o sub_first_motion_misfit.o sub_fmp_print_params.o calerr.o caperror.o sub_initSearch.o uv2lune.o
SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o sub_tt2cmt.o sub_first_motion_misfit.o sub_fmp_print_params.o calerr.o sub_initSearch.o uv2lune.o

all: $(CAP)

# Write binary output. Uncomment and then compile with the usual commands.
# There is a better way to do this... (ongoing)
#cap cap_dir: %:%.o $(SUBS) cap_sub.o
#	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio -L$(NRLIB) -lcrecipes -DWRITECAPBIN

# (Default cap) do not write binary output
cap cap_dir: %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio -L$(NRLIB) -lcrecipes


#cap cap_dir: %:%.o $(SUBS) cap_sub.o
#	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio -L$(NRLIB) -lnr

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

clean:
	rm -f $(CAP) *.o
