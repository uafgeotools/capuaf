FFLAGS = -O
FC = gfortran
#CFLAGS = ${FFLAGS} -I$(NR)/NRC_includes
#NRLIB = $(MY_BIN)
NRLIB = $(NRHOME_C)
CFLAGS = ${FFLAGS} -I$(NRLIB)/includes

CAP  = cap cap_dir

SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o sub_tt2cmt.o

all: $(CAP)

cap cap_dir: %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio -L$(NRLIB) -lcrecipes
#cap cap_dir: %:%.o $(SUBS) cap_sub.o
#	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio -L$(NRLIB) -lnr

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

clean:
	rm -f $(CAP) *.o
