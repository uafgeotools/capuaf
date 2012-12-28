FFLAGS = -O
CFLAGS = ${FFLAGS} -I$(NR)/NRC_includes
NRLIB = $(MY_BIN)

CAP  = cap cap_dir

SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o

all: $(CAP)

cap cap_dir: %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio -L$(NRLIB) -lnr

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

clean:
	rm -f $(CAP) *.o
