HTSDIR=../htslib
CFLAGS=-I$(HTSDIR) -Wall -g -O2
LDFLAGS=-L$(HTSDIR)
LIBS= #-Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -lm

all: foldreads mark5mC scanbp

foldreads: foldreads.o fold.o fit_lognorm.o kmath.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lz -Wl,-Bdynamic -lm

mark5mC: mark5mC.o fold.o aux.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread -lm

scanbp: scanbp.o fold.o aux.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread -lm

clean:
	rm -f *.o foldreads mark5mC scanbp
