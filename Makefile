HTSDIR=../htslib
CFLAGS=-I$(HTSDIR) -Wall -g -O2
LDFLAGS=-L$(HTSDIR)
LDLIBS= #-Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -lm

all: foldreads mark5mC scanbp simhbs

foldreads: foldreads.o fold.o fit_lognorm.o kmath.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lz -lm

mark5mC: mark5mC.o fold.o aux.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lhts -lm

scanbp: scanbp.o fold.o aux.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lhts  -lm

simhbs: simhbs.o fold.o kmath.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lm

test_fit: test_fit.o fold.o kmath.o fit_lognorm.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lm

clean:
	rm -f *.o foldreads mark5mC scanbp simhbs test_fit
