HTSDIR=../htslib
CFLAGS=-I$(HTSDIR) -Wall -g -O2
LDFLAGS=-L$(HTSDIR)
#HTS_LDLIBS=-Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread
HTS_LDLIBS=-lhts


all: foldreads mark5mC scanbp simhbs qualprofile

foldreads: foldreads.o fold.o fit_lognorm.o kmath.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lz -lm

mark5mC: mark5mC.o fold.o aux.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(HTS_LDLIBS) -lm

scanbp: scanbp.o fold.o aux.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(HTS_LDLIBS) -lm

simhbs: simhbs.o fold.o kmath.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lm

test_fit: test_fit.o fold.o kmath.o fit_lognorm.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lm

qualprofile: qualprofile.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lz

clean:
	rm -f *.o foldreads mark5mC scanbp simhbs test_fit qualprofile
