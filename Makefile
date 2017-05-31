HTSDIR=../../../htslib
CFLAGS=-I$(HTSDIR) -Wall -g -O2
LDFLAGS=-L$(HTSDIR)
LIBS= #-Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -lm

all: foldreads2 mark_5mC scan_pairs

foldreads2: foldreads2.o fold.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lz -Wl,-Bdynamic -lm

foldreads: foldreads.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lz -Wl,-Bdynamic -lm

mark_5mC: mark_5mC.o fold.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread -lm

scan_pairs: scan_pairs.o fold.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread -lm

clean:
	rm -f *.o foldreads foldreads2 mark_5mC scan_pairs
