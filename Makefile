HTSDIR=../../../htslib
CFLAGS=-I$(HTSDIR) -Wall -g -O2
LDFLAGS=-L$(HTSDIR)
LIBS= #-Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -lm

all: foldreads mark_5mC scan_pairs

foldreads: foldreads.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -Wl,-Bstatic -lz -Wl,-Bdynamic -lm

mark_5mC: mark_5mC.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread

scan_pairs: scan_pairs.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -Wl,-Bstatic -lhts -lz -Wl,-Bdynamic -pthread

clean:
	rm -f *.o foldreads mark_5mC scan_pairs
