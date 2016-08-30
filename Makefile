TARGET=foldreads
OBJS=foldreads.o
CFLAGS=-Wall -g -O2 -funroll-loops
LDFLAGS=
LIBS=-Wl,-Bstatic -lz -Wl,-Bdynamic

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $< -o $(TARGET) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
