CC      = gcc
CFLAGS  = -Wall -Wextra -pedantic -std=c99
LDLIBS  = -lm

.PHONY: all debug valgrind clean

# Vanlig build
all: sunscreen

sunscreen: main.c
	$(CC) $(CFLAGS) main.c $(LDLIBS) -o $@

# Debug-build: egen binær med -g -O0
sunscreen_debug: main.c
	$(CC) $(CFLAGS) -g -O0 main.c $(LDLIBS) -o $@

debug: sunscreen_debug

valgrind: sunscreen_debug
	valgrind --leak-check=yes ./sunscreen_debug

clean:
	rm -f sunscreen sunscreen_debug
