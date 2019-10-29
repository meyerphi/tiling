all: bin/tiling

LIBDIR=lib
BINDIR=bin
BIN=tiling
CC=g++
CCOPTS=-O2 -Wall -std=c++14 -o $@ -I$(LIBDIR)

SRC=src/tiling.cpp src/transducer.hpp src/tiles.hpp src/visualize.hpp lib/simple_svg.hpp

$(BINDIR):
	mkdir -p $@

$(BINDIR)/$(BIN): $(SRC) | $(BINDIR)
	$(CC) $(CCOPTS) -o $@ $<

clean:
	rm -rf $(BINDIR)

.PHONY: all clean
