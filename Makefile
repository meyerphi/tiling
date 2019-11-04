LIBDIR=lib
BINDIR=bin
SRCDIR=src
BIN=tiling
CC=g++
CCOPTS=-O2 -Wall -std=c++14 -I$(LIBDIR)

SRC=$(SRCDIR)/tiling.cpp $(SRCDIR)/transducer.hpp $(SRCDIR)/tiles.hpp $(SRCDIR)/visualize.hpp $(LIBDIR)/simple_svg.hpp

all: $(BINDIR)/$(BIN)

$(BINDIR):
	mkdir -p $@

$(BINDIR)/$(BIN): $(SRC) | $(BINDIR)
	$(CC) $(CCOPTS) -o $@ $<

clean:
	rm -rf $(BINDIR)

.PHONY: all clean
