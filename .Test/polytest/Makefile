CXX = g++
CFLAGS = -Wall -O2
VPATH = ./:$(LIBDIR)
LIBDIR = /Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs
LIB = $(LIBDIR)/*.o
INCLUDE_DIRS = -I. -I$(LIBDIR)
LDFLAGS = -lm

.PHONY: all
all: main

main: main.o
	$(MAKE) -C $(LIBDIR)
	$(CXX) -std=c++17 $^ $(LDFLAGS) $(LIB) -o $@
	@rm -f main.o

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDE_DIRS) $^

.PHONY: clean
clean:
	@rm -f *.o
	@echo make clean: done
