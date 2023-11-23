LIB = ../Libs/*.o

CXX = g++
CFLAGS = -Wall -O2
LIBDIR = ../Libs
VPATH = ./:$(LIBDIR)
INCLUDE_DIRS = -I. -I$(LIBDIR)
LDFLAGS = -lm

.PHONY: all
all: main

main: main.o
	$(MAKE) -C $(LIBDIR)
	$(CXX) $^ $(LDFLAGS) $(LIB) -o $@

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDE_DIRS) $^

.PHONY: clean
clean:
	@rm -f *.o
	@echo make clean: done
