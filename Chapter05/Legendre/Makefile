CXX = g++
CFLAGS = -Wall -O2
VPATH = ./:$(LIBDIR)
LIBDIR = /Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/
INCDIR = /Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include
OBJDIR = /Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/obj
LIB = $(OBJDIR)/*.o
INCLUDE_DIRS = -I. -I$(INCDIR)
LDFLAGS = -lm

.PHONY: all
all: main

main: main.o
	$(MAKE) -C $(LIBDIR)
	$(CXX) -std=c++23 $^ $(LDFLAGS) $(LIB) -o $@
	@rm -f main.o

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDE_DIRS) $^

.PHONY: clean
clean:
	@rm -f *.o
	@echo make clean: done

.PHONY: cleanlib
cleanlib:
	@$(MAKE) clean -C $(LIBDIR)
