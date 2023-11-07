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
all: work_grav

work_grav: work_grav.o
	$(MAKE) -C $(LIBDIR)
	$(CXX) -std=c++17 $^ $(LDFLAGS) $(LIB) -o $@
	@rm -f work_grav.o

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDE_DIRS) $^

.PHONY: clean
clean:
	@rm -f *.o
	@echo make clean: done

.PHONY: cleanlib
cleanlib:
	@$(MAKE) clean -C $(LIBDIR)
