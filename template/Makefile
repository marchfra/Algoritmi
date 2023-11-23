# Compiler stuff
CXX = g++
CFLAGS = -g -Wall -std=c++17
PYTHON = python


# Get all necessary files
SRCDIR = src
SRCS = $(wildcard $(SRCDIR)/*.cpp)

OBJDIR = obj
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

BINDIR = bin
BIN = $(BINDIR)/main

LIBDIR = /Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs
LIBNAME = algo


.PHONY: all
all: $(BIN)

# Build binaries
$(BIN): $(OBJDIR) $(OBJS) $(BINDIR)
	$(MAKE) -C $(LIBDIR)
	$(CXX) $(OBJS) -L$(LIBDIR)/lib -l$(LIBNAME) -o $@


# Build objects from sources
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@


.PHONY: plot
plot:
	$(PYTHON) plot.py


# Create directories
$(OBJDIR):
	mkdir $@

$(BINDIR):
	mkdir $@


.PHONY: lib
lib:
	$(MAKE) -C $(LIBDIR)


.PHONY: test
test:
	$(MAKE) test -C $(LIBDIR)


.PHONY: clean
clean:
	$(RM) -r $(OBJDIR) $(BINDIR) $(SRCDIR)/marchisotti_$(CLASSWORK).cpp
	@echo make clean: done


.PHONY: cleanlib
cleanlib:
	@$(MAKE) clean -C $(LIBDIR)
