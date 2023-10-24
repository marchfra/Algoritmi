# Library makefile

# Specify extensions of files to delete when cleaning
CLEANEXT = o a

# Specify the output file and install directory
OUTPUTFILE = libmylibrary.a
INSTALLDIR = ../../lib

# Source and object files
SRC = my_library.cpp
OBJ = $(SRC:.cpp = .o)

# Default target
.PHONY: all
all: $(OUTPUTFILE)

# Build libmylibrary.a from object dependencies
$(OUTPUTFILE): $(OBJ)
# 	ar rcs $(INSTALLDIR)/$@ $^
	ar ru $@ $^
	ranlib $@

# No rule to build .o from .cpp files is required; this is
# handled by make's database of implicit rules

# # Compiler and flags
# CC = g++
# CFLAGS = -Wall -O2

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: install
install: 
	mkdir -p $(INSTALLDIR)
	cp -p $(OUTPUTFILE) $(INSTALLDIR)

.PHONY: clean
clean:
	@rm -f $(INSTALLDIR)/$(OUTPUTFILE)
	@echo make clean: done

# Indicate dependencies of .ccp files on .hpp files
my_library.o: my_library.cpp
