# Main project makefile

# Compiler and flags
CC = g++
CFLAGS = -Wall -O2

# Main program source file
MAIN_SRC = src/main.cpp

# Library name and output directory
LIBRARY = my_library
LIB_DIR = lib/my_library

# Targets
all: $(LIBRARY) main

$(LIBRARY): 
	$(MAKE) -C $(LIB_DIR)

main: $(LIBRARY)
	$(CC) $(CFLAGS) $(MAIN_SRC) -o $@ -L$(LIB_DIR) -lmylibrary

clean:
	$(MAKE) -C $(LIB_DIR) clean
	rm -f $(LIBRARY)
