# Library makefile

# Compiler and flags
CC = g++
CFLAGS = -Wall -O2

# Source and object files
SRC = test.cpp
OBJ = $(SRC:.cpp = .o)

# Library name and output directory
LIBRARY = libmylibrary.a
LIB_DIR = ../../lib

# Targets
all : $(LIBRARY)

$(LIBRARY): $(OBJ)
	ar rcs $(LIB_DIR)/$@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

	clean:
		rm -f $(OBJ) $(LIB_DIR)/$(LIBRARY)
