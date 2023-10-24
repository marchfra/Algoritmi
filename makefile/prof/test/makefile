OBJ = main.o swap.o

CXX = g++
CFLAGS = -c
SRC = ../Libs
VPATH = ./:$(SRC)
INCLUDE_DIRS = -I. -I$(SRC)
LDFLAGS = -lm

main: $(OBJ)
	$(CXX) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $^

.PHONY: clean

clean:
	@rm -f *.o
	@echo make clean: done


# g++ -c -I. -I../Libs main.cpp
# make: *** No rule to make target `swap.o', needed by `main'.  Stop.
