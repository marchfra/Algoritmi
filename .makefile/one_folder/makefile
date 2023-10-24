CXX = g++
OBJ = main.o swap.o

main: $(OBJ)
	$(CXX) $^ -lm -o $@

%.o: %.cpp
	$(CXX) -c -I. $<
