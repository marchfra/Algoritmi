DEPS = ./obj/*.o

.PHONY: all
all: $(DEPS)
	@echo > /dev/null

CXX = g++
CFLAGS = -Wall -O2

./obj/%.o: ./src/%.cpp
	$(CXX) -std=c++17 -c $(CFLAGS) $(INCLUDE_DIRS) $^
	@mv *.o obj

.PHONY: clean
clean:
	@rm -f ./obj/*.o
	@echo make clean lib: done
