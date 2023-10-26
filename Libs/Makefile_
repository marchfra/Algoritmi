DEPS = ./*.o

.PHONY: all
all: $(DEPS)
	@echo > /dev/null

CXX = g++
CFLAGS = -Wall -O2

%.o: %.cpp
	$(CXX) -std=c++17 -c $(CFLAGS) $(INCLUDE_DIRS) $^

.PHONY: clean
clean:
	@rm -f *.o
	@echo make clean: done
