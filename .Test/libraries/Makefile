CXX = g++
CFLAGS = -g -Wall
BINS = librarytest libmycode.dylib runtime_librarytest libmystaticcode.a static_librarytest

all: $(BINS)

libmycode.o: libmycode.cpp mycode.hpp
	$(CXX) $(CFLAGS) -c $<

libmycode.dylib: libmycode.cpp mycode.hpp
	$(CXX) $(CFLAGS) -fPIC -shared -o $@ $< -lc

librarytest: librarytest.cpp libmycode.o
	$(CXX) $(CFLAGS) -o $@ $^

runtime_librarytest: librarytest.cpp
	# -L.  look in the . directory for libraries
	# -lmycode  link the program with libmycode
	$(CXX) $(CFLAGS) -o $@ $^ -L. -lmycode

libmystaticcode.a: libmycode.o
	ar rcs $@ $^

static_librarytest: librarytest.cpp
	# -L.  look in the . directory for libraries
	# -lmycode  link the program with libmycode
	$(CXX) $(CFLAGS) -o $@ $^ -L. -lmystaticcode


clean:
	rm *o $(BINS)
