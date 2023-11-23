# Compiler stuff
CXX = g++
CFLAGS = -g -Wall -std=c++17

# Get all necessary files
SRCDIR = src
SRCS = $(wildcard $(SRCDIR)/*.cpp)

OBJDIR = obj
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

HDRDIR = include
HDRS = $(wildcard $(HDRDIR)/*.hpp)

LIBDIR = lib
LIB = $(LIBDIR)/libmycode.a

TESTDIR = tests
TESTS = $(wildcard $(TESTDIR)/*.cpp)
TESTBINDIR = $(TESTDIR)/bin
TESTBINS = $(patsubst $(TESTDIR)/%.cpp, $(TESTDIR)/bin/%, $(TESTS))

all: $(LIB)

$(LIB): $(OBJDIR) $(OBJS) $(LIBDIR)
	ar rcs $@ $(OBJS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HDRDIR)/%.hpp
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

$(TESTDIR)/bin/%: $(TESTDIR)/%.cpp
	$(CXX) $(CFLAGS) $< -L$(LIBDIR) -lmycode -o $@

test: $(LIB) $(TESTBINDIR) $(TESTBINS)
	@for test in $(TESTBINS) ; do ./$${test} && echo "" ; done


# Create directories
$(OBJDIR):
	mkdir $@

$(LIBDIR):
	mkdir $@

$(TESTBINDIR):
	mkdir $@


clean:
	$(RM) -r $(LIBDIR) $(OBJDIR) $(TESTDIR)/bin
