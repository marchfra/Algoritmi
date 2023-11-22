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

TESTDIR = tests
TESTS = $(wildcard $(TESTDIR)/*.cpp)
TESTBINDIR = $(TESTDIR)/bin
TESTBINS = $(patsubst $(TESTDIR)/%.cpp, $(TESTDIR)/bin/%, $(TESTS))

DOCSDIR = Docs

LIBDIR = lib
LIBNAME = algo
LIB = $(LIBDIR)/lib$(LIBNAME).a


.PHONY: all
all: $(LIB)

$(LIB): $(LIBDIR) $(OBJDIR) $(OBJS)
	ar rcs $@ $(OBJS)


# Build object files from sources
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HDRDIR)/%.hpp
	$(CXX) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@


# Build test binaries from sources
$(TESTDIR)/bin/%: $(TESTDIR)/%.cpp
	$(CXX) $(CFLAGS) $< -L$(LIBDIR) -l$(LIBNAME) -o $@

test: $(LIB) $(TESTBINDIR) $(TESTBINS)
	exit 1
# 	@for test in $(TESTBINS) ; do \
# 	echo "\033[1;31;34mTesting \c" ; basename $${test%test}.hpp ; \
# 	./$${test} ; \
# 	exit 1; done

# 	@for test in $(TESTBINS) ; do \
# 	./$${test} ; \
# 	echo $$?; \
# 	done

# 	exit_codes=(15)
# 	@for test in $(TESTBINS) ; do \
# 	echo "\033[1;31;34mTesting \c" && \
# 	basename $${test%test}.hpp && \
# 	./$${test} && \
# 	exit_codes+=($$?) \
# 	echo "exit code $$?\n" ; \
# 	done
# 	echo $$exit_codes


# Create directories
$(OBJDIR):
	mkdir $@

$(LIBDIR):
	mkdir $@

$(TESTBINDIR):
	mkdir $@


# Documentation targets
.PHONY: dox
dox:
	$(MAKE) -C $(DOCSDIR)

.PHONY: cleandox
cleandox:
	$(MAKE) clean -C $(DOCSDIR)


.PHONY: clean
clean:
	$(RM) -r $(LIBDIR) $(OBJDIR) $(TESTDIR)/bin
	@echo make clean lib: done
