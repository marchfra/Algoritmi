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

$(LIB): $(LIBDIR) $(OBJDIR) $(OBJS) $(HDRS)
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
	@for test in $(TESTBINS) ; do \
	echo "\n\033[1;34mTesting \c" ; basename $${test%test}.hpp ; \
	./$${test} ; \
	done


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


# clang-format: https://gist.github.com/dtoma/61468552bbc7c0114b2e700f9247a310
.PHONY: style
style:
	@for src in $(SRCS) ; do \
		echo "Formatting $${src}..." ; \
		clang-format -i "$${src}" ; \
	done
	@echo "Done"

.PHONY: check-style
check-style:
	@for src in $(SRCS) ; do \
		var=`clang-format "$${src}" | diff "$${src}" - | wc -l` ; \
		if [ $${var} -ne 0 ] ; then \
			echo "$${src} does not respect the coding style (diff: $${var} lines)" ; \
			exit 1 ; \
		fi ; \
	done
	@echo "Style check passed"


.PHONY: clean
clean:
	$(RM) -r $(LIBDIR) $(OBJDIR) $(TESTDIR)/bin
	@echo make clean lib: done
