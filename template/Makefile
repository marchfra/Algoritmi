# Compiler stuff
CXX = g++
CFLAGS = -g -Wall -std=c++17
PYTHON = python


# Get all necessary files
SRCDIR = src
SRCS = $(wildcard $(SRCDIR)/*.cpp)

OBJDIR = obj
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

BINDIR = bin
BIN = $(BINDIR)/main

LIBDIR = /Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs
LIBNAME = algo


.PHONY: all
all: $(BIN)

# Build binaries
$(BIN): $(OBJDIR) $(OBJS) $(BINDIR)
	$(MAKE) -C $(LIBDIR)
	$(CXX) $(CFLAGS) $(OBJS) -L$(LIBDIR)/lib -l$(LIBNAME) -o $@


# Build objects from sources
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@


.PHONY: plot
plot:
	@echo "Making plot"
	@$(PYTHON) $(SRCDIR)/plot.py


# Create directories
$(OBJDIR):
	mkdir $@

$(BINDIR):
	mkdir $@


.PHONY: lib
lib:
	$(MAKE) -C $(LIBDIR)


.PHONY: test
test:
	$(MAKE) test -C $(LIBDIR)


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
	$(RM) -r $(OBJDIR) $(BINDIR)
	@echo make clean: done


.PHONY: cleanlib
cleanlib:
	@$(MAKE) clean -C $(LIBDIR)
