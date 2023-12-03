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


# clang-format
.PHONY: style
style:
	@for src in $(SRCS) ; do \
		echo "Formatting $${src}..." ; \
		clang-format -i "$${src}" ; \
	done
	@echo "Done"
# 		clang-tidy -checks='-*,readability-identifier-naming' \
# 		    -config="{CheckOptions: [ \
# 		    { key: readability-identifier-naming.NamespaceCase, value: lower_case },\
# 		    { key: readability-identifier-naming.ClassCase, value: CamelCase  },\
# 		    { key: readability-identifier-naming.StructCase, value: CamelCase  },\
# 		    { key: readability-identifier-naming.FunctionCase, value: camelBack },\
# 		    { key: readability-identifier-naming.VariableCase, value: lower_case },\
# 		    { key: readability-identifier-naming.GlobalConstantCase, value: UPPER_CASE }\
# 		    ]}" "$${src}" ; \

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
