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

.PHONY: all
all: $(BIN)

# Build release version without debug stuff
.PHONY: release
release: CFLAGS = -Wall -O2 -DNDEBUG -std=c++17
release: clean
release: $(BIN)

# Build binaries
$(BIN): $(OBJS) $(OBJDIR) $(BINDIR)
	$(CXX) $^ -o $@

# Build objects from sources
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

# Clean
.PHONY: clean
clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*

# Run src/plot.py
.PHONY: plot
plot:
	$(PYTHON) $(SRCDIR)/plot.py


# Create directories if they don't already exist
$(OBJDIR):
	mkdir $@

$(BINDIR):
	mkdir $@
