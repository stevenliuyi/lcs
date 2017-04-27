CXX = g++

# change compiler to g++-6 for macOS since the default g++ in
# macOS is a symlink to clang which doesn't support OpenMP
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	CXX = g++-6
endif

CXXFLAGS = -g -Wall -O2 -std=c++14 -fopenmp

INCLUDES = -I./include

SRCS = demo/double_gyre/double_gyre.cpp

BINDIR = ./bin

TARGET = $(BINDIR)/demo_double_gyre

all: $(TARGET)

$(TARGET): | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS) -o $(TARGET)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(TARGET)
