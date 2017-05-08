CXX = g++

# change compiler to g++-6 for macOS since the default g++ in
# macOS is a symlink to clang which doesn't support OpenMP
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	CXX = g++-6
endif

CXXFLAGS = -g -Wall -O2 -std=c++14 -fopenmp

INCLUDES = -I./include

SRC1 = demo/double_gyre/continuous_double_gyre.cpp
SRC2 = demo/double_gyre/discrete_double_gyre.cpp

BINDIR = ./bin

TARGET1 = $(BINDIR)/demo_continuous_double_gyre
TARGET2 = $(BINDIR)/demo_discrete_double_gyre

all: $(TARGET1) $(TARGET2)

$(TARGET1): | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC1) -o $(TARGET1)

$(TARGET2): | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC2) -o $(TARGET2)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(TARGET1) $(TARGET2)
