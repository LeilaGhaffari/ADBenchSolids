# Variables for Rust library
RUST_LIB_DIR ?= /home/leila/git/benchmark-ad4solids/src/ad-tools/enzyme-rust/target/release
RUST_LIB = $(RUST_LIB_DIR)/libenzyme_rust.so

# Variables for Enzyme, ADOL-C, and Tapenade paths
ENZYME_LIB ?=
ADOLC_INCLUDE ?=
ADOLC_LIB ?=

# Compilers
CC = clang-18
CXX = g++
FC = gfortran

# Flags
CFLAGS = $(OPT) -Wall -Wextra -Wunused-variable -Wunused-function -Iinclude
CXXFLAGS = $(OPT) -std=c++11 -Wall -Wextra -Wunused-variable -Wunused-function \
            -Wno-unused-parameter $(patsubst %,-I%,$(INCDIR) $(ADOLC_INCLUDE))
ifneq ($(ADOLC_LIB),)
    LDFLAGS += -L$(ADOLC_LIB) -Wl,-rpath,$(ADOLC_LIB)
endif
LDFLAGS += -L$(RUST_LIB_DIR) -Wl,-rpath,$(RUST_LIB_DIR)
LDLIBS = -ladolc -lm -lenzyme_rust

# Add Enzyme-specific flags if ENZYME_LIB is defined
ifneq ($(ENZYME_LIB),)
    CFLAGS += -Xclang -load -Xclang $(ENZYME_LIB)
endif

# Directories
SRCDIR = src
ADTOOLSDIR = $(SRCDIR)/ad-tools
INCDIR = include
BUILDDIR = build
BUILDTOOLSDIR = $(BUILDDIR)/ad-tools

# Source files
SOURCES_CXX = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(ADTOOLSDIR)/*.cpp)
SOURCES_C = $(wildcard $(ADTOOLSDIR)/*.c)

# Object files
OBJ_CXX = $(SOURCES_CXX:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)
OBJ_C = $(SOURCES_C:$(ADTOOLSDIR)/%.c=$(BUILDTOOLSDIR)/%.o)

# All object files
OBJ = $(OBJ_CXX) $(OBJ_C)

# Executable name
TARGET = $(BUILDDIR)/elasticity-exec

# Default target
all: $(TARGET)

# Build the Rust library
$(RUST_LIB):
	cd $(dir $@) && cargo build --release

# Link object files to create the single executable
$(TARGET): $(OBJ) $(RUST_LIB) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $@

# Compile C++ source files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDTOOLSDIR)/%.o: $(ADTOOLSDIR)/%.cpp | $(BUILDTOOLSDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile C source files
$(BUILDTOOLSDIR)/%.o: $(ADTOOLSDIR)/%.c | $(BUILDTOOLSDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Ensure necessary directories exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDTOOLSDIR):
	mkdir -p $(BUILDTOOLSDIR)

# Clean up build artifacts
clean:
	rm -f $(BUILDDIR)/*.o $(BUILDTOOLSDIR)/*.o $(TARGET)
	rm -f $(RUST_LIB_DIR)/*.so

print-% :
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [        flavor]: $(flavor $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info )
	@true

.PHONY: all clean
