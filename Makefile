# Variables for Rust library
RUST_PROFILE ?= release # or debug
RUST_LIB_DIR ?= $(abspath target/$(RUST_PROFILE))
RUST_LIB = $(RUST_LIB_DIR)/libenzyme_rust.a

# Variables for Enzyme, ADOL-C, and Tapenade paths
ENZYME_LIB ?=
ADOLC_INCLUDE ?=
ADOLC_LIB ?=

RUSTC = $(realpath $(shell rustup +enzyme which rustc))
RUSTC_DIR = $(abspath $(dir $(RUSTC))../..)
RUSTC_LLVM = $(RUSTC_DIR)/llvm/bin
ENZYME_LIB = $(wildcard $(RUSTC_DIR)/enzyme/lib/ClangEnzyme-*.so)

# Compilers
CC = $(RUSTC_LLVM)/clang
CXX = g++
FC = gfortran

# Flags
CFLAGS = $(OPT) -Wall -Wextra -Wunused-variable -Wunused-function -Iinclude
CXXFLAGS = $(OPT) -std=c++11 -Wall -Wextra -Wunused-variable -Wunused-function \
            -Wno-unused-parameter $(patsubst %,-I%,$(INCDIR) $(ADOLC_INCLUDE))
# OpenMP Flags
CFLAGS += -fopenmp-simd
CXXFLAGS += -fopenmp
LDFLAGS += -fopenmp
LDLIBS += -lomp

ifneq ($(ADOLC_LIB),)
    LDFLAGS += -L$(ADOLC_LIB) -Wl,-rpath,$(ADOLC_LIB)
endif
LDFLAGS += -L$(RUST_LIB_DIR)
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
$(RUST_LIB): Cargo.toml src/ad-tools/enzyme-rust/Cargo.toml $(wildcard src/ad-tools/enzyme-rust/src/*.rs)
	RUSTFLAGS='-C target-cpu=native -Z autodiff=Enable' cargo +enzyme build --profile $(RUST_PROFILE)

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
	cargo +enzyme clean

print-% :
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [        flavor]: $(flavor $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info )
	@true

.PHONY: all clean
