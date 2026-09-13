CXX      ?= g++ # e.g., make CXX=clang++
MODE     ?= dynamic
DEBUG    ?= no
NATIVE   ?= no

PROGRAM  = gdiff
BDIR = build

OS   := $(shell uname -s)
ARCH := $(shell uname -m)

COMPILER_ID := $(shell $(CXX) -v 2>&1 | grep -qi clang && echo clang || echo gcc)

CXXFLAGS  = -std=c++17 -O3 -funroll-loops -flto
CPPFLAGS  = -I vendor -DBOOST_MATH_STANDALONE
WFLAGS    = -Wall -Wextra -Wno-unused-result -Wno-unused-parameter
LDFLAGS   =
LDLIBS    = -lm -lz -lpthread

# Compiler-specific warning suppressions
ifeq ($(COMPILER_ID),clang)
  WFLAGS += -Wno-unused-command-line-argument -Wno-undefined-inline
endif

# Test whether the compiler accepts a given flag
try-cflag = $(shell echo 'int main(){return 0;}' | $(CXX) $(1) -x c++ -o /dev/null - 2>/dev/null && echo $(1))

ifeq ($(filter $(ARCH),x86_64 i386),$(ARCH))
  # x86: enable baseline SIMD extensions
  CXXFLAGS += $(call try-cflag,-msse4.2)
  CXXFLAGS += $(call try-cflag,-mpopcnt)
  CXXFLAGS += $(call try-cflag,-mbmi2)
  # optional: AVX2 (set AVX2=yes to enable)
  ifeq ($(AVX2),yes)
    CXXFLAGS += $(call try-cflag,-mavx2)
    CXXFLAGS += $(call try-cflag,-mfma)
  endif
  # optional: AVX-512 (set AVX512=yes to enable)
  ifeq ($(AVX512),yes)
    CXXFLAGS += $(call try-cflag,-mavx512f)
  endif
else ifeq ($(ARCH),arm64)
  # Apple Silicon / aarch64 — NEON is on by default; nothing extra needed.
else ifeq ($(ARCH),aarch64)
  # Linux aarch64
endif

# -march=native: optional for local-only builds, not portable binaries
ifeq ($(NATIVE),yes)
  CXXFLAGS += -march=native
endif

# OpenMP SIMD works for both GCC and Clang
CXXFLAGS += $(call try-cflag,-fopenmp-simd)

# Debug mode
ifeq ($(DEBUG),yes)
  CXXFLAGS := $(filter-out -O3 -flto,$(CXXFLAGS))
  # CXXFLAGS += -pg
  CXXFLAGS += -O0 -g -ggdb3 -fno-omit-frame-pointer
  # Address sanitizer (incompatible with LTO)
  # CXXFLAGS += $(call try-cflag,-fsanitize=address)
  # LDFLAGS  += $(call try-cflag,-fsanitize=address)
endif

# Linking mode
ifeq ($(MODE),static)
  ifneq ($(OS),Darwin)
    LDFLAGS += -static -static-libgcc
  endif
  LDFLAGS += -static-libstdc++
endif

ifneq ($(OS),Darwin)
  ifeq ($(COMPILER_ID),gcc)
    LDLIBS += -lstdc++
    LDLIBS += $(shell echo 'int main(){return 0;}' | $(CXX) -lstdc++fs -x c++ -o /dev/null - 2>/dev/null && echo -lstdc++fs)
  endif
endif

# curl (optional, dynamic only)
LCURL := 0
ifneq ($(MODE),static)
  CURL_OK := $(shell echo 'int main(){return 0;}' | $(CXX) -lcurl -x c++ -o /dev/null - 2>/dev/null && echo yes)
  ifeq ($(CURL_OK),yes)
    LDLIBS += -lcurl
    LCURL  := 1
  endif
endif
CPPFLAGS += -D_LCURL=$(LCURL)

SOURCES  = $(wildcard src/*.cpp)
OBJECTS  = $(patsubst src/%.cpp,$(BDIR)/%.o,$(SOURCES))
DEPENDS  = $(OBJECTS:.o=.d)

$(info --- $(PROGRAM): $(OS)/$(ARCH), $(COMPILER_ID), mode=$(MODE), debug=$(DEBUG), native=$(NATIVE) ---)

# Test configuration
TDIR     = test/unit
TBDIR    = $(BDIR)/test
TSOURCES = $(wildcard $(TDIR)/*.cpp)
TOBJECTS = $(patsubst $(TDIR)/%.cpp,$(TBDIR)/%.o,$(TSOURCES))
TDEPENDS = $(TOBJECTS:.o=.d)
# All project .o files except gdiff.o (which contains main())
LIB_OBJECTS = $(filter-out $(BDIR)/gdiff.o,$(OBJECTS))
TEST_BIN = $(BDIR)/test_gdiff

# Rules
.PHONY: all dynamic static debug clang tidy tidy-fix clean test-unit test-regression test

all: $(PROGRAM)

dynamic:
	$(MAKE) MODE=dynamic $(PROGRAM)

static:
	$(MAKE) MODE=static $(PROGRAM)

debug:
	$(MAKE) DEBUG=yes $(PROGRAM)

clang:
	$(MAKE) CXX=clang++ $(PROGRAM)

# Clang-tidy configuration
CLANG_TIDY     ?= clang-tidy
CLANG_TIDY_OUT = $(BDIR)/clang-tidy.out
CLANG_TIDY_FLAGS = --

# macOS system include detection
ifeq ($(OS),Darwin)
  SYSROOT := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
  ifneq ($(SYSROOT),)
    CLANG_TIDY_FLAGS += -isystem $(SYSROOT)/usr/include/c++/v1
    CLANG_TIDY_FLAGS += -isystem $(SYSROOT)/usr/include
  endif
endif

CLANG_TIDY_FLAGS += -std=c++17 -I vendor -DBOOST_MATH_STANDALONE -pthread

# Run clang-tidy and save output
tidy: $(SOURCES) .clang-tidy | $(BDIR)
	@echo "Running clang-tidy..."
	$(CLANG_TIDY) $(SOURCES) $(CLANG_TIDY_FLAGS) 2>&1 | tee $(CLANG_TIDY_OUT)
	@echo "Output saved to $(CLANG_TIDY_OUT)"

# Auto-fix clang-tidy issues
tidy-fix: $(SOURCES) .clang-tidy | $(BDIR)
	@echo "Running clang-tidy with auto-fix..."
	$(CLANG_TIDY) $(SOURCES) -fix $(CLANG_TIDY_FLAGS) 2>&1 | tee $(CLANG_TIDY_OUT)
	@echo "Fixes applied. Review changes with git diff."

$(BDIR)/%.o: src/%.cpp | $(BDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(WFLAGS) -MMD -MP -c $< -o $@

$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(BDIR):
	@mkdir -p $@

clean:
	rm -rf $(BDIR) $(PROGRAM)
	@echo "Clean."

-include $(DEPENDS)

# ── Test targets ──────────────────────────────────────────────────────────────

$(TBDIR):
	@mkdir -p $@

$(TBDIR)/%.o: $(TDIR)/%.cpp | $(TBDIR)
	$(CXX) $(CPPFLAGS) -I src $(CXXFLAGS) $(WFLAGS) -MMD -MP -c $< -o $@

$(TEST_BIN): $(TOBJECTS) $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

test-unit: $(TEST_BIN)
	@echo "Running unit tests..."
	./$(TEST_BIN) --duration

test-regression: $(PROGRAM)
	@echo "Running regression tests..."
	bash test/test_regression.sh

test: test-unit test-regression

-include $(TDEPENDS)