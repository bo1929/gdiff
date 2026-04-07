CXX      ?= g++ # e.g., make CXX=clang++
MODE     ?= dynamic
DEBUG    ?= no
NATIVE   ?= no

PROGRAM  = gidiff
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

# Rules
.PHONY: all dynamic static debug clean

all: $(PROGRAM)

dynamic:
	$(MAKE) MODE=dynamic $(PROGRAM)

static:
	$(MAKE) MODE=static $(PROGRAM)

debug:
	$(MAKE) DEBUG=yes $(PROGRAM)

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
