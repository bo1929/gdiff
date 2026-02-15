# compiler options
#--------------------------------------------
COMPILER ?= g++
mode ?= dynamic  # Default to dynamic linking

# TODO: remove -g
CXXFLAGS += -g -std=c++17 -O3 \
						-fopenmp-simd \
						-funroll-loops \
						-ftree-vectorize \
						-flto
						# -fno-trapping-math \
						# -ffast-math
# Keep it fun?
# -Wall
WFLAGS += -Wno-unused-result -Wno-unused-command-line-argument -Wno-unknown-pragmas -Wno-undefined-inline

INC = -I simde

# project files
#--------------------------------------------
PROGRAM = gidiff
OBJECTS =	build/random.o build/enc.o \
					build/MurmurHash3.o build/lshf.o \
					build/hm.o	build/rqseq.o \
					build/sketch.o build/map.o \
					build/gidiff.o

# rules
#--------------------------------------------
.PHONY: all dynamic static clean

all:
	$(MAKE) mode=dynamic $(PROGRAM)

dynamic:
	$(MAKE) mode=dynamic $(PROGRAM)

static:
	$(MAKE) mode=static $(PROGRAM)

# Check for -lcurl
CURL_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) -lcurl -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)

OS := $(shell uname -s)

$(info ===== Build mode: $(mode) =====)
ifeq ($(mode),dynamic)
	LDLIBS = -lm -lz
else ifeq ($(mode),static)
	LDLIBS = -lm -lz
	ifneq ($(OS),Darwin)
		LDLIBS += --static -static-libgcc
	endif
	LDLIBS += -static-libstdc++ 
	CURL_SUPPORTED = no
else
	LDLIBS = -lm -lz
endif

ifneq ($(OS),Darwin)
	LDLIBS += -lstdc++ -lstdc++fs
endif

L_CURL = 0
ifneq ($(CURL_SUPPORTED),no)
	ifneq ($(mode),static)
		LDLIBS += -lcurl
		L_CURL = 1
	endif
endif

VARDEF= -D _L_CURL=$(L_CURL)

ARCH := $(shell uname -m)
# Check for -mbmi2
BMI2_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) -mbmi2 -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)
ifeq ($(filter $(ARCH),x86_64 i386),$(ARCH))
	ifneq ($(BMI2_SUPPORTED),no)
		CXXFLAGS += -mbmi2
	endif
endif

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) -c src/$*.cpp -o build/$*.o $(WFLAGS) $(CXXFLAGS) $(LDLIBS) $(VARDEF) $(INC)

$(PROGRAM): $(OBJECTS)
	$(COMPILER) -o $@ $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(VARDEF) $(LDFLAGS) $(INC)

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
