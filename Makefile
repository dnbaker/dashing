.PHONY=all tests clean obj
CXX?=g++
CC?=gcc

GMATCH=$(findstring g++,$(CXX))

CLHASH_CHECKOUT = "&& git checkout master"
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -Wunused-variable -Wno-attributes -Wno-pedantic \
		# -Wduplicated-branches -Wdangling-else  # -Wsuggest-attribute=malloc   # -Wconversion
EXTRA?=
INCPLUS?=
EXTRA_LD?=
DBG:=
OS:=$(shell uname)
FLAGS=

OPT_MINUS_OPENMP= -O3 -funroll-loops\
	  -pipe -fno-strict-aliasing -march=native -mpclmul -DUSE_PDQSORT \
	-DNOT_THREADSAFE -DENABLE_COMPUTED_GOTO \
	$(FLAGS) $(EXTRA)
OPT=$(OPT_MINUS_OPENMP) -fopenmp
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++14 $(WARNINGS)
CXXFLAGS_MINUS_OPENMP=$(OPT_MINUS_OPENMP) $(XXFLAGS) -std=c++1z $(WARNINGS) -Wno-cast-align -Wno-gnu-zero-variadic-macro-arguments
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS)
LIB=-lz
LD=-L. $(EXTRA_LD) -Lbonsai/zlib

ifneq (,$(findstring g++,$(CXX)))
	ifeq ($(shell uname),Darwin)
		ifeq (,$(findstring clang,$(CXX)))
			POPCNT_CXX:=clang
		else
			POPCNT_CXX:=$(CXX)
		endif
	else
		POPCNT_CXX:=$(CXX)
	endif
endif

OBJS=$(patsubst %.c,%.o,$(wildcard src/*.c) bonsai/klib/kthread.o) $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) bonsai/klib/kstring.o clhash.o
DOBJS=$(patsubst %.c,%.do,$(wildcard src/*.c) bonsai/klib/kthread.o) $(patsubst %.cpp,%.do,$(wildcard src/*.cpp)) bonsai/klib/kstring.o clhash.o

ZSTD_INCLUDE_DIRS=bonsai/zstd/zlibWrapper bonsai/zstd/lib/common bonsai/zstd/lib
ZSTD_INCLUDE=$(patsubst %,-I%,$(ZSTD_INCLUDE_DIRS))
ZFLAGS=-DZWRAP_USE_ZSTD=1
ZCOMPILE_FLAGS= $(ZFLAGS) -lzstd
ZW_OBJS=$(patsubst %.c,%.o,$(wildcard bonsai/zstd/zlibWrapper/*.c)) libzstd.a
ALL_ZOBJS=$(ZOBJS) $(ZW_OBJS) bonsai/bonsai/clhash.o bonsai/klib/kthread.o
INCLUDE=-Ibonsai/clhash/include -I.  -Ibonsai/zlib -Ibonsai/libpopcnt -Iinclude -Ibonsai/circularqueue $(ZSTD_INCLUDE) $(INCPLUS) -Ibonsai/hll -Ibonsai/hll/vec -Ibonsai -Ibonsai/bonsai/include/

EX=$(patsubst src/%.cpp,%,$(wildcard src/*.cpp))
D_EX=$(patsubst src/%.cpp,%_d,$(wildcard src/*.cpp))


all: update dashing

d: $(D_EX)

update:
	+git submodule update --init --remote --recursive . && cd bonsai && git checkout master && git pull && make update && \
    cd linear && git checkout master && git pull && cd .. && cd .. && cd distmat && git checkout master && git pull && cd ..

libzstd.a:
	+cd bonsai/bonsai && make libzstd.a && cp libzstd.a ../../

bonsai/klib/kstring.o:
	+cd bonsai/bonsai && make ../klib/kstring.o && cd ../.. && \
	cd bonsai/bonsai && make ../klib/kthread.o && cd ../..

bonsai/bonsai/clhash.o:
	+cd bonsai/bonsai && make clhash.o && cd ../..

OBJ=bonsai/klib/kstring.o bonsai/klib/kthread.o bonsai/bonsai/clhash.o

DEPS=bonsai/hll/cbf.h bonsai/hll/bf.h bonsai/hll/hll.h

test/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

test/%.zo: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(ZCOMPILE_FLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

%.do: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

bonsai/zlib/libz.a:
	cd bonsai/zlib && ./configure && make libz.a

zobj: $(ALL_ZOBJS)

STATIC_GOMP?=$(shell $(CXX) --print-file-name=libgomp.a)

%: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_d: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -fsanitize=address -fsanitize=undefined -fsanitize=leak $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_256: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mno-avx512dq -mno-avx512vl -mno-avx512bw -mavx2 -msse4.1 -msse2 -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_512: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mavx512dq -mavx512vl -mavx2 -msse4.1 -msse2 -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_512bw: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mavx512dq -mavx512vl -mavx512bw -mavx2 -msse4.1 -msse2 -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_128: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mno-avx512dq -mno-avx512vl -mno-avx512bw -mno-avx2 -msse2 -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
%_s: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	ln -sf $(STATIC_GOMP) && \
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -static-libstdc++ -static-libgcc bonsai/zlib/libz.a  -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
%_s128: src/%.cpp $(ALL_ZOBJS) $(DEPS) bonsai/zlib/libz.a
	ln -sf $(STATIC_GOMP) && \
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mno-avx2 -msse2 -static-libstdc++ -static-libgcc bonsai/zlib/libz.a  -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
%_s512: src/%.cpp $(ALL_ZOBJS) $(DEPS) bonsai/zlib/libz.a
	ln -sf $(STATIC_GOMP) && \
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mavx512dq -mavx512vl -mavx512bw -static-libstdc++ -static-libgcc bonsai/zlib/libz.a  -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
%_di: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -fno-inline -O1 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
%_pg: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -pg -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
%_pgi: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -pg -fno-inline -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
clean:
	rm -f $(EX) $(D_EX) libzstd.a bonsai/bonsai/clhash.o clhash.o bonsai/klib/kthread.o bonsai/klib/kstring.o && cd bonsai/zstd && make clean && cd ../zlib && make clean && cd ../..
mostlyclean: clean
