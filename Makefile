.PHONY=all tests clean obj
CXX?=g++
CC?=gcc

GMATCH=$(findstring g++,$(CXX))

WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -Wunused-variable -Wno-attributes -Wno-pedantic  -Wno-ignored-attributes \
        -Wno-missing-braces -Wno-unknown-pragmas
		# -Wduplicated-branches -Wdangling-else  # -Wsuggest-attribute=malloc   # -Wconversion
EXTRA?=
INCPLUS?=
EXTRA_LD?=
DBG:=
OS:=$(shell uname)
FLAGS=
GIT_VERSION := $(shell git describe --abbrev=4 --always)

OPT_MINUS_OPENMP= -O3 -funroll-loops\
	  -pipe -fno-strict-aliasing -DUSE_PDQSORT \
	-DNOT_THREADSAFE -mpopcnt \
	$(FLAGS) $(EXTRA) #-flto

OPT=$(OPT_MINUS_OPENMP) # -lgomp /* sometimes needed */-lomp /* for clang */
ifneq (,$(findstring clang++,$(CXX)))
	OPT+=-fopenmp -lomp
else
	OPT+=-fopenmp
endif
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++14 $(WARNINGS) \
	 -DDASHING_VERSION=\"$(GIT_VERSION)\"  -fdiagnostics-color=always
CXXFLAGS_MINUS_OPENMP=$(OPT_MINUS_OPENMP) $(XXFLAGS) -std=c++1z $(WARNINGS) -Wno-cast-align -Wno-gnu-zero-variadic-macro-arguments \
	 -DDASHING_VERSION=\"$(GIT_VERSION)\"
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS)
#LIB=
LIB=-ldl
LD=-L. $(EXTRA_LD)


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
ZW_OBJS=$(patsubst %.c,%.o,bonsai/zstd/zlibWrapper/gzclose.c  bonsai/zstd/zlibWrapper/gzlib.c  bonsai/zstd/zlibWrapper/gzread.c  bonsai/zstd/zlibWrapper/gzwrite.c  bonsai/zstd/zlibWrapper/zstd_zlibwrapper.c) libzstd.a
ALL_ZOBJS=$(ZOBJS) $(ZW_OBJS) bonsai/clhash.o bonsai/klib/kthread.o bonsai/zlib/libz.a
INCLUDE=-Ibonsai/clhash/include -I.  -Ibonsai/zlib -Ibonsai/hll/libpopcnt -Iinclude -Ibonsai/circularqueue $(ZSTD_INCLUDE) $(INCPLUS) -Ibonsai/hll/vec -Ibonsai -Ibonsai/include/ \
    -Ibonsai/hll/include -Ibonsai/hll -Ibonsai/hll/include/sketch -Ibonsai/hll/vec

EX=$(patsubst src/%.cpp,%,$(wildcard src/*.cpp))
D_EX=$(patsubst src/%.cpp,%_d,$(wildcard src/*.cpp))


all: dashing cardcmp

d: $(D_EX)

libzstd.a:
	+ls libzstd.a >/dev/null || (cd bonsai && $(MAKE) libzstd.a && cp libzstd.a ../)

bonsai/klib/kstring.o:
	+cd bonsai && $(MAKE) klib/kstring.o && cd .. && \
	cd bonsai && $(MAKE) klib/kthread.o && cd ..

bonsai/clhash.o:
	+cd bonsai && $(MAKE) clhash.o && cd ..

OBJ=bonsai/klib/kstring.o bonsai/klib/kthread.o bonsai/clhash.o

DEPS=bonsai/hll/include/sketch/cbf.h bonsai/hll/include/sketch/bf.h bonsai/hll/include/sketch/hll.h bonsai/hll/include/sketch/hk.h bonsai/hll/include/sketch/ccm.h bonsai/hll/include/sketch/bbmh.h

test/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

test/%.zo: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(ZCOMPILE_FLAGS)

%.o: %.cpp libzstd.a libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) libz.a -DNDEBUG -c $< -o $@ $(LIB) -march=native

%.do: %.cpp libzstd.a libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) libz.a -c $< -o $@ $(LIB) -march=native

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB) -march=native

%.do: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB) -march=native


zobj: $(ALL_ZOBJS)

STATIC_GOMP?=$(shell $(CXX) --print-file-name=libgomp.a)

bonsai/zstd/zlibWrapper/%.c:
	cd bonsai && $(MAKE) libzstd.a

dashing.a: src/dashing.o  libzstd.a bonsai/klib/kthread.o bonsai/clhash.o $(ALL_ZOBJS)
	ar r dashing.a src/dashing.o  libzstd.a $(ALL_ZOBJS) bonsai/klib/kthread.o bonsai/clhash.o

BACKUPOBJ=src/main.o src/union.o src/hllmain.o src/mkdistmain.o src/finalizers.o src/cardests.o src/distmain.o src/construct.o src/flatten_all.o \
        $(patsubst %.cpp,%.o,$(wildcard src/sketchcmp*.cpp) $(wildcard src/sketchcore*.cpp)) src/background.o src/panel.o src/cardmain.o src/distbyseq.cpp
CARDCMPO=src/cardmain.o src/finalizers.o src/dashing.o
DASHINGSRC=src/main.cpp src/union.cpp src/hllmain.cpp src/mkdistmain.cpp src/finalizers.cpp src/cardests.cpp src/distmain.cpp src/construct.cpp src/flatten_all.cpp \
        $(wildcard src/sketchcmp*.cpp) $(wildcard src/sketchcore*.cpp) src/background.cpp src/panel.cpp src/cardmain.cpp src/distbyseq.cpp


bonsai/zstd/zlibWrapper/%.o: bonsai/zstd/zlibWrapper/%.c
	cd bonsai && $(MAKE) libzstd.a && cd zstd && $(MAKE) lib  && cd zlibWrapper && $(MAKE) $(notdir $@)

%: src/%.cpp $(ALL_ZOBJS) $(DEPS) libzstd.a $(DASHING_OBJ) $(wildcard src/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS)  -O3 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -march=native -DNDEBUG

dashing-ar: src/main.o dashing.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) dashing.a -O3 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -DNDEBUG

%: src/%.o dashing.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) dashing.a -O3 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -march=native -DNDEBUG

cardcmp: src/cardcmp.o $(ALL_ZOBJS) $(DEPS)  libzstd.a $(CARDCMPO)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -O3 $< -o $@ $(CARDCMPO) $(ZCOMPILE_FLAGS) $(LIB) -march=native -DNDEBUG -lz

dashing: src/dashing.o $(ALL_ZOBJS) $(DEPS)  libzstd.a $(BACKUPOBJ)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) $(BACKUPOBJ)  -O3 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -march=native -DNDEBUG -lz

dashing_d: $(ALL_ZOBJS) $(DEPS) libzstd.a $(DASHINGSRC)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) $(DASHINGSRC)  -O1 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -march=native src/dashing.cpp # -DNDEBUG

%0: src/%.o $(ALL_ZOBJS) $(DEPS)  libzstd.a src/main.o
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) src/main.o  -O0 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -march=native

src/%.o: src/%.cpp $(DEPS)  libzstd.a $(wildcard src/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c -O3 $<  -o $@ $(ZCOMPILE_FLAGS) $(LIB) -DNDEBUG -march=native

sparse%: src/%.cpp $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -DUSE_SPARSE -O3 $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -DNDEBUG -march=native

%_d: src/%.cpp $(ALL_ZOBJS) $(DEPS) src/main.o
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g src/main.o \
    $< -o $@ $(ZCOMPILE_FLAGS) $(LIB) -O1 -lz # -fsanitize=undefined -fsanitize=address

dashing_256: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DASHINGSRC) $(ALL_ZOBJS) -march=native -mno-avx512dq -mno-avx512vl -mno-avx512bw -mavx -mavx2 -msse4.1 -msse2 -DNDEBUG src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB)

dashing_512: $(DASHINSRC) $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DASHINGSRC) $(ALL_ZOBJS) -mavx512dq -mavx512vl -mavx2 -msse4.1 -msse2 -DNDEBUG src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB)

dashing_512bw: $(ALL_ZOBJS) $(DEPS) $(DASHINGSRC)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DASHINGSRC) $(ALL_ZOBJS) -mavx512dq -mavx512vl -mavx512bw -mavx2 -msse4.1 -msse2 -DNDEBUG src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB)

dashing_128: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DASHINGSRC) $(ALL_ZOBJS) -mno-avx512dq -mno-avx512vl -mno-avx -mno-avx512bw -mno-avx2 -msse2 -msse4.1 -DNDEBUG src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB)

libgomp.a:
	ln -sf $(STATIC_GOMP)

dashing_s: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS) libgomp.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -static-libstdc++ -static-libgcc -DNDEBUG $(DASHINGSRC) src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB) -flto
dashing_s128: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)  libgomp.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mno-avx512dq -mno-avx512vl -mno-avx512bw -mno-avx -mno-avx2 -msse2 -msse4.1 -static-libstdc++ -static-libgcc -flto \
    libgomp.a libzstd.a \
		-DNDEBUG $(DASHINGSRC) src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB) -ldl -lz -lzstd
dashing_s256: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)  libgomp.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mno-avx512dq -mno-avx512vl -mno-avx512bw -mavx2 -msse2 -static-libstdc++ -static-libgcc -flto \
    libgomp.a libzstd.a \
	-DNDEBUG $(DASHINGSRC) src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB) -ldl -lz -lzstd
dashing_s512: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)  libgomp.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -mavx512dq -mavx512vl -mavx512bw -static-libstdc++ -static-libgcc -flto  \
    libgomp.a libzstd.a \
		-DNDEBUG $(DASHINGSRC) src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB) -ldl -lz -lzstd
dashing_di: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -fno-inline -O1 $(DASHINGSRC) src/dashing.cpp -o $@ $(ZCOMPILE_FLAGS) $(LIB)
dashing_pg: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -pg -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)
dashing_pgi: $(DASHINGSRC) $(ALL_ZOBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -g -pg -fno-inline -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

bonsai/zlib/libz.a:
	+cd bonsai && $(MAKE) zlib/libz.a && cd ..

PREFIX?=/usr/local

install: dashing
	install -d $(DESTDIR)$(PREFIX)/bin/ && \
	install dashing $(DESTDIR)$(PREFIX)/bin/

release_bundle.tgz: dashing_s256 dashing_s512 dashing_s128
	rm -fr release_bundle release_bundle.tgz && mkdir release_bundle && cp dashing_s256 dashing_s512 dashing_s128 release_bundle \
	&& tar -xvzf release_bundle.tgz release_bundle

ALLSIMDS=dashing_s128 dashing_s256 dashing_s512

linux_release:
	+rm -f $(ALLSIMDS) rm release/linux/dashing_s* && \
		$(MAKE) $(ALLSIMDS) && \
		cp $(ALLSIMDS) release/linux && \
		cd release/linux && gzip -f9 $(ALLSIMDS)
osx_release:
	+rm -f $(ALLSIMDS) && \
		$(MAKE) $(ALLSIMDS) && \
        rm -f release/osx/dashing_s* && \
		cp $(ALLSIMDS) release/osx && \
		cd release/osx && \
		gzip -f9 $(ALLSIMDS)
clean:
	rm -f $(EX) $(D_EX) libzstd.a bonsai/clhash.o clhash.o \
	bonsai/klib/kthread.o bonsai/klib/kstring.o libgomp.a \
	&& cd bonsai/zstd && $(MAKE) clean && cd ../zlib && $(MAKE) clean && cd ../.. \
	&& rm -f libz.* && rm -f dashing.a
mostlyclean: clean
sparse: readfilt sparsereadfilt
