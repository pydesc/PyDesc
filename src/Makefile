LIB_PATH=$(HOME)/local/lib
INCLUDE_PATH=./include 

#flags for different builds
BASE_CFLAGS:=-Wall -Wno-unknown-pragmas -std=gnu99 -mmmx -msse -msse3 -D_GNU_SOURCE $(CUSTOM_FLAGS) -I$(INCLUDE_PATH) -I$(HOME)/local/include -I/usr/local/cuda/include/ -Iinclude/ -m64 -fopenmp # -fnested-functions#  -Wextra
BASE_LDFLAGS:= -L$(LIB_PATH)# -fopenmp
BASE_LDLIBS:=

ifeq ($(shell uname), Darwin)
	BASE_CFLAGS+= -I/opt/local/include
	BASE_LDFLAGS+= -L/opt/local/lib
	GCC=gcc-mp-4.8
	GPP=g++-mp-4.8
	LIBTOOL=glibtool
endif
ifeq ($(shell uname), Linux)
	GCC=gcc
	GPP=g++
	LIBTOOL=libtool
endif

DEBUG_CFLAGS:= -g -DSINGLE -pg #-finstrument-functions
DEBUG_LDFLAGS:= -pg -g
DEBUG_LDLIBS:=
FINAL_CFLAGS:=-ffast-math  -O3# -pg -g  #-mcpu=pentium3
FINAL_LDFLAGS:=# -pg
FINAL_LDLIBS:=

CFLAGS=$(BASE_CFLAGS) $(MODES) $(FINAL_CFLAGS) -g  -pg -DCOMPDESC_DEBUG
LDFLAGS=$(BASE_LDFLAGS) $(FINAL_LDFLAGS) -g -pg
LDLIBS=$(BASE_LDLIBS) $(FINAL_LDLIBS)


UTIL_SOURCES=util/arrays.c util/bitops.c util/fib.c util/gen_avl.c util/maps.c util/random.c util/timeout.c

CYDESC_SOURCES=cstructures.c
OVERFIT_SOURCES=overfit.c
FITDESC_SOURCES=fitdesc.c
COMPDESC_SOURCES=$(patsubst %,compdesc/%, compdesc.c rmsd.c alignment.c 01_prelim.c 02_components.c 03_build.c 04_prune.c)

SOURCES= $(CYDESC_SOURCES) $(OVERFIT_SOURCES) $(FITDESC_SOURCES) $(COMPDESC_SOURCES) $(UTIL_SOURCES)

%.d: %.c
	$(GCC) -MM $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o \1.lo $@ :,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.lo: %.c
	$(LIBTOOL) --mode=compile --tag=CC $(GCC) -c $(CFLAGS) -o $@ $<


all: ctags libcydesc.la liboverfit.la libfitdesc.la libcompdesc.la

dep: $(SOURCES:.c=.d)


libcydesc.la: $(CYDESC_SOURCES:.c=.lo) $(UTIL_SOURCES:.c=.lo)
	$(LIBTOOL) --mode=link --tag=CC $(GCC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -rpath $(LIB_PATH)

liboverfit.la: $(OVERFIT_SOURCES:.c=.lo)
	$(LIBTOOL) --mode=link --tag=CC $(GCC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -rpath $(LIB_PATH)

libfitdesc.la: $(FITDESC_SOURCES:.c=.lo) $(OVERFIT_SOURCES:.c=.lo)
	$(LIBTOOL) --mode=link --tag=CC $(GCC) $(LDFLAGS) -o $@ $^ $(LDLIBS)  -rpath $(LIB_PATH)

libcompdesc.la: $(COMPDESC_SOURCES:.c=.lo) $(UTIL_SOURCES:.c=.lo) $(OVERFIT_SOURCES:.c=.lo)
	$(LIBTOOL) --mode=link --tag=CC $(GCC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -rpath $(LIB_PATH)


clean:
	rm -fr *.o *.d .libs *.lo *.la compdesc/*.d compdesc/*.o compdesc/*.lo util/*.d util/*.o util/*.lo  tags ctags

install-dev:
	install -p cstructures.h $(INCLUDE_PATH)/cstructures.h

install: install-dev all
	$(LIBTOOL) --mode=install install libcydesc.la $(LIB_PATH)/libcydesc.la

ctags:
	ctags $(SOURCES) $(CXX_SOURCES) $(CUDA_SOURCES)


include $(SOURCES:.c=.d)


