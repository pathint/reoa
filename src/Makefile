ifeq ($(OS),Windows_NT)
    CCFLAGS += -D WIN32
    ifeq ($(PROCESSOR_ARCHITEW6432),AMD64)
        CCFLAGS += -D AMD64
    else
        ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
            CCFLAGS += -D AMD64
        endif
        ifeq ($(PROCESSOR_ARCHITECTURE),x86)
            CCFLAGS += -D IA32
        endif
    endif
else
    UNAME_S := $(shell uname -s)
    prefix := /usr/local
    ifeq ($(UNAME_S),Linux)
        CCFLAGS += -D LINUX
        CC = gcc
        CFLAGS = -O2 -Wall -fopenmp
        INCLUDES = -I. -I/usr/lib64/gcc/x86_64-pc-linux-gnu/4.9.3/include 
    endif
    ifeq ($(UNAME_S),Darwin)
        CCFLAGS += -D OSX
        CC =  gcc-mp-5 
        CFLAGS = -O2 -Wall -fopenmp
        INCLUDES = -I. -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.2.0/include
    endif
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)
        CCFLAGS += -D AMD64
    endif
    ifneq ($(filter %86,$(UNAME_P)),)
        CCFLAGS += -D IA32
    endif
    ifneq ($(filter arm%,$(UNAME_P)),)
        CCFLAGS += -D ARM
    endif
endif


LIBS = -lm -lgomp

SRCS = pcg_basic.c stat.c binomial.c fisher.c hypergeom.c random_sample.c \
       reoa.c main.c


LFLAGS =

OBJS = $(SRCS:.c=.o)

TARGETS := fisher hypergeom binomial random_sample \
           reoa

	
.PHONY: depend clean install

all: $(TARGETS)
	@echo Building all executables

binomial: pcg_basic.o stat.o binomial.o 
	$(CC) $(LFLAGS) $(LIBS) -o $@ $^ 

fisher: pcg_basic.o stat.o fisher.o 
	$(CC) $(LFLAGS) $(LIBS) -o $@ $^ 

hypergeom: pcg_basic.o stat.o hypergeom.o 
	$(CC) $(LFLAGS) $(LIBS) -o $@ $^ 

random_sample: pcg_basic.o stat.o random_sample.o 
	$(CC) $(LFLAGS) $(LIBS) -o $@ $^ 

reoa: pcg_basic.o stat.o reoa.o main.o
	$(CC) $(LFLAGS) $(LIBS) -o $@ $^ 

install: reoa 
	test -d $(prefix) || mkdir $(prefix)
	test -d $(prefix)/bin || mkdir $(prefix)/bin
	install -m 0755 reoa $(prefix)/bin
	install -m 0755 cellcomp $(prefix)/bin
	install -m 0755 onecomp $(prefix)/bin

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(TARGETS)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
