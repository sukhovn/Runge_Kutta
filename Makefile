.PHONY: all clean

PROGRAM = runge_test
LIBDIR = Runge

LIBSRC := $(wildcard $(LIBDIR)/*.cpp)
LIB_BINS := $(patsubst $(LIBDIR)/%.cpp, obj/libs/%.o, $(LIBSRC))

CC = g++ # compiler to use
CFLAGS = -c -O3 -I$(LIBDIR) # flags to use at the compliation
LINKERFLAG = -lm

all: prog

prog: $(LIB_BINS) obj/runge.o
	${CC} -o ${PROGRAM} $(LIB_BINS) obj/runge.o ${LINKERFLAG}

obj:
	@echo "Making obj directory..."
	mkdir obj
	mkdir obj/libs

obj/libs/%.o: $(LIBDIR)/%.cpp $(LIBDIR)/%.h obj
	@echo "Compiling library binary..."
	${CC} ${CFLAGS} -o $@ $<

obj/runge.o: runge.cpp
	@echo "Compiling runge.cpp..."
	${CC} ${CFLAGS} -o obj/runge.o runge.cpp

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -rvf Data
	rm -vf ${PROGRAM}