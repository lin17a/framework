# Usage:
# make        # compile all binary
# # make clean  # remove ALL binaries and objects
#
.PHONY = all clean

CC = g++                        # compiler to use

ROOTFLAGS=$(shell root-config --cflags --evelibs)
FLAGS = ${ROOTFLAGS} -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I /nfs/dust/cms/user/meyerlin/ba/framework/plugins/ -I /nfs/dust/cms/user/meyerlin/ba/framework/src/ -ggdb 
# -L/afs/desy.de/user/m/meyerlin/ba/framework/lib/boost_1_74_0 
FILENAME=example_gen_ttbar_mit_mir

# SRCS := $(wildcard *.c)
SRCS=${FILENAME}'.c'
BINS := $(SRCS:%.c=%)

# @echo ${BINS}

all: ${FILENAME}

%: %.cc
# @echo ${ROOTFLAGS}
	${CC} ${FLAGS} $< -o $@

# %.o: %.cc
# 	${CC} -c $<

clean:
	rm -rvf *.o ${BINS}
