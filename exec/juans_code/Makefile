.PHONY: default clean all

ROOTFLAGS=$(shell root-config --cflags --evelibs)
CXXFLAGS+=-std=c++17 #-Wall -Werror -Wextra -Wpedantic -O3
LDFLAGS+=-lm
CXX=g++

default: all

all: htt_code

htt_code: htt_code.o
	$(CXX) $(LDFLAGS) $(ROOTFLAGS) "$^" -o "$@"

%.o: %.c
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c "$^" -o "$@"
	
clean:
	$(RM) htt_code htt_code.o
