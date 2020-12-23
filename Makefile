### Makefile for compilation of analysis macros
### to be copied into the analysis directory or used with make -f <Makefile location>
### and perhaps edited to correct the FWK_BASE_DIR e.g. make example_gen_ttbar FWK_BASE_DIR=<fwk base dir> -f <Makefile location>

FWK_BASE_DIR = ..
CXXFLAGS = $(shell root-config --cflags --evelibs) -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare

%: %.cc $(FWK_BASE_DIR)
	@ g++ $(CXXFLAGS) -I $(FWK_BASE_DIR)/plugins/ -I $(FWK_BASE_DIR)/src/ $< -o $@
