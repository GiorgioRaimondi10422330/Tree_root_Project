# Including generic Makefile.inc
include ../../../Makefile.inc

LIB= $(TIME_ROOT_LIBRARY)
LIBDIR= $(TIME_ROOT_LIBRARY_DIR)
INC_DIR= $(TIME_ROOT_INCLUDE)

.PHONY= help all main run clean distclean
INC_DIR= $(TIME_ROOT_INCLUDE)
.DEFAUL_GOAL: all

all: main
	@echo 
	@echo "1_LinearSingleBranch test case compiled"
	@echo

main.o: main.cpp
	$(CXX) $(STDFLAG) -D=M3D1D_VERBOSE_ $(INC_DIR) $(OPTFLAGS) $(WARN)  -o main.o -c main.cpp

main: main.o
	$(CXX) $(LIBDIR) -o main main.o  $(LIB)

run: main
	./main input.param

clean:
	$(RM) main.o *~

distclean: clean
	$(RM) main

fastclean: 
	$(RM) main


