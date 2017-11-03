# Including generic Makefile.inc
include ../../../Makefile.inc

LIB= $(TIME_ROOT_LIBRARY)
LIBDIR= $(TIME_ROOT_LIBRARY_DIR)
INC_DIR= $(TIME_ROOT_INCLUDE)


.PHONY= help all main run clean distclean fastclean
.DEFAUL_GOAL: all

all: main
	@echo 
	@echo "2_Field_test_case compiled"
	@echo

help:
	@echo "make help: ------- Prints this help"
	@echo "make all: -------- Compile and prints help"
	@echo "make main: ------- Compile the program"
	@echo "make run: -------- Run the executable with input.param file"
	@echo "make clean: ------ Cleans all object files"
	@echo "make distclean: -- Cleans all object files and remove executable"

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


