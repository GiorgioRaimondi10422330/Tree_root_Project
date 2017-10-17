# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2014-2015
#
#                    Copyright D. Notaro 2015
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Domenico Notaro <domenico.not@gmail.com>
#   DATE        : November 2015
# ====================================================================

.PHONY: all doc clean distclean library help run fast run3 fast3

all: library
	$(MAKE) -C src/1_LinearSingleBranch
	$(MAKE) -C src/2_Field_test_case
	$(MAKE) -C src/3_Complex_Network

help:
	@echo "make all: build the library and complile the test case"
	@echo "make help: prints help comands"
	@echo "make library: build the static library"
	@echo "make clean: clean all objects"
	@echo "make distclean: clean all object, libraries, doc file and executables"



library: 
	$(MAKE) -C include

run: all
	$(MAKE) -C src/2_Field_test_case run

run3: all
	$(MAKE) -C src/3_Complex_Network run
	
clean:
	$(RM) -r *~ *.log
	$(MAKE) -C include clean
	$(MAKE) -C src/1_LinearSingleBranch clean
	$(MAKE) -C src/2_Field_test_case clean
	$(MAKE) -C src/3_Complex_Network clean

distclean: clean
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_LinearSingleBranch distclean
	$(MAKE) -C src/2_Field_test_case distclean
	$(MAKE) -C src/3_Complex_Network distclean
fast:	
	$(MAKE) -C include fastclean
	$(MAKE) -C src/2_Field_test_case fastclean
	$(MAKE) -C include
	$(MAKE) -C src/2_Field_test_case
	$(MAKE) -C src/2_Field_test_case run
	
fast3:
	$(MAKE) -C include fastclean
	$(MAKE) -C src/3_Complex_Network fastclean
	$(MAKE) -C include
	$(MAKE) -C src/3_Complex_Network
	$(MAKE) -C src/3_Complex_Network run
