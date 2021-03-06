# ====================================================================
#         "Computational methods for soil-root water uptake"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2016-2017
#
#                    Copyright G. Raimondi 2017
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Giorgio Raimondi <giorgio.raimondi.93@hotmail.com>
#   DATE        : November 2015
# ====================================================================

.PHONY: all doc clean distclean library help doc pdf run1 run2 run3 run4 run 5 run6

all: library
	$(MAKE) -C src/1_LinearSingleBranch
	$(MAKE) -C src/2_Field_test_case
	$(MAKE) -C src/3_Initial_P0
	$(MAKE) -C src/4_Init_No_Root
	$(MAKE) -C src/5_Draining_Root
	$(MAKE) -C src/6_Competing

help:
	@echo "make all: build the library and complile the test case"
	@echo "make help: prints help comands"
	@echo "make library: build the static library"
	@echo "make clean: clean all objects"
	@echo "make distclean: clean all object, libraries, doc file and executables"



library: 
	$(MAKE) -C include

clean:
	$(RM) -r *~ *.log
	$(MAKE) -C include clean
	$(MAKE) -C src/1_LinearSingleBranch clean
	$(MAKE) -C src/2_Field_test_case clean
	$(MAKE) -C src/3_Initial_P0 clean
	$(MAKE) -C src/4_Init_No_Root clean
	$(MAKE) -C src/5_Draining_Root clean
	$(MAKE) -C src/6_Competing clean

distclean: clean
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_LinearSingleBranch distclean
	$(MAKE) -C src/2_Field_test_case distclean
	$(MAKE) -C src/3_Initial_P0 distclean
	$(MAKE) -C src/4_Init_No_Root distclean
	$(MAKE) -C src/5_Draining_Root distclean
	$(MAKE) -C src/6_Competing distclean
	$(RM) -r ./doc

doc:
	install -d doc
	doxygen Doxyfile

pdf: doc
	$(MAKE) -C doc/latex pdf

run1:
	$(MAKE) -C src/1_LinearSingleBranch run
	
run2: 
	$(MAKE) -C src/2_Field_test_case run
	
run3: 
	$(MAKE) -C src/3_Initial_P0 run

run4:
	$(MAKE) -C src/4_Init_No_Root run

run5:
	$(MAKE) -C src/5_Draining_Root run

run6:
	$(MAKE) -C src/6_Competing run

	
	

