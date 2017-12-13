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

.PHONY: all clean distclean install help doc pdf run1 run2 run3 run4 run 5 run6

all: help

help:
	@echo "make install: build the library and complile the test case"
	@echo "make help: prints help comands"
	@echo "make library: build the static library"
	@echo "make doc: build documentation in .tex"
	@echo "make pdf: build documentation in .pdf"
	@echo "make clean: clean all objects"
	@echo "make distclean: clean all object, libraries, doc file and executables"
	@echo "make run1: run example in src/1_LinearSingleBranch"
	@echo "make run2: run example in src/2_Field_test_case"
	@echo "make run3: run example in src/3_Initial_P0"
	@echo "make run4: run example in src/4_Init_No_Root"
	@echo "make run5: run example in src/5_Draining_Root"
	@echo "make run6: run example in src/6_Competing"


install: 
	$(MAKE) -C TreeRoot_Time_solver

clean:
	$(MAKE) -C TreeRoot_Time_solver clean

distclean: 
	$(MAKE) -C TreeRoot_Time_solver distclean

doc:
	$(MAKE) -C TreeRoot_Time_solver doc


pdf: 
	$(MAKE) -C TreeRoot_Time_solver pdf


run1:
	$(MAKE) -C TreeRoot_Time_solver run1

run2: 
	$(MAKE) -C TreeRoot_Time_solver run2

run3:
	$(MAKE) -C TreeRoot_Time_solver run3

run4:
	$(MAKE) -C TreeRoot_Time_solver run4

run5:
	$(MAKE) -C TreeRoot_Time_solver run5

run6:
	$(MAKE) -C TreeRoot_Time_solver run6

	
	

