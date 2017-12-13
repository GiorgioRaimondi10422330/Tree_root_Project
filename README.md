# Computational methods for soil-root water uptake
#### *Politecnico di Milano* (ITALY)

**Author** : Giorgio Raimondi

**Mailto** : <giorgio.raimondi.93@hotmail.com>

**Date**   : December 2017

-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE
- `doc/`     : Code documentation

- `include/` : General include files

- `lib/`     : Main library (to be generated)

- `src/`     : Example sources
	- TreeRoot_Time_solver/src/1_LinearSingleBranch: Test Case for single Branch.
	- TreeRoot_Time_solver/src/2_Field_test_case: Test Case for steady uncoupled soil root problem.
	- TreeRoot_Time_solver/src/3_Initial_P0: Test Case for unsteady uncoupled soil root problem.
	- TreeRoot_Time_solver/src/4_Init_No_Root:Test Case for unsteady coupled soil root problem.
	- TreeRoot_Time_solver/src/5_Draining_Root: Test Case for unsteady coupled soil root problem in draining test.
	- TreeRoot_Time_solver/src/6_Competing: Test Case for unsteady coupled soil root problem with root competition.

- `Makefile.inc`: Specify the variable GETFEM_PREFIX for GetFEM++ library

- `Doxyfile` : Instruction to build the code documentation

- `Makefile` : Instruction to install the project (see INSTALLATION)

## INSTALLATION
### Prerequisites

You need the open source finite element library "GetFEM++"

See <http://download.gna.org/getfem/html/homepage>

Version >= 4.2 is preferible

You must modify the path to the GetFEM library in `Makefile.inc`:
``` 
GETFEM_PREFIX=/home/.../path/to/.../getfem
``` 

BEWARE: 
Recall to add the library path to LD_LIBRARY_PATH. Example:
```
$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib
```

======================

### Installation
Build the whole project with:
``` 
$ make install
``` 
It first build the (static) library "libproblem3d1d" by calling
the Makefile in `include/`:
``` 
$ make -C TreeRoot_Time_solver/include/
``` 
Then, it calls the inner makefiles provided for all examples.

It is also possible to build a single example, e.g. "1_LinearSingleBranch", with:
``` 
$(MAKE) -C TreeRoot_Time_solver library

$(MAKE) -C TreeRoot_Time_solver/src/1_LinearSingleBranch
``` 

BEWARE: 
If you want non-optimized program type:
``` 
$ make DEBUG=yes 
``` 
By defaul DEBUG=no.

The following macro are defined and exported
``` 
CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include

CXXFLAGS=-std=c++11 -D=M3D1D_VERBOSE_

OPTFLAGS=-O3 -DNDEBUG -march=native

LDFLAGS=-L$(GETFEM_PREFIX)/lib

LIBRARIES=-lgetfem
``` 
Recall that any macro may be overrlued by specifying it when calling 
make. Example: 
``` 
$ make CXXFLAGS+=-DSOMETHING OPTFLAGS=-g
``` 

======================

### Documentation
The documentation is produced by doxygen. The file Doxyfile contains 
the common doxygen configuration for all examples.
Build the code documentation with:
``` 
$ make pdf
``` 
It first fills doc/ with code documentation ($ make doc) and then compile
the .tex files to produce a portable file ($ pdflatex doc/latex/refman.tex).
You can visualize the documentation with
``` 
$ evince TreeRoot_Time_solver/doc/latex/refman.pdf
``` 

## MAKE OPTIONS
All examples are provided with a Makefile which accepts the following
options:
-  all       : print help
-  install   : install the library and compiles the example
-  clean     : as it says
-  distclean : clean and also deletes temporary file and local doc directory
Being "all" the first target of the makefile, to compile the examples is
sufficient to type make. 
In addition the external Makefile (./Makefile) has the following options:
-  doc       : produces the documentation (html, tex)
-  pdf       : produces a guide in portable format
-  run1      : run the first example
-  run2      : run the second example
-  run3      : run the third example
-  run4      : run the fourth example
-  run5      : run the fifth example
-  run6      : run the sixth example

## RUN EXAMPLES
To run a specif example simply execute
``` 
$ make run3
``` 
where the number indicates the number of the example.
Alternatively go to the related subdirectory
``` 
$ cd TreeRoot_Time_solver/src/3_Initial_P0
``` 
Build the program
``` 
$ make
``` 
Execute the program with specific input
``` 
$ ./M3D1D input.param
``` 
Each program contains the file input.param that specifies 

- Some flags to identify the particular example
  -  TEST_PARAM = 1  # import parameters in dimensionless form
  -  VTK_EXPORT = 1  # export results in vtk format
  -  ...

- The mesh
  - For the 3D mesh you can either provide instruction to build a simple
  regular mesh (TEST_GEOMETRY = 1) or the absolute path to import a mesh
  pre-built with Gmsh (.gmsh)
  - For the 1D mesh specify the path to the file of points (.pts). All
  examples come with a possible pts file

- GetFEM++ descriptors (FEM, ...)

- Problem parameters (dimensional or dimensionless)

- Boundary conditions. You can choose an arbitrary combination of
  Dirichlet-type conditions on pt and/or Robin-type conditions
  on the flux, namely:

  % Faces:   x=0  x=L  y=0  y=L  z=0  z=L

  % BC labels (DIR / MIX)

  BClabel = 'DIR  DIR  DIR  DIR  DIR  DIR'

  % BC values

  BCvalue = '0.0  0.0  0.0  0.0  0.0  0.0'
  
  
BEWARE: All paths in file param must be ABSOLUTE

##  DEV ENVIRONMENT
OS         : Ubuntu 14.04 LTS 64-bit

Processor  : Intel® Core™ i5-2410M CPU @ 2.30GHz × 4 

Compiler   : g++-4.8

GetFEM lib : 5.0

