# Project Path (this file path)
#PROJECT_DIR = /u/raimondi/Tree_root_Project
PROJECT_DIR = /media/More_SPACE/Tree_root_Project
# Path to GetFEM library
#GETFEM_LIB_DIR     = $(mkGetfemLib)
GETFEM_LIB_DIR     = /opt/getfem/lib
# Path to GetFEM includes
#GETFEM_INC_DIR     = $(mkGetfemInc)
GETFEM_INC_DIR     = /opt/getfem/include



### NOW DO NOT EDIT ANYTHING ###




# Path to time root include
TIME_ROOT_INC_DIR= $(PROJECT_DIR)/TreeRoot_Time_solver/include
# Path to time root library
TIME_ROOT_LIB_DIR= $(PROJECT_DIR)/TreeRoot_Time_solver/lib
# Path to Curved M3D1D include
TIME_ROOT_INC_DIR = $(PROJECT_DIR)/TreeRoot_Time_solver/include
# Path to Curved M3D1D library
TIME_ROOT_LIB_DIR = $(PROJECT_DIR)/TreeRoot_Time_solver/lib


# Specific Curved M3D1D Include list TIME
TIME_ROOT_INCLUDE =-I$(TIME_ROOT_INC_DIR) -I$(GETFEM_INC_DIR) 
# Specific Curved M3D1D Library list TIME
TIME_ROOT_LIBRARY =   -lroot_problem3d1d -lgetfem
# Specific Curved M3D1D Library directory TIME
TIME_ROOT_LIBRARY_DIR =  -L$(TIME_ROOT_LIB_DIR) -L$(GETFEM_LIB_DIR)



# Specific Curved M3D1D Include list
TIME_ROOT_INCLUDE =-I$(TIME_ROOT_INC_DIR) -I$(GETFEM_INC_DIR) 
# Specific Curved M3D1D Library list
TIME_ROOT_LIBRARY =   -lroot_problem3d1d -lgetfem
# Specific Curved M3D1D Library directory
TIME_ROOT_LIBRARY_DIR =  -L$(TIME_ROOT_LIB_DIR) -L$(GETFEM_LIB_DIR)

# Compiler and standard for the compiler
CXX = g++
STDFLAG = -std=c++11
WARNS = -Wall

ifeq ($(DEBUG),yes)
  OPTFLAG =-g
else
  OPTFLAG =-O3 -march=native -DNDEBUG
endif
