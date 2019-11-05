# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?=/home/ThirdEye87/amrex-master

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE

COMP      = gnu
DIM       = 3

USE_WAVE_EQ ?= FALSE
USE_MAXWELL ?= FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(USE_WAVE_EQ), TRUE)
  DEFINES += -DWAVEEQ
else ifeq ($(USE_MAXWELL), TRUE)
  DEFINES += -DMAXWELL
endif

include ../Source/Make.package

VPATH_LOCATIONS   += ../Source
INCLUDE_LOCATIONS += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
