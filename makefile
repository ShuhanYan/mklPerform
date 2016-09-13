#===============================================================================
# Copyright 2001-2016 Intel Corporation All Rights Reserved.
#
# The source code,  information  and material  ("Material") contained  herein is
# owned by Intel Corporation or its  suppliers or licensors,  and  title to such
# Material remains with Intel  Corporation or its  suppliers or  licensors.  The
# Material  contains  proprietary  information  of  Intel or  its suppliers  and
# licensors.  The Material is protected by  worldwide copyright  laws and treaty
# provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
# modified, published,  uploaded, posted, transmitted,  distributed or disclosed
# in any way without Intel's prior express written permission.  No license under
# any patent,  copyright or other  intellectual property rights  in the Material
# is granted to  or  conferred  upon  you,  either   expressly,  by implication,
# inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
# property rights must be express and approved by Intel in writing.
#
# Unless otherwise agreed by Intel in writing,  you may not remove or alter this
# notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
# suppliers or licensors in any way.
#===============================================================================

##  Content:
##      Intel(R) Math Kernel Library VML examples creation and run
##******************************************************************************

help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo "[parallel=parallel_name] [omp=omp_name] [gnu_path=gnu_lib_path]"
	@echo "name           - function name. Please see vml.lst file."
	@echo "compiler_name  - can be gnu, pgi, or intel. Default value is intel."
	@echo "                 Intel(R) C Compiler as default."
	@echo "                 If compiler=gnu then GNU C compiler will be used."
	@echo "                 If compiler=pgi then PGI C compiler will be used."
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "threading_name - can be parallel or sequential. Default value is parallel."
	@echo "parallel_name  - can be intel, pgi (only if compiler=pgi) or gnu (only if compiler=gnu). Default value is intel."
	@echo "omp_name       - can be iomp5 if parallel=intel or"
	@echo "                 can be iomp5 or gomp if parallel=gnu or"
	@echo "                 can be pgmp if parallel=pgi."
	@echo "                 Default value is iomp5."
	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
	@echo "                 Default value is /usr/lib."

##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=vssin - build by Intel(R) C compiler (as default) and
##                             run VSSIN example for 32-bit applications,
##                             static linking
##
## make soia32 compiler=gnu - build by GNU C compiler and run
##                          all examples of MKL for 32-bit applications, dynamic
##                          linking
##
## make libintel64 compiler=gnu - build by GNU C compiler and run all
##                           examples of MKL for Intel(R) 64 processor family
##                           applications, static linking
##
## make sointel64 - build by Intel(R) C compiler (as default) and
##             run all examples of MKL for Intel(R) 64 processor
##             family applications, dynamic linking
##------------------------------------------------------------------------------

include vml.lst

ifndef function
function = $(VML)
endif

ifeq (,$(filter gnu pgi,$(compiler)))
   override compiler=intel
   override parallel=intel
endif

ifneq ($(interface),ilp64)
   override interface=lp64
endif

ifneq ($(threading),sequential)
   override threading=parallel
endif

ifeq (,$(filter gnu pgi,$(parallel)))
   override parallel=intel
   override omp=iomp5
else
   ifeq ($(parallel),gnu)
      ifneq ($(omp),gomp)
         override omp=iomp5
      endif
   else
      override omp=pgmp
   endif
endif

ifndef $(gnu_path)
   gnu_path=/usr/lib
endif

RES = $(addsuffix .res ,$(function))

ifndef MKLROOT
   MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

LOPTS = 
COPTS = 

ifeq ($(interface),ilp64)
   IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART)_ilp64.$(EXT)
   COPTS += -DMKL_ILP64
else
   IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART)_lp64.$(EXT)
endif

ifeq ($(compiler),gnu)
   CC = gcc
   COPTS += -w
   COPTS += $(if $(filter ia32, $(_IA)), -m32, -m64)
   ifeq ($(RES_EXT),so)
      LOPTS = -Wl,--no-as-needed
   endif
endif

ifeq ($(compiler),intel)
   CC=icc
   COPTS += -w
endif

ifeq ($(compiler),pgi)
   CC=pgcc
   COPTS += -Minform=severe -Mnokeepobj
endif

IFACE_COMP_PART=intel

ifeq ($(parallel),gnu)
   IFACE_THREADING_PART=gnu
   GFORTRAN_LIB=-L$(gnu_path) -lgfortran
endif

ifeq ($(parallel),intel)
   IFACE_THREADING_PART=intel
   GFORTRAN_LIB=
endif

ifeq ($(parallel),pgi)
   IFACE_THREADING_PART=pgi
   GFORTRAN_LIB=
endif

ifeq ($(_IA),ia32)
   ifeq ($(compiler),intel)
       SPEC_OPT=-mia32
#This option is required by Intel(R) 11.0 compiler to produce workable binaries for Pentium(R) III.
#If you don't need it, you can remove this option.
   endif
   IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART).$(EXT)
endif

ifeq ($(threading),sequential)
   THREADING_LIB=$(MKL_PATH)/libmkl_sequential.$(EXT)
   OMP_LIB =
   GFORTRAN_LIB=
else
   THREADING_LIB=$(MKL_PATH)/libmkl_$(IFACE_THREADING_PART)_thread.$(EXT)
   ifeq ($(parallel),pgi)
      COPTS += -mp -pgf90libs
   endif
   ifeq ($(omp),gomp)
      OMP_LIB = -L$(gnu_path) -l$(omp)
   else
      ifeq ($(omp),pgmp)
         OMP_LIB =
      else
         OMP_LIB = -L$(CMPLR_PATH) -l$(omp)
      endif
   endif
endif

CORE_LIB=$(MKL_PATH)/libmkl_core.$(EXT)

MKL_LIBS = $(LOPTS) -Wl,--start-group $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)

ifeq ($(_IA),ia32)
   ifeq ($(threading),parallel)
      RES_DIR=_results/$(compiler)_$(threading)_$(parallel)_$(omp)_$(_IA)_$(RES_EXT)$Z
   else
      RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)$Z
   endif
else
   ifeq ($(threading),parallel)
      RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(parallel)_$(omp)_$(_IA)_$(RES_EXT)$Z
   else
      RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)$Z
   endif
endif

libia32 lib32:
	$(MAKE) $(RES) EXT=a _IA=ia32 RES_EXT=lib
soia32 so32:
	$(MAKE) $(RES) EXT=so _IA=ia32 RES_EXT=so
libintel64 libem64t:
	$(MAKE) $(RES) EXT=a _IA=intel64 RES_EXT=lib
sointel64 soem64t:
	$(MAKE) $(RES) EXT=so _IA=intel64 RES_EXT=so

#-------------------------------------------------------------------------------

%.res: ./source/%.c
	mkdir -p ./$(RES_DIR)
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)"; $(CC) $(SPEC_OPT) $(COPTS) -I$(MKLROOT)/include $< $(MKL_LIBS) -lm -ldl -lpthread $(GFORTRAN_LIB) -o $(RES_DIR)/$*.out
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)":$(CMPLR_PATH); $(RES_DIR)/$*.out >$(RES_DIR)/$@
#-------------------------------------------------------------------------------
