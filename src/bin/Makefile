# $@ : name of target
# $< : first dependence
# $? : all dependencies vencidas
# $^ : all dependencias (no duplicates)
# $+ : all dependencies (incl duplicates)
# $* : targets with patterns (only matching part)

###########################################################
# Usage:
# make target (for targets in Makefile or Makefile.tests, 
#                     or Makefile.user.USER, if included)
#
# make WITHMKL=1 target (using mkl libraries, not in the cluster)
# make DISK=1 target (uses temporary disk storage for MPSs)
# maxe DOUBF=1 foldedI tries the new implementation of foldedIsing
# make INTELC=1 target compiles everything with intel compiler (cluster only)
# In the cluster, make targets the parallel version of intel mkl. To
# compile with the sequential version, call
# make NOOPENMP=1 

###########################################################
# Editing this file: (from 31.1.2013)
# 
# This is the main Makefile. It should contain the general
# features regarding architectures, but the instructions for
# particular executables are now moved to Makefile.user.USER,
# which is (or can be) called from this file.
# There is a Makefile.tests always included, which can be used 
# for tests.
###########################################################

###########################################################
# User name, to determine the particular Makefile to be 
# included
export USER:=$(shell whoami)
USER_SPECIFIC_MAKEFILE=Makefile.user.$(USER)
#
###########################################################

###########################################################
# Compiler options

# 1) Compiler itself (may be changed by Makefile.arch, 
# depending on the system)

export CXX=g++

# 2) Needed option (to ensure the new program version is 
# implemented). TODO: remove this from everything but some tests.

export DEFINES:=$(DEFINES) -DNEWIMPLEM=1 

# 3) Default now is using PRIMME for the eigensolver, but if the
# system dos not support it, this can be cancelled by the Makefile.user 
# file, just by adding 
# DEFINES+= -UUSING_PRIMME

export DEFINES:=$(DEFINES) -DUSING_PRIMME=1 

# 4) Optimized code. Only when compiling for production, 
# uncomment the export instruction at the end of this comment, 
# or better add  
# CXXLABELS+= -O
# to the Makefile.user file
###  CXXLABELS+= -O3 -fsanitize=address # works for OSX (4.2.2021)


# 5) The following debugging options can be included in 
# your Makefile.user file.
# Include checking of all dims in mwArray operations
# DEFINES += -DCHECKDIMS=1
# GDB support
# CXXLABELS += -ggdb 
# More details in valgrind output
# CXXLABELS += -ggdb -fno-inline
# For Valgrind cachegrind
# CXXLABELS += -g -fsanitize=address -O

###########################################################

###########################################################
# Linker

export LD=g++ 

# Any of the following debugging/optimizing options can 
# be added to the Makefile.user file
# GDB debugging
# LDLABELS += -ggdb 
# Valgrind cachegrind
# LDLABELS += -g -fsanitize=address -O
# Production code
# LDLABELS += -O3 -fsanitize=address

###########################################################

###########################################################
#  System dependent instructions (libraries, etc)
#  
include Makefile.arch
###########################################################

###########################################################
# For the version that uses disk storage, call make as: 
# make DISK=1
ifeq ($(DISK),)
$(warning "Without disk storage")
else
$(warning "With disk storage")
DEFINES += -DDISKSTRG=1
endif
###########################################################

###########################################################
# For the version that uses disk storage, call make as: 
# make DOUBF=1
ifeq ($(DOUBF),)
$(warning "Without disk storage")
else
$(warning "Without disk storage")
DEFINES += -DTESTINGDOUB=1 -DCHECKDIMS=1
endif
###########################################################

###########################################################
# For a version that prints some more debugging info, call make as: 
# make TEST=1
ifneq ($(TEST),)
$(warning "With extra debugging functions")
DEFINES += -DTESTFUN=1
endif
###########################################################

###########################################################
# Directions to search: Where the (common) sources are
export DIRECTORIES = core/headers models/headers wrapper wrapper/multipliers core/cpp models/cpp tests programs programs/ising programs/heisenberg
# \ fortran

# What are the (really basic ones) files that I will need to compile
WRAPPER = complex.cpp wrapper.cpp svd.cpp Indices.cpp mwArray.cpp \
		eig.cpp BasicMultiplier.cpp \
	eigs.cpp lsd.cpp lss.cpp lu.cpp qr.cpp misc.cpp \
	TensorMultiplier.cpp TensorMultiplierProj.cpp \
	TensorJoinedMultiplier.cpp \
	TensorMultiplierProjMulti.cpp TensorMultiplierHermitian.cpp \
	TensorMultiplierOffset.cpp TensorMultiplierHermitianOffset.cpp \
	LinearCombinationMultiplier.cpp Transposer.cpp 

WRAPPEROBJ = $(filter %.o,$(WRAPPER:.cpp=.o))  

SRCFILES= Properties.cpp Site.cpp FileSite.cpp MPS.cpp Operator.cpp MPO.cpp \
	JoinedOperator.cpp TmpManagerExc.cpp \
	TmpManager.cpp Contractor.cpp \
	FoldedOperator.cpp DoubleOperator.cpp VariableOperator.cpp 

# Additional source files can be added with an instruction
# SRCFILES+=SomeFile.cpp 
# in the Makefile.user.{username} file
# And the same for include directories, etc






SRCOBJ = $(filter %.o,$(SRCFILES:.cpp=.o))  

SRCDIRS = $(addprefix ../,$(DIRECTORIES))

export VPATH=$(subst :, ,$(SRCDIRS))

INCLUDES+=$(addprefix -I,$(SRCDIRS)) 

OBJFILES = $(WRAPPEROBJ) $(SRCOBJ)


# Location for the headers for PRIMME eigensolver.
# The default one now is this, for v1, but the user Makefile
# can define a _v2 directory, instead
export PRIMMEH_DIR=../wrapper/primme_headers

include Makefile.tests

include $(USER_SPECIFIC_MAKEFILE)

INCLUDES+= -I. 
INCLUDES+= -I$(PRIMMEH_DIR)


###########################################################
# Actual rules to build targets
###########################################################

.phony: clean cleanc cleanall 


%.o: %.cpp 
	@(echo "  --> Compiling $@ -- ")
	$(CXX) $(CXXLABELS) $(INCLUDES) $(DEFINES) -c  $< 


$(EXEC):
	@(echo "  --> Linking $@ -- ")
	$(LD) $(LDLABELS) $(OBJFILES) $(LINKDIRS) $(LIBS) -o $@

###########################################################
# Cleanup rules
cleanall: cleanc cleantests
	rm -f Makefile.d
	rm -f $(EXEC) *.tgz

clean: cleanc

cleanc:
	rm -f *.o
