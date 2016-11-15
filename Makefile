# ---------------------------------------------------------------------------
#
#	Makefile for PSI
#
#	Devon Powell 
#	(with some useful bits from Jonathan Zrake)
#
#	Do not modify this file!
#	All user-set options live in Makefile.in
#
#	usage: make
#
# ---------------------------------------------------------------------------



# if there is no Makefile.in then use the template
ifneq ($(strip $(MAKEFILE_IN)),)
# use value of MAKEFILE_IN if provided on the command line
else ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = Makefile.in
else
MAKEFILE_IN = Makefile.in.template
endif
include $(MAKEFILE_IN)

# Source files
CC_SRC = driver.c #utils.c

COMMON = HDF_IO.hh Makefile utils.h cur3d.h
OBJ = $(NV_SRC:.cu=.o) $(CC_SRC:.c=.o)
EXE = cur3d

# base libraries and include dirs
INC = -I./ 
LIB =
LDFLAGS = -lm 

# Set up CUDA dependencies
ifdef CUDA_HOME
COMMON += cur3d.h
NV_SRC = cur3d.cu
INC += -I$(CUDA_HOME)/include 
LIB +=  -L$(CUDA_HOME)/lib64
LDFLAGS += -lcuda -lcudart -lcublas -lcufft
CFLAGS += -DCUDA
#CC = $(NVCC)
endif

# Set up HDF5 dependencies
#INC += -I$(HDF5_HOME)/include
#LIB += -L$(HDF5_HOME)/lib
#LDFLAGS += -lhdf5

# R3D
INC += -I$(R3D_HOME)
LIB += -L$(R3D_HOME)
LDFLAGS += -lr3d


all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LIB) $(OBJ) -o $@ $(LDFLAGS) $(CFLAGS)

.c.o: $(COMMON)
	$(CC) $(INC) $(DEF) $(CFLAGS) -c $< -o $@

# WHY DOESN'T THIS RULE WORK??
#.cu.o: 
	#$(NVCC) $(INC) $(DEF) $(NV_FLAGS) -c $< -o $@

cur3d.o: cur3d.cu
	$(NVCC) $(INC) $(DEF) $(NV_FLAGS) -c $< -o $@
	
clean:
	rm -f $(OBJ) *~ core $(EXE)
