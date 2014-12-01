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
CC_SRC = driver.cpp
NV_SRC = curas.cu

COMMON = HDF_IO.hh
OBJ = $(NV_SRC:.cu=.o) $(CC_SRC:.cpp=.o)
EXE = raster

# base libraries and include dirs
INC = -I./
LIB =
LDFLAGS = 

# Set up CUDA dependencies
INC += -I$(CUDA_HOME)/include 
LIB += -L$(CUDA_HOME)/lib64
LDFLAGS += -lcuda -lcudart -lcublas -lcufft 

# Set up HDF5 dependencies
INC += -I$(HDF5_HOME)/include
LIB += -L$(HDF5_HOME)/lib
LDFLAGS += -lhdf5


all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LIB) $(OBJ) -o $@ $(LDFLAGS) $(CFLAGS)

.cpp.o: $(COMMON)
	$(CC) $(INC) $(DEF) $(CC_FLAGS) -c $< -o $@

# WHY DOESN'T THIS RULE WORK??
#.cu.o: 
	#$(NVCC) $(INC) $(DEF) $(NV_FLAGS) -c $< -o $@

curas.o: curas.cu
	$(NVCC) $(INC) $(DEF) $(NV_FLAGS) -c $< -o $@
	
clean:
	rm -f $(OBJ) *~ core $(EXE)
