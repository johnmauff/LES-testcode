# Makefile
#
# defaults
USE_OPENACC_L = $(shell echo $(USE_OPENACC) | tr A-Z a-z)
USE_ISO_L     = $(shell echo $(USE_ISO) | tr A-Z a-z)
USE_MPI_L     = $(shell echo $(USE_MPI) | tr A-Z a-z)


# ------------------------------
# NVHPC compiler
# ------------------------------
ifeq (${FC},nvfortran)
  OPTS +=  -r8 -O3 -gpu=cc70,cc80 -acc -cuda -cudalib=cufft 
  ifeq (${USE_ISO_L},true)
    OPTS += -stdpar=gpu
  endif
endif
# ------------------------------
# Cray compiler
# ------------------------------
ifeq (${FC},ftn)
  OPTS +=  -s real64 -O2 -Ovector2 -Oscalar2 -Othread2 -h acc
endif
# ------------------------------
# Intel compiler
# ------------------------------
ifeq (${FC},ifort)
  #  module load oneapi/2023.0.0 
  #  module load fftw
  #OPTS += -r8 -assume byterecl -I${NCAR_INC_FFTW}/include/fftw -L${NCAR_LDFLAGS_FFTW}/lib -lfftw3
  OPTS += -r8 -O3 -assume byterecl -lfftw3
endif


ifeq (${USE_MPI_L},true)
    DM           += -DMPI
    override FC  = mpif90
endif


$(info $$OPTS is [${OPTS}])
$(info $$FC is [${FC}])

SRC = test_fftw_mpi.f

OBJS = $(addsuffix .o, $(basename $(SRC)))

.SUFFIXES:
.SUFFIXES: .f .o

all : testx

testx: $(OBJS)
	$(FC) $(OPTS) $(OBJS) -o testx

$(OBJS): %.o: %.f
	$(FC) $(OPTS) -c $<

clean:
	rm -f testx *.o *.a *.mod



