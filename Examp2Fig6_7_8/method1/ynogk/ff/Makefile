
#FOBJECTS = $(patsubst %.F90,%.o,$(patsubst %.f90,%.o,$(wildcard *90)))
#defaut: disk

#FC = mpifort
RM = rm -f

#OBJECTS = $(FOBJECTS)
AOBJECTS = ConstantsUnit.o \
	RootsUnit.o \
	ellfunctionsUnit.o \
	BlcoordinateUnit.o \
#	Pemfind.o	\
	ObseverEmitterUnit.o 
#	PemrootsUnit.o	\
#	PmFindingUnit.o	\
#	Pemfind.o	\
# include the dependencies
#include makefile.dep
 
# we should make all the objects and the archive them with AR
$(LIB_DIR)/libynogk.a: $(AOBJECTS) 
	$(AR) $(ARFLAGS) $@ $^

#$(LIB_DIR)/libynogk.so: $(AOBJECTS)
#	gcc -shared  $^ -o $@
  
%.mod: %.o %.f90 
%.o: %.f90
	$(FC) $(FFLAGS) -c -O2 -CB $<

#parameters.o: $(SRC_DIR)/parameters.f90
#	$(FC) $(FFLAGS) -c -O2 -CB $<

#disk: $(AOBJECTS)
#	$(FC) $^ -o $@ 

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~ *.a *.so
